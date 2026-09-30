# ------------------------------------------------------ #
# Helpers for reporting on DESeq2 outputs produced by
# summarize_PigX_RNAseq (column schema: log2FoldChange.<contrast>,
# padj.<contrast>).
# ------------------------------------------------------ #

classify_DE = function(
    de_table,
    contrast_name,
    lfc_cut  = 0.75,
    padj_cut = 0.05
){
    suppressPackageStartupMessages({
        library(dplyr)
        library(rlang)
    })

    lfc_col  = paste0("log2FoldChange.", contrast_name)
    padj_col = paste0("padj.",           contrast_name)
    diff_col = paste0("diff.",           contrast_name)

    stopifnot(all(c(lfc_col, padj_col) %in% colnames(de_table)))

    de_table %>%
        mutate(!!diff_col := case_when(
            .data[[lfc_col]] >  lfc_cut & .data[[padj_col]] < padj_cut ~ "Up",
            .data[[lfc_col]] < -lfc_cut & .data[[padj_col]] < padj_cut ~ "Down",
            TRUE                                                      ~ "No"
        ))
}


# ------------------------------------------------------ #
# Pseudobulk DE for an snRNA-seq Seurat object.
# Produces a DE table whose column schema matches the
# project bulk DE outputs:
#   log2FoldChange.<contrast> / padj.<contrast> / diff.<contrast>
# so it can be drop-in joined with Riboseq_Differential_Expression().
#
# group_col   : metadata column holding the per-cell label (e.g. cell_type
#               or seurat_clusters). MTJ cluster vs everything-else is the
#               target contrast.
# sample_col  : per-replicate identifier inside the Seurat object — used
#               as the unit of replication for DESeq2. The Petrany / Kim
#               atlases each have it as `sample_id` / `orig.ident`.
# group_a/b   : levels of group_col to contrast (group_a is the numerator).
# ------------------------------------------------------ #
snRNA_pseudobulk_DE = function(
    seurat_obj,
    group_col,
    group_a,
    group_b      = NULL,
    sample_col   = "sample_id",
    block_col    = NULL,          # optional blocking covariate for a paired design (e.g. mouse)
    min_cells    = 30,
    lfc_cut      = 0.75,
    padj_cut     = 0.05,
    contrast_name = NULL
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(DESeq2)
        library(dplyr)
        library(Matrix)
    })

    if(!group_col %in% colnames(seurat_obj@meta.data))
        stop("group_col not in metadata: ", group_col)
    if(!sample_col %in% colnames(seurat_obj@meta.data))
        stop("sample_col not in metadata: ", sample_col)
    if(!is.null(block_col) && !block_col %in% colnames(seurat_obj@meta.data))
        stop("block_col not in metadata: ", block_col)

    md = seurat_obj@meta.data
    md$.gid = if(is.null(group_b))
                  ifelse(md[[group_col]] == group_a, group_a, "rest")
              else
                  ifelse(md[[group_col]] %in% c(group_a, group_b),
                         as.character(md[[group_col]]),
                         NA_character_)
    keep = !is.na(md$.gid)
    md   = md[keep, , drop = FALSE]
    seurat_obj = seurat_obj[, keep]

    cdata_full = data.frame(gid = md$.gid, sample = md[[sample_col]],
                            stringsAsFactors = FALSE)
    if(!is.null(block_col))
        cdata_full$block = as.character(md[[block_col]])
    keys = factor(paste(cdata_full$gid, cdata_full$sample, sep = "__"))

    counts = SeuratObject::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")  # v5: 'slot' is defunct
    ind = sparse.model.matrix(~ 0 + keys)
    colnames(ind) = levels(keys)
    pb = as.matrix(counts %*% ind)

    n_cells = as.integer(table(keys)[colnames(pb)])
    pb       = pb[,       n_cells >= min_cells, drop = FALSE]
    keep_lev = colnames(pb)
    if(ncol(pb) < 2)
        stop("not enough pseudobulk samples after min_cells filter")

    keep_cols = c("gid","sample", if(!is.null(block_col)) "block")
    cdata = unique(cdata_full[, keep_cols, drop = FALSE])
    rownames(cdata) = paste(cdata$gid, cdata$sample, sep = "__")
    cdata = cdata[keep_lev, , drop = FALSE]
    cdata$gid = factor(cdata$gid)
    if(is.null(group_b)){
        cdata$gid = relevel(cdata$gid, ref = "rest")
        ref_level = "rest"
    } else {
        cdata$gid = relevel(cdata$gid, ref = group_b)
        ref_level = group_b
    }

    # Paired design: block on `block_col` (e.g. mouse). Put block first so `gid` is the last
    # term and results(contrast=c("gid",...)) tests the condition effect within block.
    if(!is.null(block_col)){
        cdata$block = factor(cdata$block)
        if(nlevels(cdata$block) < 2)
            stop("block_col has <2 levels after filtering — cannot fit a paired design")
        design = ~ block + gid
    } else {
        design = ~ gid
    }

    dds = DESeqDataSetFromMatrix(countData = pb, colData = cdata, design = design)
    dds = DESeq(dds, quiet = TRUE)
    res = results(dds, contrast = c("gid", group_a, ref_level))

    if(is.null(contrast_name))
        contrast_name = paste0(group_a, "_vs_", ref_level)
    out = as.data.frame(res) %>%
        tibble::rownames_to_column("gene_name") %>%
        dplyr::transmute(
            gene_name,
            !!paste0("baseMean.",       contrast_name) := baseMean,
            !!paste0("log2FoldChange.", contrast_name) := log2FoldChange,
            !!paste0("lfcSE.",          contrast_name) := lfcSE,
            !!paste0("stat.",           contrast_name) := stat,
            !!paste0("pvalue.",         contrast_name) := pvalue,
            !!paste0("padj.",           contrast_name) := padj
        ) %>%
        classify_DE(contrast_name, lfc_cut = lfc_cut, padj_cut = padj_cut)
    out
}


# ------------------------------------------------------ #
# Compute modality overlap between a bulk Ribo-seq DE and
# an snRNA-seq pseudobulk DE for the same cell-type contrast.
# Returns the 2x2 contingency, hypergeometric p-value, and a
# per-gene quadrant classification.
# ------------------------------------------------------ #
compute_modality_overlap = function(
    bulk_de,
    sc_de,
    bulk_diff_col,
    sc_diff_col,
    bulk_gene_col = "gene_name",
    sc_gene_col   = "gene_name",
    direction     = c("Up","Down","any")
){
    suppressPackageStartupMessages({ library(dplyr) })
    direction = match.arg(direction)

    bulk = bulk_de %>%
        dplyr::select(all_of(c(bulk_gene_col, bulk_diff_col))) %>%
        dplyr::rename(gene = !!bulk_gene_col, bulk_class = !!bulk_diff_col)
    sc   = sc_de %>%
        dplyr::select(all_of(c(sc_gene_col, sc_diff_col))) %>%
        dplyr::rename(gene = !!sc_gene_col, sc_class = !!sc_diff_col)

    joined = inner_join(bulk, sc, by = "gene")

    # When direction is "Up" or "Down", a gene with the *opposite* sign in
    # either modality is treated as "discordant" — kept as a separate count
    # so the bulk-only metric is not inflated with disagreeing-direction
    # genes. The hypergeometric is computed on the in-direction subset only.
    is_hit  = function(x) if(direction == "any") x != "No" else x == direction
    is_anti = function(x) if(direction == "any") FALSE
                          else x != "No" & x != direction

    joined = joined %>% mutate(
        bulk_hit  = is_hit(bulk_class),
        sc_hit    = is_hit(sc_class),
        bulk_anti = is_anti(bulk_class),
        sc_anti   = is_anti(sc_class),
        quadrant = case_when(
             bulk_hit  &  sc_hit                ~ "both",
             bulk_hit  & !sc_hit  & !sc_anti    ~ "bulk_only",
            !bulk_hit  & !bulk_anti & sc_hit    ~ "sc_only",
             bulk_anti |  sc_anti                ~ "discordant",
            TRUE                                 ~ "neither"
        )
    )

    universe = joined %>% filter(!bulk_anti, !sc_anti)
    n11 = sum( universe$bulk_hit &  universe$sc_hit)
    n10 = sum( universe$bulk_hit & !universe$sc_hit)
    n01 = sum(!universe$bulk_hit &  universe$sc_hit)
    n00 = sum(!universe$bulk_hit & !universe$sc_hit)

    pval = phyper(n11 - 1,
                  m = sum(universe$sc_hit),
                  n = sum(!universe$sc_hit),
                  k = sum(universe$bulk_hit),
                  lower.tail = FALSE)

    list(
        joined  = joined,
        contingency = matrix(c(n11, n10, n01, n00), nrow = 2,
                             dimnames = list(c("bulk_hit","bulk_no"),
                                             c("sc_hit","sc_no"))),
        n_discordant = sum(joined$quadrant == "discordant"),
        hyper_p = pval,
        n_total = nrow(joined),
        n_universe = nrow(universe),
        direction = direction
    )
}

