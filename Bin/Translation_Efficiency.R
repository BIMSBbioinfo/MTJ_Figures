# ------------------------------------------------------ #
# Translation-efficiency (TE) helper: Ribo-seq LFC
# normalized by matched RNA-seq LFC, as a DESeq2
# interaction term (~ assay + condition + assay:condition).
# ------------------------------------------------------ #

# ------------------------------------------------------ #
# Canonical TE via DESeq2 interaction.
# counts : gene x sample integer matrix
# cdata  : data.frame with columns `sample`, `assay` (riboseq/rnaseq),
#          `condition` (e.g. MTJ/CK8).
# ------------------------------------------------------ #
translation_efficiency_deseq2 = function(
    counts,
    cdata,
    lfc_cut  = 0.75,
    padj_cut = 0.05
){
    suppressPackageStartupMessages({
        library(DESeq2)
        library(dplyr)
    })
    stopifnot(all(c("assay","condition") %in% colnames(cdata)))
    cdata$assay     = factor(cdata$assay)
    cdata$condition = factor(cdata$condition)

    dds = DESeqDataSetFromMatrix(
        countData = counts,
        colData   = cdata,
        design    = ~ assay + condition + assay:condition
    )
    dds = DESeq(dds)
    res = results(dds, name = grep("assay.+condition", resultsNames(dds), value = TRUE)[1])
    res = as.data.frame(res) %>%
        tibble::rownames_to_column("gene_id") %>%
        mutate(diff_TE = case_when(
            log2FoldChange >  lfc_cut & padj < padj_cut ~ "Up",
            log2FoldChange < -lfc_cut & padj < padj_cut ~ "Down",
            TRUE                                       ~ "No"
        ))
    # Return the fitted dds too: callers (build_ribo_rna_TE) need size-factor-normalized
    # counts for the per-cell-type Ribo/RNA TE ratios, not just the interaction table.
    list(res = res, dds = dds)
}


# ------------------------------------------------------ #
# Assemble the cell-type-matched Ribo/RNA translation-efficiency dataset for
# report 02 (canvas 002). For each cell type (MTJ, CK8) TE = Ribo / RNA on
# DESeq2 size-factor-normalized counts (NOT RPKM); MTJ-vs-CK8 differential TE
# significance is the `assay:condition` interaction (translation_efficiency_deseq2).
#
# Ribo  = HA-Rpl22 pulldown only. The whole-tissue Input fraction is NOT used —
#         it contains non-myofiber cells and is not a valid cell-type translation
#         control; cell-type total RNA is the bulk Tigd4/CK8 RNA-seq instead.
# Pairing: MTJ Ribo = MTJ_HA.Rpl22 ; CK8 Ribo = CK8_HA.Rpl22 ;
#          MTJ RNA  = bulk Tigd4    ; CK8 RNA  = bulk CK8.
# ------------------------------------------------------ #
build_ribo_rna_TE = function(
    sample_sheet,
    path_gtf            = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    rna_drop            = "Tigd4_br2",
    min_mean_count      = 5,
    pseudocount         = 1,
    lfc_cut             = 0.75,
    padj_cut            = 0.05,
    protein_coding_only = TRUE
){
    suppressPackageStartupMessages({
        library(DESeq2); library(SummarizedExperiment); library(dplyr)
    })

    ribo_se = load_salmon_counts(sample_sheet, "riboseq_sedentary", path_gtf = path_gtf)
    rna_se  = load_salmon_counts(sample_sheet, "rnaseq",            path_gtf = path_gtf)

    rcd = as.data.frame(colData(ribo_se))
    ribo_se = ribo_se[, rcd$tag == "HA.Rpl22" & rcd$fraction == "pulldown"]
    rna_se  = rna_se[, !colnames(rna_se) %in% rna_drop]

    genes  = intersect(rownames(ribo_se), rownames(rna_se))
    # Restrict to protein-coding genes BEFORE the DESeq2 interaction test, so the
    # multiple-testing correction and dispersion estimation consider only coding
    # genes (drops Gm*/predicted/non-coding biotypes from the tested universe).
    if(protein_coding_only){
        annot  = read_Gene_Annotation(path_gtf)
        pc_ids = annot$gene_id[annot$gene_biotype == "protein_coding"]
        genes  = intersect(genes, pc_ids)
    }
    ribo_c = round(as.matrix(assay(ribo_se, "counts"))[genes, , drop = FALSE])
    rna_c  = round(as.matrix(assay(rna_se,  "counts"))[genes, , drop = FALSE])
    counts = cbind(ribo_c, rna_c)

    cdata = data.frame(
        sample    = colnames(counts),
        assay     = rep(c("ribo","rna"), c(ncol(ribo_c), ncol(rna_c))),
        # Tigd4-sorted bulk RNA is the MTJ total-RNA layer; CK8 stays CK8.
        condition = c(ifelse(colData(ribo_se)$cell_type == "MTJ", "MTJ", "CK8"),
                      ifelse(colData(rna_se)$cell_type  == "CK8", "CK8", "MTJ")),
        row.names = colnames(counts),
        stringsAsFactors = FALSE
    )

    fit  = translation_efficiency_deseq2(counts, cdata, lfc_cut = lfc_cut, padj_cut = padj_cut)
    dds  = fit$dds
    norm = counts(dds, normalized = TRUE)

    grp = function(a, c) rowMeans(norm[, cdata$assay == a & cdata$condition == c, drop = FALSE])
    te  = data.frame(
        gene_id  = rownames(norm),
        ribo_MTJ = grp("ribo","MTJ"), rna_MTJ = grp("rna","MTJ"),
        ribo_CK8 = grp("ribo","CK8"), rna_CK8 = grp("rna","CK8"),
        row.names = NULL, stringsAsFactors = FALSE
    ) %>%
        mutate(
            log2TE_MTJ = log2((ribo_MTJ + pseudocount) / (rna_MTJ + pseudocount)),
            log2TE_CK8 = log2((ribo_CK8 + pseudocount) / (rna_CK8 + pseudocount)),
            delta_TE   = log2TE_MTJ - log2TE_CK8,
            detectable = (ribo_MTJ + ribo_CK8) / 2 >= min_mean_count &
                         (rna_MTJ  + rna_CK8 ) / 2 >= min_mean_count
        ) %>%
        left_join(fit$res %>% dplyr::select(gene_id, int_log2FC = log2FoldChange,
                                            int_pvalue = pvalue, int_padj = padj),
                  by = "gene_id") %>%
        # Significance from the interaction padj (tests differential TE); DIRECTION from
        # the directly-computed log2TE delta (robust to DESeq2 factor-reference sign).
        mutate(
            significant = detectable & !is.na(int_padj) & int_padj < padj_cut &
                          abs(delta_TE) > lfc_cut,
            te_class = case_when(
                !significant   ~ "ns_or_concordant",
                delta_TE > 0   ~ "translation_specific_MTJ",
                TRUE           ~ "translation_specific_CK8"
            )
        )

    list(dds = dds, norm = norm, cdata = cdata, te = te)
}


# ------------------------------------------------------ #
# Exercise translation efficiency for report 02 (canvas 011/020; was report 06): MTJ Ribo-seq of each
# condition normalized by a matched MTJ snRNA-seq pseudobulk (dilbaz 2026, which has
# BOTH conditions), then sedentary vs trained compared per gene.
#
# Same TE core as build_ribo_rna_TE (joint DESeq2 size factors +
# translation_efficiency_deseq2 interaction); here the two-level factor is
# `condition` ∈ {sedentary, trained} rather than cell type, and the RNA layer is
# snRNA MTJ-pseudobulk instead of bulk. Ribo (Ensembl) is collapsed to gene symbols
# to match the symbol-indexed snRNA. Significance = assay:condition interaction
# (TE changed with training); direction from the log2-TE delta (sign-robust).
# ------------------------------------------------------ #
build_exercise_TE = function(
    sample_sheet,
    path_gtf       = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    min_mean_count = 5,
    pseudocount    = 1,
    lfc_cut        = 0.75,
    padj_cut       = 0.05
){
    suppressPackageStartupMessages({
        library(DESeq2); library(SummarizedExperiment); library(Matrix)
        library(Seurat); library(SeuratObject); library(dplyr)
    })

    annot  = read_Gene_Annotation(path_gtf)
    id2sym = setNames(annot$gene_name, annot$gene_id)

    # --- Ribo (MTJ HA pulldown), Ensembl IDs collapsed to gene symbols ---
    ribo_symbol = function(cohort, tag){
        se = load_salmon_counts(sample_sheet, cohort, path_gtf = path_gtf)
        cd = as.data.frame(colData(se))
        se = se[, cd$cell_type == "MTJ" & cd$tag == tag & cd$fraction == "pulldown"]
        m   = round(as.matrix(assay(se, "counts")))
        sym = id2sym[rownames(m)]
        keep = !is.na(sym) & nzchar(sym)
        rowsum(m[keep, , drop = FALSE], group = sym[keep])   # sum duplicate symbols
    }
    ribo_sed = ribo_symbol("riboseq_sedentary", "HA.Rpl22")   # symbols x 3 (br2 absent)
    ribo_tr  = ribo_symbol("riboseq_exercise",  "HA")         # symbols x 3

    # --- dilbaz MTJ pseudobulk per condition x replicate (deseq2=="MTJ" label) ---
    s   = load_dilbaz_2026()
    md  = s@meta.data
    mtj = which(md$deseq2 == "MTJ")
    if(length(mtj) == 0) stop("no deseq2=='MTJ' nuclei in dilbaz")
    cnt = SeuratObject::GetAssayData(s, assay = "RNA", layer = "counts")[, mtj]
    grp = factor(paste(md$condition[mtj], md$replicate[mtj], sep = "_"))
    ind = Matrix::sparse.model.matrix(~ 0 + grp); colnames(ind) = levels(grp)
    pb  = as.matrix(cnt %*% ind)                              # symbols x (Sedentary_*, Trained_*)

    # --- MTJ-label validation (oracle): Tigd4/Col22a1 enrichment vs other classes ---
    mk = intersect(c("Tigd4","Col22a1"), rownames(s))
    dat = as.matrix(SeuratObject::GetAssayData(s, assay = "RNA", layer = "data")[mk, , drop = FALSE])
    mtj_validation = data.frame(class = ifelse(md$deseq2 == "MTJ", "MTJ", "other"),
                                t(dat), check.names = FALSE) %>%
        group_by(class) %>%
        summarize(across(all_of(mk), \(x) round(mean(x), 3)),
                  n_nuclei = dplyr::n(), .groups = "drop")

    # --- combine on shared symbols ---
    genes  = Reduce(intersect, list(rownames(ribo_sed), rownames(ribo_tr), rownames(pb)))
    counts = cbind(ribo_sed[genes, , drop = FALSE],
                   ribo_tr [genes, , drop = FALSE],
                   pb      [genes, , drop = FALSE])
    cdata  = data.frame(
        sample    = colnames(counts),
        assay     = c(rep("ribo", ncol(ribo_sed) + ncol(ribo_tr)), rep("rna", ncol(pb))),
        condition = c(rep("sedentary", ncol(ribo_sed)), rep("trained", ncol(ribo_tr)),
                      ifelse(grepl("^Sedentary", colnames(pb)), "sedentary", "trained")),
        row.names = colnames(counts), stringsAsFactors = FALSE
    )

    # --- joint DESeq2 TE (shared core: ~ assay + condition + assay:condition) ---
    fit  = translation_efficiency_deseq2(counts, cdata, lfc_cut = lfc_cut, padj_cut = padj_cut)
    dds  = fit$dds
    norm = counts(dds, normalized = TRUE)

    gm = function(a, c) rowMeans(norm[, cdata$assay == a & cdata$condition == c, drop = FALSE])
    te = data.frame(
        gene_name = rownames(norm),
        ribo_sed = gm("ribo","sedentary"), rna_sed = gm("rna","sedentary"),
        ribo_tr  = gm("ribo","trained"),   rna_tr  = gm("rna","trained"),
        row.names = NULL, stringsAsFactors = FALSE
    ) %>%
        mutate(
            log2TE_sedentary = log2((ribo_sed + pseudocount) / (rna_sed + pseudocount)),
            log2TE_trained   = log2((ribo_tr  + pseudocount) / (rna_tr  + pseudocount)),
            delta_TE         = log2TE_trained - log2TE_sedentary,
            detectable       = (ribo_sed + ribo_tr) / 2 >= min_mean_count &
                               (rna_sed  + rna_tr ) / 2 >= min_mean_count
        ) %>%
        # significance from interaction padj; direction from the log2-TE delta
        left_join(fit$res %>% dplyr::select(gene_name = gene_id,
                                            int_log2FC = log2FoldChange, int_padj = padj),
                  by = "gene_name") %>%
        mutate(
            significant = detectable & !is.na(int_padj) & int_padj < padj_cut &
                          abs(delta_TE) > lfc_cut,
            te_class = case_when(
                !significant ~ "ns",
                delta_TE > 0 ~ "TE_up_trained",
                TRUE         ~ "TE_down_trained"
            )
        )

    list(dds = dds, norm = norm, cdata = cdata, te = te, mtj_validation = mtj_validation)
}
