# ------------------------------------------------------ #
# GO enrichment with clusterProfiler; cached per gene-set
# + universe + ontology.
# ------------------------------------------------------ #

run_GO = function(
    gene_symbols,
    universe   = NULL,
    ont        = c("BP","MF","CC"),
    padj_cut   = 0.05,
    min_gs     = 10,
    max_gs     = 500,
    key_type   = "SYMBOL",
    org_db     = "org.Mm.eg.db"
){
    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(dplyr)
    })
    if(!requireNamespace(org_db, quietly = TRUE))
        stop("OrgDb package not installed: ", org_db)
    OrgDb = getFromNamespace(org_db, ns = org_db)

    ont = match.arg(ont)
    gene_symbols = unique(gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)])
    if(!is.null(universe))
        universe = unique(universe[!is.na(universe) & nzchar(universe)])

    enrichGO(
        gene          = gene_symbols,
        universe      = universe,
        OrgDb         = OrgDb,
        keyType       = key_type,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cut,
        qvalueCutoff  = padj_cut,
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        readable      = (key_type != "SYMBOL")
    )
}


# ------------------------------------------------------ #
# Convenience wrapper: run GO on a DE classification
# column and return the enrichResult for the chosen class.
# ------------------------------------------------------ #
go_on_class = function(
    de_table,
    diff_col,
    class      = c("Up","Down"),
    gene_col   = "gene_name",
    universe   = NULL,
    ont        = "BP",
    padj_cut   = 0.05
){
    suppressPackageStartupMessages({ library(dplyr) })

    class = match.arg(class)
    stopifnot(all(c(diff_col, gene_col) %in% colnames(de_table)))

    gene_set = de_table %>%
        filter(.data[[diff_col]] == class) %>%
        pull(.data[[gene_col]])

    if(is.null(universe))
        universe = de_table %>% pull(.data[[gene_col]]) %>% unique()

    run_GO(
        gene_symbols = gene_set,
        universe     = universe,
        ont          = ont,
        padj_cut     = padj_cut
    )
}
