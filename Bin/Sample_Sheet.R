# ------------------------------------------------------ #
# Unified sample sheet across the three cohorts
# (riboseq_sedentary, riboseq_exercise, rnaseq).
# ------------------------------------------------------ #

load_sample_sheet = function(
    path = file.path(.scripts.path, "Documentation/sample_sheet.tsv"),
    validate = TRUE
){
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
    })

    if(!file.exists(path))
        stop("sample sheet not found: ", path)

    ss = fread(path, header = TRUE, sep = "\t", na.strings = c("NA","")) %>%
        as.data.frame()

    expected = c("sample_id","cohort","cell_type","tag","fraction",
                 "condition","replicate","bam_path","salmon_path")
    missing = setdiff(expected, colnames(ss))
    if(length(missing))
        stop("sample sheet missing columns: ", paste(missing, collapse=", "))

    if(validate){
        missing_bam    = ss$bam_path   [!file.exists(ss$bam_path)]
        missing_salmon = ss$salmon_path[!file.exists(ss$salmon_path)]
        if(length(missing_bam))
            warning("missing BAMs: ",    length(missing_bam),    " / ", nrow(ss))
        if(length(missing_salmon))
            warning("missing salmon: ", length(missing_salmon), " / ", nrow(ss))
    }

    ss
}


# ------------------------------------------------------ #
load_salmon_counts = function(
    sample_sheet,
    cohort,
    path_gtf = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf"
){
    suppressPackageStartupMessages({
        library(tximport)
        library(SummarizedExperiment)
        library(GenomicFeatures)
        library(dplyr)
    })

    ss = sample_sheet %>% filter(cohort == !!cohort)
    if(nrow(ss) == 0)
        stop("no samples for cohort: ", cohort)

    txdb  = txdbmaker::makeTxDbFromGFF(path_gtf, format = "gtf")  # Bioc >=3.19 split makeTxDbFromGFF out of GenomicFeatures
    tx2gene = AnnotationDbi::select(txdb,
        keys    = keys(txdb, "TXNAME"),
        columns = c("TXNAME","GENEID"),
        keytype = "TXNAME")

    txi = tximport(files    = setNames(ss$salmon_path, ss$sample_id),
                   type     = "salmon",
                   tx2gene  = tx2gene,
                   ignoreTxVersion = TRUE,   # salmon IDs are versioned (ENSMUST...N), GTF tx2gene is not
                   countsFromAbundance = "no")

    se = SummarizedExperiment(
        assays  = list(counts = txi$counts, abundance = txi$abundance, length = txi$length),
        colData = DataFrame(ss, row.names = ss$sample_id)
    )
    se
}
