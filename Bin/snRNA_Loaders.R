# ------------------------------------------------------ #
# Loaders for the three external mouse muscle snRNA-seq
# atlases pulled on 2026-04-24. All assemblies are mm10;
# integration with project Ribo/RNA-seq (mm39) is done
# on gene symbols.
#
# External data root:
#   Data/external/snRNA/
#     ├── kim_birchmeier_2020_E-MTAB-8623/
#     │     └── mm10.Seurat.RDS              (preprocessed Seurat)
#     ├── petrany_2020_GSE147127/
#     │     └── GSM44189{91..96}_*.h5         (CellRanger h5 × 6)
#     └── demicheli_2020_GSE143437/
#           ├── GSE143437_DeMicheli_MuSCatlas_metadata.txt.gz
#           └── GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz
# ------------------------------------------------------ #

.external_snrna_root = function(){
    file.path(Sys.getenv("LOCAL_PATH"), "MKim_MTJ", "Data", "external", "snRNA")
}


# ------------------------------------------------------ #
load_kim_birchmeier_2020 = function(
    path = file.path(.external_snrna_root(),
                     "kim_birchmeier_2020_E-MTAB-8623",
                     "mm10.Seurat.RDS")
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject)
        library(AnnotationDbi); library(org.Mm.eg.db)
    })
    if(!file.exists(path))
        stop("Kim 2020 Seurat object not found: ", path)
    # Preprocessed object was saved under an older Seurat (v3/v4); migrate it to
    # the v5 class structure or v5 generics fail (e.g. missing 'images' slot).
    obj = suppressWarnings(UpdateSeuratObject(readRDS(path)))

    # The object is indexed by Ensembl gene IDs (ENSMUSG...); the rest of the
    # project (and the marker panels in reports 05/07) works in gene symbols.
    # Relabel the RNA assay to unique symbols so Tigd4/Col22a1/... resolve.
    # Cells, metadata (clusters/celltype) and reductions (UMAP) are preserved.
    ens = rownames(obj[["RNA"]])
    if(all(grepl("^ENSMUSG", head(ens, 50)))){
        sym  = AnnotationDbi::mapIds(org.Mm.eg.db, keys = ens, column = "SYMBOL",
                                     keytype = "ENSEMBL", multiVals = "first")
        keep = !is.na(sym) & nzchar(sym) & !duplicated(sym)
        cnt  = SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")[keep, ]
        rownames(cnt) = sym[keep]
        new  = CreateSeuratObject(counts = cnt, meta.data = obj[[]],
                                  min.cells = 0, min.features = 0)
        # The hosted object ships without clustering or a UMAP; cluster de novo
        # (same params as load_petrany_2020) so reports 05/07 have seurat_clusters
        # and a umap reduction to work with.
        new  = NormalizeData(new, verbose = FALSE)
        new  = FindVariableFeatures(new, verbose = FALSE, nfeatures = 3000)
        new  = ScaleData(new, verbose = FALSE)
        new  = RunPCA(new, verbose = FALSE)
        new  = FindNeighbors(new, dims = 1:30, verbose = FALSE)
        new  = FindClusters(new, resolution = 0.6, verbose = FALSE)
        new  = RunUMAP(new, dims = 1:30, verbose = FALSE)
        obj  = new
    }
    obj
}


# ------------------------------------------------------ #
load_petrany_2020 = function(
    dir = file.path(.external_snrna_root(), "petrany_2020_GSE147127",
                    "synapse_cellranger"),
    mito_pct_max = 10,
    nfeature_min = 200,
    nfeature_max = 7000
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SeuratObject)
        library(dplyr)
    })

    # NOTE: the GEO SoupX .h5 (GSM44189*) are dense matrices with NO gene names
    # (unusable). The real CellRanger filtered_feature_bc_matrix.h5 (raw integer
    # counts + gene symbols) live on Synapse syn21676145 and are downloaded to
    # synapse_cellranger/ (5 TA samples; soleus has no CellRanger folder). See README.
    h5s = list.files(dir, pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
    if(length(h5s) == 0)
        stop("no CellRanger .h5 under ", dir,
             " (download from Synapse syn21676145; see Data README)")

    label = function(path)
        sub("_filtered_feature_bc_matrix\\.h5$", "", basename(path))

    lst = lapply(h5s, function(p){
        mat = Read10X_h5(p)
        s   = CreateSeuratObject(counts = mat, project = "Petrany2020")
        s[["sample_id"]] = label(p)
        s[["percent.mt"]] = PercentageFeatureSet(s, pattern = "^mt-")
        s = subset(s,
                   subset = nFeature_RNA > nfeature_min &
                            nFeature_RNA < nfeature_max &
                            percent.mt   < mito_pct_max)
        s
    })
    names(lst) = vapply(h5s, label, character(1))

    merged = merge(lst[[1]], y = lst[-1], add.cell.ids = names(lst))
    merged = JoinLayers(merged)   # Seurat v5: collapse per-sample count layers before processing
    merged = NormalizeData(merged, verbose = FALSE)
    merged = FindVariableFeatures(merged, verbose = FALSE, nfeatures = 3000)
    merged = ScaleData(merged, verbose = FALSE)
    merged = RunPCA(merged, verbose = FALSE)
    merged = FindNeighbors(merged, dims = 1:30, verbose = FALSE)
    merged = FindClusters(merged, resolution = 0.6, verbose = FALSE)
    merged = RunUMAP(merged, dims = 1:30, verbose = FALSE)
    merged
}


# ------------------------------------------------------ #
load_demicheli_2020 = function(
    dir = file.path(.external_snrna_root(), "demicheli_2020_GSE143437")
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(data.table)
    })

    meta_path = file.path(dir, "GSE143437_DeMicheli_MuSCatlas_metadata.txt.gz")
    raw_path  = file.path(dir, "GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz")

    meta = fread(meta_path, header = TRUE)
    raw  = fread(raw_path,  header = TRUE)

    gene_col   = names(raw)[1]
    gene_names = raw[[gene_col]]
    mat        = as.matrix(raw[, -gene_col, with = FALSE])
    rownames(mat)  = gene_names

    md = as.data.frame(meta)
    rownames(md)  = md[[1]]
    common = intersect(colnames(mat), rownames(md))
    mat    = mat[, common, drop = FALSE]
    md     = md[common, , drop = FALSE]

    s = CreateSeuratObject(counts = mat, meta.data = md, project = "DeMicheli2020")
    s = NormalizeData(s, verbose = FALSE)
    s
}


# ------------------------------------------------------ #
# iScience 2021 mouse muscle/tendon single-cell atlas (GEO GSE162307).
# 5 samples (H2B / VEH4 / VEH24 / TAM4 / TAM24) shipped as .loom (HDF5), read with
# hdf5r. The loom `matrix` is stored cells x genes here; it is oriented to genes x
# cells (rows = row_attrs/Gene, cols = col_attrs/CellID) before building the Seurat
# object. Author cluster ids (col_attrs/Clusters) are carried as `loom_cluster`;
# the object is re-clustered de novo (same params as the other loaders).
load_iscience_2021 = function(
    dir = file.path(.external_snrna_root(), "iscience_2021_GSE162307")
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(Matrix); library(hdf5r)
    })

    looms = list.files(dir, pattern = "\\.loom$", full.names = TRUE)
    if(length(looms) == 0){
        gz = list.files(dir, pattern = "\\.loom\\.gz$", full.names = TRUE)
        if(length(gz) == 0)
            stop("no .loom under ", dir, " (download GSE162307 per-GSM suppl; see README)")
        for(g in gz){ out = sub("\\.gz$", "", g); if(!file.exists(out)) system2("gunzip", c("-kf", shQuote(g))) }
        looms = list.files(dir, pattern = "\\.loom$", full.names = TRUE)
    }
    label = function(p) sub("\\.loom$", "", sub("^GSM[0-9]+_", "", basename(p)))

    read_loom = function(p){
        f = H5File$new(p, mode = "r")
        on.exit(f$close_all())
        genes = as.character(f[["row_attrs/Gene"]]$read())
        cells = as.character(f[["col_attrs/CellID"]]$read())
        m = f[["matrix"]]$read()
        if(nrow(m) == length(cells) && ncol(m) == length(genes)) m = t(m)   # -> genes x cells
        rownames(m) = make.unique(genes)
        colnames(m) = cells
        clusters = if("Clusters" %in% names(f[["col_attrs"]])) f[["col_attrs/Clusters"]]$read() else NA
        list(mat = Matrix::Matrix(m, sparse = TRUE), clusters = clusters)
    }

    lst = lapply(looms, function(p){
        d = read_loom(p)
        s = CreateSeuratObject(counts = d$mat, project = "iScience2021")
        s[["sample_id"]]    = label(p)
        s[["loom_cluster"]] = d$clusters
        s
    })
    names(lst) = vapply(looms, label, character(1))

    merged = if(length(lst) > 1) merge(lst[[1]], y = lst[-1], add.cell.ids = names(lst)) else lst[[1]]
    merged = JoinLayers(merged)
    merged = NormalizeData(merged, verbose = FALSE)
    merged = FindVariableFeatures(merged, verbose = FALSE, nfeatures = 3000)
    merged = ScaleData(merged, verbose = FALSE)
    merged = RunPCA(merged, verbose = FALSE)
    merged = FindNeighbors(merged, dims = 1:30, verbose = FALSE)
    merged = FindClusters(merged, resolution = 0.6, verbose = FALSE)
    merged = RunUMAP(merged, dims = 1:30, verbose = FALSE)
    merged
}


# ------------------------------------------------------ #
# Dilbaz 2026 proprietary exercise snRNA-seq (sedentary + trained mouse muscle).
# 10x feature-barcode matrix (symbol-indexed) + meta_data.csv. The metadata
# already carries a cell classification column `deseq2` with an explicit MTJ
# label (no de-novo clustering needed); `condition` ∈ {Sedentary, Trained},
# `replicate` ∈ {A..D}, `time` = PRE. PROPRIETARY — never commit the data.
# Used by report 02 (canvas 011/020; was report 06) for exercise Ribo/RNA normalization.
# ------------------------------------------------------ #
load_dilbaz_2026 = function(
    dir = file.path(.external_snrna_root(), "dilbaz_2026_proprietary",
                    "feature_count_matrix_mtx")
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(data.table)
    })
    if(!dir.exists(dir))
        stop("dilbaz 2026 matrix dir not found: ", dir)

    mat  = Read10X(dir, gene.column = 2)   # features.tsv col2 = gene symbol (mm39-compatible)
    meta = as.data.frame(data.table::fread(file.path(dir, "meta_data.csv.gz")))
    rownames(meta) = meta[[1]]; meta[[1]] = NULL

    common = intersect(colnames(mat), rownames(meta))
    if(length(common) == 0)
        stop("dilbaz: 0 shared cells between matrix barcodes and meta_data")
    s = CreateSeuratObject(counts = mat[, common], meta.data = meta[common, , drop = FALSE])
    s = NormalizeData(s, verbose = FALSE)   # for the Tigd4/Col22a1 MTJ-label validation
    s
}


# ------------------------------------------------------ #
# Canvas:   Prompts/canvases/078_klf2dn-snrna-analysis.md
# KLF2-DN snRNA loader — 8 libraries (4 mice x 2 legs), STARsolo GeneFull output.
# Reads each library's Solo.out/GeneFull/filtered/ matrix (10x mtx format via Read10X,
# NOT Read10X_h5 — STARsolo, not CellRanger), attaches the sample sheet (mouse/leg/
# condition), applies the project snRNA QC cuts (same as load_petrany_2020), then a
# PLAIN merge (no batch integration) + log-normalize + cluster at resolution 0.6.
# The custom reference (build_klf2dn_reference.sh) adds `Cre` and `KLF2-DN` as genes, so
# both appear in rownames() for transduction scoring / MTJ calling.
#
# Pre-filter per-cell QC (nCount/nFeature/percent.mt/sample_id) is stashed in
# Misc(obj, "pre_filter_qc") so the report can draw the QC panel BEFORE the cuts.
#
# sample_sheet: CSV/data.frame with columns `library` (= counts subdir name),
#   `mouse`, `leg`, `condition` (KLF2-DN | Cre).
# ------------------------------------------------------ #
load_klf2dn = function(
    counts_root  = "/data/local/Projects/MKim_MTJ/Results/klf2dn_starsolo_counts",
    sample_sheet = file.path(.external_snrna_root(), "klf2dn_genomeeast", "sample_sheet.csv"),
    solo_subpath = "Solo.out/GeneFull/filtered",
    mito_pct_max = 10,
    nfeature_min = 200,
    nfeature_max = 7000
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(dplyr)
    })

    ss = if(is.character(sample_sheet)) {
             if(!file.exists(sample_sheet)) stop("sample sheet not found: ", sample_sheet)
             read.csv(sample_sheet, stringsAsFactors = FALSE)
         } else sample_sheet
    req = c("library","mouse","leg","condition")
    if(!all(req %in% colnames(ss)))
        stop("sample sheet needs columns: ", paste(req, collapse = ", "))

    pre_qc = list()
    lst = lapply(seq_len(nrow(ss)), function(i){
        row = ss[i, ]
        dir = file.path(counts_root, row$library, solo_subpath)
        if(!dir.exists(dir))
            stop("STARsolo matrix dir not found for library '", row$library, "': ", dir,
                 "  (run Scripts/Bin/starsolo_count.sh first)")
        # STARsolo writes PLAIN (ungzipped) matrix.mtx/barcodes.tsv/features.tsv; Read10X on a
        # v3-style features.tsv dir insists on the .gz variants, so read explicitly with ReadMtx
        # (transparently handles .gz or plain). features.tsv col 2 = gene symbol.
        pick = function(base){ g = file.path(dir, paste0(base, ".gz"))
                               if(file.exists(g)) g else file.path(dir, base) }
        mat = ReadMtx(mtx = pick("matrix.mtx"), cells = pick("barcodes.tsv"),
                      features = pick("features.tsv"), feature.column = 2)
        s = CreateSeuratObject(counts = mat, project = "KLF2DN")
        s[["sample_id"]] = row$library
        s[["mouse"]]     = as.character(row$mouse)
        s[["leg"]]       = as.character(row$leg)
        s[["condition"]] = as.character(row$condition)
        s[["percent.mt"]] = PercentageFeatureSet(s, pattern = "^mt-")
        # stash pre-filter QC for the report's before-cuts panel
        pre_qc[[row$library]] <<- data.frame(
            sample_id  = row$library, condition = row$condition, mouse = row$mouse,
            nCount_RNA = s$nCount_RNA, nFeature_RNA = s$nFeature_RNA, percent.mt = s$percent.mt)
        subset(s, subset = nFeature_RNA > nfeature_min &
                           nFeature_RNA < nfeature_max &
                           percent.mt   < mito_pct_max)
    })
    names(lst) = ss$library

    merged = merge(lst[[1]], y = lst[-1], add.cell.ids = names(lst))   # plain merge, no integration
    merged = JoinLayers(merged)
    merged = NormalizeData(merged, verbose = FALSE)
    merged = FindVariableFeatures(merged, verbose = FALSE, nfeatures = 3000)
    merged = ScaleData(merged, verbose = FALSE)
    merged = RunPCA(merged, verbose = FALSE)
    merged = FindNeighbors(merged, dims = 1:30, verbose = FALSE)
    merged = FindClusters(merged, resolution = 0.6, verbose = FALSE)
    merged = RunUMAP(merged, dims = 1:30, verbose = FALSE)
    Misc(merged, "pre_filter_qc") = dplyr::bind_rows(pre_qc)
    merged
}


# ------------------------------------------------------ #
# Per-sample QC-and-marker processing for the KLF2-DN libraries (canvas 078,
# 2026-07-27). For EACH library independently: read the STARsolo GeneFull filtered
# matrix, apply the project QC cuts, run the standard Seurat pipeline (Normalize /
# HVG / Scale / PCA / UMAP / FindClusters res 0.6), then attach three per-cell QC
# metrics used for the per-sample UMAP/violin panels in report 08:
#   - doublet_score / doublet_class : scDblFinder (cluster-aware), on the raw counts
#   - ambient_frac                  : SoupX per-cell soup fraction — the share of a
#                                     cell's UMIs attributed to ambient RNA, from the
#                                     STARsolo RAW (unfiltered) droplet matrix. rho
#                                     (global contamination) is stashed per object in
#                                     Misc(obj,"soupx_rho").
# MTJ is labelled as the top-`pick_marker` (Tigd4) cluster within that sample.
# NB: this VISUALISES QC only — nothing is removed here.
# Returns a named list (by library) of per-sample Seurat objects.
# ------------------------------------------------------ #
process_klf2dn_per_sample = function(
    counts_root  = "/data/local/Projects/MKim_MTJ/Results/klf2dn_starsolo_counts",
    sample_sheet = file.path(.external_snrna_root(), "klf2dn_genomeeast", "sample_sheet.csv"),
    solo_filtered = "Solo.out/GeneFull/filtered",
    solo_raw      = "Solo.out/GeneFull/raw",
    pick_marker  = "Tigd4",
    mito_pct_max = 10,
    nfeature_min = 200,
    nfeature_max = 7000,
    resolution   = 0.6,
    dims         = 1:30,
    nfeatures    = 3000,
    do_doublets  = TRUE,
    do_ambient   = TRUE
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(Matrix); library(dplyr)
    })
    ss = if(is.character(sample_sheet)) read.csv(sample_sheet, stringsAsFactors = FALSE) else sample_sheet

    pick_paths = function(dir){
        f = function(base){ g = file.path(dir, paste0(base, ".gz"))
                            if(file.exists(g)) g else file.path(dir, base) }
        list(mtx = f("matrix.mtx"), cells = f("barcodes.tsv"), features = f("features.tsv"))
    }

    out = lapply(seq_len(nrow(ss)), function(i){
        row = ss[i, ]
        fdir = file.path(counts_root, row$library, solo_filtered)
        pf = pick_paths(fdir)
        mat = ReadMtx(mtx = pf$mtx, cells = pf$cells, features = pf$features, feature.column = 2)
        s = CreateSeuratObject(counts = mat, project = "KLF2DN")
        s[["sample_id"]]  = row$library
        s[["mouse"]]      = as.character(row$mouse)
        s[["condition"]]  = as.character(row$condition)
        s[["percent.mt"]] = PercentageFeatureSet(s, pattern = "^mt-")
        s = subset(s, subset = nFeature_RNA > nfeature_min &
                               nFeature_RNA < nfeature_max &
                               percent.mt   < mito_pct_max)
        # Unified transduction readout (condition-specific references, 2026-07-28): each library carries
        # only its own transgene contig — Cre for control legs, KLF2-DN for DN legs (a single combined
        # reference multimaps the shared 3' backbone and zeros both). `transduction` = per-cell UMI of
        # whichever transgene is present; `transgene` records which. NB DN-leg KLF2-DN captures only the
        # 3'UTR-proximal reads (~1/3 of the Cre-leg capture), so `transduction` is NOT cross-condition
        # comparable in magnitude — use it for detection/localisation, not for a Cre-vs-DN level test.
        tg = intersect(c("Cre","KLF2-DN"), rownames(s))
        s[["transgene"]]    = if(length(tg)) paste(tg, collapse = "+") else "none"
        s[["transduction"]] = if(length(tg))
            Matrix::colSums(SeuratObject::GetAssayData(s, assay = "RNA", layer = "counts")[tg, , drop = FALSE])
            else 0
        s = NormalizeData(s, verbose = FALSE)
        s = FindVariableFeatures(s, verbose = FALSE, nfeatures = nfeatures)
        s = ScaleData(s, verbose = FALSE)
        s = RunPCA(s, verbose = FALSE, npcs = min(50, ncol(s) - 1))
        s = FindNeighbors(s, dims = dims, verbose = FALSE)
        s = FindClusters(s, resolution = resolution, verbose = FALSE)
        s = RunUMAP(s, dims = dims, verbose = FALSE)

        # --- doublets: scDblFinder (cluster-aware) on the raw counts ------------------
        if(do_doublets){
            suppressPackageStartupMessages({
                library(scDblFinder); library(SingleCellExperiment); library(BiocParallel) })
            set.seed(42)
            sce = SingleCellExperiment(assays = list(
                      counts = SeuratObject::GetAssayData(s, assay = "RNA", layer = "counts")))
            # SerialParam: run scDblFinder in-process — avoids exporting its large kNN index to
            # future/BiocParallel workers (trips future.globals.maxSize under knit_all's plan).
            sce = scDblFinder(sce, clusters = as.character(s$seurat_clusters),
                              BPPARAM = BiocParallel::SerialParam())
            s[["doublet_score"]] = sce$scDblFinder.score
            s[["doublet_class"]] = as.character(sce$scDblFinder.class)
        }

        # --- ambient: SoupX per-cell soup fraction (needs the RAW droplet matrix) -----
        if(do_ambient){
            suppressPackageStartupMessages({ library(SoupX) })
            rdir = file.path(counts_root, row$library, solo_raw)
            pr = pick_paths(rdir)
            tod = ReadMtx(mtx = pr$mtx, cells = pr$cells, features = pr$features, feature.column = 2)
            toc = SeuratObject::GetAssayData(s, assay = "RNA", layer = "counts")
            common = intersect(rownames(tod), rownames(toc))
            sc = SoupChannel(tod[common, ], toc[common, ], calcSoupProfile = TRUE)
            sc = setClusters(sc, setNames(as.character(s$seurat_clusters), colnames(toc)))
            sc = tryCatch(autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE),
                          error = function(e) setContaminationFraction(sc, 0.1))
            adj = adjustCounts(sc, roundToInt = FALSE)
            amb = 1 - Matrix::colSums(adj[common, colnames(toc)]) / Matrix::colSums(toc[common, ])
            s[["ambient_frac"]] = pmax(0, pmin(1, amb[colnames(s)]))
            Misc(s, "soupx_rho") = unique(sc$metaData$rho)[1]
        }

        # --- MTJ = top-pick_marker (Tigd4) cluster per sample -------------------------
        if(pick_marker %in% rownames(s)){
            sc_expr = FetchData(s, vars = c("seurat_clusters", pick_marker))
            mtj_cl = sc_expr %>% group_by(seurat_clusters) %>%
                summarize(p = mean(.data[[pick_marker]] > 0), .groups = "drop") %>%
                arrange(desc(p)) %>% dplyr::slice(1) %>% pull(seurat_clusters) %>% as.character()
            s[["mtj"]] = as.character(s$seurat_clusters) == mtj_cl
        }
        s
    })
    names(out) = ss$library
    out
}


# ------------------------------------------------------ #
# Kim 2020 MTJ signature module score + INTEGRATED (Harmony) MTJ definition.
# These provide a cluster-robust alternative to the per-sample top-Tigd4 MTJ call
# (which is sensitive to per-sample Louvain boundaries).
#
# kim_mtj_module_score(): AddModuleScore of the Kim MTJ marker set onto one object.
# integrate_klf2dn(): merge the 8 per-sample objects, integrate across samples with
#   Seurat-native RPCA (no external dep — harmony/batchelor are not in the project renv),
#   cluster ONCE jointly, score the Kim MTJ signature, and call the integrated MTJ
#   cluster (highest mean module score). Returns the integrated Seurat object with
#   `kimMTJ1` (score), joint `seurat_clusters`, and `mtj_int` (integrated MTJ membership).
# ------------------------------------------------------ #
kim_mtj_module_score = function(obj, sig = NULL, name = "kimMTJ", ctrl = 50){
    suppressPackageStartupMessages(library(Seurat))
    if(is.null(sig)) sig = kim_mtj_signature()$gene_name
    sig = intersect(sig, rownames(obj))
    AddModuleScore(obj, features = list(sig), name = name, ctrl = min(ctrl, length(sig)))
}

integrate_klf2dn = function(
    resolution = 0.6,
    dims       = 1:30,
    nfeatures  = 3000,
    sig_n_top  = 200
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(dplyr)
    })
    # IntegrateLayers (RPCA) parallelises anchor-finding via future and tries to export the
    # merged object (~9 GB) to workers, tripping future.globals.maxSize under knit_all's plan.
    # Force sequential (no export) for the duration, same rationale as scDblFinder's SerialParam.
    oplan = future::plan("sequential"); on.exit(future::plan(oplan), add = TRUE)
    oopt  = options(future.globals.maxSize = 32 * 1024^3); on.exit(options(oopt), add = TRUE)

    ps  = process_klf2dn_per_sample()
    sig = kim_mtj_signature(n_top = sig_n_top)$gene_name

    # merge -> unify -> split RNA into per-sample layers for Seurat v5 integration
    m = merge(ps[[1]], y = ps[-1], add.cell.ids = names(ps))
    m = JoinLayers(m)
    m[["RNA"]] = split(m[["RNA"]], f = m$sample_id)

    m = NormalizeData(m, verbose = FALSE)
    m = FindVariableFeatures(m, verbose = FALSE, nfeatures = nfeatures)
    m = ScaleData(m, verbose = FALSE)
    m = RunPCA(m, verbose = FALSE, npcs = 50)
    # Seurat-native RPCA integration across the per-sample layers (no harmony dep in renv)
    m = IntegrateLayers(m, method = RPCAIntegration, orig.reduction = "pca",
                        new.reduction = "integrated.rpca", verbose = FALSE)
    m = FindNeighbors(m, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
    m = FindClusters(m, resolution = resolution, verbose = FALSE)
    m = RunUMAP(m, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
    m = JoinLayers(m)

    # Kim MTJ module score; integrated MTJ = cluster with highest mean score
    m = kim_mtj_module_score(m, sig = sig, name = "kimMTJ")
    sc = FetchData(m, vars = c("seurat_clusters", "kimMTJ1"))
    mtj_cl = sc %>% group_by(seurat_clusters) %>%
        summarize(s = mean(kimMTJ1), .groups = "drop") %>%
        arrange(desc(s)) %>% dplyr::slice(1) %>% pull(seurat_clusters) %>% as.character()
    m[["mtj_int"]] = as.character(m$seurat_clusters) == mtj_cl
    Misc(m, "mtj_int_cluster") = mtj_cl
    m
}


# ------------------------------------------------------ #
# PER-MOUSE MTJ definition: integrate each mouse's two contralateral legs on their own
# (RPCA), cluster jointly, and call MTJ = highest Kim-module-score cluster — so both legs
# of a mouse share ONE MTJ definition, the cleanest basis for the paired Cre-vs-KLF2-DN
# %MTJ comparison (no cross-mouse batch, no per-leg clustering wobble). Returns a NAMED
# LIST of per-mouse Seurat objects (RPCA UMAP + joint `seurat_clusters` + `kimMTJ1` +
# `mtj_mouse` logical; `Misc(o,"mtj_cluster")` = the MTJ cluster id). Derive the per-leg
# %MTJ summary from `mtj_mouse` in the caller.
# ------------------------------------------------------ #
mtj_per_mouse_integrated = function(
    resolution = 0.6,
    dims       = 1:30,
    nfeatures  = 3000,
    sig_n_top  = 200
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(dplyr)
    })
    oplan = future::plan("sequential"); on.exit(future::plan(oplan), add = TRUE)
    oopt  = options(future.globals.maxSize = 32 * 1024^3); on.exit(options(oopt), add = TRUE)

    ps   = process_klf2dn_per_sample()
    sig  = kim_mtj_signature(n_top = sig_n_top)$gene_name
    meta = data.frame(sample    = names(ps),
                      mouse     = sapply(ps, function(o) unique(as.character(o$mouse))),
                      condition = sapply(ps, function(o) unique(as.character(o$condition))),
                      stringsAsFactors = FALSE)
    mice = sort(unique(meta$mouse))

    out = lapply(mice, function(mm){
        legs = meta$sample[meta$mouse == mm]
        m = merge(ps[[legs[1]]], y = ps[legs[-1]], add.cell.ids = legs)
        m = JoinLayers(m)
        m[["RNA"]] = split(m[["RNA"]], f = m$sample_id)
        m = NormalizeData(m, verbose = FALSE)
        m = FindVariableFeatures(m, verbose = FALSE, nfeatures = nfeatures)
        m = ScaleData(m, verbose = FALSE)
        m = RunPCA(m, verbose = FALSE, npcs = 50)
        m = IntegrateLayers(m, method = RPCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.rpca", verbose = FALSE)
        m = FindNeighbors(m, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
        m = FindClusters(m, resolution = resolution, verbose = FALSE)
        m = RunUMAP(m, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
        m = JoinLayers(m)
        m = kim_mtj_module_score(m, sig = sig, name = "kimMTJ")
        sc = FetchData(m, vars = c("seurat_clusters", "kimMTJ1"))
        mtj_cl = sc %>% group_by(seurat_clusters) %>%
            summarize(s = mean(kimMTJ1), .groups = "drop") %>%
            arrange(desc(s)) %>% dplyr::slice(1) %>% pull(seurat_clusters) %>% as.character()
        m$mtj_mouse = as.character(m$seurat_clusters) == mtj_cl
        Misc(m, "mtj_cluster") = mtj_cl
        m
    })
    names(out) = mice
    out
}


# ------------------------------------------------------ #
# Transgene⁺ vs transgene⁻ per-cell DE, pooled within condition (Cre legs or DN legs).
# `which="cre"` pools the 4 Cre legs, `"dn"` the 4 DN legs; transgene⁺ = transduction>0.
# Wilcoxon FindMarkers on the log-norm RNA. (Previously cacheFile-memoized; now recomputed.)
# not recomputed on every report render. Returns list(de, n_pos, n_neg, med_nf_pos, med_nf_neg).
# ------------------------------------------------------ #
transgene_pos_markers = function(
    which   = "cre",
    logfc   = 0.25,
    min_pct = 0.1
){
    suppressPackageStartupMessages({
        library(Seurat); library(SeuratObject); library(dplyr)
    })
    oplan = future::plan("sequential"); on.exit(future::plan(oplan), add = TRUE)
    oopt  = options(future.globals.maxSize = 32 * 1024^3); on.exit(options(oopt), add = TRUE)

    ps   = process_klf2dn_per_sample()
    cond = if(which == "cre") "Cre" else "KLF2-DN"
    nms  = names(ps)[sapply(ps, function(o) unique(as.character(o$condition)) == cond)]
    m = merge(ps[[nms[1]]], y = ps[nms[-1]]); m = JoinLayers(m)
    m$tgpos = ifelse(m$transduction > 0, "pos", "neg")
    Idents(m) = "tgpos"
    de = FindMarkers(m, ident.1 = "pos", ident.2 = "neg",
                     logfc.threshold = logfc, min.pct = min_pct)
    de$gene = rownames(de)
    list(de = de, n_pos = sum(m$tgpos == "pos"), n_neg = sum(m$tgpos == "neg"),
         med_nf_pos = median(m$nFeature_RNA[m$tgpos == "pos"]),
         med_nf_neg = median(m$nFeature_RNA[m$tgpos == "neg"]))
}


# ------------------------------------------------------ #
# Published MTJ transcriptome signature from the Kim/Birchmeier 2020 atlas
# (used by report 02 to compare our translatome against a published MTJ
# transcriptome). The MTJ cluster is identified as the cluster with the
# highest fraction of Col22a1+ nuclei (same rule as report 07's Petrany MTJ
# call); its positive markers (FindMarkers, MTJ-cluster vs rest) are the
# signature. Returns gene_name + avg_log2FC + p_val_adj.
# ------------------------------------------------------ #
kim_mtj_signature = function(
    n_top    = 200,
    padj_cut = 0.05
){
    suppressPackageStartupMessages({
        library(Seurat); library(dplyr)
    })
    s = load_kim_birchmeier_2020()
    if(!"Col22a1" %in% rownames(s))
        stop("Col22a1 not found in Kim object; cannot identify the MTJ cluster")

    col = FetchData(s, vars = c("Col22a1", "seurat_clusters"))
    mtj = col %>%
        group_by(seurat_clusters) %>%
        summarize(pct_col22a1 = mean(Col22a1 > 0), .groups = "drop") %>%
        arrange(desc(pct_col22a1)) %>%
        dplyr::slice(1) %>%
        pull(seurat_clusters) %>% as.character()

    Idents(s) = "seurat_clusters"
    mk = FindMarkers(s, ident.1 = mtj, only.pos = TRUE, logfc.threshold = 0.25)
    mk %>%
        tibble::rownames_to_column("gene_name") %>%
        filter(p_val_adj < padj_cut) %>%
        arrange(desc(avg_log2FC)) %>%
        head(n_top) %>%
        dplyr::select(gene_name, avg_log2FC, p_val_adj)
}


# ------------------------------------------------------ #
# Petrany 2020 MTJ signature — identical method to kim_mtj_signature(),
# on the Petrany atlas. The MTJ cluster is the one with the highest fraction
# of Col22a1+ nuclei (auto-picks cluster 26, the established Petrany MTJ);
# its FindMarkers positive markers are the signature.
# Returns gene_name + avg_log2FC + p_val_adj.
# ------------------------------------------------------ #
petrany_mtj_signature = function(
    n_top    = 200,
    padj_cut = 0.05
){
    suppressPackageStartupMessages({
        library(Seurat); library(dplyr)
    })
    s = load_petrany_2020()
    if(!"Col22a1" %in% rownames(s))
        stop("Col22a1 not found in Petrany object; cannot identify the MTJ cluster")

    col = FetchData(s, vars = c("Col22a1", "seurat_clusters"))
    mtj = col %>%
        group_by(seurat_clusters) %>%
        summarize(pct_col22a1 = mean(Col22a1 > 0), .groups = "drop") %>%
        arrange(desc(pct_col22a1)) %>%
        dplyr::slice(1) %>%
        pull(seurat_clusters) %>% as.character()

    Idents(s) = "seurat_clusters"
    mk = FindMarkers(s, ident.1 = mtj, only.pos = TRUE, logfc.threshold = 0.25)
    mk %>%
        tibble::rownames_to_column("gene_name") %>%
        filter(p_val_adj < padj_cut) %>%
        arrange(desc(avg_log2FC)) %>%
        head(n_top) %>%
        dplyr::select(gene_name, avg_log2FC, p_val_adj)
}


# ------------------------------------------------------ #
# Cluster the highest fraction of a marker⁺ nuclei express (e.g. the MTJ cluster
# by Col22a1). Returns the cluster label as a string.
# ------------------------------------------------------ #
top_marker_cluster = function(seurat_obj, gene, group_col = "seurat_clusters"){
    suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
    if(!gene %in% rownames(seurat_obj)) return(NA_character_)
    FetchData(seurat_obj, vars = c(group_col, gene)) %>%
        group_by(.data[[group_col]]) %>%
        summarize(p = mean(.data[[gene]] > 0), .groups = "drop") %>%
        arrange(desc(p)) %>% dplyr::slice(1) %>% pull(1) %>% as.character()
}


# ------------------------------------------------------ #
# Per-cluster mean log-normalized expression (gene x cluster) — for within-dataset
# cluster-vs-cluster scatterplots.
# ------------------------------------------------------ #
cluster_mean_expr = function(seurat_obj, group_col = "seurat_clusters"){
    suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })
    d  = SeuratObject::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    cl = as.character(seurat_obj@meta.data[[group_col]])
    levs = sort(unique(cl))
    m = vapply(levs, function(c) Matrix::rowMeans(d[, cl == c, drop = FALSE]),
               numeric(nrow(d)))
    colnames(m) = levs
    m
}


# ------------------------------------------------------ #
# Cluster-vs-rest log2 fold-change (Seurat::FoldChange) + mean expression A, for MA
# plots and cross-dataset concordance. Returns gene_name / avg_log2FC / pct.1 / pct.2 / A.
# ------------------------------------------------------ #
cluster_vs_rest_lfc = function(seurat_obj, ident, group_col = "seurat_clusters"){
    suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })
    Idents(seurat_obj) = group_col
    fc = Seurat::FoldChange(seurat_obj, ident.1 = as.character(ident))
    A  = Matrix::rowMeans(SeuratObject::GetAssayData(seurat_obj, assay = "RNA", layer = "data"))
    data.frame(gene_name = rownames(fc), avg_log2FC = fc$avg_log2FC,
               pct.1 = fc$pct.1, pct.2 = fc$pct.2, A = A[rownames(fc)],
               row.names = NULL)
}


# ------------------------------------------------------ #
# Tigd4 (+ control markers) summary per cluster.
# dataset_tag is used to prefix the output so multiple
# datasets can be bound together for a joint table.
# ------------------------------------------------------ #
marker_expression_by_cluster = function(
    seurat_obj,
    cluster_col = "seurat_clusters",
    genes       = c("Tigd4","Col22a1","Ckm","Pdgfra","Tnmd","Scx","Mkx"),
    dataset_tag = NULL
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(tidyr)
    })

    genes = intersect(genes, rownames(seurat_obj))
    if(length(genes) == 0)
        stop("none of the requested genes found in ", dataset_tag %||% "object")

    expr = FetchData(seurat_obj, vars = c(cluster_col, genes))
    colnames(expr)[1] = "cluster"

    long = expr %>%
        tidyr::pivot_longer(-cluster, names_to = "gene", values_to = "expr") %>%
        group_by(cluster, gene) %>%
        summarize(pct_expressing = mean(expr > 0) * 100,
                  mean_expr      = mean(expr),
                  n_cells        = n(),
                  .groups        = "drop")

    if(!is.null(dataset_tag))
        long$dataset = dataset_tag

    long
}
