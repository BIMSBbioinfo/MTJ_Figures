# MKim_MTJ — analysis code archive

Code and processed supporting tables for the myotendinous-junction (MTJ)
translatome / single-nucleus study (Kim & Franke). This is the archival code
deposit that accompanies the manuscript's *Data and Code Availability* statement.

> **AI-assistance disclosure.** This analysis code was written with the help of
> AI (Anthropic's Claude, via Claude Code). All analyses, figures, and
> interpretations were reviewed and validated by the authors. The AI assistance
> covered code authoring, refactoring, and documentation; it did not generate
> data or determine scientific conclusions.

## What this archive is (and is not)

This is a **code + provenance snapshot**, not a turnkey pipeline. The reports
resolve file paths through a project-root global and read raw sequencing data
from absolute host paths; the raw data and the large third-party reference
datasets are **not** included (size and independent provenance — see
`datasets.yaml`). To actually re-run the analysis you need the full project
layout plus those external datasets. What *is* included is every piece of
analysis code and the small in-house tables the reports consume.

## Directory map

```
Scripts/Zenodo/
├── Bin/            Analysis helper functions (R) + provenance shell scripts
│   ├── Sample_Sheet.R              sample-sheet + salmon count loaders
│   ├── Differential_Expression.R   bulk DE (DESeq2)
│   ├── Translation_Efficiency.R    Ribo-Tag / RNA translation-efficiency models
│   ├── Exercise_Classification.R   exercise response classes
│   ├── GO_Enrichment.R             GO enrichment (clusterProfiler)
│   ├── DE_Reporting.R              snRNA pseudobulk DE + modality overlap
│   ├── snRNA_Loaders.R             loaders for the published snRNA atlases
│   ├── Read_Documentation.R        loaders for RBP/POSTAR/mass-spec tables
│   │                               (not used by reports 00–08; see datasets.yaml)
│   ├── build_klf2dn_reference.sh   custom STARsolo reference (mm39 + transgenes)
│   └── starsolo_count.sh           STARsolo counting of the KLF2-DN libraries
├── Reports/        Numbered R Markdown reports (00–08) + knit_all.R
├── Data/           Small processed tables the reports read (~1 MB; see below)
├── renv.lock       Exact package versions (R 4.5.0 / Bioconductor 3.21)
├── SESSION.md      Human-readable environment summary
├── datasets.yaml   Manifest of every dataset the analysis needs (bundled or not)
├── LICENSE         MIT License
├── run_all.sh      Documented run recipe (archival; not self-executing)
└── README.md       This file
```

### Reports

| Report | Question |
|--------|----------|
| `00_Report` | Integration / index — per-reviewer-question deliverable map |
| `01_sample_sheet_QC` | Sample-sheet QC |
| `02_sedentary_DE` | Sedentary Ribo-Tag vs RNA translation efficiency (MTJ vs CK8) |
| `03_bulk_downregulated` | Bulk RNA down-regulated genes + GO |
| `04_exercise_DE` | Exercise contrasts + response classes |
| `05_tigd4_specificity_snRNA` | Tigd4 / MTJ-marker specificity across published snRNA atlases |
| `07_ribotag_vs_snRNA` | Ribo-Tag vs snRNA marker concordance |
| `08_klf2dn_snRNA` | KLF2-DN effect on mouse muscle (paired snRNA) |

The legacy exploratory drafts (`differential_expression.rmd`, `tigd_rip.rmd`)
are intentionally excluded — they are superseded and not part of the pipeline.

### Bundled data (`Data/`)

Only the small, in-house processed tables the numbered reports read: the MTJ
marker list and bioRxiv supplement (xlsx), the KLF2-DN construct FASTAs, the two
sample sheets, the mm10→mm39 liftOver chain, and the original Documentation
README. Everything else the analysis needs — raw sequencing, published snRNA
atlases, and large third-party reference data (e.g. POSTAR3, ~684 MB) — is listed
in `datasets.yaml` with its public source and expected path.

## Reproducing

1. Install R 4.5.0 (Bioconductor 3.21) and restore the library:
   `R --vanilla -e 'renv::restore()'` from the project root.
2. Obtain every dataset marked `bundled: false` in `datasets.yaml` and place it
   at its documented path.
3. Render the reports:
   `R --vanilla -e 'source("Scripts/Bin/knit_all.R"); knit_all()'`

See `run_all.sh` for the annotated command sequence and `SESSION.md` for the
environment details.

## License

Released under the MIT License — see `LICENSE`.

## Citation

Please cite the accompanying manuscript (Kim & Franke) and this archive's Zenodo
DOI. Published datasets re-used here retain their original accessions
(see `datasets.yaml`).
