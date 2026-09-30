#!/usr/bin/env bash
# ============================================================================
# run_all.sh — DOCUMENTED RUN RECIPE (archival, not a standalone runner)
# ============================================================================
# This Zenodo bundle is a CODE + PROVENANCE snapshot. It is NOT self-executing:
# the reports resolve paths through a project-root global and read raw data from
# absolute /data/local/ paths (see datasets.yaml). To reproduce the analysis you
# need the full project layout and the external datasets. The steps below are the
# recipe, kept as commands for reference rather than run blindly.
#
# Prerequisites
#   - R 4.5.0 with Bioconductor 3.21 (see SESSION.md / renv.lock)
#   - the renv library restored from renv.lock
#   - every dataset in datasets.yaml obtained and placed at its documented path
#
# 1. Restore the environment (from the project root, where renv.lock lives):
#        R --vanilla -e 'renv::restore()'
#    (On the original host R 4.5 is invoked as `R45`; --vanilla skips the
#     personal .Rprofile. Adjust to your R 4.5 launcher.)
#
# 2. Render every numbered report into Results/<report>/<yymmdd>/:
#        R --vanilla -e 'source("Scripts/Bin/knit_all.R"); knit_all()'
#    (In this bundle knit_all.R is copied under Reports/ for reference; the
#     canonical copy lives at Scripts/Bin/knit_all.R in the full project tree,
#     which is the path knit_all() expects. Run from the project root.)
#
# 3. To render a single report, e.g. report 02:
#        R --vanilla -e 'source("Scripts/Bin/knit_all.R"); knit_all(include="02_sedentary_DE")'
#
# See README.md for the directory map, dataset acquisition, and provenance.
# ============================================================================

echo "run_all.sh is a documented recipe, not an executable runner."
echo "Open this file and follow the commented steps. See README.md and datasets.yaml."
exit 0
