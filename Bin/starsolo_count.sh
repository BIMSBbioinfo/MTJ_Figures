#!/usr/bin/env bash
# Canvas:     Prompts/canvases/078_klf2dn-snrna-analysis.md
# Purpose:    STARsolo per-library counting for the KLF2-DN 10x 3' snRNA data (GeneFull / intron-inclusive)
# Inputs:     per-library R1/R2 FASTQs; STAR index from build_klf2dn_reference.sh; 10x barcode whitelist
# Outputs:    <OUT>/<sample>/Solo.out/GeneFull/filtered/{matrix.mtx,barcodes.tsv,features.tsv}
# Upstream:   Scripts/Bin/build_klf2dn_reference.sh (index/)
# Downstream: Scripts/Bin/snRNA_Loaders.R::load_klf2dn(); report 08
#
# 10x Chromium 3' GEX. THIS RUN = GEM-X 3' v4 (IGBMC report S26064: kit PN-1000691, R1=28/R2=85).
# CB=16bp always; UMI len auto-detected from R1 length. IMPORTANT: R1 length does NOT determine the
# whitelist — GEM-X v4 and 3' v3.1 are BOTH 28bp (16+12) but use DIFFERENT whitelists, so --whitelist
# is a required arg, not inferred:
#   GEM-X 3' v4 : R1 28bp, UMI 12, whitelist 3M-3pgex-may-2023.txt   (THIS DATA)
#   3' v3/v3.1  : R1 28bp, UMI 12, whitelist 3M-february-2018.txt
#   3' v2       : R1 26bp, UMI 10, whitelist 737K-august-2016.txt
# --soloFeatures GeneFull = pre-mRNA/intron-inclusive counting, the correct mode for single-NUCLEUS
# data (= CellRanger --include-introns).
#
# Usage:
#   starsolo_count.sh --sample S1 --r1 "s1_*_R1_*.fastq.gz" --r2 "s1_*_R2_*.fastq.gz" \
#       --index DIR --whitelist FILE [--umi-len N] [--out DIR] [--threads N]
#   R1/R2 accept comma-separated lists (multi-lane) matching STAR's readFilesIn order.
set -euo pipefail

INDEX="/data/local/vfranke/Annotation/mm39/klf2dn_starsolo/index"
OUT="/data/local/Projects/MKim_MTJ/Results/klf2dn_starsolo_counts"
THREADS=8
UMI_LEN=""            # empty -> auto from R1 length
CELLFILTER="CellRanger2.2 10000 0.99 10"   # 2.7.8a has no EmptyDrops_CR; CR2.2 knee, ~10k expected nuclei
MULTIMAP="Unique"    # --soloMultiMappers mode. Unique (default) = drop multimappers. EM/Uniform/Rescue
                     # rescue reads that multimap; used for DN legs where KLF2-DN CDS multimaps with
                     # endogenous Klf2. EM writes an extra GeneFull/*/UniqueAndMult-EM.mtx matrix.
WHITELIST=""
SAMPLE=""; R1=""; R2=""
STAR="${STAR:-STAR}"

while [ $# -gt 0 ]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --r1) R1="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --index) INDEX="$2"; shift 2;;
    --whitelist) WHITELIST="$2"; shift 2;;
    --umi-len) UMI_LEN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --multimappers) MULTIMAP="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$SAMPLE" ] && [ -n "$R1" ] && [ -n "$R2" ] && [ -n "$WHITELIST" ] || {
  echo "required: --sample --r1 --r2 --whitelist" >&2; exit 2; }

# --- auto-detect chemistry from R1 read length (first fastq, first read) ------------------
# NB: `awk ...exit` closes the pipe early -> zcat gets SIGPIPE; guard so pipefail+errexit
# don't abort the script on that expected 141.
first_r1="${R1%%,*}"
set +o pipefail
R1LEN=$(zcat -f "$first_r1" 2>/dev/null | awk 'NR==2{print length($0); exit}')
set -o pipefail
echo "[count] $SAMPLE: R1 length = $R1LEN"
if [ -z "$UMI_LEN" ]; then
  case "$R1LEN" in
    28) UMI_LEN=12; echo "[count] R1=28 -> UMI 12 (GEM-X v4 or 3' v3.1). Whitelist is NOT inferred — using --whitelist=$WHITELIST (GEM-X v4 = 3M-3pgex-may-2023.txt)";;
    26) UMI_LEN=10; echo "[count] R1=26 -> UMI 10 (3' v2); expect whitelist 737K-august-2016.txt";;
    *)  echo "[count] WARNING: unexpected R1 length $R1LEN; set --umi-len explicitly" >&2;
        UMI_LEN=12;;
  esac
fi

mkdir -p "$OUT/$SAMPLE"
# NB: do NOT cd into the output dir — that would break relative R1/R2 paths. Use --outFileNamePrefix.

# STARsolo: R2 (cDNA) is the first readFilesIn arg, R1 (CB+UMI) the second (STAR convention).
"$STAR" --runMode alignReads \
     --genomeDir "$INDEX" \
     --readFilesIn "$R2" "$R1" \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist "$WHITELIST" \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen "$UMI_LEN" \
     --soloFeatures GeneFull \
     --soloMultiMappers "$MULTIMAP" \
     --soloCellFilter $CELLFILTER \
     --outSAMtype None \
     --runThreadN "$THREADS" \
     --outFileNamePrefix "$OUT/$SAMPLE/"
echo "[count] $SAMPLE done -> $OUT/$SAMPLE/Solo.out/GeneFull/filtered/"
