#!/usr/bin/env bash
# Canvas:     Prompts/canvases/078_klf2dn-snrna-analysis.md
# Purpose:    Build the custom STARsolo reference = mm39 (Ensembl 113) + KLF2-DN + Cre transgene contigs
# Inputs:     Documentation/KLF2_Construct/klf2.fa (KLF2-DN), paav-ck8-cre-nls.fasta (Cre);
#             mm39 Ensembl-113 genome DNA FASTA (downloaded if absent) + GTF
# Outputs:    <REF_DIR>/klf2dn.genome.fa, klf2dn.annotation.gtf, index/ (STAR genome index)
# Upstream:   canvas 078 (transgene FASTAs)
# Downstream: Scripts/Bin/starsolo_count.sh (uses index/); report 08
#
# STARsolo, not CellRanger (decided 2026-07-26; CellRanger not installed, STAR 2.7.8a via guix).
# The two transgenes are appended as single-exon contigs so transduced nuclei can be scored by
# KLF2-DN / Cre reads. NOTE the Cre contig is the FULL pAAV plasmid (4692 bp) annotated as one gene
# "Cre" — sufficient for transduction detection / pseudobulk; refine to the Cre ORF only if per-base
# specificity is ever needed.
#
# Usage:
#   build_klf2dn_reference.sh [--augment-only] \
#       [--ref-dir DIR] [--gtf FILE] [--genome-fa FILE] [--threads N] [--sjdb-overhang N]
#   --augment-only : build the augmented FASTA + GTF and stop (skip the STAR index).
set -euo pipefail

# --- defaults ---------------------------------------------------------------
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ANNOT_DIR="/data/local/vfranke/Annotation/mm39/EnsemblGenes"
REF_DIR="/data/local/vfranke/Annotation/mm39/klf2dn_starsolo"
GTF="$ANNOT_DIR/Mus_musculus.GRCm39.113.chr.gtf"
GENOME_FA="$ANNOT_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GENOME_URL="https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
KLF2_FA="$PROJ/Documentation/KLF2_Construct/klf2.fa"
CRE_FA="$PROJ/Documentation/KLF2_Construct/paav-ck8-cre-nls.fasta"
# KLF2-DN 3'UTR (user-provided 2026-07-28). The original klf2.fa is CDS-only, so 3' snRNA reads
# — which land at the transcript's 3' end — had nothing to map to and KLF2-DN was undetectable.
# Appended to the KLF2-DN contig (kept out of Documentation/ per project rules) to rescue detection.
KLF2_3UTR="CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAAGACAGCAAGGGGGAGGATTGGGAAGAGAATAGCAGGCATGCTGGGGACTGATTTTGTAGGTAACCACGTGCGGACCGAGCGGCCGCAGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGCTGCCTGCAGG"
THREADS=8
SJDB_OVERHANG=84           # readlen-1; this GEM-X v4 run has R2=85bp -> 84 (IGBMC report S26064)
AUGMENT_ONLY=0
TRANSGENE=both             # both|cre|dn — which transgene contig(s) to append.
STAR="${STAR:-STAR}"       # expects STAR on PATH (guix: ~/.guix-profile/bin/STAR)

# CONDITION-SPECIFIC references (decided 2026-07-28): the KLF2-DN 3'UTR is shared with the Cre pAAV
# backbone, so putting BOTH contigs in one reference makes shared-backbone 3' reads multimap between
# them and get discarded (KLF2-DN stays 0). Building a Cre-ONLY reference (for Cre/control legs) and a
# DN-ONLY reference (for KLF2-DN legs) lets those reads map uniquely to the single transgene contig.
# In DN legs both plasmids are co-delivered and share the 3' end, so the DN-ref "KLF2-DN" count is the
# total transduction readout (= Cre + KLF2-DN transcripts) — the user's "Cre = DN-construct readout".

while [ $# -gt 0 ]; do
  case "$1" in
    --augment-only) AUGMENT_ONLY=1; shift;;
    --transgene)    TRANSGENE="$2"; shift 2;;
    --ref-dir)      REF_DIR="$2"; shift 2;;
    --gtf)          GTF="$2"; shift 2;;
    --genome-fa)    GENOME_FA="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
case "$TRANSGENE" in both|cre|dn) ;; *) echo "--transgene must be both|cre|dn" >&2; exit 2;; esac
add_dn()  { [ "$TRANSGENE" = both ] || [ "$TRANSGENE" = dn  ]; }
add_cre() { [ "$TRANSGENE" = both ] || [ "$TRANSGENE" = cre ]; }

mkdir -p "$REF_DIR"
AUG_FA="$REF_DIR/klf2dn.genome.fa"
AUG_GTF="$REF_DIR/klf2dn.annotation.gtf"

# --- 1. transgene contigs: rename headers to clean single-token contig names -------------
# klf2.fa header is ">KLF2-DN"; paav header is ">pAAV CK8-Cre-NLS" (space -> must rename to one token).
tg_fa() {  # $1 = fasta, $2 = contig name, [$3 = extra seq appended (e.g. 3'UTR)]
  awk -v n="$2" 'BEGIN{print ">"n} !/^>/{print toupper($0)}' "$1"
  if [ -n "${3:-}" ]; then echo "${3^^}"; fi   # if/then (not &&) so an empty $3 doesn't return non-zero under set -e
}
seqlen() { awk '!/^>/{L+=length($0)} END{print L}' "$1"; }
KLF2_LEN=$(( $(seqlen "$KLF2_FA") + ${#KLF2_3UTR} )); CRE_LEN=$(seqlen "$CRE_FA")
echo "[build] transgene lengths: KLF2-DN=$KLF2_LEN (CDS $(seqlen "$KLF2_FA") + 3'UTR ${#KLF2_3UTR})  Cre=$CRE_LEN"

# --- 2. GTF records: one gene/transcript/exon spanning each contig, + strand -------------
gtf_records() {  # $1 = contig, $2 = length
  local c="$1" L="$2"
  printf '%s\tcustom\tgene\t1\t%s\t.\t+\t.\tgene_id "%s"; gene_name "%s"; gene_biotype "transgene";\n' "$c" "$L" "$c" "$c"
  printf '%s\tcustom\ttranscript\t1\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s-1"; gene_name "%s"; gene_biotype "transgene";\n' "$c" "$L" "$c" "$c" "$c"
  printf '%s\tcustom\texon\t1\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s-1"; exon_number "1"; gene_name "%s"; gene_biotype "transgene";\n' "$c" "$L" "$c" "$c" "$c"
}

# --- 3. download genome if missing -------------------------------------------------------
if [ "$AUGMENT_ONLY" -eq 0 ] && [ ! -s "$GENOME_FA" ]; then
  echo "[build] genome DNA FASTA absent -> downloading $GENOME_URL"
  curl -fSL "$GENOME_URL" -o "$GENOME_FA.gz"
  gunzip -f "$GENOME_FA.gz"
fi

# --- 4. assemble augmented FASTA + GTF ---------------------------------------------------
echo "[build] transgene set = $TRANSGENE"
echo "[build] writing augmented GTF -> $AUG_GTF"
{
  cat "$GTF"
  if add_dn;  then gtf_records "KLF2-DN" "$KLF2_LEN"; fi
  if add_cre; then gtf_records "Cre"     "$CRE_LEN";  fi
} > "$AUG_GTF"

if [ "$AUGMENT_ONLY" -eq 1 ]; then
  echo "[build] --augment-only: writing transgene-only FASTA (no genome) for inspection"
  {
    if add_dn;  then tg_fa "$KLF2_FA" "KLF2-DN" "$KLF2_3UTR"; fi
    if add_cre; then tg_fa "$CRE_FA" "Cre"; fi
  } > "$REF_DIR/transgenes.fa"
  echo "[build] done (augment-only). transgenes.fa + $AUG_GTF written."
  exit 0
fi

echo "[build] writing augmented genome FASTA -> $AUG_FA"
{
  cat "$GENOME_FA"
  if add_dn;  then tg_fa "$KLF2_FA" "KLF2-DN" "$KLF2_3UTR"; fi
  if add_cre; then tg_fa "$CRE_FA" "Cre"; fi
} > "$AUG_FA"

# --- 5. STAR genome index ---------------------------------------------------------------
echo "[build] STAR genomeGenerate (threads=$THREADS, sjdbOverhang=$SJDB_OVERHANG)"
mkdir -p "$REF_DIR/index"
"$STAR" --runMode genomeGenerate \
     --genomeDir "$REF_DIR/index" \
     --genomeFastaFiles "$AUG_FA" \
     --sjdbGTFfile "$AUG_GTF" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --genomeSAindexNbases 14 \
     --runThreadN "$THREADS"
echo "[build] done. index at $REF_DIR/index"
