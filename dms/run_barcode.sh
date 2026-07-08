#!/usr/bin/env bash
# End-to-end per-barcode DMS analysis:
#   trim reads to the DNAP window -> count variants (vsearch derep) ->
#   analyze (CSV + plots) -> nearest-reference mismatch distribution.
#
#   usage: run_barcode.sh SORTED_BAM SAMPLE OUTDIR
#     SORTED_BAM  full coordinate-sorted BAM aligned to pFR494_pRT300_rham_wt
#     SAMPLE      label, e.g. bc2026
#     OUTDIR      results directory (created)
#
#   env (with defaults):
#     REFS       variant_refs.fasta                              (REQUIRED)
#     WIN_START  1790     WIN_END 3448     THREADS 8
set -euo pipefail

BAM="${1:?usage: run_barcode.sh SORTED_BAM SAMPLE OUTDIR}"
SAMPLE="${2:?need SAMPLE label}"
OUTDIR="${3:?need OUTDIR}"
REFS="${REFS:?set REFS=/path/variant_refs.fasta}"
WS="${WIN_START:-1790}"; WE="${WIN_END:-3448}"; THREADS="${THREADS:-8}"

DMS_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DMS_DIR/.." && pwd)"        # for `uv run`
mkdir -p "$OUTDIR/bam"

# contig name + length from the BAM header
CTG=$(samtools view -H "$BAM" | awk '/^@SQ/{sub("SN:","",$2); print $2; exit}')
CTGLEN=$(samtools view -H "$BAM" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^LN:/){sub("LN:","",$i); print $i}; exit}')
echo "[INFO] $SAMPLE: contig=$CTG len=$CTGLEN  window=${WS}..${WE}"

# 1) trim to the window: keep reads overlapping it, hard-clip the flanks off each read
FLANKS="$OUTDIR/flanks.bed"
printf "%s\t0\t%d\n%s\t%d\t%d\n" "$CTG" $((WS-1)) "$CTG" "$WE" "$CTGLEN" > "$FLANKS"
TRIM="$OUTDIR/bam/${SAMPLE}.DNAP_${WS}_${WE}.bam"
echo "[INFO] trimming -> $TRIM"
samtools view -b "$BAM" "${CTG}:${WS}-${WE}" 2>/dev/null \
  | samtools ampliconclip --both-ends --hard-clip -b "$FLANKS" -o - - 2>/dev/null \
  | samtools sort -@ "$THREADS" -o "$TRIM" - 2>/dev/null
samtools index "$TRIM"
echo "[INFO] trimmed reads: $(samtools view -c "$TRIM")"

# 2) count variants (derep, both strands)
THREADS="$THREADS" bash "$DMS_DIR/count_variants.sh" "$TRIM" "$REFS" "$OUTDIR/variant_counts"

# 3) analyze: CSV + distribution/per-segment plots + read-fate breakdown
CN="$OUTDIR/variant_counts"
( cd "$REPO" && uv run python dms/analyze_counts.py \
    --counts "$CN/countsA_derep.tsv" --derep "$CN/derep.fasta" \
    --refs "$REFS" --outdir "$CN/analysis" --sample "$SAMPLE" )

# 4) nearest-reference mismatch distribution for full-length 'other' reads
( cd "$REPO" && uv run python dms/other_mismatch_dist.py \
    --derep "$CN/derep.fasta" --refs "$REFS" \
    --outdir "$CN/analysis" --sample "$SAMPLE" )

echo "[DONE] $SAMPLE -> $CN/analysis"
