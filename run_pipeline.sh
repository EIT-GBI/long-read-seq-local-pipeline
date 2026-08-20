#!/usr/bin/env bash
# Multi-sample long-read pipeline (filter -> align -> sort/index -> QC -> call/filter -> export -> consensus)
# Supports Nanopore (FASTQ) and PacBio HiFi (unaligned BAM or FASTQ).
# - Parallelizes ACROSS samples with GNU parallel
# - Uses THREADS_PER_SAMPLE WITHIN each sample for minimap2/samtools/chopper
# - Organizes outputs into qc/, alignment/, variants/, consensus/, logs/, tmp/
#
# Reads are single-end (no R1/R2). uBAM input is streamed through samtools fastq
# into the aligner without writing an intermediate FASTQ.

set -euo pipefail

# Load config (default: config.txt, or pass path as first argument)
CONFIG_FILE="${1:-config.txt}"
source "$CONFIG_FILE"

# Any args after the config path are treated as explicit input files to run,
# e.g.:  run_pipeline.sh config.txt /path/a.bam /path/b.bam
# (set -u safe expansion of a possibly-empty list)
CLI_FILES=("${@:2}")

# Platform -> minimap2 preset + bcftools mpileup profile
case "$PLATFORM" in
  nanopore)
    MM2_PRESET="map-ont"
    MPILEUP_X="ont"
    ;;
  pacbio-hifi)
    MM2_PRESET="map-hifi"
    MPILEUP_X="pacbio-ccs"
    ;;
  *)
    echo "ERROR: PLATFORM must be 'nanopore' or 'pacbio-hifi' (got '$PLATFORM')"
    exit 1
    ;;
esac

# Output directory layout
QC_DIR="$OUTPUT_DIR/qc"
NANOPLOT_DIR="$QC_DIR/nanoplot"
FLAGSTAT_DIR="$QC_DIR/flagstat"

ALIGN_DIR="$OUTPUT_DIR/alignment"
BAM_DIR="$ALIGN_DIR/bam"

VARIANT_DIR="$OUTPUT_DIR/variants"
BCF_DIR="$VARIANT_DIR/bcf"
VCF_DIR="$VARIANT_DIR/vcf"
CSV_DIR="$VARIANT_DIR/csv"

CONSENSUS_DIR="$OUTPUT_DIR/consensus"

BW_DIR="$OUTPUT_DIR/bigwig"

TMP_DIR="$OUTPUT_DIR/tmp"
LOG_DIR="$OUTPUT_DIR/logs"

mkdir -p \
  "$NANOPLOT_DIR" "$FLAGSTAT_DIR" \
  "$BAM_DIR" \
  "$BCF_DIR" "$VCF_DIR" "$CSV_DIR" \
  "$CONSENSUS_DIR" \
  "$BW_DIR" \
  "$TMP_DIR" "$LOG_DIR"

# Reference indexing (once)
# minimap2 index is preset-specific, so name it accordingly.
MMI="${REFERENCE_GENOME}.${MM2_PRESET}.mmi"
if [[ ! -f "$MMI" ]]; then
  echo "[INFO] Building minimap2 index ($MM2_PRESET): $MMI"
  minimap2 -x "$MM2_PRESET" -d "$MMI" "$REFERENCE_GENOME"
fi

if [[ ! -f "${REFERENCE_GENOME}.fai" ]]; then
  echo "[INFO] Indexing reference for samtools: $REFERENCE_GENOME"
  samtools faidx "$REFERENCE_GENOME"
fi

CHROM_SIZES="${REFERENCE_GENOME}.chrom.sizes"
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "[INFO] Generating chrom sizes: $CHROM_SIZES"
  cut -f1,2 "${REFERENCE_GENOME}.fai" > "$CHROM_SIZES"
fi

BEDGRAPHTOBIGWIG="bedGraphToBigWig"

# Helper: derive sample name by stripping read/archive extensions.
# PacBio demux files end in .bcNNNN -> use just that barcode as the sample name
# (e.g. m84334_..._s4.hifi_reads.bc2033.bam -> bc2033). Otherwise keep the full stem.
get_sample_name() {
  local b
  b="$(basename "$1")"
  b="${b%.gz}"
  b="${b%.fastq}"
  b="${b%.fq}"
  b="${b%.bam}"
  if [[ "$b" =~ (bc[0-9]+)$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "$b"
  fi
}

# Process one sample function
process_one() {
  local INPUT_FILE="$1"
  local SAMPLENAME LOG
  SAMPLENAME="$(get_sample_name "$INPUT_FILE")"
  LOG="$LOG_DIR/${SAMPLENAME}.log"

  echo "[INFO] ===== Sample: $SAMPLENAME ====="
  echo "[INFO] Input: $INPUT_FILE"
  echo "[INFO] Log: $LOG"

  # Build a reads-producing command depending on input type.
  # uBAM -> samtools fastq (streamed); FASTQ -> cat. Either way piped to chopper.
  local -a READS_CMD
  case "$INPUT_FILE" in
    *.bam)
      echo "[INFO] uBAM input; extracting reads with samtools fastq (no intermediate FASTQ)"
      READS_CMD=(samtools fastq -@ "$THREADS_PER_SAMPLE" "$INPUT_FILE")
      ;;
    *.fastq|*.fq|*.fastq.gz|*.fq.gz)
      READS_CMD=(cat "$INPUT_FILE")
      ;;
    *)
      echo "[WARN] Unrecognized input type for '$SAMPLENAME': $INPUT_FILE (skipping)"
      return 0
      ;;
  esac

  # Basic gzip integrity check for gzipped FASTQ
  if [[ "$INPUT_FILE" == *.gz ]]; then gzip -t "$INPUT_FILE" >>"$LOG" 2>&1; fi

  # Filter -> align -> sort, fully streamed (no intermediate FASTQ on disk)
  echo "[INFO] Filter + align + sort: $SAMPLENAME (platform=$PLATFORM, preset=$MM2_PRESET)"
  local BAM_SORTED="$BAM_DIR/${SAMPLENAME}.sorted.bam"

  "${READS_CMD[@]}" 2>>"$LOG" \
  | chopper -q "$MIN_READ_QUAL" -l "$MIN_READ_LEN" --threads "$THREADS_PER_SAMPLE" 2>>"$LOG" \
  | minimap2 -a -x "$MM2_PRESET" -t "$THREADS_PER_SAMPLE" "$MMI" - 2>>"$LOG" \
  | samtools sort -@ "$THREADS_PER_SAMPLE" -m "$SORT_MEM" \
      -T "$TMP_DIR/${SAMPLENAME}" -o "$BAM_SORTED" - >>"$LOG" 2>&1

  samtools index -@ "$THREADS_PER_SAMPLE" "$BAM_SORTED" >>"$LOG" 2>&1

  # Long-read QC (alignment-based, so no FASTQ needed)
  if [[ "$RUN_NANOPLOT" == "1" ]]; then
    echo "[INFO] NanoPlot: $SAMPLENAME"
    NanoPlot -t "$THREADS_PER_SAMPLE" --bam "$BAM_SORTED" \
      -o "$NANOPLOT_DIR/$SAMPLENAME" -p "${SAMPLENAME}_" >>"$LOG" 2>&1
  fi

  # BigWig via bedgraph
  echo "[INFO] BigWig: $SAMPLENAME"
  local BG_TMP="$TMP_DIR/${SAMPLENAME}.bedgraph"
  bedtools genomecov -ibam "$BAM_SORTED" -bg \
    | sort -k1,1 -k2,2n > "$BG_TMP"
  "$BEDGRAPHTOBIGWIG" "$BG_TMP" "$CHROM_SIZES" "$BW_DIR/${SAMPLENAME}.bw" >>"$LOG" 2>&1
  rm -f "$BG_TMP"

  # Post-alignment QC
  samtools flagstat "$BAM_SORTED" > "$FLAGSTAT_DIR/${SAMPLENAME}.flagstat.txt"

  # Variant calling + filtering (haploid).
  # NOTE: bcftools' model is short-read oriented; -X tunes indel/error params per
  # platform, but expect noisier indels than a dedicated long-read caller (Clair3/Medaka).
  echo "[INFO] Call/filter: $SAMPLENAME"
  local BCF_OUT="$BCF_DIR/${SAMPLENAME}.calls.bcf"

  bcftools mpileup \
    -X "$MPILEUP_X" \
    -Ou \
    -f "$REFERENCE_GENOME" \
    -q "$MIN_MAPQ" \
    -a AD,DP \
    "$BAM_SORTED" 2>>"$LOG" \
  | bcftools call \
      --ploidy 1 \
      -mv \
      -Ou 2>>"$LOG" \
  | bcftools filter \
      -e "FMT/DP<$MIN_DEPTH || QUAL<$MIN_QUAL" \
      -Ob \
      -o "$BCF_OUT" 2>>"$LOG"

  bcftools index "$BCF_OUT" >>"$LOG" 2>&1

  # IGV-friendly VCF.gz (+tabix)
  echo "[INFO] Export VCF.gz: $SAMPLENAME"
  local VCFGZ="$VCF_DIR/${SAMPLENAME}.calls.vcf.gz"
  bcftools view -Oz -o "$VCFGZ" "$BCF_OUT" >>"$LOG" 2>&1
  tabix -p vcf "$VCFGZ" >>"$LOG" 2>&1

  # Legacy-friendly CSV
  echo "[INFO] Export CSV: $SAMPLENAME"
  local CSV_OUT="$CSV_DIR/${SAMPLENAME}.calls.csv"
  bcftools query \
    -f '%CHROM,%POS,%REF,%ALT,%DP,[%AD],%QUAL,%FILTER\n' \
    "$BCF_OUT" > "$CSV_OUT"

  # macOS sed -i compatibility: create .bak and remove
  sed -i.bak '1i\
CHROM,POS,REF,ALT,DP,AD(ref,alt),QUAL,FILTER
' "$CSV_OUT" && rm -f "${CSV_OUT}.bak"

  # consensus FASTA
  echo "[INFO] Consensus FASTA: $SAMPLENAME"
  bcftools consensus -f "$REFERENCE_GENOME" "$BCF_OUT" > "$CONSENSUS_DIR/${SAMPLENAME}.consensus.fna"

  echo "[INFO] Done: $SAMPLENAME"
}

# export functions and vars for parallel to find them
export -f process_one get_sample_name
export INPUT_DIR OUTPUT_DIR REFERENCE_GENOME THREADS_PER_SAMPLE SAMPLES_PARALLEL MIN_MAPQ MIN_DEPTH MIN_QUAL RUN_NANOPLOT SORT_MEM
export PLATFORM MM2_PRESET MPILEUP_X MMI MIN_READ_LEN MIN_READ_QUAL
export QC_DIR NANOPLOT_DIR FLAGSTAT_DIR ALIGN_DIR BAM_DIR VARIANT_DIR BCF_DIR VCF_DIR CSV_DIR CONSENSUS_DIR BW_DIR TMP_DIR LOG_DIR
export CHROM_SIZES BEDGRAPHTOBIGWIG

# Decide which files to process. Precedence:
#   1) files given on the command line   (run_pipeline.sh config.txt a.bam b.bam)
#   2) INPUT_FILES=( ... ) set in config  (explicit list, absolute or relative)
#   3) otherwise scan INPUT_DIR with INPUT_PATTERN
# (set -u safe expansion for the possibly-unset INPUT_FILES array)
shopt -s nullglob
INPUT_FILES=(${INPUT_FILES[@]+"${INPUT_FILES[@]}"})

INPUT_LIST=()
if [[ "${#CLI_FILES[@]}" -gt 0 ]]; then
  INPUT_LIST=("${CLI_FILES[@]}")
  echo "[INFO] Using ${#INPUT_LIST[@]} file(s) from the command line"
elif [[ "${#INPUT_FILES[@]}" -gt 0 ]]; then
  INPUT_LIST=("${INPUT_FILES[@]}")
  echo "[INFO] Using ${#INPUT_LIST[@]} file(s) from INPUT_FILES in $CONFIG_FILE"
else
  for pat in "${INPUT_PATTERN[@]}"; do
    for f in "$INPUT_DIR"/$pat; do
      INPUT_LIST+=("$f")
    done
  done
  echo "[INFO] Scanned $INPUT_DIR: found ${#INPUT_LIST[@]} input file(s)"
fi

# Validate: non-empty, and every explicit file actually exists.
if [[ "${#INPUT_LIST[@]}" -eq 0 ]]; then
  echo "ERROR: No input reads. Give files on the command line, set INPUT_FILES,"
  echo "       or ensure '$INPUT_DIR' contains files matching:"
  printf '  - %s\n' "${INPUT_PATTERN[@]}"
  exit 1
fi

MISSING=0
for f in "${INPUT_LIST[@]}"; do
  if [[ ! -f "$f" ]]; then echo "ERROR: input file not found: $f"; MISSING=1; fi
done
[[ "$MISSING" -eq 0 ]] || exit 1
echo "[INFO] Platform: PLATFORM=$PLATFORM (preset=$MM2_PRESET, mpileup -X $MPILEUP_X)"
echo "[INFO] Parallel: SAMPLES_PARALLEL=$SAMPLES_PARALLEL, THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE"
echo "[INFO] Read filter: MIN_READ_LEN=$MIN_READ_LEN, MIN_READ_QUAL=$MIN_READ_QUAL"
echo "[INFO] Variant filters: MIN_MAPQ=$MIN_MAPQ, MIN_DEPTH=$MIN_DEPTH, MIN_QUAL=$MIN_QUAL"
echo "[INFO] NanoPlot: RUN_NANOPLOT=$RUN_NANOPLOT"
echo "[INFO] Output layout:"
echo "  QC:        $QC_DIR"
echo "  Alignment: $ALIGN_DIR"
echo "  Variants:  $VARIANT_DIR"
echo "  Consensus: $CONSENSUS_DIR"
echo "  BigWig:    $BW_DIR"
echo "  Logs:      $LOG_DIR"
echo "  Tmp:       $TMP_DIR"
echo

# Run across samples in parallel
printf "%s\n" "${INPUT_LIST[@]}" \
    | parallel -j "$SAMPLES_PARALLEL" --linebuffer process_one {}

echo "[INFO] All samples completed."
