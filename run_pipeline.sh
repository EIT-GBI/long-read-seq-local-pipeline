#!/usr/bin/env bash
# Multi-sample Illumina pipeline (QC -> align -> sort/index -> call/filter -> export -> consensus)
# - Parallelizes ACROSS samples with xargs -P
# - Uses THREADS_PER_SAMPLE WITHIN each sample for bwa/samtools/fastqc
# - Organizes outputs into qc/, alignment/, variants/, consensus/, logs/, tmp/
#

set -euo pipefail

# Load config
source config.txt

# Output directory layout
QC_DIR="$OUTPUT_DIR/qc"
FASTQC_DIR="$QC_DIR/fastqc"
FLAGSTAT_DIR="$QC_DIR/flagstat"
TRIM_DIR="$OUTPUT_DIR/trimmed"

ALIGN_DIR="$OUTPUT_DIR/alignment"
BAM_DIR="$ALIGN_DIR/bam"

VARIANT_DIR="$OUTPUT_DIR/variants"
BCF_DIR="$VARIANT_DIR/bcf"
VCF_DIR="$VARIANT_DIR/vcf"
CSV_DIR="$VARIANT_DIR/csv"


CONSENSUS_DIR="$OUTPUT_DIR/consensus"

TMP_DIR="$OUTPUT_DIR/tmp"
LOG_DIR="$OUTPUT_DIR/logs"

mkdir -p \
  "$FASTQC_DIR" "$FLAGSTAT_DIR" \
  "$BAM_DIR" \
  "$TRIM_DIR" \
  "$BCF_DIR" "$VCF_DIR" "$CSV_DIR" \
  "$CONSENSUS_DIR" \
  "$TMP_DIR" "$LOG_DIR"


# Reference indexing (once)
if [[ ! -f "${REFERENCE_GENOME}.bwt" ]]; then
  echo "[INFO] Indexing reference for bwa: $REFERENCE_GENOME"
  bwa index "$REFERENCE_GENOME"
fi

if [[ ! -f "${REFERENCE_GENOME}.fai" ]]; then
  echo "[INFO] Indexing reference for samtools: $REFERENCE_GENOME"
  samtools faidx "$REFERENCE_GENOME"
fi

# Helpers: sample name and R2 path
get_sample_name() {
  local r1="$1"
  local b
  b="$(basename "$r1")"
  b="${b%%_R1*}"
  b="${b%%_1.fastq*}"
  b="${b%%_1.fq*}"
  echo "$b"
}

get_r2_path() {
  local r1="$1"
  local r2="$r1"
  if [[ "$r1" == *"_R1"* ]]; then
    r2="${r2/_R1/_R2}"
  elif [[ "$r1" == *"_1.fastq.gz"* ]]; then
    r2="${r2/_1.fastq.gz/_2.fastq.gz}"
  elif [[ "$r1" == *"_1.fastq"* ]]; then
    r2="${r2/_1.fastq/_2.fastq}"
  elif [[ "$r1" == *"_1.fq.gz"* ]]; then
    r2="${r2/_1.fq.gz/_2.fq.gz}"
  elif [[ "$r1" == *"_1.fq"* ]]; then
    r2="${r2/_1.fq/_2.fq}"
  fi
  echo "$r2"
}

# Process one sample function
process_one() {
  local R1_FILE="$1"
  local SAMPLENAME R2_FILE LOG
  SAMPLENAME="$(get_sample_name "$R1_FILE")"
  R2_FILE="$(get_r2_path "$R1_FILE")"
  LOG="$LOG_DIR/${SAMPLENAME}.log"

  if [[ ! -f "$R2_FILE" ]]; then
    echo "[WARN] Missing R2 for '$SAMPLENAME': $R2_FILE (skipping)"
    return 0
  fi

  echo "[INFO] ===== Sample: $SAMPLENAME ====="
  echo "[INFO] R1: $R1_FILE / R2: $R2_FILE"
  echo "[INFO] Log: $LOG"

  # Basic gzip integrity checks
  if [[ "$R1_FILE" == *.gz ]]; then gzip -t "$R1_FILE" >>"$LOG" 2>&1; fi
  if [[ "$R2_FILE" == *.gz ]]; then gzip -t "$R2_FILE" >>"$LOG" 2>&1; fi

  # trimming first
  local R1_TRIMMED="$TRIM_DIR/${SAMPLENAME}_R1.trimmed.fastq"
  local R2_TRIMMED="$TRIM_DIR/${SAMPLENAME}_R2.trimmed.fastq"
  echo "[INFO] Trimming: $SAMPLENAME"
  fastp \
    -w "$THREADS_PER_SAMPLE" \
    -i "$R1_FILE" -I "$R2_FILE" \
    -o "${R1_TRIMMED}" -O "${R2_TRIMMED}" \
    -h "$TRIM_DIR/${SAMPLENAME}.fastp.html" \
    -j "$TRIM_DIR/${SAMPLENAME}.fastp.json" >>"$LOG" 2>&1

  # Pre-alignment QC
  if [[ "$RUN_FASTQC" == "1" ]]; then
    echo "[INFO] FastQC: $SAMPLENAME"
    fastqc -t "$THREADS_PER_SAMPLE" -o "$FASTQC_DIR" "$R1_TRIMMED" "$R2_TRIMMED" >>"$LOG" 2>&1
  fi

  # Align -> sorted BAM 
  echo "[INFO] Align + sorting: $SAMPLENAME"
  local BAM_SORTED="$BAM_DIR/${SAMPLENAME}.sorted.bam"

  bwa mem -t "$THREADS_PER_SAMPLE" \
    "$REFERENCE_GENOME" "$R1_TRIMMED" "$R2_TRIMMED" 2>>"$LOG" \
  | samtools view -@ "$THREADS_PER_SAMPLE" -bS - 2>>"$LOG" \
  | samtools sort -@ "$THREADS_PER_SAMPLE" -m "$SORT_MEM" \
      -T "$TMP_DIR/${SAMPLENAME}" -o "$BAM_SORTED" - >>"$LOG" 2>&1

  # remove trimmed FASTQ to save space
  rm -f "$R1_TRIMMED" "$R2_TRIMMED"

  samtools index -@ "$THREADS_PER_SAMPLE" "$BAM_SORTED" >>"$LOG" 2>&1

  # Post-alignment QC
  samtools flagstat "$BAM_SORTED" > "$FLAGSTAT_DIR/${SAMPLENAME}.flagstat.txt"

  # Variant calling + filtering (haploid)
  echo "[INFO] Call/filter: $SAMPLENAME"
  local BCF_OUT="$BCF_DIR/${SAMPLENAME}.calls.bcf"

  bcftools mpileup \
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

# export functions and vars for xargs to find them
export -f process_one get_sample_name get_r2_path
export INPUT_DIR OUTPUT_DIR REFERENCE_GENOME THREADS_PER_SAMPLE SAMPLES_PARALLEL MIN_MAPQ MIN_DEPTH MIN_QUAL RUN_FASTQC SORT_MEM
export QC_DIR FASTQC_DIR FLAGSTAT_DIR ALIGN_DIR BAM_DIR VARIANT_DIR BCF_DIR VCF_DIR CSV_DIR CONSENSUS_DIR TMP_DIR LOG_DIR TRIM_DIR

# find R1 files. if glob matches nothing, return error
shopt -s nullglob

R1_LIST=()
for pat in "${FASTQ_PATTERN_R1[@]}"; do
  for f in "$INPUT_DIR"/$pat; do
    R1_LIST+=("$f")
  done
done

if [[ "${#R1_LIST[@]}" -eq 0 ]]; then
  echo "ERROR: No R1 FASTQ files found in '$INPUT_DIR' with patterns:"
  printf '  - %s\n' "${FASTQ_PATTERN_R1[@]}"
  exit 1
fi

echo "[INFO] Found ${#R1_LIST[@]} R1 files in $INPUT_DIR"
echo "[INFO] Parallel: SAMPLES_PARALLEL=$SAMPLES_PARALLEL, THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE"
echo "[INFO] Filters: MIN_MAPQ=$MIN_MAPQ, MIN_DEPTH=$MIN_DEPTH, MIN_QUAL=$MIN_QUAL"
echo "[INFO] FastQC: RUN_FASTQC=$RUN_FASTQC"
echo "[INFO] Output layout:"
echo "  QC:        $QC_DIR"
echo "  Trimmed:  $TRIM_DIR"
echo "  Alignment: $ALIGN_DIR"
echo "  Variants:  $VARIANT_DIR"
echo "  Consensus: $CONSENSUS_DIR"
echo "  Logs:      $LOG_DIR"
echo "  Tmp:       $TMP_DIR"



# Run across samples in parallel
#printf "%s\n" "${R1_LIST[@]}" \
#  | xargs -P "$SAMPLES_PARALLEL" -n 1 -I {} bash -lc 'process_one "$@"' _ {}

# Let's try GNU parallel instead
printf "%s\n" "${R1_LIST[@]}" \
    | parallel -j "$SAMPLES_PARALLEL" --linebuffer process_one {}

echo "[INFO] All samples completed."
