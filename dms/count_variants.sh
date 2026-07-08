#!/usr/bin/env bash
# Count reads supporting each DNAP variant with vsearch --derep_fulllength.
#
#   usage: count_variants.sh BAM refs.fasta OUTDIR [--with-global]
#
#   BAM         trimmed+sorted BAM, reads hard-clipped to the DNAP window 1790..3448
#   refs.fasta  per-variant references from build_variant_refs.py (includes >WT)
#   OUTDIR      output directory
#   --with-global  ALSO run the (slow) usearch_global best-hit mapping. Off by default:
#                  against ~10k near-identical refs it is exhaustive and impractical
#                  at full scale (see the counting notes).
#
# METHOD (derep, exact): concatenate refs FIRST then reads and derep_fulllength.
#   vsearch keeps the first member of each identical group as the seed, so refs
#   become the seeds and --sizeout counts identical reads. reads_for_variant = size-1.
#
#   This counts ONLY byte-identical, full-length, error-free reads. Partial reads
#   (~37% here) and any read with a single HiFi error are not counted. It is an exact
#   LOWER BOUND, biased toward the cleanest reads. We report the uncounted fraction so
#   the bias is explicit.
set -euo pipefail

BAM="${1:?usage: count_variants.sh BAM refs.fasta OUTDIR [--with-global]}"
REFS="${2:?need variant_refs.fasta}"
OUT="${3:?need output dir}"
WITH_GLOBAL=0
[[ "${4:-}" == "--with-global" ]] && WITH_GLOBAL=1
THREADS="${THREADS:-8}"
mkdir -p "$OUT"

# 0) reads -> FASTA.
#    NB: samtools fasta reverse-complements reverse-strand reads back to their
#    ORIGINAL sequencing orientation, so ~half the reads come out as the minus
#    (coding) strand. That is why the derep/search steps below use `--strand both`
#    -- with `--strand plus` every reverse-oriented read is missed (it lands in
#    "other" as an exact RC of a reference, which is what we initially saw).
samtools fasta -@ "$THREADS" "$BAM" > "$OUT/reads.fasta" 2>/dev/null
N_READS=$(grep -c '^>' "$OUT/reads.fasta")
echo "[INFO] reads: $N_READS"

# full list of variant/WT names (to emit zero-count dropouts too)
grep '^>' "$REFS" | sed 's/^>//; s/;.*//' | sort > "$OUT/.all_names"
N_REFS=$(wc -l < "$OUT/.all_names" | tr -d ' ')

# 1) derep exact, refs first
cat "$REFS" "$OUT/reads.fasta" \
| vsearch --derep_fulllength - --strand both --sizeout --minseqlength 1 \
    --output "$OUT/derep.fasta" 2>"$OUT/derep.log"

# 2) per-variant read counts (size-1) for ref-seeded clusters
grep '^>' "$OUT/derep.fasta" \
  | sed -E 's/^>([^;]+);size=([0-9]+).*/\1\t\2/' \
  | awk -F'\t' 'NF==2 && $1 ~ /^(DNAP_seg|WT)/ {print $1"\t"($2-1)}' \
  > "$OUT/.counted"

# 3) merge with the full name list so every variant appears (0 if unseen), sort desc
awk -F'\t' 'NR==FNR{c[$1]=$2; next}{print $1"\t"(($1 in c)?c[$1]:0)}' \
    "$OUT/.counted" "$OUT/.all_names" \
  | sort -t$'\t' -k2,2nr > "$OUT/countsA_derep.tsv"

# 4) summary
COUNTED=$(awk -F'\t' '$1!="WT"{s+=$2} END{print s+0}' "$OUT/countsA_derep.tsv")
WT=$(awk -F'\t' '$1=="WT"{print $2+0}' "$OUT/countsA_derep.tsv")
NONZERO=$(awk -F'\t' '$1!="WT" && $2>0' "$OUT/countsA_derep.tsv" | wc -l | tr -d ' ')
ZERO=$(awk -F'\t' '$1!="WT" && $2==0' "$OUT/countsA_derep.tsv" | wc -l | tr -d ' ')
{
  echo "reads_total            $N_READS"
  echo "reads_counted_variant  $COUNTED"
  echo "reads_counted_WT       $WT"
  echo "reads_uncounted        $((N_READS - COUNTED - WT))   (partial or with >=1 error)"
  echo "variants_total         $((N_REFS-1))"
  echo "variants_with_reads    $NONZERO"
  echo "variants_zero_dropout  $ZERO"
} | tee "$OUT/summaryA.txt"
echo "[INFO] per-variant counts -> $OUT/countsA_derep.tsv"
rm -f "$OUT/.all_names" "$OUT/.counted"

# 5) optional slow best-hit mapping
if [[ "$WITH_GLOBAL" == "1" ]]; then
  echo "[INFO] running usearch_global (slow; exhaustive over near-identical refs)…"
  vsearch --usearch_global "$OUT/reads.fasta" --db "$REFS" \
      --id 0.98 --strand both --maxaccepts 0 --maxrejects 0 --top_hits_only \
      --threads "$THREADS" --userout "$OUT/hits.tsv" --userfields query+target+id 2>"$OUT/global.log"
  cut -f1 "$OUT/hits.tsv" | sort | uniq -c > "$OUT/.qh"
  awk '$1==1{print $2}' "$OUT/.qh" > "$OUT/.uq"
  LC_ALL=C sort -k1,1 "$OUT/hits.tsv" > "$OUT/.hs"
  LC_ALL=C join <(sort "$OUT/.uq") "$OUT/.hs" | awk '{print $2}' | sort | uniq -c \
    | sort -rn | awk '{print $2"\t"$1}' > "$OUT/countsB_global.tsv"
  echo "[INFO] best-hit counts -> $OUT/countsB_global.tsv (ambiguous: $(awk '$1>1' "$OUT/.qh" | wc -l | tr -d ' '))"
  rm -f "$OUT/.qh" "$OUT/.uq" "$OUT/.hs"
fi
echo "[DONE]"
