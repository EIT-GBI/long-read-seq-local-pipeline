#!/usr/bin/env bash
# Run the DMS analysis locally for every aligned BAM in a folder.
# For each *.sorted.bam it calls run_barcode.sh (trim -> count -> analyze -> mismatch),
# writing per-barcode results to OUT_BASE/<bc>_result/.
#
#   usage: dms/run_local_batch.sh BAM_DIR OUT_BASE
#     BAM_DIR   folder with the downloaded bcXXXX.sorted.bam files
#     OUT_BASE  where to put the <bc>_result/ folders
#
#   env:
#     REFS   variant_refs.fasta (default: local mutational_scanning path)
#     WIN_START/WIN_END/THREADS  passed through to run_barcode.sh
#
#   example:
#     dms/run_local_batch.sh ~/Data/KhoaChung/r84334_20260723/bams \
#                            ~/Data/KhoaChung/r84334_20260723/dms
set -euo pipefail

BAM_DIR="${1:?usage: run_local_batch.sh BAM_DIR OUT_BASE}"
OUT_BASE="${2:?need OUT_BASE}"
REFS="${REFS:-/Users/cristian.soitu/Data/KhoaChung/mutational_scanning/variant_refs.fasta}"
HERE="$(cd "$(dirname "$0")" && pwd)"

[[ -f "$REFS" ]] || { echo "ERROR: refs not found: $REFS" >&2; exit 1; }

shopt -s nullglob
bams=("$BAM_DIR"/*.sorted.bam)
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: no *.sorted.bam in $BAM_DIR" >&2; exit 1; }
echo "[INFO] ${#bams[@]} BAM(s) in $BAM_DIR | refs=$REFS"

ok=(); failed=()
for bam in "${bams[@]}"; do
  name="$(basename "$bam")"; name="${name%.sorted.bam}"; name="${name%.bam}"
  if [[ "$name" =~ (bc[0-9]+) ]]; then bc="${BASH_REMATCH[1]}"; else bc="$name"; fi
  echo "==== $bc ===="
  # guard against partial/corrupt files (e.g. a download still in progress)
  if ! samtools quickcheck "$bam" 2>/dev/null; then
    echo "[WARN] $bc: BAM failed quickcheck (incomplete/corrupt?) — skipping" >&2
    failed+=("$bc"); continue
  fi
  if REFS="$REFS" bash "$HERE/run_barcode.sh" "$bam" "$bc" "$OUT_BASE/${bc}_result"; then
    ok+=("$bc")
  else
    echo "[WARN] $bc failed" >&2
    failed+=("$bc")
  fi
done

echo "======================================================================"
echo "[DONE] ${#ok[@]} ok: ${ok[*]:-none}"
[[ ${#failed[@]} -gt 0 ]] && echo "[FAIL] ${#failed[@]}: ${failed[*]}"
echo "results under $OUT_BASE/<bc>_result/variant_counts/analysis/"
