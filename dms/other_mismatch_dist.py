#!/usr/bin/env python3
"""Distance-to-nearest-reference for the full-length 'other' reads.

For every full-length (1659 bp) read that did NOT exactly match a variant or WT,
compute the minimum number of mismatches to the closest reference (any designed
variant OR WT), and report how many are 1, 2, 3, or >3 mismatches away.

Exact shortcut (so we don't compare 60k reads x 10k refs base-by-base):
each reference is WT with at most one codon changed, therefore
    dist(S, variant_v) = dist(S, WT) - hd(S_c, WT_c) + hd(S_c, alt_c)
where c is v's mutated codon. So we only inspect codons where S differs from WT.
Both read orientations are considered (reads come off in either direction).

Run:
  uv run python dms/other_mismatch_dist.py \
      --derep  .../variant_counts/derep.fasta \
      --refs   .../mutational_scanning/variant_refs.fasta \
      --outdir .../variant_counts/analysis --sample bc2027
"""
from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from dms_common import WIN_LEN, iter_fasta, load_codon_model, nearest_ref_dist, revcomp


def main(derep: Path, refs: Path, outdir: Path, sample: str):
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"=== {sample}: nearest-reference distance for full-length 'other' reads ===")

    wt, alts = load_codon_model(refs)
    print(f"model built: WT + {sum(len(v) for v in alts.values()):,} variant codons "
          f"across {len(alts)} positions")

    size_re = re.compile(r";size=(\d+)")
    by_reads: Counter[int] = Counter()   # weighted by read count
    by_seq: Counter[int] = Counter()     # distinct sequences
    n_seq = n_reads = 0
    for header, seq in iter_fasta(derep):
        seed = header.split(";")[0]
        if seed.startswith("DNAP_seg") or seed == "WT" or len(seq) != WIN_LEN:
            continue                      # counted, or partial -> skip
        size = int(m.group(1)) if (m := size_re.search(header)) else 1
        d = min(nearest_ref_dist(seq, wt, alts),
                nearest_ref_dist(revcomp(seq), wt, alts))
        by_reads[d] += size
        by_seq[d] += 1
        n_seq += 1
        n_reads += size

    # tabulate 1,2,3,>3
    def bucket(counter):
        out = {k: counter.get(k, 0) for k in (1, 2, 3)}
        out[">3"] = sum(v for k, v in counter.items() if k > 3)
        return out
    b_reads, b_seq = bucket(by_reads), bucket(by_seq)

    df = pd.DataFrame({"mismatches": ["1", "2", "3", ">3"],
                       "reads": [b_reads[k] for k in (1, 2, 3, ">3")],
                       "reads_pct": [round(100*b_reads[k]/max(n_reads, 1), 1) for k in (1, 2, 3, ">3")],
                       "distinct_seqs": [b_seq[k] for k in (1, 2, 3, ">3")]})
    print(f"\n{sample}: full-length 'other' reads by mismatches to nearest ref (variant or WT)")
    print(df.to_string(index=False))
    print(f"total: reads={n_reads:,}  distinct_seqs={n_seq:,}")
    df.to_csv(outdir / f"{sample}_other_fulllen_mismatch_dist.csv", index=False)

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.barplot(data=df, x="mismatches", y="reads", ax=ax, color="#DD8452")
    for i, r in df.iterrows():
        ax.text(i, r["reads"], f'{int(r["reads"]):,}', ha="center", va="bottom")
    ax.set(title=f"{sample}: full-length 'other' reads vs nearest reference",
           xlabel="mismatches to nearest variant/WT", ylabel="reads")
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_other_fulllen_mismatch_dist.png", dpi=150); plt.close(fig)
    print(f"done -> {outdir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--derep", type=Path, required=True)
    p.add_argument("--refs", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample", required=True)
    a = p.parse_args()
    main(a.derep, a.refs, a.outdir, a.sample)
