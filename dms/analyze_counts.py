#!/usr/bin/env python3
"""Analyze DMS variant counts.

Given the outputs of count_variants.sh for one barcode, this:
  1. exports a tidy per-variant CSV (name parsed into segment / codon / WT / mut),
  2. plots the count distribution and per-segment breakdowns, and
  3. characterises the "other" reads (those NOT counted by exact derep) to work out
     what they are — partial reads, single-error near-variants, or something foreign.

Run:
  uv run python dms/analyze_counts.py \
      --counts   .../variant_counts/countsA_derep.tsv \
      --derep    .../variant_counts/derep.fasta \
      --refs     .../mutational_scanning/variant_refs.fasta \
      --outdir   .../variant_counts/analysis \
      --sample   bc2027
"""
from __future__ import annotations

import argparse
import heapq
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from dms_common import WIN_LEN, VARIANT_RE, iter_fasta, read_fasta


# --------------------------------------------------------------------------- IO
def parse_counts(tsv: Path) -> tuple[pd.DataFrame, int]:
    """Load the name<TAB>count table; split names into columns. Returns
    (variant dataframe, WT read count)."""
    df = pd.read_csv(tsv, sep="\t", names=["name", "count"])
    wt_count = int(df.loc[df["name"] == "WT", "count"].sum())
    df = df[df["name"] != "WT"].copy()
    meta = df["name"].str.extract(VARIANT_RE)
    meta.columns = ["seg_num", "codon", "wt_aa", "mut_aa"]
    df["segment"] = "seg" + meta["seg_num"].str.zfill(2)
    df["codon"] = meta["codon"].astype(int)
    df["wt_aa"] = meta["wt_aa"]
    df["mut_aa"] = meta["mut_aa"]
    return df.sort_values("count", ascending=False), wt_count


# ------------------------------------------------------------------ derep study
def analyze_derep(derep: Path, top_k: int = 50) -> dict:
    """One streaming pass over derep.fasta.

    Clusters seeded by a reference (name starts with DNAP_seg / WT) are the counted
    reads; clusters seeded by a read id are the "other" (uncounted) reads. For the
    "other" set we record read totals, lengths and cluster sizes, and keep the
    `top_k` most abundant full-length sequences for a nearest-variant comparison.
    """
    size_re = re.compile(r";size=(\d+)")
    variant_reads = wt_reads = 0
    other_rows = []  # (size, length) per distinct other-sequence
    heap: list[tuple[int, str]] = []  # (size, seq) capped at top_k, full-length only

    for header, seq in iter_fasta(derep):
        m = size_re.search(header)
        size = int(m.group(1)) if m else 1
        seed = header.split(";")[0]
        if seed.startswith("DNAP_seg"):
            variant_reads += size - 1          # minus the reference itself
        elif seed == "WT":
            wt_reads += size - 1
        else:                                   # read-seeded => "other"
            other_rows.append((size, len(seq)))
            if len(seq) == WIN_LEN:
                if len(heap) < top_k:
                    heapq.heappush(heap, (size, seq))
                elif size > heap[0][0]:
                    heapq.heapreplace(heap, (size, seq))

    other = pd.DataFrame(other_rows, columns=["size", "length"])
    return {
        "variant_reads": variant_reads,
        "wt_reads": wt_reads,
        "other": other,
        "top_full_other": sorted(heap, reverse=True),  # [(size, seq), ...]
    }


def nearest_variant(query_seqs: list[str], refs: dict[str, str]) -> pd.DataFrame:
    """For each query, find the reference with the fewest mismatches (numpy)."""
    names = [n for n, s in refs.items() if len(s) == WIN_LEN]
    ref_arr = np.frombuffer("".join(refs[n] for n in names).encode(), dtype=np.uint8)
    ref_arr = ref_arr.reshape(len(names), WIN_LEN)
    wt = np.frombuffer(refs["WT"].encode(), dtype=np.uint8) if "WT" in refs else None

    rows = []
    for seq in query_seqs:
        q = np.frombuffer(seq.encode(), dtype=np.uint8)
        d = (ref_arr != q).sum(axis=1)
        j = int(d.argmin())
        rows.append({
            "nearest_variant": names[j],
            "mismatches_to_nearest": int(d[j]),
            "mismatches_to_WT": int((wt != q).sum()) if wt is not None else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------- plotting
def plot_all(df: pd.DataFrame, wt_count: int, study: dict, outdir: Path, sample: str):
    sns.set_theme(style="whitegrid")
    total_reads = study["variant_reads"] + study["wt_reads"] + int(study["other"]["size"].sum())

    # 1) per-variant count distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(df["count"], bins=60, ax=ax, color="#4C72B0")
    med = df["count"].median()
    ax.axvline(med, color="crimson", ls="--", label=f"median = {med:.0f}")
    ax.set(title=f"{sample}: reads per variant (n={len(df):,})",
           xlabel="exact reads per variant", ylabel="number of variants")
    ax.legend()
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_dist_per_variant.png", dpi=150); plt.close(fig)

    # 2) per-segment total counts
    seg_tot = df.groupby("segment")["count"].sum().reset_index()
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.barplot(data=seg_tot, x="segment", y="count", ax=ax, color="#55A868")
    ax.set(title=f"{sample}: total exact reads per segment", xlabel="", ylabel="reads")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_segment_totals.png", dpi=150); plt.close(fig)

    # 3) per-segment distribution of per-variant counts
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.boxplot(data=df, x="segment", y="count", ax=ax, color="#C7CCE0", fliersize=1)
    ax.set(title=f"{sample}: per-variant count spread by segment", xlabel="", ylabel="reads per variant")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_segment_boxplot.png", dpi=150); plt.close(fig)

    # 4) read-category breakdown (where did ALL reads go?)
    other = study["other"]
    other_full = int(other.loc[other["length"] == WIN_LEN, "size"].sum())
    other_partial = int(other.loc[other["length"] != WIN_LEN, "size"].sum())
    cats = pd.DataFrame({
        "category": ["exact variant", "exact WT", "other: full-length\n(has error)", "other: partial length"],
        "reads": [study["variant_reads"], study["wt_reads"], other_full, other_partial],
    })
    cats["pct"] = 100 * cats["reads"] / total_reads
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=cats, x="reads", y="category", ax=ax, color="#8172B3")
    for i, r in cats.iterrows():
        ax.text(r["reads"], i, f'  {r["reads"]:,} ({r["pct"]:.0f}%)', va="center")
    ax.set(title=f"{sample}: where all {total_reads:,} reads went", xlabel="reads", ylabel="")
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_read_categories.png", dpi=150); plt.close(fig)

    # 5) 'other' reads: length spectrum + cluster-size structure
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.histplot(other, x="length", weights="size", bins=80, ax=axes[0], color="#CCB974")
    axes[0].axvline(WIN_LEN, color="crimson", ls="--", label=f"full length {WIN_LEN}")
    axes[0].set(title=f"{sample}: 'other' reads by length", xlabel="read length (bp)", ylabel="reads"); axes[0].legend()
    recurrent = other.assign(kind=np.where(other["size"] > 1, "recurrent (size>1)", "singleton"))
    grp = recurrent.groupby("kind")["size"].sum().reset_index()
    sns.barplot(data=grp, x="kind", y="size", ax=axes[1], color="#CCB974")
    axes[1].set(title=f"{sample}: 'other' reads by sequence recurrence", xlabel="", ylabel="reads")
    fig.tight_layout(); fig.savefig(outdir / f"{sample}_other_reads.png", dpi=150); plt.close(fig)

    return cats, total_reads


# -------------------------------------------------------------------------- main
def main(counts: Path, derep: Path, refs: Path, outdir: Path, sample: str, top_k: int = 50):
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"=== DMS count analysis — {sample} ===")

    df, wt_count = parse_counts(counts)
    df.to_csv(outdir / f"{sample}_variant_counts.csv", index=False)
    print(f"variants: {len(df):,} | WT exact reads: {wt_count:,} | CSV -> {sample}_variant_counts.csv")

    print("streaming derep.fasta (categorising counted vs other reads)…")
    study = analyze_derep(derep, top_k=top_k)

    cats, total_reads = plot_all(df, wt_count, study, outdir, sample)

    # nearest-variant trace for the most abundant full-length 'other' sequences
    ref_map = read_fasta(refs)
    top_seqs = [s for _, s in study["top_full_other"]]
    if top_seqs:
        near = nearest_variant(top_seqs, ref_map)
        near.insert(0, "reads", [sz for sz, _ in study["top_full_other"]])
        near.to_csv(outdir / f"{sample}_other_top{top_k}.csv", index=False)

    # ---- summaries ----
    print(f"\n{sample}: read fate")
    print(cats[["category", "reads", "pct"]].assign(
        category=lambda d: d["category"].str.replace("\n", " ")).to_string(index=False))
    print(f"total: {total_reads:,}")

    seg = df.groupby("segment")["count"].agg(["sum", "median", "min", "max"]).reset_index()
    seg.to_csv(outdir / f"{sample}_segment_summary.csv", index=False)
    print(f"\n{sample}: per-segment counts")
    print(seg.to_string(index=False))

    if top_seqs:
        print(f"\nTop {top_k} 'other' full-length sequences — nearest designed variant:")
        for d, n in near["mismatches_to_nearest"].value_counts().sort_index().items():
            print(f"  {n} of top {top_k} are {d} mismatch(es) from a designed variant")

    print(f"done — plots + CSVs in {outdir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", type=Path, required=True, help="countsA_derep.tsv")
    p.add_argument("--derep", type=Path, required=True, help="derep.fasta from the count run")
    p.add_argument("--refs", type=Path, required=True, help="variant_refs.fasta")
    p.add_argument("--outdir", type=Path, required=True, help="output directory for CSV + plots")
    p.add_argument("--sample", required=True, help="label, e.g. bc2027")
    p.add_argument("--top-k", dest="top_k", type=int, default=50,
                   help="# of most-abundant 'other' seqs to trace to nearest variant")
    a = p.parse_args()
    main(a.counts, a.derep, a.refs, a.outdir, a.sample, a.top_k)
