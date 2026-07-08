#!/usr/bin/env python3
"""Partition error into sequencing vs template (PCR / synthesis) using base QVs,
and use the vector backbone as a no-synthesis control.

A PacBio HiFi per-base quality is the CCS *consensus confidence*, so a substitution
(read vs reference) at LOW QV is a sequencing/CCS error, while one at HIGH QV is a
base read confidently but differing from the reference -> a TEMPLATE error (already
present before sequencing: PCR or oligo synthesis). We bin every aligned base by QV
and compute the substitution-error rate per QV bin; the rate at the top QV bin is the
QV-independent template floor.

Regions ('tracks'):
  * DNAP window (1790..3448): the Twist-synthesised insert. Designed variants are
    excluded (a read's intended codon is detected and not counted as error).
  * backbone (any coords outside the window): NOT Twist-synthesised and carries no
    designed variants, so every mismatch is error. This isolates sequencing+PCR
    error WITHOUT synthesis error. Then:
        DNAP floor - backbone floor  (same sample)  ~ Twist synthesis error
        backbone floor: PCR vs RE-only              ~ PCR polymerase error

Errors are called by comparing each read base to the reference FASTA (no MD tag needed).

Run:
  uv run python dms/error_qv_profile.py \
      --track bc2027:dnap:/…/bc2027.DNAP_1790_3448.bam \
      --track bc2027:3600-5258:/…/bc2027…sorted.bam \
      --track bc2026:dnap:/…/bc2026.DNAP_1790_3448.bam \
      --track bc2026:3600-5258:/…/bc2026…sorted.bam \
      --ref /…/pFR494_pRT300_rham_wt.fa \
      --variant-refs /…/variant_refs.fasta \
      --outdir /…/qv_error --max-reads 150000
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns

from dms_common import CONTIG, WIN0, WINDOW, load_codon_model


def process_track(label, region, bam_path, fasta, alts, max_reads) -> pd.DataFrame:
    is_dnap = region == "dnap"
    if is_dnap:
        start0, end0, region_lab = WIN0, WINDOW[1], "DNAP"
    else:
        a, b = region.split("-")
        start0, end0, region_lab = int(a) - 1, int(b), f"backbone({region})"
    ref_seq = fasta.fetch(CONTIG, start0, end0).upper()

    total = defaultdict(int)
    errors = defaultdict(int)
    bam = pysam.AlignmentFile(bam_path, "rb")
    n = 0
    for read in bam.fetch(CONTIG, start0, end0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        quals, seq = read.query_qualities, read.query_sequence
        if quals is None or seq is None:
            continue
        bases = {}  # rel_pos -> (base, qv, is_mismatch)
        for qpos, rpos in read.get_aligned_pairs():
            if qpos is None or rpos is None or not (start0 <= rpos < end0):
                continue
            rel = rpos - start0
            rb = ref_seq[rel]
            if rb == "N":
                continue
            base, qv = seq[qpos].upper(), quals[qpos]
            bases[rel] = (base, qv, base != rb)
            total[qv] += 1
        if not bases:
            continue
        mism = [p for p, (_, _, mm) in bases.items() if mm]
        designed = set()
        if is_dnap:  # exclude the read's intended variant codon
            for ci in {p // 3 for p in mism}:
                trip = (3 * ci, 3 * ci + 1, 3 * ci + 2)
                if all(x in bases for x in trip):
                    codon = "".join(bases[x][0] for x in trip)
                    if codon in alts.get(ci, ()):
                        designed.update(trip)
        for p in mism:
            if p not in designed:
                errors[bases[p][1]] += 1
        n += 1
        if n >= max_reads:
            break
    bam.close()
    qvs = sorted(total)
    df = pd.DataFrame({"sample": label, "region": region_lab, "qv": qvs,
                       "total": [total[q] for q in qvs],
                       "errors": [errors.get(q, 0) for q in qvs]})
    df["error_rate"] = df["errors"] / df["total"]
    df["reads_used"] = n
    return df


def summarise(df) -> dict:
    top = df["qv"].max()
    hi, lo = df[df["qv"] == top], df[df["qv"] < top]
    err, tot = df["errors"].sum(), df["total"].sum()
    return {
        "reads": int(df["reads_used"].iloc[0]),
        "overall": err / tot,
        "template_floor": hi["errors"].sum() / hi["total"].sum(),
        "seq_lowqv": (lo["errors"].sum() / lo["total"].sum()) if lo["total"].sum() else 0.0,
        "pct_err_topqv": 100 * hi["errors"].sum() / err if err else 0.0,
    }


def main(track: list[str], ref: Path, variant_refs: Path, outdir: Path, max_reads: int):
    outdir.mkdir(parents=True, exist_ok=True)
    print("=== QV-stratified error profile (sequencing vs template) + backbone control ===")
    _, alts = load_codon_model(variant_refs)
    fasta = pysam.FastaFile(str(ref))

    frames = []
    for spec in track:
        label, region, path = spec.split(":", 2)
        print(f"processing {label} / {region} …")
        frames.append(process_track(label, region, path, fasta, alts, max_reads))
    data = pd.concat(frames, ignore_index=True)
    data.to_csv(outdir / "qv_error_profile.csv", index=False)

    summ = []
    for (label, region), df in data.groupby(["sample", "region"]):
        s = summarise(df); s.update(sample=label, region=region); summ.append(s)
    summ = pd.DataFrame(summ)[["sample", "region", "reads", "overall",
                               "template_floor", "seq_lowqv", "pct_err_topqv"]]
    print("\nError partition by base quality (per track):")
    with pd.option_context("display.float_format", lambda v: f"{v:.2e}"):
        print(summ.to_string(index=False))
    summ.to_csv(outdir / "qv_error_summary.csv", index=False)

    sns.set_theme(style="whitegrid")
    data["track"] = data["sample"] + " " + data["region"]

    # 1) error rate vs QV (log y)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.lineplot(data=data, x="qv", y="error_rate", hue="track", marker="o", ax=ax)
    ax.set_yscale("log")
    ax.set(title="Substitution-error rate vs base quality",
           xlabel="base quality (QV)", ylabel="error rate per base (log)")
    ax.legend(title="", fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "error_rate_vs_qv.png", dpi=150); plt.close(fig)

    # 2) template floor (top-QV) per track, grouped by sample, colored by region
    fig, ax = plt.subplots(figsize=(7.5, 5))
    sns.barplot(data=summ, x="sample", y="template_floor", hue="region", ax=ax)
    for c in ax.containers:
        ax.bar_label(c, fmt="%.1e", fontsize=8)
    ax.set(title="Template error floor (top-QV mismatch rate)\nDNAP−backbone ≈ synthesis; backbone PCR−RE ≈ PCR",
           xlabel="", ylabel="error/base at max QV")
    ax.legend(title="")
    fig.tight_layout(); fig.savefig(outdir / "template_floor_comparison.png", dpi=150); plt.close(fig)

    print(f"done -> {outdir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--track", action="append", required=True,
                   help="LABEL:REGION:BAM  (REGION='dnap' or 'START-END'); repeatable")
    p.add_argument("--ref", type=Path, required=True, help="plasmid FASTA (indexed .fai)")
    p.add_argument("--variant-refs", dest="variant_refs", type=Path, required=True,
                   help="variant_refs.fasta (for designed-codon exclusion)")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--max-reads", dest="max_reads", type=int, default=150_000)
    a = p.parse_args()
    main(a.track, a.ref, a.variant_refs, a.outdir, a.max_reads)
