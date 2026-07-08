"""Shared helpers for the DMS analysis scripts (count / plot / error profiling).

Holds the things that would otherwise be copied across scripts: FASTA IO, the
DNAP-window geometry, the variant-name regex, and the per-variant reference model
plus the exact nearest-reference distance used to classify reads.
"""
from __future__ import annotations

import re
from collections import defaultdict

# DNAP window on the plasmid: 1-based inclusive, plus derived 0-based constants.
CONTIG = "pFR494_pRT300_rham_wt"
WINDOW = (1790, 3448)                    # == 553 codons
WIN0 = WINDOW[0] - 1                      # 0-based start
WIN_LEN = WINDOW[1] - WINDOW[0] + 1       # 1659

# Variant name, e.g. DNAP_seg01_DNAP_002_P_to_A -> (seg_num, codon, wt_aa, mut_aa)
VARIANT_RE = re.compile(r"^DNAP_seg(\d+)_DNAP_(\d+)_([A-Z*])_to_([A-Z*])$")

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def iter_fasta(path):
    """Yield (name, sequence) per record, joining wrapped lines."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def read_fasta(path) -> dict:
    """Whole FASTA as {name: sequence}."""
    return dict(iter_fasta(path))


def load_codon_model(variant_refs) -> tuple[str, dict]:
    """From variant_refs.fasta return (wt_window, {codon_index: {designed alt 3-mers}}).

    Codons tile the window every 3 bp, so codon i occupies window[3i:3i+3].
    """
    seqs = read_fasta(variant_refs)
    wt = seqs["WT"]
    alts: dict[int, set] = defaultdict(set)
    for name, s in seqs.items():
        if name != "WT" and len(s) == WIN_LEN:
            for ci in {j // 3 for j in range(WIN_LEN) if s[j] != wt[j]}:
                alts[ci].add(s[3 * ci:3 * ci + 3])
    return wt, alts


def _hd3(a: str, b: str) -> int:
    return (a[0] != b[0]) + (a[1] != b[1]) + (a[2] != b[2])


def nearest_ref_dist(seq: str, wt: str, alts: dict) -> int:
    """Fewest mismatches from a full-window `seq` to any reference (WT or a variant).

    Shortcut: every reference is WT with at most one codon changed, so
        dist(seq, variant) = dist(seq, WT) - hd(seq_c, WT_c) + hd(seq_c, alt_c),
    and only codons where seq differs from WT can beat the plain WT distance.
    """
    diff_codons = {j // 3 for j in range(WIN_LEN) if seq[j] != wt[j]}
    d_wt = sum(_hd3(seq[3 * c:3 * c + 3], wt[3 * c:3 * c + 3]) for c in diff_codons)
    best = d_wt
    for c in diff_codons & alts.keys():
        sc = seq[3 * c:3 * c + 3]
        best = min(best, d_wt - _hd3(sc, wt[3 * c:3 * c + 3]) + min(_hd3(sc, a) for a in alts[c]))
    return best
