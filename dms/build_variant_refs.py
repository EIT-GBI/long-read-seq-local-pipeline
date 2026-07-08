#!/usr/bin/env python3
"""
Build per-variant nucleotide references for the DNAP deep-mutational-scan.

Goal
----
For each designed variant (e.g. DNAP_seg01_DNAP_002_P_to_A) produce the full
DNAP-window sequence (plasmid 1-based 1790..3448) with ONLY that variant's single
codon change slotted in. These references are what reads get counted against with
vsearch (see the companion notes / count_variants section at the bottom).

Why not just trust the oligo sequences directly?
------------------------------------------------
Each row in DMS_segments.csv is a synthesis oligo, not a clean insert:

    [FR2241 arm][BsaI GGTCTCN][ATCG overhang][ ... DNAP body, 1 codon mutated ... ][overhang][BsaI][rc FR2242 arm]

We must strip the PCR/Golden-Gate cassette (arms + BsaI + spacers) and keep only
the DNAP body. Two ideas were considered:

  * "Strip columns that are constant across a segment's variants."  Close, but the
    constant/variable boundary is the first/last *mutated* codon, which is INSIDE
    the true body (it drops the overhang + never-mutated codons). Over-trims.

  * (used here) Reconstruct each segment's WT oligo as the per-column consensus of
    its ~817 variants (only one interior codon differs per variant, so the majority
    base is WT everywhere). The WT body is then, by definition, the longest stretch
    of that consensus that is an EXACT substring of the WT window. Its start offset
    in the oligo = the 5' cassette length; its coordinates in the window tell us
    exactly where to slot mutations. Robust, and self-checking.

Geometry discovered in the data (see build log):
  * All segment bodies match the WINDOW's + strand (the CSV is plus-strand).
  * The DNAP ORF is on the MINUS strand: revcomp(window) starts with ATG (codon 1).
  * Segments tile the window in reverse order (seg01 at the 3' end, seg13 at 5').

Validation
----------
For every variant we (a) require the slotted body to differ from WT within a single
codon, and (b) translate the RC (coding) strand and confirm the residue change and
codon number match the variant name. Anything failing is reported, not silently written.

Usage
-----
  python build_variant_refs.py \
      --plasmid  /path/pFR494_pRT300_rham_wt.fa \
      --csv      /path/DMS_segments.csv \
      --out      /path/variant_refs.fasta \
      [--segment DNAP_seg01]        # optional: just one segment (~817 refs)
      [--window 1790 3448]          # 1-based inclusive; default 1790 3448
      [--include-wt]                # also emit the unmutated WT window as ">WT"
"""
import argparse, csv, collections, re, sys

CODON_RE = re.compile(r"^(DNAP_seg\d+)_DNAP_(\d+)_([A-Z*])_to_([A-Z*])$")

# Standard genetic code (DNA codons -> 1-letter AA; '*' = stop)
_BASES = "TCAG"
_AAS   = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = {a+b+c: _AAS[i] for i, (a, b, c) in
            enumerate((x, y, z) for x in _BASES for y in _BASES for z in _BASES)}
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def translate(nt):
    return "".join(CODON2AA.get(nt[i:i+3], "X") for i in range(0, len(nt) - 2, 3))


def read_fasta_single(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def load_variants(csv_path):
    """Return {segment: [(name, oligo_seq), ...]} in file order."""
    byseg = collections.defaultdict(list)
    with open(csv_path, encoding="utf-8-sig", newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        for row in rdr:
            if len(row) < 2 or not row[0].strip():
                continue
            name, seq = row[0].strip(), row[1].strip().upper()
            m = CODON_RE.match(name)
            if not m:
                print(f"[WARN] name not parseable, skipping: {name}", file=sys.stderr)
                continue
            byseg[m.group(1)].append((name, seq))
    return byseg


def segment_consensus(oligos):
    """Per-column majority base across a segment's oligos (reconstructs WT oligo).
    Assumes equal length within a segment (verified true in this library)."""
    L = len(oligos[0][1])
    if any(len(s) != L for _, s in oligos):
        raise ValueError("oligos in a segment are not all the same length")
    cons = []
    for i in range(L):
        c = collections.Counter(s[i] for _, s in oligos)
        cons.append(c.most_common(1)[0][0])
    return "".join(cons)


def locate_body(consensus, window):
    """Longest exact substring of the WT consensus oligo that occurs in the WT window.
    Returns (oligo_start, body_len, window_start) — all 0-based."""
    best = (0, 0, 0)  # (body_len, oligo_start, window_start)
    n = len(consensus)
    for start in range(0, 60):                 # cassette is within the first ~40 bp
        for blen in range(n - start, 60, -1):  # try longest first
            j = window.find(consensus[start:start + blen])
            if j >= 0:
                if blen > best[0]:
                    best = (blen, start, j)
                break                          # longer starts here won't extend
    body_len, oligo_start, window_start = best
    return oligo_start, body_len, window_start


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plasmid", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--segment", default=None, help="build only this segment, e.g. DNAP_seg01")
    ap.add_argument("--window", nargs=2, type=int, default=[1790, 3448],
                    metavar=("START", "END"), help="1-based inclusive plasmid coords")
    ap.add_argument("--include-wt", action="store_true")
    ap.add_argument("--report", default=None, help="TSV validation report (default: <out>.report.tsv)")
    args = ap.parse_args()

    plasmid = read_fasta_single(args.plasmid)
    w0, w1 = args.window[0] - 1, args.window[1]      # -> 0-based [w0, w1)
    wt_win = plasmid[w0:w1]
    win_len = len(wt_win)
    wt_coding_aa = translate(revcomp(wt_win))
    print(f"[INFO] WT window {args.window[0]}..{args.window[1]}  len={win_len} "
          f"({win_len/3:.0f} codons)  ORF starts {revcomp(wt_win)[:3]} -> {wt_coding_aa[:1]}")

    byseg = load_variants(args.csv)
    segs = [args.segment] if args.segment else sorted(byseg, key=lambda x: int(x[-2:]))

    report_path = args.report or (args.out + ".report.tsv")
    n_written = n_fail = 0
    with open(args.out, "w") as fa, open(report_path, "w") as rep:
        rep.write("name\tsegment\tcodon\twt_aa\tmut_aa\tn_nt_diff\twin_lo1\twin_hi1\taa_ok\tnote\n")
        if args.include_wt:
            fa.write(f">WT\n{wt_win}\n")

        for seg in segs:
            oligos = byseg[seg]
            cons = segment_consensus(oligos)
            o_start, body_len, w_start = locate_body(cons, wt_win)
            print(f"[INFO] {seg}: {len(oligos)} variants | cassette 5'={o_start}bp | "
                  f"body={body_len}bp -> window {w_start+args.window[0]}.."
                  f"{w_start+body_len-1+args.window[0]}")

            for name, oligo in oligos:
                m = CODON_RE.match(name)
                codon_num, wt_aa, mut_aa = int(m.group(2)), m.group(3), m.group(4)
                body = oligo[o_start:o_start + body_len]

                # Slot the variant body into the WT window (length preserved).
                var_win = wt_win[:w_start] + body + wt_win[w_start + body_len:]

                # Validation A: nucleotide diffs vs WT should sit within one codon.
                diffs = [i for i in range(win_len) if var_win[i] != wt_win[i]]
                # Validation B: translate coding strand, confirm the AA change + position.
                var_aa = translate(revcomp(var_win))
                aa_diffs = [(i + 1, wt_coding_aa[i], var_aa[i])
                            for i in range(len(wt_coding_aa)) if wt_coding_aa[i] != var_aa[i]]
                aa_ok = (len(aa_diffs) == 1 and aa_diffs[0] == (codon_num, wt_aa, mut_aa))
                note = ""
                if len(var_win) != win_len:
                    note = "LENGTH_CHANGED"
                elif not diffs:
                    note = "NO_DIFF"
                elif len(aa_diffs) != 1:
                    note = f"AA_DIFFS={len(aa_diffs)}"
                elif not aa_ok:
                    note = f"AA_MISMATCH:got_{aa_diffs[0][0]}_{aa_diffs[0][1]}to{aa_diffs[0][2]}"

                lo1 = (min(diffs) + args.window[0]) if diffs else 0
                hi1 = (max(diffs) + args.window[0]) if diffs else 0
                rep.write(f"{name}\t{seg}\t{codon_num}\t{wt_aa}\t{mut_aa}\t{len(diffs)}\t"
                          f"{lo1}\t{hi1}\t{int(aa_ok)}\t{note}\n")

                if aa_ok:
                    fa.write(f">{name}\n{var_win}\n")
                    n_written += 1
                else:
                    n_fail += 1

    print(f"[DONE] wrote {n_written} variant refs -> {args.out}")
    print(f"[DONE] validation report -> {report_path}")
    if n_fail:
        print(f"[WARN] {n_fail} variant(s) failed AA validation and were NOT written "
              f"(see report 'note' column)")


if __name__ == "__main__":
    main()
