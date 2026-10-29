#!/usr/bin/env python3
"""
Collect the evaluation files written by swisstree/reproduce.sh into one Markdown report.

  python summarize.py --eval out/eval --out out/summary.md

Reads   eval/ldo_<data>_<repair>_<truth>.tsv   per_family.tsv of evaluate_orthologs.py
        eval/of_<variant>_<truth>.tsv          --per-family output of evaluate_orthofinder.py
        eval/clusters_<data>_<repair>.txt      cluster_analysis.py output
<data> is 'base' (conversion st) or 'frag' (conversion stf: gene-name matching, fragments removed).
<truth> is 'curated' (pairs of the curated trees; ST006, whose tree has no event tags, is left out) or 'st006' (ST006
only, against the curators' pair list gene_relationships.txt).
Totals pool the gene pairs of the families (micro average), as evaluate_orthologs.py does.
"""
import argparse
import csv
import os
import re

from ldo_utils import prf

EXCLUDE = {"curated": {"ST006"}, "st006": set()}

# (row label, LDO key without the data prefix or OrthoFinder key), in report order
LDO_RUNS = [("LDO, no repair", "none"), ("LDO, count", "count"), ("LDO, count+sparse", "count+sparse"),
            ("LDO, count+sparse + age test (--count-age 1.0)", "count+sparse+age")]
OF_RUNS = [("OrthoFinder, pooled, true species tree", "of_pooled"),
           ("OrthoFinder, pooled, patched normalisation", "of_pooled_fixnorm"),
           ("OrthoFinder, per family, true species tree", "of_perfam"),
           ("OrthoFinder, per family, patched normalisation", "of_perfam_fixnorm")]


# Per-family counts of an LDO evaluation.
#   path   per_family.tsv of evaluate_orthologs.py
# Returns {family: (truth pairs, predicted pairs, correct pairs)}.
def read_ldo(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["tree_id"]] = (int(r["truth_pairs"]), int(r["pred_pairs"]), int(r["tp"]))
    return out


# Per-family counts of an OrthoFinder evaluation.
#   path   --per-family file of evaluate_orthofinder.py
# Returns {family: (truth pairs, predicted pairs, correct pairs)}.
def read_of(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["family"]] = (int(r["truth"]), int(r["pred"]), int(r["tp"]))
    return out


# Per-family counts of one evaluation file, if it exists.
#   E       eval folder
#   key     file key without the truth suffix (ldo_base_count, of_pooled, ...)
#   truth   "curated" or "st006"
# Returns {family: (truth, predicted, correct)}, or None if the file is missing.
def read_eval(E, key, truth):
    path = os.path.join(E, f"{key}_{truth}.tsv")
    if not os.path.exists(path):
        return None
    return read_ldo(path) if key.startswith("ldo") else read_of(path)


# Micro-averaged precision, recall, F1 of per-family counts.
#   d       {family: (truth pairs, predicted pairs, correct pairs)} (read_ldo() / read_of())
#   truth   "curated" or "st006": the families in EXCLUDE[truth] are left out
# Returns (precision, recall, F1).
def total(d, truth):
    t = p = tp = 0
    for f, (a, b, c) in d.items():
        if f in EXCLUDE.get(truth, set()):
            continue
        t += a; p += b; tp += c
    return prf(t, p, tp)


# "P / R / F1" of per-family counts, or "-".
#   d       per-family counts, or None
#   truth   "curated" or "st006"
# Returns the table cell.
def cell(d, truth):
    return " / ".join(f"{x:.3f}" for x in total(d, truth)) if d else "-"


# Cluster precision and recall of a cluster_analysis.py report.
#   path   text file with the printed output of cluster_analysis.py
# Returns the string "precision / recall", or "-" if the file is missing or not recognised.
def clusters(path):
    if not os.path.exists(path):
        return "-"
    txt = open(path).read()
    p = re.search(r"cluster precision \(exact \+ polytomy\)\s*:\s*([0-9.]+)", txt)
    r = re.search(r"cluster recall:\s*([0-9.]+)", txt)
    return f"{float(p.group(1)):.3f} / {float(r.group(1)):.3f}" if p and r else "-"


# Read all the evaluation files and write the Markdown report.
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    E = a.eval

    # Main table: curated truth, 8 families (conversion st; stf gives the same numbers, it differs only in ST006).
    L = ["# Results", "",
         "Ortholog pairs, precision / recall / F1, gene pairs pooled over the families (micro average).",
         "Truth: the curated SwissTree trees, 8 families (ST006 is scored separately below: its tree file has no "
         "event tags).",
         "Clusters: cluster precision / recall of the reconstructed trees against the curated trees, all 9 families.",
         "", "| Run | Curated, 8 families | Clusters P / R |", "|---|---|---|"]
    for label, r in LDO_RUNS:
        d = read_eval(E, "ldo_base_" + r, "curated")
        if d:
            L.append(f"| {label} | {cell(d, 'curated')} | {clusters(os.path.join(E, f'clusters_base_{r}.txt'))} |")
    for label, k in OF_RUNS:
        d = read_eval(E, k, "curated")
        if d:
            L.append(f"| {label} | {cell(d, 'curated')} | - |")

    # ST006 table: every method against the curators' pair list, on both conversions (OrthoFinder: st only).
    rows = []
    for label, r in LDO_RUNS:
        b, f = read_eval(E, "ldo_base_" + r, "st006"), read_eval(E, "ldo_frag_" + r, "st006")
        if b or f:
            rows.append(f"| {label} | {cell(b, 'st006')} | {cell(f, 'st006')} |")
    for label, k in OF_RUNS:
        d = read_eval(E, k, "st006")
        if d:
            rows.append(f"| {label} | {cell(d, 'st006')} | - |")
    if rows:
        L += ["", "## ST006, curators' pair list", "",
              "ST006 against the ortholog pairs of the curators' gene_relationships.txt. st: plain conversion; stf: "
              "leaves also matched by gene name and fragments removed (OrthoFinder was run on st only).", "",
              "| Run | st | stf |", "|---|---|---|"] + rows

    # Per-tree table, curated truth: precision / recall of every family.
    base = [("LDO none", "ldo_base_none"), ("LDO count", "ldo_base_count"),
            ("LDO count+sparse", "ldo_base_count+sparse"), ("LDO count+sparse+age", "ldo_base_count+sparse+age"),
            ("OF unpatched", "of_pooled"), ("OF patched", "of_pooled_fixnorm")]
    cols = [(c, read_eval(E, k, "curated")) for c, k in base]
    cols = [(c, d) for c, d in cols if d]
    if cols:
        fams = sorted(set().union(*[set(d) for _, d in cols]) - EXCLUDE["curated"])
        L += ["", "## Per tree, curated truth", "", "Precision / recall.", "",
              "| Tree | " + " | ".join(c for c, _ in cols) + " |", "|---|" + "---|" * len(cols)]
        for f in fams:
            row = []
            for _, d in cols:
                if f in d:
                    P, R, _ = prf(*d[f])
                    row.append(f"{P:.3f} / {R:.3f}")
                else:
                    row.append("-")
            L.append(f"| {f} | " + " | ".join(row) + " |")
        L.append("| **All** | " + " | ".join(
            "{:.3f} / {:.3f}".format(*total(d, "curated")[:2]) for _, d in cols) + " |")
    open(a.out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
