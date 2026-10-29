#!/usr/bin/env python3
"""
Collect the evaluations written by ensembl_run.sh into OUT/summary.md.

  python ensembl_summarize.py --out OUT

Reads   OUT/prepare_report.tsv               families kept, genes, species, events (ensembl_prepare.py)
        OUT/eval/ldo_<repair>_<dub>.tsv       per_family.tsv of evaluate_orthologs.py
        OUT/eval/clusters_<repair>_<dub>.tsv  cluster_analysis.py --out
        OUT/eval/of_<variant>_<dub>.tsv       evaluate_orthofinder.py --per-family (variant: pooled, pooled_fixnorm)
<dub> is the reading of the Ensembl dubious nodes in the truth: dup or spec.
Ortholog pairs: precision / recall / F1 with the gene pairs pooled over the families (micro average), and the
per-family average precision / recall. Three truths: dubious nodes as Dup, as Spec, and the families without dubious
nodes only (where both readings give the same truth).
"""
import argparse
import csv
import os

from ldo_utils import prf

ORDER = ["none", "count", "count+sparse", "count+sparse+age", "count+sparse+mutual", "count+sparse+mutual+age"]
NAMES = {"none": "LDO, no repair", "count": "LDO, count", "count+sparse": "LDO, count+sparse",
         "count+sparse+age": "LDO, count+sparse, age test (--count-age 1.0)",
         "count+sparse+mutual": "LDO, count+sparse+mutual", "count+sparse+mutual+age": "LDO, count+sparse+mutual, age test",
         "of_pooled": "OrthoFinder, pooled, true species tree", "of_pooled_fixnorm": "OrthoFinder, pooled, patched "
         "normalisation"}
OF_ORDER = ["of_pooled", "of_pooled_fixnorm"]
SIZE_BINS = [(2, 10), (11, 25), (26, 50), (51, 100), (101, 200), (201, None)]   # genes per family (None = no limit)
BIN_METHODS = ["count+sparse", "count+sparse+age", "count+sparse+mutual+age", "of_pooled"]          # methods shown by family size


# Rows of a TSV file with a header line.
#   path   file
# Returns a list of dicts (column name -> value).
def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


# Per-family pair counts of an evaluation file, in the column names of evaluate_orthologs.py.
#   path   per_family.tsv of evaluate_orthologs.py, or --per-family file of evaluate_orthofinder.py (columns family,
#          truth, pred, tp)
# Returns a list of dicts with tree_id, truth_pairs, pred_pairs, tp.
def read_counts(path):
    rows = read_tsv(path)
    if rows and "family" in rows[0]:
        rows = [dict(tree_id=r["family"], truth_pairs=r["truth"], pred_pairs=r["pred"], tp=r["tp"]) for r in rows]
    return rows


# Pooled precision / recall / F1 and per-family average precision / recall of some families.
#   rows   rows of a per_family.tsv of evaluate_orthologs.py
#   keep   set of the families to count
# Returns ("P / R / F1", "P / R"); the average leaves out families without true pairs (a family without predicted
# pair counts as precision 1).
def scores(rows, keep):
    sel = [r for r in rows if r["tree_id"] in keep]
    t = sum(int(r["truth_pairs"]) for r in sel)
    p = sum(int(r["pred_pairs"]) for r in sel)
    tp = sum(int(r["tp"]) for r in sel)
    pooled = "{:.3f} / {:.3f} / {:.3f}".format(*prf(t, p, tp))
    avg = [prf(int(r["truth_pairs"]), int(r["pred_pairs"]), int(r["tp"]))[:2] for r in sel if int(r["truth_pairs"])]
    mean = (f"{sum(x[0] for x in avg) / len(avg):.3f} / {sum(x[1] for x in avg) / len(avg):.3f}" if avg else "-")
    return pooled, mean


# Pooled cluster precision (exact + polytomy) / recall (reference clusters other than nested duplications).
#   path   cluster_analysis.py --out file
#   keep   set of the families to count
# Returns "P / R", or "-" if the file is missing.
def cluster_scores(path, keep):
    if not os.path.exists(path):
        return "-"
    sel = [r for r in read_tsv(path) if r["tree_id"] in keep]
    P = sum(int(r["predicted"]) for r in sel)
    C = sum(int(r["exact"]) + int(r["polytomy"]) for r in sel)
    R = sum(int(r["reference"]) - int(r["nested_dup"]) for r in sel)
    F = sum(int(r["found"]) for r in sel)
    return f"{C / P if P else 1:.3f} / {F / R if R else 1:.3f}"


# Markdown table of the scores by family size: one line per size bin, one column per method; the gene pairs are
# pooled within each bin.
#   data   {(method, dub): rows of read_counts()}
#   rep    kept rows of prepare_report.tsv (tree_id, genes)
#   dub    reading of the dubious nodes, "dup" or "spec"
#   runs   methods with evaluation files
# Returns the list of Markdown lines.
def size_table(data, rep, dub, runs):
    methods = [m for m in BIN_METHODS if m in runs]
    L = ["| Genes per family | families | true pairs | " + " | ".join(NAMES[m] for m in methods) + " |",
         "|---|---|---|" + "---|" * len(methods)]
    for lo, hi in SIZE_BINS:
        keep = {r["tree_id"] for r in rep if int(r["genes"]) >= lo and (hi is None or int(r["genes"]) <= hi)}
        if not keep:
            continue
        n_true = sum(int(r["truth_pairs"]) for r in data[(methods[0], dub)] if r["tree_id"] in keep)
        label = f"{lo}-{hi}" if hi else f"> {lo - 1}"
        L.append(f"| {label} | {len(keep)} | {n_true} | " +
                 " | ".join(scores(data[(m, dub)], keep)[0] for m in methods) + " |")
    return L


# Write OUT/summary.md.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    E = os.path.join(a.out, "eval")

    rep = [r for r in read_tsv(os.path.join(a.out, "prepare_report.tsv")) if r["status"] == "kept"]
    fams = {r["tree_id"] for r in rep}
    clean = {r["tree_id"] for r in rep if int(r["dubious"]) == 0}
    S = {k: sum(int(r[k]) for r in rep) for k in ("genes", "removed_leaves", "spec", "dup", "dubious", "other",
                                                  "spec_conflict", "dubious_spec_ok")}
    truths = [("dup", fams, f"dubious = Dup, {len(fams)} families"),
              ("spec", fams, f"dubious = Spec, {len(fams)} families"),
              ("dup", clean, f"no dubious node, {len(clean)} families")]

    L = ["# Ensembl Compara: LDO reconstruction", "",
         f"{len(fams)} families, {S['genes']} genes ({S['removed_leaves']} leaves removed: no aligned sequence or "
         f"species not in the species tree). Internal nodes of the Ensembl trees: {S['spec']} speciation "
         f"({S['spec_conflict']} that the species tree cannot explain), {S['dup']} duplication, {S['dubious']} "
         f"dubious ({S['dubious_spec_ok']} of them would be valid speciations), {S['other']} other. "
         f"{len(clean)} families have no dubious node.", "",
         "Truth: the Ensembl trees. Dubious nodes read as duplications or as speciations; on the families without "
         "dubious nodes both readings agree.", "",
         "## Ortholog pairs, pooled (precision / recall / F1)", "",
         "| Method | " + " | ".join(t[2] for t in truths) + " |", "|---|" + "---|" * len(truths)]
    # methods with evaluation files: LDO repair settings (key = setting) and OrthoFinder variants (key = of_<variant>)
    files = {r: f"ldo_{r}" for r in ORDER} | {k: k for k in OF_ORDER}
    runs = [r for r in ORDER + OF_ORDER if os.path.exists(os.path.join(E, f"{files[r]}_dup.tsv"))]
    data = {(r, d): read_counts(os.path.join(E, f"{files[r]}_{d}.tsv")) for r in runs for d in ("dup", "spec")}
    for r in runs:
        L.append(f"| {NAMES[r]} | " + " | ".join(scores(data[(r, d)], k)[0] for d, k, _ in truths) + " |")
    L += ["", "## Ortholog pairs, average per family (precision / recall)", "",
          "Families without any true ortholog pair are left out; a family without predicted pair counts as "
          "precision 1.", "",
          "| Method | " + " | ".join(t[2] for t in truths) + " |", "|---|" + "---|" * len(truths)]
    for r in runs:
        L.append(f"| {NAMES[r]} | " + " | ".join(scores(data[(r, d)], k)[1] for d, k, _ in truths) + " |")
    for dub, title in (("dup", "dubious = Dup"), ("spec", "dubious = Spec")):
        L += ["", f"## Ortholog pairs by family size, {title} (precision / recall / F1)", "",
              "Gene pairs pooled within each size bin (genes per family), so that no single family dominates.", ""]
        L += size_table(data, rep, dub, runs)
    L += ["", "## Clusters (precision / recall)", "",
          "As cluster_analysis.py, pooled; nested duplications depend on the reading of the dubious nodes.", "",
          "| Repair | " + " | ".join(t[2] for t in truths) + " |", "|---|" + "---|" * len(truths)]
    for r in [r for r in runs if r in ORDER]:
        L.append(f"| {r} | " + " | ".join(cluster_scores(os.path.join(E, f"clusters_{r}_{d}.tsv"), k)
                                           for d, k, _ in truths) + " |")
    out = os.path.join(a.out, "summary.md")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
