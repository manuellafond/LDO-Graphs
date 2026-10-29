#!/usr/bin/env python3
"""
Collect the per-replicate evaluations written by simphy_run.sh into OUT/summary.md and OUT/all_families.tsv.

  python simphy_summarize.py --out OUT

Reads OUT/<setting>/<seed>/eval/ldo_<repair>.tsv     (evaluate_orthologs.py per_family.tsv)
                               eval/clusters_<repair>.tsv (cluster_analysis.py --out)
                               eval/of_<mode>.tsv         (evaluate_orthofinder.py --per-family)
      and prepare_report.tsv (genes, species, Spec / Dup nodes of each labelled gene tree).
Precision / recall pool the gene pairs (or clusters) over all families of a setting (micro average).
"""
import argparse
import csv
import glob
import os
import re

from ldo_utils import prf

LDO_ORDER = ["none", "count", "count+sparse", "count+sparse+age"]
OF_ORDER = ["pooled", "perfam", "perfam_fixnorm"]
NAMES = {"ldo_none": "LDO, no repair", "ldo_count": "LDO, count", "ldo_count+sparse": "LDO, count+sparse",
         "ldo_count+sparse+age": "LDO, count+sparse, age test (--count-age 1.0)",
         "of_pooled": "OrthoFinder, pooled per replicate", "of_perfam": "OrthoFinder, per family",
         "of_perfam_fixnorm": "OrthoFinder, per family, patched normalisation"}


# Short name of a simphy setting folder.
#   s   folder name (e.g. simpphy_p10_dl2_t0_s300)
# Returns the "dl<N>" part (dl2), or s itself if there is none.
def setting_label(s):
    m = re.search(r"_(dl\d+)_", s)
    return m.group(1) if m else s


# Read every replicate's evaluation files and write summary.md and all_families.tsv.
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []        # setting, seed, family, method, truth, pred, tp
    clus = []        # setting, seed, family, repair, predicted, correct, ref_scored, found
    info = {}        # (setting, seed, family) -> (genes, species, spec, dup)
    for ev in sorted(glob.glob(os.path.join(a.out, "*", "*", "eval"))):
        rep = os.path.dirname(ev)
        seed = os.path.basename(rep)
        setting = setting_label(os.path.basename(os.path.dirname(rep)))
        pr = os.path.join(rep, "prepare_report.tsv")
        if os.path.exists(pr):
            for r in csv.DictReader(open(pr), delimiter="\t"):
                info[(setting, seed, r["family"])] = (int(r["genes"]), int(r["species"]),
                                                      int(r["spec_nodes"]), int(r["dup_nodes"]))
        for path in sorted(glob.glob(os.path.join(ev, "ldo_*.tsv"))):
            m = os.path.basename(path)[:-4]
            for r in csv.DictReader(open(path), delimiter="\t"):
                rows.append((setting, seed, r["tree_id"], m, int(r["truth_pairs"]), int(r["pred_pairs"]), int(r["tp"])))
        for path in sorted(glob.glob(os.path.join(ev, "of_*.tsv"))):
            m = os.path.basename(path)[:-4]
            for r in csv.DictReader(open(path), delimiter="\t"):
                rows.append((setting, seed, r["family"], m, int(r["truth"]), int(r["pred"]), int(r["tp"])))
        for path in sorted(glob.glob(os.path.join(ev, "clusters_*.tsv"))):
            rp = os.path.basename(path)[len("clusters_"):-4]
            for r in csv.DictReader(open(path), delimiter="\t"):
                clus.append((setting, seed, r["tree_id"], rp, int(r["predicted"]),
                             int(r["exact"]) + int(r["polytomy"]),
                             int(r["reference"]) - int(r["nested_dup"]), int(r["found"])))
    if not rows:
        raise SystemExit(f"no evaluation files under {a.out}")

    with open(os.path.join(a.out, "all_families.tsv"), "w") as f:
        f.write("setting\tseed\tfamily\tgenes\tspecies\tspec_nodes\tdup_nodes\tmethod\ttruth_pairs\tpred_pairs\ttp\n")
        for s, sd, fam, m, t, p, tp in rows:
            g = info.get((s, sd, fam), ("", "", "", ""))
            f.write("\t".join(map(str, (s, sd, fam) + g + (m, t, p, tp))) + "\n")

    methods = [m for m in ["ldo_" + x for x in LDO_ORDER] + ["of_" + x for x in OF_ORDER] if any(r[3] == m for r in rows)]
    settings = sorted({r[0] for r in rows}, key=lambda s: (int(re.sub(r"\D", "", s) or 0), s))

    # Markdown table of pooled precision / recall / F1, one line per method, one column per setting.
    #   keep    function of a row of `rows` -> True if the row is counted
    #   title   section title
    #   note    sentence printed under the title
    # Returns the list of Markdown lines.
    def table(keep, title, note):
        L = [f"## {title}", "", note, ""]
        L.append("| Method | " + " | ".join(settings) + " | All |")
        L.append("|---|" + "---|" * (len(settings) + 1))
        for m in methods:
            cells = []
            for s in settings + [None]:
                t = p = tp = 0
                for r in rows:
                    if r[3] == m and (s is None or r[0] == s) and keep(r):
                        t += r[4]; p += r[5]; tp += r[6]
                cells.append("{:.3f} / {:.3f} / {:.3f}".format(*prf(t, p, tp)) if t or p else "-")
            L.append(f"| {NAMES.get(m, m)} | " + " | ".join(cells) + " |")
        return L

    L = ["# simphy: LDO reconstruction vs OrthoFinder", "",
         "Truth: the simphy gene trees labelled with simphy's own events (simphy_prepare.py --labels simphy). "
         "Ortholog pairs, precision / recall / F1, gene pairs pooled over the families of each setting.", ""]
    n_rep = len({(r[0], r[1]) for r in rows})
    L.append(f"{n_rep} replicates, {len({(r[0], r[1], r[2]) for r in rows})} families.")
    L.append("")
    L.append("| Setting | replicates | families | genes / family | Dup nodes / family |")
    L.append("|---|---|---|---|---|")
    for s in settings:
        fam = [v for k, v in info.items() if k[0] == s]
        reps = len({k[1] for k in info if k[0] == s})
        if fam:
            L.append(f"| {s} | {reps} | {len(fam)} | {sum(x[0] for x in fam) / len(fam):.1f} | "
                     f"{sum(x[3] for x in fam) / len(fam):.2f} |")
    L.append("")
    L += table(lambda r: True, "Ortholog pairs, all families",
               "Each method on all the families it was run on (per-family OrthoFinder may cover a subset only).")
    # Markdown table of per-family average precision / recall, one line per method, one column per setting.
    #   keep    function of a row of `rows` -> True if the row is counted
    #   title   section title
    #   note    sentence printed under the title
    # Returns the list of Markdown lines.
    def macro_table(keep, title, note):
        L = [f"## {title}", "", note, ""]
        L.append("| Method | families | " + " | ".join(settings) + " | All |")
        L.append("|---|---|" + "---|" * (len(settings) + 1))
        for m in methods:
            cells, n = [], 0
            for s in settings + [None]:
                v = [(r[6] / r[5] if r[5] else 1.0, r[6] / r[4]) for r in rows
                     if r[3] == m and (s is None or r[0] == s) and r[4] > 0 and keep(r)]
                if s is None:
                    n = len(v)
                cells.append(f"{sum(x[0] for x in v) / len(v):.3f} / {sum(x[1] for x in v) / len(v):.3f}" if v else "-")
            L.append(f"| {NAMES.get(m, m)} | {n} | " + " | ".join(cells) + " |")
        return L

    pf = {(r[0], r[1], r[2]) for r in rows if r[3].startswith("of_perfam")}
    allf = {(r[0], r[1], r[2]) for r in rows}
    if pf and pf != allf:
        L.append("")
        L += table(lambda r: (r[0], r[1], r[2]) in pf, "Ortholog pairs, families with per-family OrthoFinder runs",
                   f"All methods restricted to the {len(pf)} families that have per-family OrthoFinder results.")
    L.append("")
    L += macro_table(lambda r: True, "Average per family (precision / recall)",
                     "Each family scored on its own, then averaged; families without any ortholog pair (all genes in "
                     "one species) are left out for every method. A family with no predicted pair counts as precision 1.")
    if pf and pf != allf:
        L.append("")
        L += macro_table(lambda r: (r[0], r[1], r[2]) in pf, "Average per family, families with per-family OrthoFinder runs",
                         "As above, restricted to the families with per-family OrthoFinder results.")
    if clus:
        L += ["", "## Clusters (LDO trees)", "",
              "Cluster precision (exact + polytomy) / recall (reference clusters other than nested duplications), "
              "as cluster_analysis.py, pooled.", "",
              "| Repair | " + " | ".join(settings) + " | All |", "|---|" + "---|" * (len(settings) + 1)]
        for rp in [x for x in LDO_ORDER if any(c[3] == x for c in clus)]:
            cells = []
            for s in settings + [None]:
                sel = [c for c in clus if c[3] == rp and (s is None or c[0] == s)]
                P = sum(c[4] for c in sel); C = sum(c[5] for c in sel)
                R = sum(c[6] for c in sel); F = sum(c[7] for c in sel)
                cells.append(f"{C / P if P else 1:.3f} / {F / R if R else 1:.3f}" if sel else "-")
            L.append(f"| {rp} | " + " | ".join(cells) + " |")
    out = os.path.join(a.out, "summary.md")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
