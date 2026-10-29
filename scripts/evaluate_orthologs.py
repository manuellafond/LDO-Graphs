#!/usr/bin/env python3
"""
Score the ortholog pairs of reconstructed gene trees against true ortholog pairs.

  python evaluate_orthologs.py --trees trees.txt --truth ldo/orthologs.tsv --out OUT_DIR

trees.txt : "tree_id<TAB>newick" lines (ldo_builder.py), internal nodes labelled Spec / Dup; two genes are predicted
            orthologs iff their last common ancestor is a Spec node.
truth     : "tree_id<TAB>gene1<TAB>gene2" (true_orthologs.py), restricted to the genes of the reconstructed tree.
Writes OUT_DIR/per_family.tsv and prints the micro-averaged precision / recall / F1.
"""
import argparse
import os
import sys

from ldo_utils import parse_newick, prf


# Leaves, ortholog pairs and numbers of Spec / Dup nodes of a reconstructed tree.
#   text   Newick string, internal nodes labelled Spec / Dup
# Returns (list of leaves, set of ortholog pairs (g, h) with g < h, number of Spec nodes, number of Dup nodes).
def orthologs_from_newick(text):
    pairs, counts = set(), {"Spec": 0, "Dup": 0}

    # Post-order walk: adds the pairs of the Spec nodes to `pairs` and counts the nodes.
    #   x   current node (dict tree)
    # Returns the leaves below x.
    def rec(x):
        if not x["children"]:
            return [x["name"]]
        groups = [rec(c) for c in x["children"]]
        if x["name"] in counts:
            counts[x["name"]] += 1
        if x["name"] == "Spec":
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    for g in groups[i]:
                        for h in groups[j]:
                            pairs.add((g, h) if g < h else (h, g))
        return [g for grp in groups for g in grp]
    leaves = rec(parse_newick(text))
    return leaves, pairs, counts["Spec"], counts["Dup"]


# Score every family and write per_family.tsv.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fam = {}
    for line in open(a.trees):
        line = line.strip()
        if line:
            tid, nwk = line.split("\t", 1)
            leaves, pairs, spec, dup = orthologs_from_newick(nwk)
            fam[tid] = dict(leaves=set(leaves), pred=pairs, spec=spec, dup=dup, truth=set())
    if not fam:
        sys.exit("No trees read from " + a.trees)
    for line in open(a.truth):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3 or p[0] not in fam:
            continue
        d = fam[p[0]]
        if p[1] in d["leaves"] and p[2] in d["leaves"]:
            d["truth"].add((p[1], p[2]) if p[1] < p[2] else (p[2], p[1]))

    tot = [0, 0, 0]
    with open(os.path.join(a.out, "per_family.tsv"), "w") as f:
        f.write("tree_id\tn_genes\tspec_nodes\tdup_nodes\ttruth_pairs\tpred_pairs\ttp\tprecision\trecall\tf1\texact\n")
        for tid in sorted(fam):
            d = fam[tid]
            tp = len(d["pred"] & d["truth"])
            p, r, f1 = prf(len(d["truth"]), len(d["pred"]), tp)
            row = (tid, len(d["leaves"]), d["spec"], d["dup"], len(d["truth"]), len(d["pred"]), tp, p, r, f1,
                   int(d["pred"] == d["truth"]))
            f.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in row) + "\n")
            tot[0] += len(d["truth"]); tot[1] += len(d["pred"]); tot[2] += tp
    p, r, f1 = prf(*tot)
    print(f"{len(fam)} families. Pairs: truth={tot[0]} predicted={tot[1]} correct={tot[2]}")
    print(f"Micro precision={p:.4f} recall={r:.4f} F1={f1:.4f}")


if __name__ == "__main__":
    main()
