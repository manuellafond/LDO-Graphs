#!/usr/bin/env python3
"""
Compare the clusters (gene sets of internal nodes) of reconstructed trees with those of the reference trees.

  python cluster_analysis.py --pred trees.txt --true-trees TREES_DIR [--out clusters.tsv] [--dubious dup|spec]

--dubious: label of the Ensembl "dubious" nodes of the reference trees (default dup); it matters for the nested
duplications. Use the same value as for true_orthologs.py.

A predicted cluster (internal node; root and leaves excluded) is
  exact     a cluster of the reference tree,
  polytomy  a union of >= 2 (not all) children of one multifurcating reference node (a binary resolution: no error),
  wrong     anything else.
Reference clusters missing from the prediction are 'nested dup' (a duplication whose parent is a duplication, which
Algorithm 3 merges into one node by design) or 'missed'.
Cluster precision = (exact + polytomy) / predicted; recall = found / (reference clusters other than nested dup).
"""
import argparse
import collections
import os

from ldo_utils import event_label, leaf_gene, load_json_tree, parse_newick


# Clusters of a reconstructed tree.
#   root   root node of the reconstructed tree (dict tree, internal nodes named Spec / Dup)
# Returns a list of (frozenset of leaves, Spec / Dup label), one per internal node, root last.
def pred_clusters(root):
    out = []

    # Post-order walk; appends the cluster of every internal node to `out`.
    #   x   current node
    # Returns the leaf set of x.
    def rec(x):
        if not x["children"]:
            return frozenset([x["name"]])
        s = frozenset().union(*[rec(c) for c in x["children"]])
        out.append((s, x["name"]))
        return s
    rec(root)
    return out


# Internal nodes of a reference tree.
#   path      JSON gene tree file
#   dubious   label of the dubious nodes, "dup" or "spec" (see ldo_utils.event_label())
# Returns a list of dict(label = Spec / Dup, parent = index of the parent node (-1 for the root), set = frozenset of
# the leaves below, ch = list of the leaf sets of the children), root first.
def true_nodes(path, dubious="dup"):
    _, root = load_json_tree(path)
    nodes = []

    # Pre-order walk; appends every internal node to `nodes`.
    #   x        current node (JSON dict)
    #   parent   index of the parent node in `nodes` (-1 for the root)
    # Returns (leaf set of x, node index of x or -1 for a leaf).
    def rec(x, parent):
        if not x.get("children"):
            return frozenset([leaf_gene(x)]), -1
        idx = len(nodes)
        nodes.append(dict(label=event_label(x, dubious), parent=parent, ch=[], set=None))
        sets = []
        for c in x["children"]:
            s, _ = rec(c, idx)
            sets.append(s)
            nodes[idx]["ch"].append(s)
        nodes[idx]["set"] = frozenset().union(*sets)
        return nodes[idx]["set"], idx
    rec(root, -1)
    return nodes


# Compare the clusters of one family.
#   pred_root   root node of the reconstructed tree (dict tree)
#   tnodes      internal nodes of the reference tree (true_nodes())
# Returns a Counter: predicted, exact, polytomy, wrong, wrong_Spec, wrong_Dup (predicted clusters, root excluded);
# reference, ref_found, ref_nested_dup, ref_missed (reference clusters, root excluded).
def analyse(pred_root, tnodes):
    leaves_all = tnodes[0]["set"] if tnodes else frozenset()
    by_set = {n["set"]: n for n in tnodes}
    order = sorted(tnodes, key=lambda n: len(n["set"]))
    res = collections.Counter()
    pc = pred_clusters(pred_root)
    root_set = pc[-1][0]
    found = set()
    for s, label in pc:
        if s == root_set:
            continue
        res["predicted"] += 1
        if s in by_set:
            res["exact"] += 1
            found.add(s)
            continue
        host = next((n for n in order if s <= n["set"]), None)       # smallest reference node containing s
        ok = False
        if host is not None and len(host["ch"]) > 2:
            inside = [c for c in host["ch"] if c <= s]
            partial = [c for c in host["ch"] if (c & s) and not c <= s]
            ok = 2 <= len(inside) < len(host["ch"]) and not partial
        if ok:
            res["polytomy"] += 1
        else:
            res["wrong"] += 1
            res["wrong_" + label] += 1
    for n in tnodes:
        if n["set"] == leaves_all or n["parent"] == -1:
            continue
        res["reference"] += 1
        if n["set"] in found:
            res["ref_found"] += 1
        elif n["label"] == "Dup" and tnodes[n["parent"]]["label"] == "Dup":
            res["ref_nested_dup"] += 1
        else:
            res["ref_missed"] += 1
    return res


# Analyse every family, print the totals and optionally write one line per family.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--true-trees", required=True)
    ap.add_argument("--out")
    ap.add_argument("--dubious", choices=["dup", "spec"], default="dup", help="label of Ensembl dubious nodes (dup)")
    a = ap.parse_args()

    tot, rows = collections.Counter(), []
    for line in open(a.pred):
        line = line.strip()
        if not line:
            continue
        tid, nwk = line.split("\t", 1)
        path = os.path.join(a.true_trees, tid + ".json")
        if not os.path.exists(path):
            continue
        r = analyse(parse_newick(nwk), true_nodes(path, a.dubious))
        tot.update(r)
        rows.append((tid, r))
    P, E, Y, W = tot["predicted"], tot["exact"], tot["polytomy"], tot["wrong"]
    R, F = tot["reference"], tot["ref_found"]
    Rn = R - tot["ref_nested_dup"]
    print(f"families: {len(rows)}")
    print(f"predicted clusters (root excluded): {P}")
    print(f"  exact (also a reference cluster)            : {E:5d}  ({100 * E / max(P, 1):.1f}%)")
    print(f"  polytomy (binary resolution of a polytomy)  : {Y:5d}  ({100 * Y / max(P, 1):.1f}%)")
    print(f"  wrong (cut through a reference cluster)     : {W:5d}  ({100 * W / max(P, 1):.1f}%)   "
          f"Spec {tot['wrong_Spec']}, Dup {tot['wrong_Dup']}")
    print(f"  => cluster precision (exact + polytomy)     : {(E + Y) / max(P, 1):.3f}   strict (exact only): "
          f"{E / max(P, 1):.3f}")
    print(f"reference clusters (root excluded): {R}; found {F}; nested duplications (not produced by design) "
          f"{tot['ref_nested_dup']}; missed {tot['ref_missed']}")
    print(f"  => cluster recall: {F / max(Rn, 1):.3f}   (reference clusters other than nested duplications: {Rn})")
    if a.out:
        with open(a.out, "w") as f:
            f.write("tree_id\tpredicted\texact\tpolytomy\twrong\treference\tfound\tnested_dup\tmissed\n")
            for tid, r in rows:
                f.write(f"{tid}\t{r['predicted']}\t{r['exact']}\t{r['polytomy']}\t{r['wrong']}\t{r['reference']}\t"
                        f"{r['ref_found']}\t{r['ref_nested_dup']}\t{r['ref_missed']}\n")


if __name__ == "__main__":
    main()
