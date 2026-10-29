#!/usr/bin/env python3
"""
True ortholog pairs of labelled gene trees (Ensembl-like JSON): two genes are orthologs iff their last common
ancestor is a speciation node.

  python true_orthologs.py --trees TREES_DIR --out OUT_DIR [--dubious dup|spec]

--dubious: label of the Ensembl "dubious" nodes (default dup); other trees have none.

Outputs in OUT_DIR:
  genes.tsv       tree_id  gene  taxon_id
  orthologs.tsv   tree_id  gene1  gene2
"""
import argparse
import glob
import os

from ldo_utils import event_label, leaf_gene, leaf_taxon, load_json_tree


# Flatten a JSON gene tree into a node list (depth first with a stack, so parents come before their children).
#   tid       tree id (used to name leaves without a gene id: <tid>_leaf<index>)
#   root      root node of the JSON gene tree
#   dubious   label of the dubious nodes, "dup" or "spec" (see ldo_utils.event_label())
# Returns a list of dict(children=[indices], label="Spec" / "Dup" / "Leaf", gene (leaves), tax (leaves)).
def tree_nodes(tid, root, dubious="dup"):
    nodes, stack = [], [(root, -1)]
    while stack:
        x, parent = stack.pop()
        idx = len(nodes)
        ch = x.get("children") or []
        if ch:
            nodes.append(dict(children=[], label=event_label(x, dubious), gene=None, tax=None))
            for c in ch:
                stack.append((c, idx))
        else:
            gene = leaf_gene(x) or f"{tid}_leaf{idx}"
            nodes.append(dict(children=[], label="Leaf", gene=gene, tax=leaf_taxon(x)))
        if parent >= 0:
            nodes[parent]["children"].append(idx)
    return nodes


# Ortholog pairs of a node list: for every speciation node, all pairs of genes taken in two different children.
#   nodes   node list of tree_nodes()
# Returns a list of (gene, gene).
def ortholog_pairs(nodes):
    pairs, below = [], [None] * len(nodes)
    for v in reversed(range(len(nodes))):          # children have larger indices than their parent
        nd = nodes[v]
        if nd["label"] == "Leaf":
            below[v] = [nd["gene"]]
            continue
        parts = [below[c] for c in nd["children"]]
        if nd["label"] == "Spec":
            for i in range(len(parts)):
                for j in range(i + 1, len(parts)):
                    pairs += [(g, h) for g in parts[i] for h in parts[j]]
        below[v] = [g for p in parts for g in p]
        for c in nd["children"]:
            below[c] = None
    return pairs


# Write genes.tsv and orthologs.tsv for all the trees of a folder.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="folder of *.json gene trees")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dubious", choices=["dup", "spec"], default="dup", help="label of Ensembl dubious nodes (dup)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "genes.tsv"), "w") as fg, open(os.path.join(a.out, "orthologs.tsv"), "w") as fo:
        for path in sorted(glob.glob(os.path.join(a.trees, "**", "*.json"), recursive=True)):
            tid, root = load_json_tree(path)
            tid = tid or os.path.splitext(os.path.basename(path))[0]
            nodes = tree_nodes(tid, root, a.dubious)
            for nd in nodes:
                if nd["label"] == "Leaf":
                    fg.write(f"{tid}\t{nd['gene']}\t{nd['tax']}\n")
            for g, h in ortholog_pairs(nodes):
                fo.write(f"{tid}\t{g}\t{h}\n")


if __name__ == "__main__":
    main()
