#!/usr/bin/env python3
"""
Predict an LDO graph from sequence distances: for every gene g and every other species s, add the arc g -> h where h
is the gene of species s closest to g (paralogs included). Ties are broken by gene name and flagged.

  python predict_ldo.py --trees TREES_DIR --out OUT_DIR [--model poisson|kimura|p|jc] [--min-overlap 50]

TREES_DIR holds Ensembl-like JSON gene trees whose leaves carry the aligned sequences (sequence.mol_seq.seq); only the
leaves are used. --model: poisson / kimura / p for proteins, jc (Jukes-Cantor) for DNA.

Outputs in OUT_DIR:
  genes.tsv       tree_id  gene  taxon
  ldo_edges.tsv   tree_id  source  source_taxon  target  target_taxon  distance  tie
"""
import argparse
import glob
import os
import sys

import numpy as np

from ldo_utils import alignment_matrix, distance_row, load_alignment

EPS = 1e-9


# Closest gene of every other species, for every gene of an aligned family (the predicted LDO arcs).
#   fam           aligned family: list of (gene, taxon, aligned sequence) (ldo_utils.load_alignment())
#   model         distance model (see ldo_utils.distance_row())
#   min_overlap   minimum number of compared columns; pairs below it are never chosen
# Returns {(i, j): (distance, tie)}: arc fam[i] -> fam[j], tie = 1 if several genes of j's species were equally close
# (the first by gene name is kept).
def predict_family(fam, model, min_overlap):
    genes = [g for g, _, _ in fam]
    tax = [t for _, t, _ in fam]
    A, V = alignment_matrix([s for _, _, s in fam])
    species = sorted(set(tax))
    idx = {s: np.array([k for k, t in enumerate(tax) if t == s]) for s in species}
    arcs = {}
    for i in range(len(fam)):
        d = distance_row(A, V, i, model, min_overlap)
        for s in species:
            if s == tax[i]:
                continue
            cand = idx[s]
            dc = d[cand]
            if not np.isfinite(dc).any():
                continue
            best = float(dc.min())
            near = sorted((int(cand[k]) for k in range(len(cand)) if dc[k] - best <= EPS), key=lambda k: genes[k])
            arcs[(i, near[0])] = (best, int(len(near) > 1))
    return arcs


# Predict the LDO graph of every family of the folder and write genes.tsv and ldo_edges.tsv.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", choices=["poisson", "kimura", "p", "jc"], default="poisson")
    ap.add_argument("--min-overlap", type=int, default=50, help="minimum number of compared columns")
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.trees, "*.json")))
    if not files:
        sys.exit(f"No .json files in {os.path.abspath(a.trees)}")
    os.makedirs(a.out, exist_ok=True)
    n_trees = n_arcs = 0
    with open(os.path.join(a.out, "genes.tsv"), "w") as fg, open(os.path.join(a.out, "ldo_edges.tsv"), "w") as fe:
        fe.write("tree_id\tsource\tsource_taxon\ttarget\ttarget_taxon\tdistance\ttie\n")
        for path in files:
            tid = os.path.basename(path)[:-5]
            fam = load_alignment(path)
            if not fam or len(fam) < 2:
                print(f"{tid}: sequences missing or not aligned, skipped", file=sys.stderr)
                continue
            arcs = predict_family(fam, a.model, a.min_overlap)
            for g, t, _ in fam:
                fg.write(f"{tid}\t{g}\t{t}\n")
            for (i, j), (d, tie) in sorted(arcs.items()):
                fe.write(f"{tid}\t{fam[i][0]}\t{fam[i][1]}\t{fam[j][0]}\t{fam[j][1]}\t{d:.6g}\t{tie}\n")
            n_trees += 1
            n_arcs += len(arcs)
    print(f"{n_trees} families, {n_arcs} arcs", file=sys.stderr)


if __name__ == "__main__":
    main()
