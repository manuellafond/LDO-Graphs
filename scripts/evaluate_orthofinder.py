#!/usr/bin/env python3
"""
Score the ortholog pairs predicted by OrthoFinder against true ortholog pairs, as evaluate_orthologs.py scores the
LDO trees.

  python evaluate_orthofinder.py --truth ldo/orthologs.tsv --genes ldo/genes.tsv --run RESULTS_DIR [RESULTS_DIR ...] \
      [--families ST001 ...] [--label NAME] [--per-family out.tsv]

--truth     true ortholog pairs: tree_id, gene1, gene2 (true_orthologs.py or truth_from_relationships.py).
--genes     genes.tsv of true_orthologs.py (tree_id, gene, taxon): the genes given to OrthoFinder, and their family.
--run       OrthoFinder Results_* folders: one (a pooled run) or several (one per family); the predicted pairs of all
            the folders given are merged before scoring.
--families  score only these families (default: all).
--label     name printed at the start of the summary line and in the first column of --per-family.
Predicted pairs are read from the files Orthologues/Orthologues_*/A__v__B.tsv that OrthoFinder writes (OrthoFinder's
own folder names). Truth and predictions are both limited to the genes of --genes; a family is scored only if it has
at least one true pair. A predicted pair counts only if both genes are in the same family; pairs joining two families
(possible only in a pooled run) cannot be true, and are reported in the summary line (cross-family=) but not scored.

Prints one line: label, pooled precision / recall / F1, numbers of true / predicted / correct pairs, cross-family pairs.
--per-family writes one line per scored family: run (label), family, truth, pred, tp, precision, recall, f1.
"""
import argparse
import collections
import glob
import os

from ldo_utils import prf, read_genes


# A gene pair in sorted order.
#   a, b   gene ids
# Returns (a, b) or (b, a), the smaller first.
def norm(a, b):
    return (a, b) if a < b else (b, a)


# Ortholog pairs listed in OrthoFinder's Orthologues/ folder of a results folder. Each file A__v__B.tsv has a header
# line, then one line per orthology relation: orthogroup, genes of species A (comma-separated), genes of species B;
# every gene of A is an ortholog of every gene of B on that line.
#   res   OrthoFinder Results_* folder
# Returns a set of sorted gene pairs (empty, with a warning, if the folder has no Orthologues/).
def read_orthologs(res):
    pairs = set()
    for f in glob.glob(os.path.join(res, "Orthologues", "Orthologues_*", "*__v__*.tsv")):
        with open(f) as fh:
            next(fh)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                for x in (g.strip() for g in p[1].split(",") if g.strip()):
                    for y in (g.strip() for g in p[2].split(",") if g.strip()):
                        pairs.add(norm(x, y))
    if not os.path.isdir(os.path.join(res, "Orthologues")):
        print(f"WARNING: {res} has no Orthologues folder (failed OrthoFinder run?); counted as no prediction")
    return pairs


# Score the runs: print the totals and optionally write one line per family.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--run", nargs="+", required=True)
    ap.add_argument("--families", nargs="*")
    ap.add_argument("--label", default="")
    ap.add_argument("--per-family")
    a = ap.parse_args()

    # Genes file: taxon and family of every gene (the genes given to OrthoFinder).
    tax, fam = {}, {}
    for f, g, t in read_genes(a.genes):
        tax[g], fam[g] = t, f

    # True ortholog pairs per family. A pair is kept if its family is requested (--families; all if not given) and
    # both genes are in the genes file under that family, so that truth and predictions are limited to the same genes
    # (a safety net: with the truth and genes files of the same run, as in the pipelines, nothing is removed).
    # Families without any true pair (e.g. a single species) do not appear in `truth`, so they are not scored.
    truth = collections.defaultdict(set)
    for f, x, y in read_genes(a.truth):
        if (not a.families or f in a.families) and fam.get(x) == f and fam.get(y) == f:
            truth[f].add(norm(x, y))

    # Predicted pairs: merged over all the --run folders. Pairs of two genes of the same species are dropped (they
    # cannot be orthologs; OrthoFinder does not list any, this is a safeguard), then pairs with a gene that is not in
    # the genes file.
    pred = set()
    for r in a.run:
        pred |= {p for p in read_orthologs(r) if tax.get(p[0]) != tax.get(p[1])}
    pred = {p for p in pred if p[0] in fam and p[1] in fam}

    # Pairs joining two families (possible only in a pooled run, when OrthoFinder puts genes of two of our families
    # in one orthogroup): they cannot be true, but belong to no family; counted and reported, not scored.
    cross = sum(1 for p in pred if fam[p[0]] != fam[p[1]])

    # The other predicted pairs, grouped by family (only the scored families).
    by_fam = collections.defaultdict(set)
    for p in pred:
        if fam[p[0]] == fam[p[1]] and fam[p[0]] in truth:
            by_fam[fam[p[0]]].add(p)

    # Per-family counts: true pairs, predicted pairs and correct pairs (predicted and true), with precision / recall
    # / F1; the counts are also summed over the families for the totals (micro average).
    rows = []
    n_true = n_pred = n_correct = 0
    for f in sorted(truth):
        true_pairs = truth[f]
        pred_pairs = by_fam.get(f, set())
        correct = len(true_pairs & pred_pairs)
        precision, recall, f1 = prf(len(true_pairs), len(pred_pairs), correct)
        rows.append((a.label, f, len(true_pairs), len(pred_pairs), correct, precision, recall, f1))
        n_true += len(true_pairs)
        n_pred += len(pred_pairs)
        n_correct += correct
    p, r, f1 = prf(n_true, n_pred, n_correct)

    # One summary line: pooled precision / recall / F1, counts, cross-family pairs.
    print(f"{a.label:<22} P={p:.3f} R={r:.3f} F1={f1:.3f}   truth={n_true} pred={n_pred} tp={n_correct}  "
          f"cross-family={cross}")

    # Optional --per-family file: one line per scored family.
    if a.per_family:
        with open(a.per_family, "w") as fh:
            fh.write("run\tfamily\ttruth\tpred\ttp\tprecision\trecall\tf1\n")
            for row in rows:
                fh.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in row) + "\n")


if __name__ == "__main__":
    main()
