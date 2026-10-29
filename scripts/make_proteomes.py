#!/usr/bin/env python3
"""
Build OrthoFinder input proteomes (one FASTA per species, named <taxon>.fa) from the raw family sequences written by
swisstree_convert.py.

  python make_proteomes.py --genes ldo/genes.tsv --seqs st/seqs --out of_in

Writes OUT/pooled/<taxon>.fa (all families together) and OUT/fam/<family>/<taxon>.fa (one folder per family).
"""
import argparse
import collections
import os

from ldo_utils import read_fasta, read_genes


# Write the input folder of one OrthoFinder run. OrthoFinder takes a folder with one FASTA file per species (a
# "proteome"); this creates the folder and writes <taxon>.fa with the genes of that taxon. The sequences are the raw
# ones (no alignment gaps): OrthoFinder does its own similarity search. Called once per family (of_in/fam/<family>/,
# for the per-family runs) and once with all the families together (of_in/pooled/, for the pooled run).
#   folder      output folder (created if needed)
#   per_taxon   {taxon: [(gene, sequence)]}
def write_proteomes(folder, per_taxon):
    os.makedirs(folder, exist_ok=True)
    for t, recs in per_taxon.items():
        with open(os.path.join(folder, t + ".fa"), "w") as f:
            f.write("".join(f">{g}\n{s}\n" for g, s in recs))


# Split the sequences of every family by taxon, per family and pooled.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="genes.tsv of true_orthologs.py (tree_id, gene, taxon)")
    ap.add_argument("--seqs", required=True, help="folder with <family>.fa")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tax, fams = {}, collections.defaultdict(list)
    for tid, g, t in read_genes(a.genes):
        tax[g] = t
        fams[tid].append(g)
    pooled = collections.defaultdict(list)
    for fam in sorted(fams):
        per = collections.defaultdict(list)
        for h, s in read_fasta(os.path.join(a.seqs, fam + ".fa")):
            g = h.split()[0]
            if g in tax:
                per[tax[g]].append((g, s))
                pooled[tax[g]].append((g, s))
        write_proteomes(os.path.join(a.out, "fam", fam), per)
    write_proteomes(os.path.join(a.out, "pooled"), pooled)
    print(f"{len(fams)} families, {len(pooled)} species, {sum(map(len, pooled.values()))} genes -> {a.out}")


if __name__ == "__main__":
    main()
