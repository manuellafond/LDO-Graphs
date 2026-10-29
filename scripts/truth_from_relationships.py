#!/usr/bin/env python3
"""
Curated ortholog pairs from the SwissTree gene_relationships.txt files ("geneA geneB ortholog|paralog|...").

  python truth_from_relationships.py --raw swisstree_raw --genes ldo/genes.tsv --out curator_orthologs.tsv

Names in these files drop the accession (HoxA9_HUMAN vs HoxA9_HUMAN_P31269); they are matched to the converted genes by
full name or by the name without its last token. Output: tree_id, gene1, gene2 (ortholog pairs). Families without a
gene_relationships.txt file are skipped.
"""
import argparse
import collections
import os

from ldo_utils import read_genes


# Lookup table of the genes of a family.
#   genes   list of the gene names of the family (e.g. HoxA9_HUMAN_P31269)
# Returns {lower-case key: gene}, the keys being the full name and the name without its last _token (hoxa9_human).
def name_keys(genes):
    key = {}
    for g in genes:
        key[g.lower()] = g
        p = g.split("_")
        if len(p) >= 3:
            key.setdefault("_".join(p[:-1]).lower(), g)
    return key


# Write the ortholog pairs of every family that has a gene_relationships.txt file.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--genes", required=True, help="genes.tsv of true_orthologs.py")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fam_genes = collections.defaultdict(list)
    for tid, g, _ in read_genes(a.genes):
        fam_genes[tid].append(g)
    with open(a.out, "w") as out:
        for f in sorted(fam_genes):
            path = os.path.join(a.raw, f, "gene_relationships.txt")
            if not os.path.exists(path):
                print(f, "no gene_relationships.txt")
                continue
            key, cnt, miss = name_keys(fam_genes[f]), collections.Counter(), set()
            for line in open(path, errors="replace"):
                p = line.split()
                if len(p) < 3:
                    continue
                x, y = key.get(p[0].lower()), key.get(p[1].lower())
                if x is None or y is None:
                    miss.update(q for q, k in ((p[0], x), (p[1], y)) if k is None)
                    continue
                cnt[p[2]] += 1
                if p[2] == "ortholog":
                    out.write(f"{f}\t{min(x, y)}\t{max(x, y)}\n")
            print(f, dict(cnt), f"names not matched: {len(miss)}")


if __name__ == "__main__":
    main()
