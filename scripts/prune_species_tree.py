#!/usr/bin/env python3
"""
Species trees for OrthoFinder's -s option, one per OrthoFinder input folder.

  python prune_species_tree.py st/species_tree.nwk of_in/pooled of_in/fam/ST001 ...

orthofinder -s needs a species tree whose leaves are exactly the species of its input folder (one <taxon>.fa file per
species, written by make_proteomes.py): a per-family folder holds only some of the species of the full tree. For every
folder given, the full species tree is pruned to the taxa of its .fa files (ldo_utils.prune: nodes left with one child
are collapsed) and written next to the folder as <folder>.sp.nwk: rooted topology, leaves = taxon ids, no branch
lengths.

OrthoFinder is given the true species tree on purpose: the LDO reconstruction also uses it, so both methods start from
the same information. Used for SwissTree only; simphy_prepare.py writes the pruned trees of the simphy runs itself.
"""
import glob
import os
import sys

from ldo_utils import first_tree_text, parse_newick, prune, write_topology


# Write the pruned species tree of every folder given on the command line.
def main():
    root = parse_newick(first_tree_text(sys.argv[1]))
    for d in sys.argv[2:]:
        d = d.rstrip("/")
        taxa = {os.path.basename(p)[:-3] for p in glob.glob(os.path.join(d, "*.fa"))}
        with open(d + ".sp.nwk", "w") as f:
            f.write(write_topology(prune(root, taxa)))
        print(f"{d}.sp.nwk: {len(taxa)} species")


if __name__ == "__main__":
    main()
