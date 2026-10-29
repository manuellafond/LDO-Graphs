#!/usr/bin/env python3
"""
Convert one simphy replicate (s_tree.trees, g_trees<NNN>.trees, dataset_<NNN>.fasta) into the inputs of the LDO
pipeline and of OrthoFinder.

  python simphy_prepare.py --run simphy/<setting>/<ssim_...>/1 --out OUT

Gene trees are labelled with simphy's own events (--labels simphy, the default): the locus-tree event (Sp / Dup) of
every gene-tree node, read from <NNN>.mapsl (locus-tree node -> event) and <NNN>l1g.maplg (gene-tree node -> locus-tree
node; the gene-tree nodes are listed in post-order). Without a .mapsl file (simphy writes none when the duplication and
loss rates are 0, as in dl0), all internal nodes are speciations. With --labels lca, the gene
trees are labelled by LCA mapping on the species tree instead: with M(v) = lca_S(species below v), an internal node v is
a duplication iff M(v) = M(c) for some child c, otherwise a speciation. LCA mapping misses the duplications followed by
complementary losses (hidden paralogy: the gene tree looks like a speciation), so its truth calls some paralogs
orthologs (18% of the true pairs over the 25 replicates of the paper) and favours methods that also label by the species
tree. A leaf "xx_yy_zz" belongs to species "xx".
Gene ids are made unique across the families of the replicate: "<family>.<leaf>" (e.g. 005.19_17_0).

Writes in OUT:
  trees/<family>.json          labelled gene tree, aligned sequences in the leaves
  species_tree.nwk             species tree topology, leaves = species labels
  of_in/fam/<family>/<sp>.fa   one OrthoFinder input folder per family (sequences without gaps)
  of_in/fam/<family>.sp.nwk    species tree pruned to the species of that family (for orthofinder -s)
  of_in/pooled/<sp>.fa         all families of the replicate together, and of_in/pooled.sp.nwk
  prepare_report.tsv           family, genes, species, Spec nodes, Dup nodes
"""
import argparse
import glob
import json
import os
import re
import sys

from ldo_utils import (common_prefix_len, first_tree_text, parse_newick, prune, read_fasta, species_paths,
                       tree_leaves, write_topology)


# Species of a simphy gene-tree leaf.
#   leaf   leaf name "xx_yy_zz" (species xx, locus yy, individual zz)
# Returns the species label "xx".
def species_of(leaf):
    return leaf.split("_")[0]


# Label every internal node of a gene tree by LCA mapping (sets node["event"] = "speciation" / "duplication").
#   groot   root node of the gene tree (dict tree)
#   paths   {species: path} of the species tree (ldo_utils.species_paths())
# Returns (number of Spec nodes, number of Dup nodes).
def label_events(groot, paths):
    counts = [0, 0]

    # Post-order walk: stores x["_lca"] = path of lca_S(species below x) and labels x.
    #   x   current node
    # Returns the set of species below x.
    def rec(x):
        if not x["children"]:
            sp = {species_of(x["name"])}
        else:
            sp = set().union(*[rec(c) for c in x["children"]])
        ps = [paths[s] for s in sp]
        x["_lca"] = tuple(ps[0][:common_prefix_len(ps)])
        if x["children"]:
            dup = any(c["_lca"] == x["_lca"] for c in x["children"])
            x["event"] = "duplication" if dup else "speciation"
            counts[dup] += 1
        return sp
    rec(groot)
    return counts


# Label every internal node of a gene tree with simphy's events (sets node["event"] = "speciation" / "duplication").
#   groot   root node of the gene tree (dict tree)
#   run     simphy replicate folder
#   fam     family number as in the file names (e.g. "005")
# Returns (number of Spec nodes, number of Dup nodes).
def label_simphy(groot, run, fam):
    nodes = []

    # Post-order list of the nodes (the order of the .maplg file).
    #   x   current node
    def post(x):
        for c in x["children"]:
            post(c)
        nodes.append(x)
    post(groot)
    mapsl, maplg = os.path.join(run, f"{fam}.mapsl"), os.path.join(run, f"{fam}l1g.maplg")
    if not os.path.exists(mapsl):          # no duplication or loss in the run (no locus-tree file): all speciations
        events = ["Sp" if x["children"] else "Leaf" for x in nodes]
    else:
        lt_event = {}
        for line in open(mapsl).read().splitlines()[1:]:              # Lt_node, Lt_paralog, Lt_event, St_node
            f = line.split("\t")
            lt_event[f[0].strip("'")] = f[2]
        rows = [line.split("\t") for line in open(maplg).read().splitlines()[1:]]   # Gt_node, Lt_node, Lt_paralog
        if len(rows) != len(nodes):
            sys.exit(f"{maplg}: {len(rows)} nodes, gene tree has {len(nodes)}")
        events = []
        for x, r in zip(nodes, rows):
            if not x["children"] and r[0].strip("'") != x["name"]:
                sys.exit(f"{maplg}: leaf {r[0]} where the gene tree has {x['name']} (order differs)")
            events.append(lt_event[r[1].strip("'")])
    counts = [0, 0]
    for x, ev in zip(nodes, events):
        if x["children"]:
            if ev not in ("Sp", "Dup"):
                sys.exit(f"{run} family {fam}: unexpected event {ev} (transfers / gene conversion are not supported)")
            x["event"] = "duplication" if ev == "Dup" else "speciation"
            counts[ev == "Dup"] += 1
    return counts


# Ensembl-like JSON of a labelled gene tree.
#   x      root node of the labelled gene tree (dict tree, after label_events())
#   fam    family number (gene ids become "<fam>.<leaf>")
#   seqs   {leaf name: aligned sequence}
# Returns the JSON node (dict).
def to_json(x, fam, seqs):
    if not x["children"]:
        return {"id": {"accession": f"{fam}.{x['name']}"},
                "taxonomy": {"id": int(species_of(x["name"]))},
                "branch_length": x["bl"] if x["bl"] is not None else 0.0,
                "sequence": {"mol_seq": {"seq": seqs[x["name"]].upper()}}}
    return {"events": {"type": x["event"]},
            "branch_length": x["bl"] if x["bl"] is not None else 0.0,
            "children": [to_json(c, fam, seqs) for c in x["children"]]}


# Write the input folder of one OrthoFinder run. OrthoFinder takes a folder with one FASTA file per species (a
# "proteome"); this creates the folder and writes <species>.fa with the genes of that species. The sequences are the
# raw ones (no alignment gaps): OrthoFinder does its own similarity search. Called once per family
# (of_in/fam/<family>/, for the per-family runs) and once with all the families together (of_in/pooled/, for the
# pooled run).
#   folder        output folder (created if needed)
#   per_species   {species: [(gene, sequence)]}
def write_proteomes(folder, per_species):
    os.makedirs(folder, exist_ok=True)
    for sp, recs in per_species.items():
        with open(os.path.join(folder, sp + ".fa"), "w") as f:
            f.write("".join(f">{g}\n{s}\n" for g, s in recs))


# Convert every family of the replicate.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="simphy replicate folder (contains s_tree.trees)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", choices=("simphy", "lca"), default="simphy",
                    help="true events: simphy's own (default) or LCA mapping on the species tree")
    a = ap.parse_args()

    sroot = parse_newick(first_tree_text(os.path.join(a.run, "s_tree.trees")))
    paths = species_paths(sroot)
    os.makedirs(os.path.join(a.out, "trees"), exist_ok=True)
    with open(os.path.join(a.out, "species_tree.nwk"), "w") as f:
        f.write(write_topology(sroot) + "\n")

    fams = sorted(re.search(r"g_trees(\d+)\.trees$", p).group(1)
                  for p in glob.glob(os.path.join(a.run, "g_trees*.trees")))
    pooled, report = {}, []
    for fam in fams:
        groot = parse_newick(first_tree_text(os.path.join(a.run, f"g_trees{fam}.trees")))
        seqs = {h.split()[0]: s for h, s in read_fasta(os.path.join(a.run, f"dataset_{fam}.fasta"))}
        lv = tree_leaves(groot)
        missing = [x["name"] for x in lv if x["name"] not in seqs]
        unknown = [x["name"] for x in lv if species_of(x["name"]) not in paths]
        if missing or unknown:
            sys.exit(f"{a.run} family {fam}: leaves without sequence {missing[:3]} / unknown species {unknown[:3]}")
        if len({len(seqs[x["name"]]) for x in lv}) != 1:
            sys.exit(f"{a.run} family {fam}: sequences are not aligned (different lengths)")
        if len(lv) < 2:
            n_spec, n_dup = 0, 0
        elif a.labels == "simphy":
            n_spec, n_dup = label_simphy(groot, a.run, fam)
        else:
            n_spec, n_dup = label_events(groot, paths)
        with open(os.path.join(a.out, "trees", fam + ".json"), "w") as f:
            json.dump({"tree": {"id": fam, "tree": to_json(groot, fam, seqs)}}, f)

        per = {}
        for x in lv:
            rec = (f"{fam}.{x['name']}", seqs[x["name"]].replace("-", "").upper())
            per.setdefault(species_of(x["name"]), []).append(rec)
            pooled.setdefault(species_of(x["name"]), []).append(rec)
        folder = os.path.join(a.out, "of_in", "fam", fam)
        write_proteomes(folder, per)
        with open(folder + ".sp.nwk", "w") as f:
            f.write(write_topology(prune(sroot, set(per))) + "\n")
        report.append((fam, len(lv), len(per), n_spec, n_dup))

    write_proteomes(os.path.join(a.out, "of_in", "pooled"), pooled)
    with open(os.path.join(a.out, "of_in", "pooled.sp.nwk"), "w") as f:
        f.write(write_topology(prune(sroot, set(pooled))) + "\n")
    with open(os.path.join(a.out, "prepare_report.tsv"), "w") as f:
        f.write("family\tgenes\tspecies\tspec_nodes\tdup_nodes\n")
        for r in report:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{a.run}: {len(report)} families, {sum(r[1] for r in report)} genes, "
          f"{sum(r[3] for r in report)} Spec / {sum(r[4] for r in report)} Dup nodes -> {a.out}")


if __name__ == "__main__":
    main()
