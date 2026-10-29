#!/usr/bin/env python3
"""
Turn downloaded Ensembl Compara gene trees (ensembl_download.py) into the inputs of the LDO pipeline: gene trees with
aligned sequences whose taxa are leaves of the species tree, and that species tree with taxon-id leaves.

  python ensembl_prepare.py --trees DATA/trees --species-newick DATA/compara_species_tree.nh --out OUT \
      [--max-genes 1500] [--min-species 2]

Species tree: the leaves of the Compara species tree are matched to the taxa of the gene trees, by taxon id (a
numeric leaf label) or by name (the leaf label against the taxonomy names of the gene-tree leaves, ignoring case and
punctuation, so 'homo_sapiens' matches 'Homo sapiens'). It is pruned to the matched taxa that occur in the kept
trees, leaves renamed to their taxon ids, unary nodes collapsed.

Gene trees: leaves without gene id or aligned sequence, and leaves whose taxon is not in the species tree, are
removed (unary nodes collapsed, branch lengths added); alignment columns with no residue left are removed (this does
not change any distance). Families with more than --max-genes genes or fewer than --min-species species are skipped.
Event types are kept as Ensembl gives them: "dubious" nodes stay dubious, and the scripts that read the true labels
(true_orthologs.py, cluster_analysis.py) choose how to label them with --dubious dup|spec.

Writes in OUT:
  trees/<tree_id>.json   the kept gene trees
  seqs/<tree_id>.fa      their sequences without gaps (input of make_proteomes.py, for OrthoFinder)
  species_tree.nwk       the species tree (leaves = taxon ids)
  prepare_report.tsv     per family: tree_id, genes, species, removed_leaves, and the internal nodes: spec, dup,
                         dubious, other (e.g. gene_split), spec_conflict (speciation nodes that the species tree
                         cannot explain: their children's species are not separated by one node of the species tree),
                         dubious_spec_ok (dubious nodes that would be valid speciations); status (kept / reason)
"""
import argparse
import collections
import glob
import json
import os
import re
import sys

from ldo_utils import (INVALID, common_prefix_len, event_type, first_tree_text, get_root, is_dubious, json_leaves,
                       leaf_gene, leaf_seq, parse_newick, prune, species_paths, tree_leaves, write_topology)


# Lower-case name without punctuation, for matching species names.
#   s   name
# Returns the normalised string ('Homo sapiens' and 'homo_sapiens' both give 'homosapiens').
def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


# Taxa of the gene trees and their names.
#   files   JSON gene tree files
# Returns {taxon id (string): set of normalised names (from the leaves' taxonomy fields)}.
def collect_taxa(files):
    taxa = {}
    for p in files:
        _, root = get_root(json.load(open(p)))
        for x in json_leaves(root) if root else []:
            tx = x.get("taxonomy")
            if isinstance(tx, dict) and "id" in tx:
                names = taxa.setdefault(str(int(tx["id"])), set())
                names.update(norm(v) for v in tx.values() if isinstance(v, str) and not v.isdigit())
    return taxa


# Rename the leaves of the Compara species tree to taxon ids.
#   root    root node of the Compara species tree (dict tree)
#   taxa    {taxon id: names} of collect_taxa()
# Returns the set of matched taxon ids; leaves that match no taxon (or a taxon already matched) get an empty name.
def match_species_tree(root, taxa):
    by_name = {}
    for t, names in sorted(taxa.items()):
        for n in names:
            by_name.setdefault(n, t)
    matched = set()
    for x in tree_leaves(root):
        lab = x["name"]
        t = lab if lab.isdigit() and lab in taxa else by_name.get(norm(lab))
        if t is None and "_" in lab:                 # labels such as 'Homo_sapiens_9606' or 'homo_sapiens.GRCh38'
            t = by_name.get(norm(lab.split(".")[0]))
        x["name"] = t if t is not None and t not in matched else ""
        if t is not None:
            matched.add(t)
    return matched


# Copy of a JSON gene tree restricted to some leaves; unary nodes are collapsed (the child keeps its own fields, with
# the branch lengths added).
#   x      JSON node
#   keep   function of a leaf node -> True if the leaf is kept
# Returns the pruned copy, or None if no leaf is kept.
def prune_json(x, keep):
    ch = x.get("children")
    if not ch:
        return dict(x) if keep(x) else None
    kids = [c for c in (prune_json(c, keep) for c in ch) if c is not None]
    if not kids:
        return None
    if len(kids) == 1:
        c = kids[0]
        c["branch_length"] = (c.get("branch_length") or 0.0) + (x.get("branch_length") or 0.0)
        return c
    out = dict(x)
    out["children"] = kids
    return out


# Whether children with these species sets can be separated by speciations only: with y = lca_S(all their species),
# the species of each child lie below one child of y, and children below the same child of y form (recursively) a
# valid group again.
#   child_sets   list of the sets of taxa below each child
#   sp           {taxon: path} of the species tree
# Returns True or False.
def spec_compatible(child_sets, sp):
    if len(child_sets) < 2:
        return True
    union = set().union(*child_sets)
    if len(union) < 2:
        return False
    k = common_prefix_len([sp[t] for t in union])
    groups = {}
    for cs in child_sets:
        side = {sp[t][k] for t in cs}
        if len(side) != 1:
            return False
        groups.setdefault(next(iter(side)), []).append(cs)
    if len(groups) < 2:
        return False
    return all(spec_compatible(g, sp) for g in groups.values() if len(g) > 1)


# Count the internal nodes of a gene tree by event, and check the speciations against the species tree.
#   root   root node of the (pruned) JSON gene tree
#   sp     {taxon: path} of the species tree
# Returns a dict: spec, dup, dubious, other, spec_conflict, dubious_spec_ok.
def count_events(root, sp):
    c = dict(spec=0, dup=0, dubious=0, other=0, spec_conflict=0, dubious_spec_ok=0)

    # Post-order walk; returns the set of taxa below x.
    #   x   current JSON node
    def rec(x):
        if not x.get("children"):
            return {str(int(x["taxonomy"]["id"]))}
        sets = [rec(ch) for ch in x["children"]]
        ev = event_type(x)
        ok = spec_compatible(sets, sp)
        if is_dubious(ev):
            c["dubious"] += 1
            c["dubious_spec_ok"] += ok
        elif "speciation" in ev and "duplication" not in ev:
            c["spec"] += 1
            c["spec_conflict"] += not ok
        elif "duplication" in ev:
            c["dup"] += 1
        else:
            c["other"] += 1
        return set().union(*sets)
    rec(root)
    return c


# Remove the alignment columns where no leaf has a residue (in place).
#   root   root node of the JSON gene tree (all leaves aligned, same length)
def drop_empty_columns(root):
    leaves = json_leaves(root)
    seqs = [leaf_seq(x).upper() for x in leaves]
    keep = [i for i in range(len(seqs[0])) if any(s[i].encode() not in INVALID.tobytes() for s in seqs)]
    for x, s in zip(leaves, seqs):
        new = "".join(s[i] for i in keep)
        seq = x["sequence"]
        if isinstance(seq.get("mol_seq"), dict):
            seq["mol_seq"]["seq"] = new
        else:
            seq["seq"] = new


# Prepare every tree, write the trees, the species tree and the report.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="folder of the downloaded JSON gene trees")
    ap.add_argument("--species-newick", required=True, help="Compara species tree (Newick)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-genes", type=int, default=1500, help="skip larger families (1500)")
    ap.add_argument("--min-species", type=int, default=2, help="skip families with fewer species (2)")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.trees, "*.json")))
    if not files:
        sys.exit(f"no .json files in {a.trees}")
    taxa = collect_taxa(files)
    sroot = parse_newick(first_tree_text(a.species_newick))
    matched = match_species_tree(sroot, taxa)
    missing = sorted(set(taxa) - matched)
    print(f"species tree: {len(tree_leaves(sroot))} leaves; gene trees: {len(taxa)} taxa, {len(matched)} matched")
    if missing:
        print(f"WARNING: {len(missing)} taxa not in the species tree (their genes are removed): {missing[:20]}",
              file=sys.stderr)
    if not matched:
        sample = [x["name"] for x in tree_leaves(parse_newick(first_tree_text(a.species_newick)))[:8]]
        sys.exit(f"no taxon matched; species-tree leaf labels look like {sample}, gene-tree taxa like "
                 f"{list(taxa.items())[:3]}")
    sp = species_paths(prune(sroot, matched))

    os.makedirs(os.path.join(a.out, "trees"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "seqs"), exist_ok=True)
    used, rows = set(), []
    for p in files:
        tid, root = get_root(json.load(open(p)))
        tid = tid or os.path.basename(p)[:-5]
        n0 = len(json_leaves(root))
        # alignment length: the most common sequence length (all leaves have it in an Ensembl answer)
        lengths = collections.Counter(len(leaf_seq(x)) for x in json_leaves(root) if leaf_seq(x))
        L = lengths.most_common(1)[0][0] if lengths else None

        # A leaf is kept if it has a gene id, an aligned sequence of the family's length, and a taxon of the tree.
        #   x   leaf node
        def keep(x):
            s = leaf_seq(x)
            return (leaf_gene(x) is not None and s is not None and len(s) == L
                    and str(int(x["taxonomy"]["id"])) in sp)
        root = prune_json(root, keep) if L else None
        genes = json_leaves(root) if root else []
        species = {str(int(x["taxonomy"]["id"])) for x in genes}
        c = count_events(root, sp) if root else dict.fromkeys(
            ["spec", "dup", "dubious", "other", "spec_conflict", "dubious_spec_ok"], 0)
        status = "kept"
        if len(genes) > a.max_genes:
            status = f"skipped: more than {a.max_genes} genes"
        elif len(species) < a.min_species:
            status = f"skipped: fewer than {a.min_species} species"
        rows.append((tid, len(genes), len(species), n0 - len(genes), c["spec"], c["dup"], c["dubious"], c["other"],
                     c["spec_conflict"], c["dubious_spec_ok"], status))
        if status != "kept":
            continue
        drop_empty_columns(root)
        used |= species
        with open(os.path.join(a.out, "trees", tid + ".json"), "w") as f:
            json.dump({"tree": {"id": tid, "tree": root}}, f)
        with open(os.path.join(a.out, "seqs", tid + ".fa"), "w") as f:
            for x in genes:
                f.write(f">{leaf_gene(x)}\n{leaf_seq(x).replace('-', '').upper()}\n")

    with open(os.path.join(a.out, "species_tree.nwk"), "w") as f:
        f.write(write_topology(prune(sroot, used)) + "\n")
    with open(os.path.join(a.out, "prepare_report.tsv"), "w") as f:
        f.write("tree_id\tgenes\tspecies\tremoved_leaves\tspec\tdup\tdubious\tother\tspec_conflict\t"
                "dubious_spec_ok\tstatus\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    kept = [r for r in rows if r[-1] == "kept"]
    tot = [sum(r[i] for r in kept) for i in range(1, 10)]
    print(f"{len(kept)}/{len(rows)} families kept: {tot[0]} genes ({tot[2]} leaves removed), {len(used)} species; "
          f"internal nodes: {tot[3]} spec ({tot[7]} in conflict with the species tree), {tot[4]} dup, "
          f"{tot[5]} dubious ({tot[8]} would be valid speciations), {tot[6]} other; "
          f"{sum(1 for r in kept if r[6] == 0)} families without dubious node")


if __name__ == "__main__":
    main()
