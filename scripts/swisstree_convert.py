#!/usr/bin/env python3
"""
Convert the SwissTree ftp release to the inputs of the LDO pipeline: Ensembl-like JSON gene trees with MAFFT
alignments in the leaves, and a species tree whose leaves are NCBI taxon ids.

  python swisstree_convert.py --raw swisstree_raw --out OUT [--mafft mafft] [--match-gene-name] [--min-length-frac F]

Per family the curated tree consensus_tree.nhx is used:
  - an internal node is a duplication if it carries Ev=duplication / D=Y / D=T (case-insensitive), else a speciation;
  - families with horizontal / endosymbiotic transfer (Ev=EGT/HGT) are skipped;
  - the species of a leaf is the token of its name that is a species code of speciestree.nhx (POP1_HUMAN_Q9HBU9 ->
    HUMAN -> 9606); leaves of other species, and leaves without a sequence in sequences.fst / sequences.txt (matched by
    name or accession; with --match-gene-name also by the GN= field), are pruned (unary nodes collapsed);
  - with --min-length-frac F, sequences shorter than F x the family median length (fragments) are pruned too.

Outputs in OUT: trees/<family>.json, species_tree.nwk, seqs/<family>.fa (raw) and seqs/<family>.aln.fa (aligned).
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from ldo_utils import first_tree_text, iter_nodes, parse_newick, prune, read_fasta, write_topology

MIN_GENES = 4               # families with fewer genes (after pruning) are skipped
MIN_SPECIES = 3             # ... or fewer species
MIN_SEQ_FRACTION = 0.6      # skipped if fewer leaves than this fraction have a sequence
MIN_KEPT_FRACTION = 0.5     # skipped if fewer leaves than this fraction remain (known species + sequence)


# Species tree of the release. Its leaves are named CODE__name__taxid; they are renamed to their taxon id (leaves
# without a taxon id get an empty name).
#   path   speciestree.nhx
# Returns (root node of the renamed tree, {species code: taxon id}).
def read_species_tree(path):
    root = parse_newick(first_tree_text(path))
    code2tax = {}
    for x in iter_nodes(root):
        if not x["children"]:
            p = x["name"].split("__")
            if len(p) >= 3 and p[-1].isdigit():
                code2tax[p[0]] = p[-1]
                x["name"] = p[-1]
            else:
                x["name"] = ""
    return root, code2tax


# Lower-case keys under which a FASTA header can be found: first word, its |-separated fields, and the first word
# without its last _token.
#   h   FASTA header (without '>')
# Returns the list of keys.
def header_keys(h):
    first = h.split()[0] if h.split() else h
    keys = [first]
    if "|" in first:
        keys += [x for x in first.split("|")[1:] if x]
    t = first.split("_")
    if len(t) >= 3:
        keys.append("_".join(t[:-1]))
    return [k.lower() for k in keys]


# Index of the sequences by header_keys() (the first sequence wins).
#   seqs   list of (header, sequence) of read_fasta()
# Returns {key: sequence}.
def build_seq_index(seqs):
    idx = {}
    for h, s in seqs:
        if s:
            for k in header_keys(h):
                idx.setdefault(k, s)
    return idx


# Index of the sequences by (gene name, species code), from UniProt headers 'db|ACC|ENTRY_SPECIES ... GN=name ...'.
#   seqs   list of (header, sequence) of read_fasta()
# Returns {(lower-case gene name, species code): [sequences]}.
def build_gene_name_index(seqs):
    out = {}
    for h, s in seqs:
        m = re.search(r"GN=(\S+)", h)
        if m and s:
            code = h.split()[0].split("|")[-1].split("_")[-1]
            out.setdefault((m.group(1).lower(), code), []).append(s)
    return out


# Sequence of a leaf, looked up by full name, name without the id, or id alone.
#   idx    index of build_seq_index()
#   leaf   leaf name GENE_SPECIES_id
# Returns the sequence, or None.
def find_seq(idx, leaf):
    t = leaf.split("_")
    for k in (leaf, "_".join(t[:-1]) if len(t) >= 3 else None, t[-1]):
        if k and k.lower() in idx:
            return idx[k.lower()]
    return None


# Fallback for leaves whose id is not in the FASTA (e.g. an Ensembl id): the sequence with the same GN= gene name and
# species code, only if exactly one sequence matches.
#   gn_idx   index of build_gene_name_index()
#   leaf     leaf name GENE_SPECIES_id
# Returns the sequence, or None.
def find_seq_by_gene_name(gn_idx, leaf):
    t = leaf.split("_")
    if len(t) < 2:
        return None
    hits = gn_idx.get((t[0].lower(), t[1]), [])
    return hits[0] if len(hits) == 1 else None


# Event of an internal node from its NHX tags.
#   tags   {key: value} NHX tags of the node (Ev or D is used)
# Returns 'dup', 'spec' or 'hgt'.
def dup_value(tags):
    e = (tags.get("Ev") or tags.get("D") or "").lower()
    if "egt" in e or "hgt" in e:
        return "hgt"
    if e in ("y", "t", "yes", "true") or "dup" in e:
        return "dup"
    return "spec"


# Align sequences with MAFFT --auto.
#   mafft   MAFFT executable
#   seqs    list of (name, sequence) (names are not passed to MAFFT)
# Returns the list of aligned sequences (upper case), in the same order.
def align(mafft, seqs):
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "in.fa")
        with open(inp, "w") as f:
            for k, (_, s) in enumerate(seqs):
                f.write(f">s{k}\n{s}\n")
        r = subprocess.run([mafft, "--auto", "--quiet", inp], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("mafft failed: " + r.stderr[:300])
    out, name, buf = {}, None, []
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return [out[f"s{k}"].upper() for k in range(len(seqs))]


# Ensembl-like JSON of a pruned gene tree (missing branch lengths become 1.0).
#   x      root node of the dict tree (internal nodes with x["event"] = 'dup' / 'spec')
#   info   {leaf: (taxon id, raw sequence)}
#   aln    {leaf: aligned sequence}
# Returns the JSON node (dict).
def to_json(x, info, aln):
    if not x["children"]:
        g = x["name"]
        return {"id": {"accession": g}, "taxonomy": {"id": int(info[g][0])},
                "branch_length": x["bl"] if x["bl"] is not None else 1.0,
                "sequence": {"mol_seq": {"seq": aln[g]}}}
    return {"events": {"type": "duplication" if x.get("event") == "dup" else "speciation"},
            "branch_length": x["bl"] if x["bl"] is not None else 1.0,
            "children": [to_json(c, info, aln) for c in x["children"]]}


# Convert one family: write trees/<fam>.json, seqs/<fam>.fa and seqs/<fam>.aln.fa in a.out.
#   fam        family name (ST001, ...)
#   d          folder of the family in the release
#   code2tax   {species code: taxon id} of read_species_tree()
#   a          parsed command-line options (out, mafft, match_gene_name, min_length_frac)
# Returns the set of taxa of the family, or None (with a message) if it is skipped.
def convert_family(fam, d, code2tax, a):
    # Input files: the curated tree and the first sequence file found (sequences.fst, else sequences.txt).
    tree_file = os.path.join(d, "consensus_tree.nhx")
    seq_files = [p for p in (os.path.join(d, "sequences.fst"), os.path.join(d, "sequences.txt")) if os.path.exists(p)]
    if not os.path.exists(tree_file) or not seq_files:
        print(f"{fam}: skipped (no consensus tree or sequence file)", file=sys.stderr)
        return None

    # Parse the curated tree and label every internal node 'dup' / 'spec' / 'hgt' from its NHX tags.
    root = parse_newick(first_tree_text(tree_file))
    nodes = list(iter_nodes(root))
    leaves = [x for x in nodes if not x["children"]]
    for x in nodes:
        if x["children"]:
            x["event"] = dup_value(x["tags"])

    # Families with a horizontal / endosymbiotic transfer are not gene trees in the sense of Algorithm 3: skip.
    if any(x.get("event") == "hgt" for x in nodes):
        print(f"{fam}: skipped (horizontal / endosymbiotic gene transfer)", file=sys.stderr)
        return None

    # Indexes of the family's sequences: by header keys, and (with --match-gene-name) by gene name + species code.
    seqs = read_fasta(seq_files[0], letters_only=True)
    idx = build_seq_index(seqs)
    gn_idx = build_gene_name_index(seqs) if a.match_gene_name else None

    # Sequence of a leaf (with the gene-name fallback if requested).
    #   name   leaf name
    # Returns the sequence, or None.
    def seq_of(name):
        s = find_seq(idx, name)
        if s is None and gn_idx is not None:
            s = find_seq_by_gene_name(gn_idx, name)
        return s

    # Leaves kept: a species code known in the species tree and a sequence; duplicate leaf names are kept once.
    # info[leaf] = (taxon id, raw sequence).
    keep, info = set(), {}
    for x in leaves:
        code = next((c for c in x["name"].split("_")[1:] if c in code2tax), None)
        seq = seq_of(x["name"])
        if code is None or seq is None or x["name"] in keep:
            continue
        keep.add(x["name"])
        info[x["name"]] = (code2tax[code], seq)

    # Skip families where too few leaves have a sequence (< MIN_SEQ_FRACTION) or too few are kept
    # (< MIN_KEPT_FRACTION): the pruned tree would not represent the curated family.
    n_seq = sum(1 for x in leaves if seq_of(x["name"]) is not None)
    if not leaves or n_seq / len(leaves) < MIN_SEQ_FRACTION:
        print(f"{fam}: skipped (only {n_seq}/{len(leaves)} leaves have a sequence)", file=sys.stderr)
        return None
    if len(keep) < MIN_KEPT_FRACTION * len(leaves):
        print(f"{fam}: skipped (only {len(keep)}/{len(leaves)} leaves left)", file=sys.stderr)
        return None

    # With --min-length-frac F: drop the fragments, sequences shorter than F x the family median length.
    if a.min_length_frac > 0 and keep:
        lens = sorted(len(info[g][1]) for g in keep)
        med = lens[(len(lens) - 1) // 2]
        frag = {g for g in keep if len(info[g][1]) < a.min_length_frac * med}
        if frag:
            print(f"{fam}: dropping {len(frag)} fragments (< {a.min_length_frac} x median {med}): "
                  + ", ".join(sorted(frag)), file=sys.stderr)
            keep -= frag
            for g in frag:
                info.pop(g)

    # Prune the tree to the kept leaves (unary nodes collapsed); skip it if too few genes or species remain.
    root = prune(root, keep)
    if root is None or len(keep) < MIN_GENES or len({v[0] for v in info.values()}) < MIN_SPECIES:
        print(f"{fam}: skipped (too few genes / species)", file=sys.stderr)
        return None

    # Align the raw sequences with MAFFT, in the leaf order of iter_nodes().
    order = [x["name"] for x in iter_nodes(root) if not x["children"]]
    raw = [info[g][1] for g in order]
    aln = align(a.mafft, list(zip(order, raw)))

    # Write seqs/<fam>.fa (raw, for OrthoFinder), seqs/<fam>.aln.fa (aligned) and trees/<fam>.json (labelled tree
    # with the aligned sequences in the leaves).
    with open(os.path.join(a.out, "seqs", fam + ".fa"), "w") as f:
        f.writelines(f">{g}\n{s}\n" for g, s in zip(order, raw))
    with open(os.path.join(a.out, "seqs", fam + ".aln.fa"), "w") as f:
        f.writelines(f">{g}\n{s}\n" for g, s in zip(order, aln))
    tree = to_json(root, info, dict(zip(order, aln)))
    json.dump({"tree": {"id": fam, "tree": tree}}, open(os.path.join(a.out, "trees", fam + ".json"), "w"))

    # Report and return the taxa of the family (main() collects them to prune the species tree).
    taxa = {v[0] for v in info.values()}
    print(f"{fam}: ok ({len(keep)}/{len(leaves)} genes, {len(taxa)} species)", file=sys.stderr)
    return taxa


# Convert all the families of the release and write the species tree of the taxa used.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True, help="SwissTree release folder (speciestree.nhx, ST001/, ...)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--match-gene-name", action="store_true",
                    help="find the sequence of a leaf by gene name (GN=) + species when its id is not in the FASTA")
    ap.add_argument("--min-length-frac", type=float, default=0.0,
                    help="drop sequences shorter than this fraction of the family median length (default 0 = keep)")
    a = ap.parse_args()
    if shutil.which(a.mafft) is None:
        sys.exit(f"'{a.mafft}' not found: install MAFFT or give --mafft PATH")
    sp_root, code2tax = read_species_tree(os.path.join(a.raw, "speciestree.nhx"))
    for sub in ("trees", "seqs"):
        os.makedirs(os.path.join(a.out, sub), exist_ok=True)
    used = set()
    fams = sorted(os.path.basename(d) for d in glob.glob(os.path.join(a.raw, "ST*")) if os.path.isdir(d))
    for fam in fams:
        taxa = convert_family(fam, os.path.join(a.raw, fam), code2tax, a)
        if taxa:
            used |= taxa
    with open(os.path.join(a.out, "species_tree.nwk"), "w") as f:
        f.write(write_topology(prune(sp_root, used)) + "\n")
    print(f"species tree with {len(used)} taxa", file=sys.stderr)


if __name__ == "__main__":
    main()
