"""
Functions shared by the scripts of the pipeline: Newick / NHX trees, species-tree paths, the Ensembl-like JSON gene
trees, FASTA and genes.tsv files, alignment distances, precision / recall.

A "dict tree" below is the result of parse_newick(): every node is a dict
  name      leaf name or internal label ('' if none)
  bl        branch length above the node (float, or None if absent)
  tags      NHX tags {key: value} of the node ({} if none)
  children  list of child nodes ([] for a leaf)
"""
import json
import re
import sys

import numpy as np

sys.setrecursionlimit(1000000)


# =============================================================================================== Newick / NHX trees

# Read the first tree of a Newick / NHX file (everything up to the first ';').
#   path   Newick / NHX file
# Returns the text of the first tree, with its ';' (or the first line if there is no ';').
def first_tree_text(path):
    txt = open(path, errors="replace").read().replace("\r", "")
    i = txt.find(";")
    return txt[:i + 1] if i >= 0 else txt.strip().split("\n")[0]


# Parse a Newick / NHX string into a dict tree (see the module description). Internal nodes may carry a label
# (e.g. Spec / Dup in the reconstructed trees), NHX tags [&&NHX:key=value:...] are stored in 'tags', other [...]
# comments are ignored, and quotes around names are removed.
#   text   Newick / NHX string of one tree
# Returns the root node.
def parse_newick(text):
    text = text.strip()
    pos, n = 0, len(text)
    stack, root = [], None

    # Read the label, branch length and [comment] that follow a node, starting at the current position `pos` of
    # `text` (pos is advanced past them).
    #   node   node dict whose name / bl / tags are filled in
    def read_tail(node):
        nonlocal pos
        j = pos
        while j < n and text[j] not in "[:,();":
            j += 1
        node["name"] = text[pos:j].strip().strip("'\"")
        pos = j
        while pos < n and text[pos] in "[:":
            if text[pos] == "[":
                k = text.index("]", pos)
                body = text[pos + 1:k]
                if body.startswith("&&NHX"):
                    for kv in body.split(":")[1:]:
                        if "=" in kv:
                            key, val = kv.split("=", 1)
                            node["tags"][key.strip()] = val.strip()
                pos = k + 1
            else:
                j = pos + 1
                while j < n and text[j] not in "[:,();":
                    j += 1
                try:
                    node["bl"] = float(text[pos + 1:j])
                except ValueError:
                    pass
                pos = j

    while pos < n:
        c = text[pos]
        if c == "(":
            node = dict(name="", bl=None, tags={}, children=[])
            if stack:
                stack[-1]["children"].append(node)
            stack.append(node)
            pos += 1
        elif c == ",":
            pos += 1
        elif c == ")":
            node = stack.pop()
            pos += 1
            read_tail(node)
            root = node
        elif c == ";":
            break
        elif c.isspace():
            pos += 1
        else:
            leaf = dict(name="", bl=None, tags={}, children=[])
            read_tail(leaf)
            if stack:
                stack[-1]["children"].append(leaf)
            else:
                root = leaf
    return root


# Iterate over all nodes of a dict tree (root first, then the children in reverse order, depth first). This order is
# also the order in which leaves are written to the sequence files, so it is kept fixed.
#   root   root node of the dict tree
# Yields the nodes.
def iter_nodes(root):
    stack = [root]
    while stack:
        x = stack.pop()
        yield x
        stack.extend(x["children"])


# Leaves of a dict tree, left to right.
#   node   root node of the (sub)tree
# Returns the list of leaf nodes.
def tree_leaves(node):
    if not node["children"]:
        return [node]
    out = []
    for c in node["children"]:
        out += tree_leaves(c)
    return out


# Copy of a dict tree restricted to some leaves; unary nodes are collapsed (the child survives, with the branch
# lengths added). The input tree is not modified.
#   node   root node of the (sub)tree
#   keep   set of the leaf names to keep
# Returns the root of the pruned copy, or None if no leaf is kept.
def prune(node, keep):
    if not node["children"]:
        return dict(node) if node["name"] in keep else None
    ch = [c for c in (prune(c, keep) for c in node["children"]) if c is not None]
    if not ch:
        return None
    if len(ch) == 1:
        c = ch[0]
        if c["bl"] is not None or node["bl"] is not None:
            c["bl"] = (c["bl"] or 0.0) + (node["bl"] or 0.0)
        return c
    out = dict(node)
    out["children"] = ch
    return out


# Newick topology of a dict tree (leaf names only, no internal labels or lengths).
#   node   root node of the dict tree
# Returns the Newick string, with the final ';'.
def write_topology(node):
    # Recursive part, without the ';'.
    #   x   current node
    def rec(x):
        if not x["children"]:
            return x["name"]
        return "(" + ",".join(rec(c) for c in x["children"]) + ")"
    return rec(node) + ";"


# ======================================================================================================= species tree

# Map every leaf (taxon) of a species tree to its path from the root: the list of child indices taken at each level.
# Two taxa share the first k entries iff their LCA is at depth >= k; entry k is the child of the LCA on the way down.
#   root   root node of the species tree (dict tree, leaves named by taxon)
# Returns {taxon: path (list of ints)}.
def species_paths(root):
    paths = {}

    # Depth-first walk carrying the current path.
    #   x      current node
    #   path   path from the root to x
    def rec(x, path):
        if not x["children"]:
            paths[x["name"]] = path
        for i, c in enumerate(x["children"]):
            rec(c, path + [i])
    rec(root, [])
    return paths


# Read a species-tree Newick file and return species_paths() of it.
#   path   Newick file of the species tree
# Returns {taxon: path}.
def read_species_paths(path):
    return species_paths(parse_newick(first_tree_text(path)))


# Length of the common prefix of several paths (= depth of the LCA of the corresponding taxa).
#   paths   list of species_paths() paths (at least one)
# Returns the length (int).
def common_prefix_len(paths):
    k = 0
    for cols in zip(*paths):
        if len(set(cols)) != 1:
            break
        k += 1
    return k


# =================================================================================================== JSON gene trees

# Find the tree id and root node of an Ensembl-like JSON gene tree (any of the usual wrappings: {"tree": {...}},
# {"tree": {"tree": {...}}}, with an "id" at any level).
#   d   parsed JSON object
# Returns (tree id or None, root node), or (None, None) if d is not a gene tree.
def get_root(d):
    if not isinstance(d, dict) or not isinstance(d.get("tree"), dict):
        return None, None
    tid = d["id"] if isinstance(d.get("id"), str) else None
    t = d["tree"]
    while "children" not in t and isinstance(t.get("tree"), dict):
        if tid is None and isinstance(t.get("id"), str):
            tid = t["id"]
        t = t["tree"]
    if tid is None and isinstance(t.get("id"), str):
        tid = t["id"]
    return tid, t


# Load a JSON gene tree file.
#   path   JSON file
# Returns (tree id or None, root node), as get_root().
def load_json_tree(path):
    return get_root(json.load(open(path)))


# Event type of an internal node of a JSON gene tree, lower case ("" if none).
#   node   internal node (JSON dict, with "events": {"type": ...})
# Returns the event type string (e.g. "speciation", "duplication", "dubious", "gene_split").
def event_type(node):
    ev = node.get("events")
    return ((ev.get("type") if isinstance(ev, dict) else None) or "").lower()


# Whether an event type is an Ensembl "dubious" node: a duplication node with a low duplication-consistency score,
# that Ensembl Compara marks as neither a sure duplication nor a speciation ("dubious", or the older
# "speciation_or_duplication").
#   ev   event type of event_type()
# Returns True or False.
def is_dubious(ev):
    return "dubious" in ev or ev == "speciation_or_duplication"


# Event of an internal node of a JSON gene tree: a plain speciation is Spec, a dubious node (Ensembl) is labelled as
# chosen, anything else (duplication, gene_split, no event) is Dup.
#   node      internal node (JSON dict, with "events": {"type": ...})
#   dubious   "dup" (default) or "spec": the label of the dubious nodes
# Returns "Spec" or "Dup".
def event_label(node, dubious="dup"):
    ev = event_type(node)
    if is_dubious(ev):
        return "Spec" if dubious == "spec" else "Dup"
    return "Spec" if ("speciation" in ev and "duplication" not in ev) else "Dup"


# Gene id of a JSON leaf: id.accession, or id if it is a string, or sequence.id[0].accession.
#   node   leaf node (JSON dict)
# Returns the gene id, or None.
def leaf_gene(node):
    i = node.get("id")
    if isinstance(i, dict) and i.get("accession"):
        return i["accession"]
    if isinstance(i, str):
        return i
    ids = (node.get("sequence") or {}).get("id")
    if isinstance(ids, list) and ids and isinstance(ids[0], dict):
        return ids[0].get("accession")
    return None


# Taxon id of a JSON leaf (taxonomy.id).
#   node   leaf node (JSON dict)
# Returns the taxon id as a string.
def leaf_taxon(node):
    return str(int(node["taxonomy"]["id"]))


# Aligned sequence of a JSON leaf (sequence.seq or sequence.mol_seq.seq).
#   node   leaf node (JSON dict)
# Returns the sequence, or None.
def leaf_seq(node):
    s = node.get("sequence")
    if not isinstance(s, dict):
        return None
    if isinstance(s.get("seq"), str):
        return s["seq"]
    m = s.get("mol_seq")
    if isinstance(m, dict) and isinstance(m.get("seq"), str):
        return m["seq"]
    return None


# Leaves of a JSON gene tree, in the order of a depth-first walk with a stack (last child first).
#   root   root node (JSON dict)
# Returns the list of leaf nodes.
def json_leaves(root):
    stack, out = [root], []
    while stack:
        x = stack.pop()
        ch = x.get("children")
        if ch:
            stack.extend(ch)
        else:
            out.append(x)
    return out


# Aligned family of a JSON gene tree file.
#   path   JSON gene tree file with aligned sequences in the leaves
# Returns a list of (gene, taxon, aligned sequence), or None if the leaves do not all carry a gene id and sequences of
# the same length.
def load_alignment(path):
    _, root = load_json_tree(path)
    if root is None:
        return None
    fam = [(leaf_gene(x), leaf_taxon(x), leaf_seq(x)) for x in json_leaves(root)]
    if not fam or not all(s and g for g, _, s in fam) or len({len(s) for _, _, s in fam}) != 1:
        return None
    return fam


# ============================================================================================================ files

# Read a FASTA file.
#   path           FASTA file
#   letters_only   if True, every character of the sequence lines that is not a letter is dropped (raw sequences);
#                  otherwise only white space is removed (alignments keep '-')
# Returns a list of (header line without '>', sequence).
def read_fasta(path, letters_only=False):
    seqs, name, buf = [], None, []
    for line in open(path, errors="replace"):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif name is not None:
            buf.append(re.sub(r"[^A-Za-z]", "", line) if letters_only else line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


# Read a three-column TSV file such as genes.tsv (tree_id, gene, taxon) or orthologs.tsv (tree_id, gene1, gene2);
# a header line starting with "tree_id" is skipped.
#   path   TSV file
# Returns a list of (column 1, column 2, column 3).
def read_genes(path):
    out = []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0] != "tree_id":
            out.append((p[0], p[1], p[2]))
    return out


# ========================================================================================================= distances

# Characters that are not residues in an alignment.
INVALID = np.frombuffer(b"-X*.?", dtype=np.uint8)


# Turn aligned sequences into matrices.
#   seqs   list of aligned sequences (all of the same length)
# Returns (A, V): A[k] = bytes of sequence k (uint8 matrix), V[k] = True where sequence k has a residue.
def alignment_matrix(seqs):
    A = np.empty((len(seqs), len(seqs[0])), dtype=np.uint8)
    for k, s in enumerate(seqs):
        A[k] = np.frombuffer(s.upper().encode(), dtype=np.uint8)
    return A, ~np.isin(A, INVALID)


# Distances from sequence i to all sequences. p = 1 - identical residues / columns where at least one of the two has a
# residue.
#   A, V          matrices of alignment_matrix()
#   i             index of the sequence
#   model         "p" (p-distance), "poisson" (-ln(1-p)), "kimura" (-ln(1-p-0.2p^2)) for proteins, "jc" (Jukes-Cantor)
#                 for DNA
#   min_overlap   pairs with fewer columns where at least one of the two has a residue get +inf
# Returns a numpy array of distances (+inf for i itself, too short overlaps and saturated pairs).
def distance_row(A, V, i, model, min_overlap):
    both = V[i] | V
    overlap = both.sum(axis=1).astype(float)
    match = ((A[i] == A) & both).sum(axis=1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 - match / overlap
        if model == "p":
            d = p
        elif model == "poisson":
            d = -np.log(1.0 - p)
        elif model == "jc":
            d = -0.75 * np.log(1.0 - 4.0 * p / 3.0)
        else:
            d = -np.log(1.0 - p - 0.2 * p * p)
    d[(overlap < min_overlap) | ~np.isfinite(d)] = np.inf
    d[i] = np.inf
    return d


# Distance function over the aligned family of a JSON gene tree file (rows computed on demand and cached).
#   path          JSON gene tree file with aligned sequences in the leaves
#   model         distance model, as distance_row()
#   min_overlap   as distance_row()
# Returns a function dist(g, h) of two gene ids, or None if the file has no usable alignment.
def family_distances(path, model, min_overlap):
    fam = load_alignment(path)
    if not fam:
        return None
    index = {g: k for k, (g, _, _) in enumerate(fam)}
    A, V = alignment_matrix([s for _, _, s in fam])
    rows = {}

    # Distance between two genes.
    #   g, h   gene ids
    # Returns the distance (float; nan if one of them is not in the alignment).
    def dist(g, h):
        if g not in index or h not in index:
            return float("nan")
        if g not in rows:
            rows[g] = distance_row(A, V, index[g], model, min_overlap)
        return float(rows[g][index[h]])
    return dist


# ============================================================================================================ scores

# Precision, recall and F1 (1.0 when nothing is expected / predicted).
#   n_true   number of true items
#   n_pred   number of predicted items
#   tp       number of predicted items that are true
# Returns (precision, recall, F1).
def prf(n_true, n_pred, tp):
    p = tp / n_pred if n_pred else 1.0
    r = tp / n_true if n_true else 1.0
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)
