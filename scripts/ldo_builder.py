#!/usr/bin/env python3
"""
Reconstruct gene trees from LDO graphs (Algorithm 3), with optional repair of spurious arcs. Deterministic.

  python ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
      --repair none|count|sparse|count+sparse|count+sparse+mutual --trees-out trees.txt \
      [--count-age 1.0 --seqs-trees TREES_DIR --model poisson|jc]

Input files
  --graph FILE      the LDO graphs of all the families (ldo_edges.tsv of predict_ldo.py). Tab-separated, one arc per
                    line; the FIRST LINE IS ALWAYS SKIPPED (header); lines with fewer than 6 columns are ignored.
                      1 tree_id        family id; the families are reconstructed independently of each other
                      2 source         gene g
                      3 source_taxon   taxon id of g
                      4 target         gene h: the gene of h's species closest to g (the LDO of g in that species)
                      5 target_taxon   taxon id of h
                      6 distance       distance g -> h (>= 0); used only by the repairs
                      7 tie            optional, ignored
                    Expected: arcs only between genes of the same family and of DIFFERENT species, at most one arc
                    from a gene to each other species, gene ids unique within a family.
  --genes FILE      optional, tree_id / gene / taxon (tab-separated, header line optional): adds the genes that have
                    no arc (they become separate leaves under a duplication).
  --species-tree    rooted Newick, leaves = the taxon ids used above (may be multifurcating); a family with a taxon
                    that is not a leaf of it is skipped with a message.
  --seqs-trees DIR  only with --count-age: DIR/<tree_id>.json, the gene tree JSON with the aligned sequences in the
                    leaves (leaf id.accession = gene id), as written by swisstree_convert.py / simphy_prepare.py.

Output (--trees-out): one line "tree_id<TAB>newick;" per family, internal nodes labelled Spec / Dup.

Algorithm 3: a disconnected graph is a duplication of its components; a connected graph is a speciation whose children
are the parts given by the children of x = lca_S(species of the graph). Before a connected graph is treated as a
speciation, heuristic fixes (repairs) are attempted, in this order; the first one that splits the graph makes it a
duplication of the groups it returns, and each group is then reconstructed recursively.

Definitions. The graph has an arc g -> h when h is the gene of h's species closest to g (the predicted LDO). The
speciation would separate the genes into parts P1, P2, ... (one per child of x; two or more). A gene h is a "strong
target" if at least --min-support genes from parts other than h's part have an arc to h; these genes are its
"choosers". The "strong target species" is the species with the most strong targets, over all the species of the
graph (ties: the first taxon id in sorted order); it must have at least two.

 count   "several targets in one species". In a true speciation, all the genes of one part are equally related to the
         genes of a species s of another part, so they should all choose the same gene of s. The count repair applies
         if a strong target species s exists, and uses only that species (no other species is tried):
         1. grouping: every strong target of s starts a group, together with its choosers; the other genes are added
            by label propagation (in sorted order, a gene takes the group most frequent among its already grouped
            neighbours through out- and in-arcs, ties to the first group; up to 10 rounds; genes never reached join
            the first group). At least two groups are needed.
         2. ratio test: over the arcs with both ends in the graph, the split is accepted if
            median(distance of the arcs between two different groups) >= --ratio x median(distance of the arcs inside
            a group) (default ratio 2: arcs between paralog groups span an older divergence).
         3. age test, only with --count-age G and only if the ratio test failed: the split is also accepted if
            median(alignment distance between the strong targets of s, over all pairs) >= G x median(distance of the
            arcs from the choosers to these targets, one per chooser). Under a true speciation the targets would
            differ only by duplications after it and be closer to each other than to their choosers; a ratio >= 1
            means their duplication is at least as old as the speciation. Needs --seqs-trees.
         If neither test accepts the split, the count repair does nothing.

 sparse  "a few arcs hold the graph together", tried only if the count repair did not split the graph. The graph is
         taken as undirected (an edge between g and h if there is an arc in either direction; weight 1).
         1. cut: with deg(g) = number of neighbours, the genes are sorted by v(g) / sqrt(deg(g)), where v is the
            eigenvector of the second largest eigenvalue of the normalised adjacency matrix D^-1/2 A D^-1/2. Every
            prefix of this order with >= --sparse-min-side genes on each side is a candidate cut; the one with the
            lowest conductance = (edges cut) / min(sum of deg on each side) is kept. It is a candidate duplication
            if its conductance <= --sparse-conductance (default 0.3).
         2. refinement: in sorted order, a gene moves to the other side if it is attached to it by more weight than
            to its own side, a mutual arc (g -> h and h -> g) counting --mutual-weight (default 2) and a one-way arc
            1; a side never drops below --sparse-min-side genes; repeated until nothing moves (at most 20 rounds).
         3. bridge check: a gene is a bridge if it has >= --sparse-bridge-mutual (default 2) mutual arcs to each
            side and median(distance of its mutual arcs to the farther side) <= --sparse-bridge-ratio (default 1.5)
            x median(distance to the nearer side), the distance of a mutual arc being the mean of its two directions.
            Such a gene looks like an ortholog of both sides (a speciation lies above the duplication); the cut is
            rejected if there are more than --sparse-max-bridge bridges (default 0; -1 = no check).
         4. ratio test: the cut is accepted if median(distance of the arcs crossing the cut) >= --sparse-ratio
            (default 1.5) x median(distance of the other arcs of the graph).
         Graphs with more than 3000 genes are not searched for a sparse cut.

 mutual  "the arcs are mostly one-way", tried only if neither count nor sparse split the graph (opt-in: --repair
         count+sparse+mutual). An arc g -> h is reciprocal if h -> g is an arc too. In a family made of many small
         ortholog groups (a long series of old duplications with gene losses), most arcs are forced: a gene whose own
         group has no gene of species s must choose a paralog of s, and such arcs are rarely reciprocal and longer than
         the arcs between orthologs.
         1. trigger: fewer than --mutual-share (default 0.5) of the arcs of the graph are reciprocal.
         2. groups: the connected components of the reciprocal arcs; then, repeatedly in sorted order, a gene without
            group joins the group of its closest target (shortest out-arc) that has one; genes never reached form
            their own group. At least two groups are needed.
         3. duplication location: the duplicated species (those with genes in two or more groups) must have the
            same LCA in the species tree as all the species of the graph. Otherwise the duplication lies lower in
            the species tree, below a speciation (e.g. vertebrate paralogs with their invertebrate ortholog); the
            graph is then left to the speciation, and the repair may apply again to the part that holds the
            duplicated species. No duplicated species: no split.
         4. ratio test: the split is accepted if median(distance of the arcs between two different groups) >=
            --mutual-ratio (default 1.5) x median(distance of the arcs inside a group).
"""
import argparse
import collections
import math
import os
import statistics
import sys

import numpy as np

from ldo_utils import common_prefix_len, family_distances, read_species_paths

SPARSE_MAX_SIZE = 3000          # graphs with more genes are not searched for a sparse cut
REFINE_ROUNDS = 20              # rounds of the side refinement after a sparse cut
PROPAGATE_ROUNDS = 10           # rounds of label propagation in the count repair
REPAIRS = ("count", "sparse", "mutual")   # the repairs, in the order they are tried


# ------------------------------------------------------------------------------------------------------ graph input

# Read the LDO graphs of all the families (see "Input files" in the module description).
#   path         --graph file (tree_id, source, source taxon, target, target taxon, distance)
#   genes_path   --genes file (tree_id, gene, taxon), or None; genes listed only there (no arcs) are added to tax
# Returns {tree_id: dict(arcs={(source, target): distance}, tax={gene: taxon})}.
def load_graph(path, genes_path=None):
    fam = collections.defaultdict(lambda: dict(arcs={}, tax={}))
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            d = fam[p[0]]
            d["arcs"][(p[1], p[3])] = float(p[5])
            d["tax"][p[1]] = p[2]
            d["tax"][p[3]] = p[4]
    if genes_path:
        with open(genes_path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3 and p[0] != "tree_id":
                    fam[p[0]]["tax"].setdefault(p[1], p[2])
    return fam


# Kinds of repair of a --repair string.
#   repair   "none", or "count", "sparse", "mutual" joined by + (e.g. "count+sparse")
# Returns the list of repair kinds (empty for "none"); exits with a message on an unknown kind.
def repair_kinds(repair):
    kinds = [p for p in repair.split("+") if p and p != "none"]
    if any(p not in REPAIRS for p in kinds):
        sys.exit(f"unknown repair '{repair}'; use none or count, sparse, mutual joined by + (e.g. count+sparse)")
    return kinds


# -------------------------------------------------------------------------------------------------- Algorithm 3

# Connected components of a gene set in the undirected version of the graph, in a fixed order.
#   V     set of genes
#   out   arcs, out[g] = {h: distance} (arcs leaving V are ignored)
# Returns a list of sets of genes (ordered by their smallest gene).
def components(V, out):
    adj = collections.defaultdict(set)
    for g in V:
        for h in out.get(g, ()):
            if h in V:
                adj[g].add(h)
                adj[h].add(g)
    seen, comps = set(), []
    for g in sorted(V):
        if g in seen:
            continue
        stack, comp = [g], {g}
        seen.add(g)
        while stack:
            u = stack.pop()
            for w in sorted(adj[u]):
                if w not in seen:
                    seen.add(w)
                    comp.add(w)
                    stack.append(w)
        comps.append(comp)
    return comps


# Newick string of a nested tuple tree (internal nodes labelled Spec / Dup).
#   t   tuple tree of Builder.build()
# Returns the Newick string, without the final ';'.
def newick(t):
    if t[0] == "leaf":
        return t[1]
    return "(" + ",".join(newick(c) for c in t[1]) + ")" + t[0]


# Algorithm 3 on the LDO graph of one family, with the repairs (count, sparse) described in the module description.
# Member variables:
#   tid      tree id (family name); the alignment for --count-age is read from <opts.seqs_trees>/<tid>.json
#   tax      {gene: taxon} of all the genes of the family
#   sp       {taxon: path} of species_paths() of the species tree
#   kinds    repairs to try, in this order: a sublist of ["count", "sparse", "mutual"] (empty = plain Algorithm 3)
#   opts     parsed command-line options: the parameters of the repairs (min_support, ratio, count_age, model,
#            min_overlap, seqs_trees, sparse_*, mutual_weight, mutual_share, mutual_ratio)
#   out      out[g] = {h: distance}: the arcs g -> h of the graph, in sorted order
#   inn      inn[h] = [g, ...]: the genes with an arc to h, in sorted order
#   _dist    distance function of the family for the age test (see distances()), or None
#   _dist_loaded   True once _dist has been computed (it is computed only if the age test is reached)
class Builder:

    # Store the family and its arcs.
    #   tid      tree id of the family
    #   d        the family as returned by load_graph(): dict(arcs={(g, h): distance}, tax={gene: taxon})
    #   sp       {taxon: path} of the species tree
    #   repair   the --repair string ("none", "count", "count+sparse", "count+sparse+mutual", ...)
    #   opts     parsed command-line options (only needed if repair is not "none")
    def __init__(self, tid, d, sp, repair="none", opts=None):
        self.tid, self.tax, self.sp, self.opts = tid, d["tax"], sp, opts
        self.kinds = repair_kinds(repair)
        self.out = collections.defaultdict(dict)
        for (g, h), x in sorted(d["arcs"].items()):
            self.out[g][h] = x
        self.inn = collections.defaultdict(list)
        for g in sorted(self.out):
            for h in self.out[g]:
                self.inn[h].append(g)
        self._dist, self._dist_loaded = None, False

    # Parts of a gene set given by the children of x = lca_S(species of the genes). Only called on connected sets of
    # two or more genes, which always have two or more species (arcs join genes of different species).
    #   V   set of genes
    # Returns a list of sets of genes, one per child of x that has genes (at least two), in the order of the species
    # tree.
    def parts(self, V):
        k = common_prefix_len([self.sp[t] for t in {self.tax[g] for g in V}])
        groups = collections.defaultdict(set)
        for g in V:
            groups[self.sp[self.tax[g]][k]].add(g)
        return [groups[i] for i in sorted(groups)]

    # Recursive construction of the tree of a gene set (Algorithm 3, with the repairs).
    #   V   genes (any iterable)
    # Returns a nested tuple tree: ('leaf', gene) or ('Spec' | 'Dup', [children]).
    def build(self, V):
        V = set(V)
        if len(V) == 1:
            return ("leaf", next(iter(V)))
        comps = components(V, self.out)
        if len(comps) > 1:
            return ("Dup", [self.build(c) for c in sorted(comps, key=min)])
        parts = self.parts(V)
        if self.kinds:
            groups = self.repair(V, parts)
            if groups is not None:
                return ("Dup", [self.build(set(x)) for x in sorted(groups, key=min)])
        return ("Spec", [self.build(p) for p in parts])

    # Try the repairs of self.kinds in order (count, then sparse) on a connected gene set.
    #   V       set of genes of the (connected) graph
    #   parts   parts of V under the speciation (parts())
    # Returns the groups of a duplication (from the first repair that splits the graph), or None.
    def repair(self, V, parts):
        if "count" in self.kinds:
            groups = self.count_repair(V, parts)
            if groups is not None:
                return groups
        if "sparse" in self.kinds:
            groups = self.sparse_cut_repair(V)
            if groups is not None:
                return groups
        if "mutual" in self.kinds:
            return self.mutual_repair(V)
        return None

    # Distance function of the family's alignment, for the age test of the count repair (computed on first use).
    # Returns dist(g, h) (ldo_utils.family_distances()), or None if the family has no usable alignment.
    def distances(self):
        if not self._dist_loaded:
            o = self.opts
            self._dist = family_distances(os.path.join(o.seqs_trees, self.tid + ".json"), o.model, o.min_overlap)
            self._dist_loaded = True
        return self._dist

    # -------------------------------------------------------------------------------------------- count repair

    # Label propagation: every unlabelled gene of V takes the most frequent label of its labelled neighbours (out and
    # in arcs; ties to the smallest label), repeated until nothing changes (at most PROPAGATE_ROUNDS rounds); genes
    # never reached get label 0.
    #   V       set of genes of the graph
    #   label   {gene: group number} of the already grouped genes; completed in place for every gene of V
    def propagate(self, V, label):
        for _ in range(PROPAGATE_ROUNDS):
            changed = False
            for g in sorted(V):
                if g in label:
                    continue
                votes = collections.Counter()
                for h in self.out.get(g, ()):
                    if h in V and h in label:
                        votes[label[h]] += 1
                for h in self.inn.get(g, ()):
                    if h in V and h in label:
                        votes[label[h]] += 1
                if votes:
                    best = max(votes.values())
                    label[g] = min(k for k, v in votes.items() if v == best)
                    changed = True
            if not changed:
                break
        for g in sorted(V):
            label.setdefault(g, 0)

    # Arcs with both ends in a gene set.
    #   V   set of genes
    # Returns the list of arcs (g, h) (in no particular order: only their medians are used).
    def arcs_in(self, V):
        return [(g, h) for g in V for h in self.out.get(g, ()) if h in V]

    # Count repair (see the module description). Parameters from self.opts: min_support (--min-support), ratio
    # (--ratio), count_age (--count-age, None = no age test).
    #   V       set of genes of the (connected) graph
    #   parts   parts of V under the speciation (parts())
    # Returns the groups of a duplication (list of sets of genes), or None.
    def count_repair(self, V, parts):
        o = self.opts
        # Part of each gene (index in parts).
        pid = {g: i for i, p in enumerate(parts) for g in p}

        # Choosers of every gene: choosers[s][h] = the genes of other parts that have an arc to h (h of species s).
        choosers = collections.defaultdict(lambda: collections.defaultdict(set))
        for g in sorted(V):
            for h in self.out.get(g, ()):
                if h in V and pid[g] != pid[h]:
                    choosers[self.tax[h]][h].add(g)

        # Strong target species: the species with the most strong targets (genes with >= --min-support choosers);
        # it needs at least two, ties go to the first species in sorted order. None: no count repair.
        best = None
        for s in sorted(choosers):
            strong = {t: c for t, c in choosers[s].items() if len(c) >= o.min_support}
            if len(strong) >= 2 and (best is None or len(strong) > len(best[1])):
                best = (s, strong)
        if best is None:
            return None
        s, strong = best

        # Initial groups: each strong target of s starts its own group, and its choosers join it (a chooser of
        # several targets stays in the first one).
        label = {t: i for i, t in enumerate(sorted(strong))}
        for t in sorted(strong):
            for g in sorted(strong[t]):
                label.setdefault(g, label[t])

        # The other genes join a group by label propagation; collect the groups. Fewer than two: no split.
        self.propagate(V, label)
        groups = collections.defaultdict(set)
        for g, i in label.items():
            groups[i].add(g)
        if len(groups) < 2:
            return None

        # Distances of the arcs of the graph, split into arcs between two groups and arcs inside a group.
        arcs = self.arcs_in(V)
        cross = [self.out[g][h] for g, h in arcs if label[g] != label[h]]
        within = [self.out[g][h] for g, h in arcs if label[g] == label[h]]
        if not cross or not within:
            return None

        # Ratio test: accept if the arcs between groups are >= --ratio times longer (median) than inside groups.
        if statistics.median(cross) >= o.ratio * statistics.median(within):
            return list(groups.values())

        # Age test (only with --count-age): accept if the strong targets are >= --count-age times farther apart
        # (median over all pairs, alignment distance) than from their choosers (median of the chooser -> target arcs).
        dist = self.distances() if o.count_age is not None else None
        if dist is not None:
            ts = sorted(strong)
            between = [dist(x, y) for i, x in enumerate(ts) for y in ts[i + 1:]]
            between = [x for x in between if math.isfinite(x)]
            toward = [self.out[g][t] for t in ts for g in strong[t]]
            if between and toward:
                mt = statistics.median(toward)
                if (statistics.median(between) / mt if mt > 0 else float("inf")) >= o.count_age:
                    return list(groups.values())

        # Neither test accepted the split.
        return None

    # ------------------------------------------------------------------------------------------- sparse repair

    # Refine a two-sided cut: a gene moves to the other side when it is attached to it by a larger weight of arcs than
    # to its own side (a mutual arc g <-> h counts --mutual-weight, a one-way arc 1); repeated until nothing moves (at
    # most REFINE_ROUNDS rounds). A gene does not leave a side that has --sparse-min-side genes or fewer.
    #   V      set of genes of the graph
    #   side   genes of one side of the cut (the other side is V - side)
    # Returns the refined side (new set).
    def refine_sides(self, V, side):
        min_side, mutual_weight = self.opts.sparse_min_side, self.opts.mutual_weight
        side = set(side)
        rest = set(V) - side
        for _ in range(REFINE_ROUNDS):
            moved = False
            for g in sorted(V):
                own, other = (side, rest) if g in side else (rest, side)
                if len(own) <= min_side:
                    continue
                w_own = w_other = 0.0
                for h in V:
                    if h == g:
                        continue
                    fwd, back = h in self.out.get(g, ()), g in self.out.get(h, ())
                    if not fwd and not back:
                        continue
                    w = mutual_weight if (fwd and back) else 1.0
                    if h in own:
                        w_own += w
                    else:
                        w_other += w
                if w_other > w_own:
                    own.discard(g)
                    other.add(g)
                    moved = True
            if not moved:
                break
        return side

    # Bridge genes of a cut: genes with >= --sparse-bridge-mutual mutual arcs to each side whose median mutual-arc
    # distance to the farther side is at most --sparse-bridge-ratio times that to the nearer side (orthologs of both
    # families).
    #   V            set of genes of the graph
    #   side, rest   the two sides of the cut
    # Returns the sorted list of bridge genes.
    def bridge_genes(self, V, side, rest):
        min_mutual, ratio = self.opts.sparse_bridge_mutual, self.opts.sparse_bridge_ratio

        # Mutual arcs between a gene and a group, the distance of a mutual arc being the mean of its two directions.
        #   g     gene
        #   grp   set of genes
        # Returns (number of mutual arcs, median distance or None).
        def stats(g, grp):
            d = [(self.out[g][h] + self.out[h][g]) / 2 for h in grp
                 if h != g and h in self.out.get(g, ()) and g in self.out.get(h, ())]
            return len(d), (statistics.median(d) if d else None)
        out = []
        for g in V:
            n1, m1 = stats(g, side)
            n2, m2 = stats(g, rest)
            if n1 >= min_mutual and n2 >= min_mutual:
                lo, hi = min(m1, m2), max(m1, m2)
                if hi <= ratio * max(lo, 1e-9):
                    out.append(g)
        return sorted(out)

    # Sparse repair (see the module description). Parameters from self.opts: sparse_min_side, sparse_conductance,
    # mutual_weight, sparse_max_bridge, sparse_bridge_mutual, sparse_bridge_ratio, sparse_ratio.
    #   V   set of genes of the (connected) graph
    # Returns the two sides of a duplication ([set, set]), or None.
    def sparse_cut_repair(self, V):
        o = self.opts
        # Too large (eigen-decomposition is cubic) or too small for two sides of --sparse-min-side genes: no cut.
        n = len(V)
        if n > SPARSE_MAX_SIZE or n < 2 * o.sparse_min_side:
            return None

        # Undirected graph on the genes (numbered 0..n-1 in sorted order): an edge if there is an arc either way.
        genes = sorted(V)
        idx = {g: i for i, g in enumerate(genes)}
        adj = [set() for _ in range(n)]
        for g in genes:
            for h in self.out.get(g, ()):
                if h in idx and h != g:
                    adj[idx[g]].add(idx[h])
                    adj[idx[h]].add(idx[g])

        # Degrees (number of neighbours); a gene without neighbours cannot be placed by the spectral order.
        deg = np.array([len(x) for x in adj], dtype=float)
        if (deg == 0).any():
            return None

        # Normalised adjacency matrix N = D^-1/2 A D^-1/2.
        dm = 1.0 / np.sqrt(deg)
        N = np.zeros((n, n))
        for i in range(n):
            for j in adj[i]:
                N[i, j] = dm[i] * dm[j]

        # Spectral order: genes sorted by v(g) / sqrt(deg(g)), v = eigenvector of the second largest eigenvalue
        # (eigh returns the eigenvalues in increasing order, so it is the column -2).
        _, vec = np.linalg.eigh(N)
        order = np.argsort(vec[:, -2] * dm, kind="stable")

        # Sweep: S grows one gene at a time along the order; keep the prefix with the lowest conductance
        # cut(S) / min(vol(S), vol(rest)), among those with >= --sparse-min-side genes on each side.
        # When gene i joins S, its edges to S stop being cut and its other edges become cut.
        vol = deg.sum()
        in_s = [False] * n
        cut = vol_s = 0.0
        best = None
        for k in range(n - 1):
            i = int(order[k])
            to_in = sum(1 for j in adj[i] if in_s[j])
            cut += (len(adj[i]) - to_in) - to_in
            vol_s += deg[i]
            in_s[i] = True
            if k + 1 >= o.sparse_min_side and n - (k + 1) >= o.sparse_min_side:
                phi = cut / min(vol_s, vol - vol_s)
                if best is None or phi < best[0] - 1e-12:
                    best = (phi, k + 1)

        # No candidate, or the best cut is not sparse enough (conductance > --sparse-conductance): no cut.
        if best is None or best[0] > o.sparse_conductance:
            return None

        # The two sides of the best cut, after moving genes that are more attached to the other side.
        side = {genes[int(order[k])] for k in range(best[1])}
        side = self.refine_sides(V, side)
        rest = set(genes) - side

        # Bridge check (unless --sparse-max-bridge -1): reject the cut if too many genes look like orthologs of
        # both sides.
        if o.sparse_max_bridge >= 0:
            if len(self.bridge_genes(V, side, rest)) > o.sparse_max_bridge:
                return None

        # Ratio test (unless --sparse-ratio 0): the arcs crossing the cut must be >= --sparse-ratio times longer
        # (median) than the other arcs of the graph.
        if o.sparse_ratio > 0:
            cut_d = [x for g in genes for h, x in self.out.get(g, {}).items()
                     if h in idx and ((g in side) != (h in side))]
            in_d = [x for g in genes for h, x in self.out.get(g, {}).items()
                    if h in idx and ((g in side) == (h in side))]
            if not cut_d or not in_d or statistics.median(cut_d) < o.sparse_ratio * statistics.median(in_d):
                return None

        # All checks passed: duplication of the two sides.
        return [side, rest]


    # ------------------------------------------------------------------------------------------- mutual repair

    # Mutual repair (see the module description). Parameters from self.opts: mutual_share (--mutual-share),
    # mutual_ratio (--mutual-ratio).
    #   V   set of genes of the (connected) graph
    # Returns the groups of a duplication (list of sets of genes), or None.
    def mutual_repair(self, V):
        o = self.opts
        # Arcs of the graph and the reciprocal ones; trigger only if few arcs are reciprocal.
        arcs = self.arcs_in(V)
        mutual = [(g, h) for g, h in arcs if g in self.out.get(h, ())]
        if not arcs or len(mutual) >= o.mutual_share * len(arcs):
            return None

        # Groups: connected components of the reciprocal arcs (each labelled by its smallest gene).
        adj = collections.defaultdict(set)
        for g, h in mutual:
            adj[g].add(h)
            adj[h].add(g)
        label = {}
        for g in sorted(V):
            if g in label or g not in adj:
                continue
            label[g], stack = g, [g]
            while stack:
                u = stack.pop()
                for w in sorted(adj[u]):
                    if w not in label:
                        label[w] = g
                        stack.append(w)

        # The other genes join the group of their closest target that has a group (ties: smallest gene), in sorted
        # order, repeated until nothing changes; genes never reached form their own group.
        for _ in range(len(V)):
            changed = False
            for g in sorted(V):
                if g in label:
                    continue
                targets = sorted((x, h) for h, x in self.out.get(g, {}).items() if h in V and h in label)
                if targets:
                    label[g] = label[targets[0][1]]
                    changed = True
            if not changed:
                break
        for g in V:
            label.setdefault(g, g)
        groups = collections.defaultdict(set)
        for g, i in label.items():
            groups[i].add(g)
        if len(groups) < 2:
            return None

        # Duplication location: the species with genes in >= 2 groups must have the same LCA as all the species of V.
        count = collections.Counter(t for grp in groups.values() for t in {self.tax[g] for g in grp})
        duplicated = [self.sp[t] for t, n in count.items() if n >= 2]
        if not duplicated or common_prefix_len(duplicated) != common_prefix_len([self.sp[t] for t in count]):
            return None

        # Ratio test: accept if the arcs between groups are >= --mutual-ratio times longer (median) than inside groups.
        cross = [self.out[g][h] for g, h in arcs if label[g] != label[h]]
        within = [self.out[g][h] for g, h in arcs if label[g] == label[h]]
        if not cross or not within or statistics.median(cross) < o.mutual_ratio * statistics.median(within):
            return None
        return list(groups.values())

# ---------------------------------------------------------------------------------------------------------- main

# Read the options and graphs, reconstruct every family, write the trees.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", required=True, help="ldo_edges.tsv of predict_ldo.py")
    ap.add_argument("--genes", help="genes.tsv of predict_ldo.py (adds genes without arcs)")
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--trees-out", required=True)
    ap.add_argument("--repair", default="count", help="none, or count / sparse / mutual joined by + (default count)")
    ap.add_argument("--min-support", type=int, default=3, help="count: genes that must choose a target (default 3)")
    ap.add_argument("--ratio", type=float, default=2.0,
                    help="count: minimal median(arcs between groups) / median(arcs inside groups) (default 2)")
    ap.add_argument("--count-age", type=float, default=None,
                    help="count: also accept the split if the strong targets are >= this many times farther apart than "
                         "the median arc choosing them (default off; needs --seqs-trees)")
    ap.add_argument("--seqs-trees", help="folder of the JSON trees with the alignments (for --count-age)")
    ap.add_argument("--model", choices=["poisson", "kimura", "p", "jc"], default="poisson",
                    help="distance model for --count-age (as predict_ldo.py)")
    ap.add_argument("--min-overlap", type=int, default=50, help="as predict_ldo.py, for --count-age")
    ap.add_argument("--sparse-conductance", type=float, default=0.3, help="sparse: maximal conductance (default 0.3)")
    ap.add_argument("--sparse-min-side", type=int, default=2, help="sparse: genes on each side (default 2)")
    ap.add_argument("--sparse-ratio", type=float, default=1.5,
                    help="sparse: minimal median(arcs crossing the cut) / median(other arcs) (default 1.5; 0 = off)")
    ap.add_argument("--mutual-weight", type=float, default=2.0, help="sparse refinement: weight of a mutual arc (2)")
    ap.add_argument("--sparse-max-bridge", type=int, default=0,
                    help="sparse: reject the cut if more bridge genes than this (default 0; -1 = no check)")
    ap.add_argument("--sparse-bridge-ratio", type=float, default=1.5, help="sparse: bridge distance ratio (1.5)")
    ap.add_argument("--sparse-bridge-mutual", type=int, default=2, help="sparse: mutual arcs to each side (2)")
    ap.add_argument("--mutual-share", type=float, default=0.5,
                    help="mutual: tried only if fewer than this share of the arcs are reciprocal (default 0.5)")
    ap.add_argument("--mutual-ratio", type=float, default=1.5,
                    help="mutual: minimal median(arcs between groups) / median(arcs inside groups) (default 1.5)")
    a = ap.parse_args()
    if a.count_age is not None and not a.seqs_trees:
        sys.exit("--count-age needs --seqs-trees (the tree json files with the alignments)")
    repair_kinds(a.repair)

    sp = read_species_paths(a.species_tree)
    fams = load_graph(a.graph, a.genes)
    with open(a.trees_out, "w") as tf:
        for tid in sorted(fams):
            d = fams[tid]
            missing = {t for t in d["tax"].values() if t not in sp}
            if missing:
                print(f"{tid}: taxa not in the species tree: {sorted(missing)[:5]}; skipped", file=sys.stderr)
                continue
            tree = Builder(tid, d, sp, a.repair, a).build(set(d["tax"]))
            tf.write(f"{tid}\t{newick(tree)};\n")


if __name__ == "__main__":
    main()
