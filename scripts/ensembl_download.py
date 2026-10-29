#!/usr/bin/env python3
"""
Download Ensembl Compara gene trees (JSON, with event labels, branch lengths and the aligned protein sequences in the
leaves) from the Ensembl REST server, and the Compara species tree.

  python ensembl_download.py --out DATA --genes HBB HBA1 MB HOXA1 TP53 [--prune-taxon 9443]
  python ensembl_download.py --out DATA --chromosome 21 --max-trees 300 --prune-taxon 9443
  python ensembl_download.py --out DATA --gene-file ensg_ids.txt

Queries (several can be combined; the trees are downloaded in this order until --max-trees distinct trees):
  --genes        gene symbols of --species (default homo_sapiens)
  --gene-file    one Ensembl gene id per line
  --chromosome   the protein-coding genes of that chromosome of --species, in random order (--seed)
Every gene gives the Compara gene tree it belongs to (one file per tree; a tree reached twice is downloaded once).
--prune-taxon T keeps only the genes of the species below NCBI taxon T (e.g. 9443 = primates, 40674 = mammals).

Writes in DATA:
  trees/<tree_id>.json             one gene tree per file ({"tree": {"id": ..., "tree": root}}); the leaves carry
                                   id.accession (gene id), taxonomy (id, scientific_name) and sequence.mol_seq.seq
                                   (aligned protein); internal nodes carry events.type (speciation, duplication,
                                   dubious, gene_split)
  compara_species_tree.nh          the Compara species tree of --compara (from the Ensembl ftp site), unless
                                   --no-species-tree; if that download fails, get the file by hand (see the message)
  download_log.tsv                 query, tree id or "no tree"
Files already in DATA/trees are kept: an interrupted download restarts with the same command.
Needs network access to rest.ensembl.org (and ftp.ensembl.org for the species tree).
"""
import argparse
import json
import os
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from ldo_utils import get_root, json_leaves, leaf_seq

SPECIES_TREE_URL = "https://ftp.ensembl.org/pub/release-116/compara/species_trees/{compara}_species-tree_Ensembl.nh"
MIN_INTERVAL = 0.085            # seconds between requests (Ensembl allows about 15 per second)
_last = [0.0]


# GET a JSON document from the REST server, with rate limiting and retries (HTTP 429 / 5xx / network errors).
#   server    REST server URL (e.g. https://rest.ensembl.org)
#   path      endpoint path (e.g. /genetree/member/id/homo_sapiens/ENSG00000139618)
#   params    {name: value} query parameters
#   retries   number of attempts
# Returns the parsed JSON, or None (HTTP 400 / 404: no tree for this gene; or all attempts failed).
def get(server, path, params=None, retries=6):
    url = server + path + ("?" + urllib.parse.urlencode(params) if params else "")
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "Accept": "application/json"})
    for attempt in range(retries):
        wait = MIN_INTERVAL - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(float(e.headers.get("Retry-After", 2)))
                continue
            if e.code in (500, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            print(f"  HTTP {e.code} for {url}", file=sys.stderr)
            return None
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            print(f"  network error ({e}), retrying", file=sys.stderr)
            time.sleep(2 ** attempt)
    return None


# Protein-coding gene ids of a chromosome (the overlap endpoint accepts windows of at most 5 Mb).
#   server     REST server URL
#   species    Ensembl species name (e.g. homo_sapiens)
#   chrom      chromosome name (e.g. "21")
# Returns the sorted list of gene ids.
def chromosome_gene_ids(server, species, chrom):
    info = get(server, f"/info/assembly/{species}/{chrom}")
    if not info:
        sys.exit(f"chromosome {chrom} of {species} not found")
    step, ids = 4_900_000, set()
    for start in range(1, info["length"] + 1, step):
        end = min(start + step - 1, info["length"])
        res = get(server, f"/overlap/region/{species}/{chrom}:{start}-{end}",
                  {"feature": "gene", "biotype": "protein_coding"})
        ids.update(g["id"] for g in (res or []))
    return sorted(ids)


# Whether every leaf of a gene tree carries an aligned sequence, all of the same length.
#   root   root node of the JSON gene tree
# Returns True or False.
def aligned(root):
    seqs = [leaf_seq(x) for x in json_leaves(root)]
    return bool(seqs) and all(seqs) and len({len(s) for s in seqs}) == 1


# Download the gene tree of one gene, with the aligned protein sequences. If the tree comes back without a complete
# alignment, it is requested again by tree id.
#   a      parsed command-line options (server, species, compara, prune_taxon)
#   kind   "symbol" or "id"
#   val    gene symbol or gene id
#   seen   ids of the trees already downloaded (no second request for their alignment)
# Returns (tree id, root node), or (None, None) if the gene has no tree.
def fetch_tree(a, kind, val, seen):
    params = {"compara": a.compara, "sequence": "protein", "aligned": 1}
    if a.prune_taxon:
        params["prune_taxon"] = a.prune_taxon
    data = get(a.server, f"/genetree/member/{kind}/{a.species}/{val}", params)
    tid, root = get_root(data) if data else (None, None)
    if root is None:
        return None, None
    tid = tid or val
    if tid not in seen and not aligned(root):
        again = get(a.server, f"/genetree/id/{tid}", params)
        _, root2 = get_root(again) if again else (None, None)
        if root2 is not None and aligned(root2):
            root = root2
        else:
            print(f"  {tid}: no complete alignment in the answer; leaves without sequence will be dropped later",
                  file=sys.stderr)
    return tid, root


# Download the Compara species tree to DATA/compara_species_tree.nh.
#   a   parsed command-line options (out, compara)
def fetch_species_tree(a):
    dst = os.path.join(a.out, "compara_species_tree.nh")
    if os.path.exists(dst):
        return
    url = SPECIES_TREE_URL.format(compara=a.compara)
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            open(dst, "wb").write(r.read())
        print(f"species tree: {url} -> {dst}")
    except (urllib.error.URLError, OSError) as e:
        print(f"WARNING: could not download the species tree ({e}).\n  Download {a.compara}_species-tree_Ensembl.nh "
              f"from the Ensembl ftp site (pub/release-N/compara/species_trees/) and save it as {dst}.", file=sys.stderr)


# Run the queries and write the trees.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", nargs="*", help="gene symbols of --species")
    ap.add_argument("--gene-file", help="file with one Ensembl gene id per line")
    ap.add_argument("--chromosome", help="sample the protein-coding genes of this chromosome of --species")
    ap.add_argument("--species", default="homo_sapiens", help="species of the queried genes (homo_sapiens)")
    ap.add_argument("--max-trees", type=int, default=100, help="stop after this many distinct trees (100)")
    ap.add_argument("--seed", type=int, default=1, help="order of the chromosome genes (1)")
    ap.add_argument("--prune-taxon", type=int, help="keep only the species below this NCBI taxon (e.g. 9443)")
    ap.add_argument("--compara", default="vertebrates", help="Compara division (vertebrates)")
    ap.add_argument("--server", default="https://rest.ensembl.org", help="REST server (e.g. a mirror)")
    ap.add_argument("--no-species-tree", action="store_true", help="do not download the Compara species tree")
    a = ap.parse_args()

    queries = [("symbol", g) for g in (a.genes or [])]
    if a.gene_file:
        queries += [("id", x.strip()) for x in open(a.gene_file) if x.strip()]
    if a.chromosome:
        ids = chromosome_gene_ids(a.server, a.species, a.chromosome)
        random.Random(a.seed).shuffle(ids)
        queries += [("id", x) for x in ids]
    if not queries:
        sys.exit("give --genes, --gene-file or --chromosome")

    treedir = os.path.join(a.out, "trees")
    os.makedirs(treedir, exist_ok=True)
    seen = {f[:-5] for f in os.listdir(treedir) if f.endswith(".json")}       # resume
    with open(os.path.join(a.out, "download_log.tsv"), "a") as log:
        for kind, val in queries:
            if len(seen) >= a.max_trees:
                break
            tid, root = fetch_tree(a, kind, val, seen)
            log.write(f"{val}\t{tid or 'no tree'}\n")
            if root is None:
                print(f"{val}: no tree", file=sys.stderr)
                continue
            if tid in seen:
                continue
            seen.add(tid)
            with open(os.path.join(treedir, tid + ".json"), "w") as f:
                json.dump({"tree": {"id": tid, "tree": root}}, f)
            print(f"{val} -> {tid}: {len(json_leaves(root))} genes  ({len(seen)}/{a.max_trees})")
    if not a.no_species_tree:
        fetch_species_tree(a)


if __name__ == "__main__":
    main()
