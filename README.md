## `ldo_builder.py`: from LDO graphs to gene trees and orthologs

`scripts/ldo_builder.py` implements Algorithm 3. Given the LDO graph of a gene family (every gene points to its
putative least diverged ortholog) and the species tree, it reconstructs the
gene tree with every internal node labelled as a speciation (`Spec`) or a duplication (`Dup`). The orthologs are then
read off the tree: two genes are orthologs if and only if their lowest common ancestor is a `Spec` node.

The reconstruction is top-down and recursive:

* if the graph is disconnected (ignoring arc directions), the root is a `Dup` node whose children are the trees of
  the connected components;
* if it is connected, the root is a `Spec` node whose children are the trees of the parts given by the species
  partition at x = lca(species of the graph) in the species tree: one part per child of x that has genes.

On the true LDO graph, this recovers all orthologs and paralogs exactly. On a predicted graph (built from sequence
distances by `predict_ldo.py`), optional *checks* (called repairs in the code: `count`, `sparse`, `mutual`) are tried
before each speciation. They detect false speciation signals and turn them into duplications. They are described in
detail in the paper and at the top of the script, and the paper's LDO+check is `--repair count+sparse+mutual --count-age 1.0`.

```bash
python scripts/ldo_builder.py --graph ldo_edges.tsv --species-tree species_tree.nwk --trees-out trees.txt \
    --repair none                                   # plain Algorithm 3
python scripts/ldo_builder.py --graph ldo_edges.tsv --genes genes.tsv --species-tree species_tree.nwk \
    --trees-out trees.txt --repair count+sparse+mutual --count-age 1.0 --seqs-trees TREES_DIR --model poisson
```

The default is `--repair count`, so give `--repair none` explicitly for plain Algorithm 3.



# Reproducing the results: LDO-graph reconstruction (Algorithm 3) vs OrthoFinder: SwissTree and simphy

Everything needed to reproduce the results of this work from the raw data:

* **SwissTree** (2014.0 ftp release, 9 usable families): LDO pipeline (`predict_ldo.py` → `ldo_builder.py`) and the
  OrthoFinder 3 baseline, scored against the curated trees (8 families) and, for ST006, the curators' pair list.
* **simphy** simulations (25 replicates × 100 gene families, 25 species, duplication/loss settings dl0–dl5): gene
  trees labelled with simphy's own events, LDO pipeline on the DNA alignments, OrthoFinder on DNA, same scoring.
* **Ensembl Compara** gene trees (downloaded from the Ensembl REST server): LDO pipeline on the protein alignments,
  and OrthoFinder, scored against the Ensembl trees, with their "dubious" nodes read as duplications and as
  speciations. The trees change with every Ensembl release; `ensembl/example_summary.md` is one run, not a reference.

## Layout

```
setup.sh                          OrthoFinder 3 (pinned commit + opt-in patch) and one Python env; writes env/env.sh
orthofinder_waterfall_norm.patch  opt-in OrthoFinder score-normalisation fallback
scripts/                          all Python scripts (one copy, shared by both pipelines; common code in ldo_utils.py)
swisstree/reproduce.sh            SwissTree pipeline: swisstree_raw -> OUT/summary.md
swisstree/expected_summary.md     report of a from-scratch run (PERFAM=1)
simphy/simphy_run.sh              simphy pipeline: simphy root folder -> OUT/summary.md, OUT/all_families.tsv
simphy/expected/                  summary.md and all_families.tsv of the full run (25 replicates)
ensembl/ensembl_run.sh            Ensembl pipeline: downloaded trees -> OUT/summary.md
ensembl/example_summary.md        report of a run on 158 primate gene trees (chromosome sample, --prune-taxon 9443)
```

## Quick start

```bash
# 0. software (needs git, Python >= 3.10 with venv, network access to GitHub and PyPI)
bash setup.sh env
source env/env.sh
#    simphy only: BLAST+ (blastn, makeblastdb) for OrthoFinder on DNA, e.g. apt-get install ncbi-blast+

# 1. SwissTree (~20 min with 2 threads, PERFAM=1 included); swisstree_raw = the SwissTree 2014.0 ftp release
PERFAM=1 bash swisstree/reproduce.sh swisstree_raw st_out
cat st_out/summary.md                       # compare with swisstree/expected_summary.md

# 2. simphy (~2 h with 2 cores for the 25 replicates, per-family OrthoFinder on families 1-10)
OF_MODES="pooled perfam perfam_fixnorm" PERFAM_FAMILIES=1-10 bash simphy/simphy_run.sh simphy sim_out
cat sim_out/summary.md                      # compare with simphy/expected/summary.md

# 3. Ensembl (needs network access to rest.ensembl.org and ftp.ensembl.org for the download)
$PY scripts/ensembl_download.py --out ens_data --chromosome 21 --max-trees 300 --prune-taxon 9443   # primates
bash ensembl/ensembl_run.sh ens_data ens_out
cat ens_out/summary.md
```

Environment variables:

| Variable | Used by | Meaning (default) |
|---|---|---|
| `PY`, `OF`, `MAFFT`, `MAFFT_BINARIES` | both | set by `env/env.sh`: Python, OrthoFinder, bundled MAFFT 7.511 |
| `THREADS` | both | threads for OrthoFinder (2) |
| `PERFAM` | SwissTree | `1` = also one OrthoFinder run per family (0) |
| `OF_MODES` | simphy | `pooled`, `perfam`, `perfam_fixnorm` (`pooled`; empty = no OrthoFinder) |
| `PERFAM_FAMILIES`, `JOBS` | simphy | families of the per-family runs, e.g. `1-10` (all); parallel runs (2) |
| `COUNT_AGE` | simphy | γ of the age test, run as `count+sparse+age` (1.0; empty = skip) |
| `ONLY` | simphy | regular expression on the replicate paths (all) |
| `SPECIES_NEWICK`, `MAX_GENES` | Ensembl | Compara species tree (`DATA/compara_species_tree.nh`); larger families skipped (1500) |
| `REPAIRS`, `COUNT_AGE` | Ensembl | as for simphy |
| `ORTHOFINDER` | Ensembl | `1` = also run OrthoFinder, pooled, standard and patched (1) |

A simphy replicate with a `done` file is skipped, so an interrupted run restarts with the same command; a replicate
that fails is listed in `OUT/failed.txt` and the run goes on.

## Software

* `setup.sh` clones OrthoFinder at commit `c85f0034679a9812959655f1cd60c74a30597cf1` (GitHub main, 13 May 2026), applies
  `orthofinder_waterfall_norm.patch` and installs it in `env/venv`. OrthoFinder brings numpy (the only Python
  dependency of the scripts) and DIAMOND, MCL, FAMSA, FastTree, FastME and MAFFT 7.511. `orthofinder --version` prints
  3.1.4.post1.dev1 or 3.1.5.post1.dev1 depending on the git tags fetched; it is the same commit.
* The patch does nothing unless `OF_NORM_LOG`, `OF_NORM_MIN`, `OF_NORM_PARAMS` are set, so unpatched and patched
  results come from the same install.
* BLAST+ 2.12.0 was used for simphy (OrthoFinder `-d -S blastn`; DIAMOND is protein only).

## Scripts

| Script | Role |
|---|---|
| `ldo_utils.py` | shared functions: Newick / NHX parsing, pruning and writing, species-tree paths, JSON gene trees, FASTA and genes.tsv reading, alignment distances, precision / recall |
| `swisstree_convert.py` | SwissTree release -> JSON trees with MAFFT alignments + species tree (`--match-gene-name`, `--min-length-frac`) |
| `truth_from_relationships.py` | curators' ortholog pairs from `gene_relationships.txt` (ST006 truth) |
| `simphy_prepare.py` | simphy replicate -> LCA-labelled JSON trees with alignments + OrthoFinder inputs |
| `true_orthologs.py` | true ortholog pairs (and genes.tsv) of labelled trees |
| `predict_ldo.py` | predicted LDO graph: closest gene of each species (`--model` poisson / kimura / p / jc) |
| `ldo_builder.py` | Algorithm 3 with the count repair (optional `--count-age`), the sparse cut and the opt-in mutual repair (`--repair count+sparse+mutual`, not used by the runners) |
| `evaluate_orthologs.py`, `cluster_analysis.py` | scoring of the reconstructed trees (pairs, clusters) |
| `make_proteomes.py`, `prune_species_tree.py`, `fit_norm_params.py`, `evaluate_orthofinder.py` | OrthoFinder inputs, normalisation parameters and scoring |
| `ensembl_download.py`, `ensembl_prepare.py` | Ensembl Compara gene trees (with alignments) and species tree -> JSON trees + species tree of the pipeline |
| `summarize.py`, `simphy_summarize.py`, `ensembl_summarize.py` | the `summary.md` reports |

Changes made in this work, compared with the handoff scripts (`ldo_builder.py` was called `ldo_repair.py`): `--match-gene-name` and `--min-length-frac` in
`swisstree_convert.py`; `--model jc` in `predict_ldo.py`; `--count-age` in `ldo_builder.py` (off by default; with it off
the trees are those of the handoff version); `orthofinder_waterfall_norm.patch` for OrthoFinder (inactive by default).
The code was then reduced to what the results need: options and code paths never used here were removed (Ensembl
download and top-k output of `predict_ldo.py`; the `cut` repair, the diagnostics and the tested-but-rejected sparse
options of `ldo_builder.py`; the true LDO graph of `compute_ldo.py`, now `true_orthologs.py`), and repeated code moved to
`ldo_utils.py`. All outputs were checked to be identical before and after.

## Expected results

Ortholog pairs, precision / recall / F1 with gene pairs pooled over families (micro average), unless stated.

**SwissTree** (`swisstree/expected_summary.md` has the per-tree table)

Curated trees, 8 families (ST001–005, ST008–010); clusters over all 9 families.

| Run | Curated, 8 families | Clusters P / R |
|---|---|---|
| LDO, no repair | 0.843 / 0.999 / 0.914 | 0.726 / 0.684 |
| LDO, count | 0.898 / 0.997 / 0.945 | 0.805 / 0.780 |
| LDO, count+sparse | 0.947 / 0.987 / 0.967 | 0.831 / 0.818 |
| LDO, count+sparse + age test | 0.970 / 0.982 / **0.976** | 0.820 / 0.799 |
| OrthoFinder, pooled | 0.991 / 0.819 / 0.897 | |
| OrthoFinder, pooled, patched | 0.996 / 0.962 / **0.979** | |
| OrthoFinder, per family | 0.990 / 0.753 / 0.855 | |
| OrthoFinder, per family, patched | 0.996 / 0.962 / 0.979 | |

ST006 against the curators' pair list (`st` = plain conversion, 64 genes; `stf` = gene-name matching and fragments
removed, 75 genes; OrthoFinder on `st` only):

| Run | `st` | `stf` |
|---|---|---|
| LDO, count+sparse | 0.623 / 0.962 / 0.756 | 0.684 / 0.960 / 0.799 |
| LDO, count+sparse + age test | 0.655 / 0.928 / 0.768 | 0.719 / 0.951 / 0.819 |
| OrthoFinder, pooled | 0.841 / 0.959 / 0.896 | |
| OrthoFinder, pooled, patched | 0.864 / 0.935 / 0.898 | |

**simphy**, pooled pairs, F1 (`simphy/expected/summary.md` has precision / recall, per-family averages, clusters)

| Method | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|
| LDO, count+sparse | 1.000 | 0.796 | 0.507 | 0.358 | 0.323 | 0.434 |
| LDO, count+sparse + age test | 1.000 | 0.829 | 0.595 | 0.572 | 0.427 | 0.578 |
| OrthoFinder, pooled | 0.789 | 0.759 | 0.663 | 0.554 | 0.459 | 0.582 |

Per-family averages (precision / recall), all settings: LDO count+sparse 0.776 / 0.972, + age test 0.800 / 0.918,
OrthoFinder pooled 0.860 / 0.702. (Truth: simphy's own events. With the LCA-mapping truth of earlier versions, the
"All" F1 were 0.505, 0.646 and 0.659.)
