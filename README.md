# LDO-graph reconstruction (Algorithm 3) vs OrthoFinder: SwissTree and simphy

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

## SwissTree pipeline (`swisstree/reproduce.sh`)

Run from the output folder; `S` = `scripts/`, `RAW` = `swisstree_raw/`.

```bash
# conversion: st = the 9 families; stf = ST006 fixed (only ST006 differs)
python $S/swisstree_convert.py --raw $RAW --out st  --mafft $MAFFT
python $S/swisstree_convert.py --raw $RAW --out stf --mafft $MAFFT --match-gene-name --min-length-frac 0.5
# truth (for d in st stf)
python $S/true_orthologs.py --trees $d/trees --out $d/ldo_cur                    # curated trees
python $S/truth_from_relationships.py --raw $RAW --genes $d/ldo_cur/genes.tsv --out $d/curator_orthologs.tsv  # ST006
# LDO (for d in st stf; r in none count count+sparse)
python $S/predict_ldo.py --trees $d/trees --out $d/pred
python $S/ldo_builder.py --graph $d/pred/ldo_edges.tsv --genes $d/pred/genes.tsv --species-tree $d/species_tree.nwk \
    --repair $r --trees-out $d/trees_$r.txt
python $S/ldo_builder.py ... --repair count+sparse --count-age 1.0 --seqs-trees $d/trees --trees-out $d/trees_count+sparse+age.txt
python $S/evaluate_orthologs.py --trees <trees> --truth $d/ldo_cur/orthologs.tsv --out <dir>
python $S/evaluate_orthologs.py --trees <ST006 line of trees> --truth $d/curator_orthologs.tsv --out <dir>
python $S/cluster_analysis.py   --pred  <trees> --true-trees $d/trees
# OrthoFinder (inputs from st)
python $S/make_proteomes.py --genes st/ldo_cur/genes.tsv --seqs st/seqs --out of_in
python $S/prune_species_tree.py st/species_tree.nwk of_in/pooled of_in/fam/ST*
orthofinder -f of_in/pooled -o of_out/pooled -s of_in/pooled.sp.nwk -X   # DIAMOND search + standard analysis
W=$(ls -d of_out/pooled/Results_*/WorkingDirectory)                         # (all 9 families, true species tree)
OF_NORM_LOG=$PWD/of_out/normfit.tsv orthofinder -b $W -s of_in/pooled.sp.nwk -n fitlog -X -og
NORM=$(python $S/fit_norm_params.py of_out/normfit.tsv --min-hits 30)                       # 0.636,-0.807
OF_NORM_MIN=30 OF_NORM_PARAMS=$NORM orthofinder -b $W -s of_in/pooled.sp.nwk -n fixnorm30 -X   # "patched"
python $S/evaluate_orthofinder.py --truth st/ldo_cur/orthologs.tsv --genes st/ldo_cur/genes.tsv --run <Results> ...
python $S/evaluate_orthofinder.py --truth st/curator_orthologs.tsv --genes st/ldo_cur/genes.tsv --run <Results> \
    --families ST006 ...
# PERFAM=1: the same per family (of_in/fam/STxxx, -s of_in/fam/STxxx.sp.nwk), the patched run reusing its search
python $S/summarize.py --eval eval --out summary.md
```

* **Truth.** ST001–005 and ST008–010: the ortholog pairs of the curated trees (`consensus_tree.nhx` event tags); this
  is the main result. ST004 has no curators' pair list; for the other seven the tree pairs contain every pair of
  the list. ST006: its tree file has no event tags, so it is scored separately against the `ortholog` pairs of the
  curators' `gene_relationships.txt` (`truth_from_relationships.py`, names matched to our genes), for every method.
  Truth pairs are limited to the genes kept by the conversion.
* **`stf`**: `--match-gene-name` finds 18 ST006 sequences of leaves named by Ensembl ids (unique `GN=` + species match);
  `--min-length-frac 0.5` drops 7 homeodomain fragments (60–100 aa, family median 282).

## simphy pipeline (`simphy/simphy_run.sh`)

Every folder under the simphy root that contains `s_tree.trees` is a replicate. Per replicate, in
`OUT/<setting>/seed<NNNN>/`:

```bash
python $S/simphy_prepare.py --run $RUN --out .          # simphy's events, trees/*.json with alignments, of_in/
python $S/true_orthologs.py --trees trees --out ldo
python $S/predict_ldo.py --trees trees --out pred --model jc
python $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
    --repair $r --trees-out ldo_trees_$r.txt             # r = none, count, count+sparse
python $S/ldo_builder.py ... --repair count+sparse --count-age 1.0 --seqs-trees trees --model jc   # count+sparse+age
python $S/evaluate_orthologs.py ...; python $S/cluster_analysis.py ... --out eval/clusters_$r.tsv
OF_NORM_LOG=$PWD/of_out/normfit.tsv orthofinder -f of_in/pooled -o of_out/pooled -d -S blastn -s of_in/pooled.sp.nwk -X
orthofinder -f of_in/fam/$f -o of_out/perfam_$f -d -S blastn -s of_in/fam/$f.sp.nwk -X -t 1 -a 1          # perfam
OF_NORM_MIN=30 OF_NORM_PARAMS=$(python $S/fit_norm_params.py of_out/normfit.tsv) orthofinder ...           # perfam_fixnorm
python $S/evaluate_orthofinder.py --truth ldo/orthologs.tsv --genes ldo/genes.tsv --run <Results> ...
# after all replicates
python $S/simphy_summarize.py --out OUT
```

* **Labels (truth)**: every gene-tree node gets simphy's own locus-tree event (Sp / Dup), read from `<NNN>.mapsl`
  and `<NNN>l1g.maplg` (no `.mapsl` when there are no duplications or losses, as in dl0: all speciations). Species of
  leaf `xx_yy_zz` = `xx`. `--labels lca` labels by LCA mapping on the species tree instead (the truth of earlier
  versions): it misses the duplications followed by complementary losses (hidden paralogy), which look like
  speciations: 1390 of the 38,617 duplication nodes, so that 266,151 of its 1,470,106 "true" pairs (18%) are
  paralogs (25 replicates). Since Algorithm 3 also labels by the species tree, that truth favoured it. ILS is
  negligible (population size 10), so LCA mapping never mistakes a speciation for a duplication. Rebuilding from the
  true LDO graphs of these labelled trees still gives 1.000 / 1.000 (the hidden paralogs have no arc between them).
* **Gene ids** are `<family>.<leaf>` (unique in a replicate), so families can be pooled in one OrthoFinder run.
* **`dataset_NNN.fasta`** is already aligned (INDELible without indels) and used as is; OrthoFinder gets it ungapped.
* **OrthoFinder per family** loses its hits (species pairs with < 2 hits are discarded; a single-copy family: recall 0).
  `perfam_fixnorm` normalises pairs with < 30 hits with the parameters fitted in the pooled run of the replicate. The
  pooled run needs no patch and is the recommended baseline. Families with < 2 species are not run per family.

## Ensembl pipeline (`ensembl/ensembl_run.sh`)

The download is a separate step, since it needs network access (the scripts use only the Python standard library):

```bash
python $S/ensembl_download.py --out DATA --genes HBB HBA1 MB HOXA1 TP53 [--prune-taxon 9443]   # by gene symbol
python $S/ensembl_download.py --out DATA --chromosome 21 --max-trees 300 --prune-taxon 9443   # sample of a chromosome
python $S/ensembl_download.py --out DATA --gene-file ensg_ids.txt                             # by gene id
```

It writes one JSON tree per Compara gene tree in `DATA/trees/` (event labels, branch lengths and the aligned protein
sequences in the leaves; with `--prune-taxon`, only the species below that NCBI taxon) and the Compara species tree
`DATA/compara_species_tree.nh` (from `ftp.ensembl.org/pub/release-116/compara/species_trees/`; if that fails, the script
says where to put the file). An interrupted download restarts with the same command. Then, in `OUT/`:

```bash
python $S/ensembl_prepare.py --trees DATA/trees --species-newick DATA/compara_species_tree.nh --out .
python $S/true_orthologs.py --trees trees --out truth_dup  --dubious dup      # dubious nodes = duplications
python $S/true_orthologs.py --trees trees --out truth_spec --dubious spec     # dubious nodes = speciations
python $S/predict_ldo.py --trees trees --out pred                            # Poisson distance, proteins
python $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
    --repair $r --trees-out trees_$r.txt                                     # r = none, count, count+sparse
python $S/ldo_builder.py ... --repair count+sparse --count-age 1.0 --seqs-trees trees    # count+sparse+age
python $S/evaluate_orthologs.py --trees trees_$r.txt --truth truth_$dub/orthologs.tsv --out ...   # dub = dup, spec
python $S/cluster_analysis.py --pred trees_$r.txt --true-trees trees --dubious $dub --out eval/clusters_${r}_$dub.tsv
# OrthoFinder (ORTHOFINDER=1), as for SwissTree: all families pooled, true species tree, standard then patched
python $S/make_proteomes.py --genes truth_dup/genes.tsv --seqs seqs --out of_in
python $S/prune_species_tree.py species_tree.nwk of_in/pooled
orthofinder -f of_in/pooled -o of_out/pooled -s of_in/pooled.sp.nwk -X              # search + standard analysis
OF_NORM_LOG=$PWD/of_out/normfit.tsv orthofinder -b $W -s of_in/pooled.sp.nwk -n fitlog -X -og     # fit log
OF_NORM_MIN=30 OF_NORM_PARAMS=$(python $S/fit_norm_params.py of_out/normfit.tsv) orthofinder -b $W ... -n fixnorm30
python $S/evaluate_orthofinder.py --truth truth_$dub/orthologs.tsv --genes truth_$dub/genes.tsv --run <Results> ...
python $S/ensembl_summarize.py --out .
```

* **Preparation** (`ensembl_prepare.py`): the Compara species tree is matched to the taxa of the trees (by taxon id or
  by species name) and pruned to them, with taxon-id leaves. Leaves without an aligned sequence or of a species not
  in the species tree are removed; families with more than `MAX_GENES` genes or fewer than 2 species are skipped.
  `prepare_report.tsv` counts per family the speciation, duplication, dubious and other nodes, the speciations the
  species tree cannot explain, and the dubious nodes that would be valid speciations.
* **Dubious nodes.** Ensembl Compara calls a duplication node "dubious" when its two subtrees share few species (low
  duplication-consistency score): it may be a speciation in a tree with errors. The prepared trees keep them as they
  are; the truth is computed both ways (`--dubious dup` / `spec` of `true_orthologs.py` and `cluster_analysis.py`;
  without the option, as in the other pipelines, dubious nodes are duplications). The predictions use only the
  sequences, so they do not depend on this choice. The report has three columns: dubious = Dup, dubious = Spec, and
  the families without dubious nodes (both readings agree there).
* **Truth.** The Ensembl trees as they are. Their speciation nodes are not checked against the species tree (the
  report counts those that the species tree cannot explain).
* **How Ensembl labels its nodes** (seen on the 158 primate trees): every speciation node and every dubious node has
  children without a common species, so reading the dubious nodes as Spec amounts to labelling by species overlap.
  With the Compara species tree, 665 of the 4391 speciation nodes (15.1%) contradict the species tree (the gene
  tree disagrees with it). Algorithm 3 cannot produce such a node, but the ortholog pairs that go through it can still
  be predicted, through a different speciation of the reconstructed tree. Measured on these trees, LDO+check recovers
  98.7% (dubious = Dup) / 92.2% (dubious = Spec) of those pairs, OrthoFinder 97.2% / 88.3%: the contradicting
  speciations do not penalise LDO+check.
* **Species tree.** Use the Ensembl Compara species tree of the division and release of the download (e.g.
  `https://ftp.ensembl.org/pub/release-116/compara/species_trees/vertebrates_species-tree_Ensembl.nh`; since the
  move to the new Ensembl platform (August 2026), `pub/current_compara/` may not hold it), given with
  `SPECIES_NEWICK=...` if it is not in `DATA/compara_species_tree.nh`. Its leaves are matched by species name.

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

## Verification

This directory was run from scratch (fresh unzip of both datasets, fresh `setup.sh`, OrthoFinder cloned again):
`swisstree/reproduce.sh` with `PERFAM=1` reproduced `swisstree/expected_summary.md` exactly, and `simphy/simphy_run.sh`
on one replicate per setting (`ONLY=seed3001`, all OrthoFinder modes) gave evaluation files byte-identical to those
of the full run in `simphy/expected/` (55 files, plus the fitted OrthoFinder parameters).

After the last code changes (`ldo_repair.py` restructured and renamed `ldo_builder.py`, the relabelled truth dropped,
comments added), `swisstree/reproduce.sh` (`PERFAM=1`) was run again: all reconstructed trees, curated-truth scores,
cluster scores and OrthoFinder scores were byte-identical to the previous run, and `expected_summary.md` is that run's
report. `ldo_builder.py` gave byte-identical trees to `ldo_repair.py` on all 25 simphy replicates (4 settings each),
and `simphy/simphy_run.sh` on `ONLY=seed3001` with all OrthoFinder modes reproduced the evaluation files of the full
run (OrthoFinder per-family counts identical).

After the OrthoFinder runs were simplified (the first run of each input folder now takes the true species tree and
is the scored standard run, instead of a search-only run without it) and `evaluate_orthofinder.py` was commented and
its counting rewritten, `swisstree/reproduce.sh` (`PERFAM=1`) was run again: all 32 evaluation files and
`summary.md` were byte-identical to the previous run, with the same fitted parameters (0.636, -0.807).
`evaluate_orthofinder.py` gave identical output to its previous version on the SwissTree and simphy OrthoFinder results
(25 comparisons).

**simphy truth from simphy's events.** `simphy_prepare.py` now labels the gene trees with simphy's own events
instead of LCA mapping. The predictions were not recomputed (they do not depend on the labels); all the evaluation
files of the 25 replicates were recomputed against the new truth. With `--labels lca`, the same rescoring reproduced
the previous evaluation files byte for byte (one replicate checked, 38 files; the OrthoFinder files differ only by the
older file format, with identical counts). `simphy/expected/` is the rescored run.

**Ensembl pipeline.** It was run on a real download (158 primate gene trees, 5353 genes; `ensembl/example_summary.md`).
Before that, since Ensembl could not be reached from the machine where this bundle was built, it was tested without
real Ensembl data: `ensembl_download.py` against a local stand-in of the REST endpoints it uses
(trees by gene symbol / id, the second request for a missing alignment, chromosome sampling, resume, genes without
tree), and the rest on the 9 SwissTree families turned into Ensembl-like trees (species names in the leaves and the
species tree, about a third of the duplications and 5% of the speciations made dubious, a gene_split node, a leaf
without sequence, a leaf of a species missing from the species tree). The reconstructed trees were byte-identical to
those of the SwissTree pipeline (all four settings, the changed family excepted), the dubious-as-Spec truth added
exactly the pairs separated by the dubious nodes, and the SwissTree and simphy truth and cluster files were unchanged
by the new `--dubious` option (default dup).

## Caveats

* All repair thresholds, the age test γ, the OrthoFinder normalisation parameters and the `stf` choices were developed
  on these data; there is no held-out set.
* The age test trades tree topology for orthology: cluster precision / recall drop (simphy 0.689 / 0.665 → 0.612 /
  0.609; SwissTree 0.831 / 0.818 → 0.820 / 0.799) because the groups built by the count repair are imprecise.
* Algorithm 3 cannot produce a curated speciation node that contradicts the species tree (OrthoFinder labels by
  species overlap and can), but the ortholog pairs through such a node can still be predicted through another
  speciation; on Ensembl, LDO+check recovers these pairs at least as well as OrthoFinder (see the Ensembl section).
  Relabelling the curated trees from the species
  tree turns 7 of the 298 curated speciation nodes of the 8 tagged families into duplications (none the other way): there the curated trees and the
  species tree disagree. (An earlier version also scored against these relabelled trees; it penalised both methods
  alike and was dropped.)
* LDO is given the family partition; pooled OrthoFinder has to find it.
* The fair SwissTree OrthoFinder baseline would use the full QfO 2011 reference proteomes (not run here).
* Ensembl: the first real download may need small fixes (the REST answers were not seen from here). The Ensembl trees
  are themselves inferred (TreeBeST reconciliation), so they are a reference, not a ground truth; their speciation
  nodes that the species tree cannot explain are counted in the report but not corrected.
