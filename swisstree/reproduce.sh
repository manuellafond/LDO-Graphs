#!/bin/bash
# Reproduce all results (LDO pipeline and OrthoFinder baseline) from the SwissTree 2014.0 ftp download.
#
#   source env/env.sh                                   # from ../setup.sh (defines PY, OF, MAFFT, MAFFT_BINARIES)
#   bash swisstree/reproduce.sh SWISSTREE_RAW [OUT]     # OUT defaults to ./out and must not exist yet
#
# Environment: THREADS (default 2), PERFAM=1 to also run OrthoFinder once per family (about 15 extra minutes).
# Time with 2 threads: about 15 minutes, almost all of it the OrthoFinder DIAMOND search.
# The final report is OUT/summary.md; every step's output is in OUT/logs/.
#
# True orthologs (two genes are orthologs iff their last common ancestor is a speciation):
#   ST001-005, ST008-010  from the curated trees (consensus_tree.nhx event tags), true_orthologs.py. Main result:
#                         every method scored on these 8 families, gene pairs pooled. ST004 has no curators' pair
#                         list; for the others the tree pairs include every pair of the list.
#   ST006                 its curated tree has no event tags: the 'ortholog' pairs of the curators' list
#                         (gene_relationships.txt, matched to our genes by name), truth_from_relationships.py.
#                         Scored separately, every method, on st and stf.
#   Truth pairs are limited to the genes kept by the conversion; st and stf differ only in ST006.
#   Not used: ST007 (gene transfer), ST011 / ST012 (species missing from the species tree).
set -euo pipefail

# Arguments and programs: RAW = the SwissTree release, OUT = the output folder, S = the scripts folder;
# PY / OF / MAFFT come from env/env.sh (plain names otherwise), T = OrthoFinder threads.
RAW=$(realpath "$1")
OUT=$(realpath -m "${2:-out}")
S=$(cd "$(dirname "$0")/.." && pwd)/scripts
PY=${PY:-python3}; OF=${OF:-orthofinder}; MAFFT=${MAFFT:-mafft}; T=${THREADS:-2}

# Start from a new, empty output folder (never overwrite a previous run); everything below runs inside it.
[ -e "$OUT" ] && { echo "$OUT exists; give a new output folder" >&2; exit 1; }
mkdir -p "$OUT/eval" "$OUT/logs"
cd "$OUT"
export PYTHONDONTWRITEBYTECODE=1

# Progress message with the time.
step() { echo "== $(date +%H:%M:%S) $*"; }

# ---------------------------------------------------------------- 1. conversion
# SwissTree release -> JSON gene trees with MAFFT alignments in the leaves, species tree, raw and aligned FASTA.
# st : the 9 usable families (ST001-006, 008-010; ST007 has a gene transfer, ST011/012 species missing)
# stf: same, plus ST006 leaves matched by gene name (18 genes named by Ensembl ids) and fragments
#      (< 0.5 x family median length, 7 ST006 genes) removed. Only ST006 differs from st.
step "convert (MAFFT alignments)"
$PY $S/swisstree_convert.py --raw "$RAW" --out st  --mafft "$MAFFT"                                         > logs/convert_st.log 2>&1
$PY $S/swisstree_convert.py --raw "$RAW" --out stf --mafft "$MAFFT" --match-gene-name --min-length-frac 0.5 > logs/convert_stf.log 2>&1

# ---------------------------------------------------------------- 2. truth
# For each conversion: the ortholog pairs (and genes.tsv) of the curated trees, in <d>/ldo_cur/, and the ortholog
# pairs of the curators' lists, in <d>/curator_orthologs.tsv (only its ST006 pairs are used; see the header).
step "truth: curated trees, curators' pair lists"
for d in st stf; do
  $PY $S/true_orthologs.py --trees $d/trees --out $d/ldo_cur > logs/true_orthologs_$d.log 2>&1
  $PY $S/truth_from_relationships.py --raw "$RAW" --genes $d/ldo_cur/genes.tsv --out $d/curator_orthologs.tsv \
      > logs/curator_pairs_$d.log 2>&1
done

# ---------------------------------------------------------------- 3. LDO prediction + reconstruction
# Score one file of reconstructed trees:
#   $1 trees file (ldo_builder.py output)   $2 key of the eval files (ldo_<data>_<repair>)   $3 conversion folder
# writes eval/<key>_curated.tsv (curated-tree pairs), eval/<key>_st006.tsv (ST006 against the curators' list) and
# eval/clusters_<data>_<repair>.txt (clusters against the curated trees).
evaluate_ldo() {
  local tr=$1 key=$2 d=$3
  $PY $S/evaluate_orthologs.py --trees $tr --truth $d/ldo_cur/orthologs.tsv --out eval/tmp > logs/eval_${key}_curated.log
  cp eval/tmp/per_family.tsv eval/${key}_curated.tsv
  grep -P "^ST006\t" $tr > eval/tmp_st006.txt
  $PY $S/evaluate_orthologs.py --trees eval/tmp_st006.txt --truth $d/curator_orthologs.tsv --out eval/tmp \
      > logs/eval_${key}_st006.log
  cp eval/tmp/per_family.tsv eval/${key}_st006.tsv
  $PY $S/cluster_analysis.py --pred $tr --true-trees $d/trees > eval/clusters_${key#ldo_}.txt 2>&1
}
for d in st stf; do
  tag=$([ $d = st ] && echo base || echo frag)
  step "LDO: predict_ldo.py + ldo_builder.py on $d"
  # Predicted LDO graph of every family: for each gene, the closest gene of every other species.
  $PY $S/predict_ldo.py --trees $d/trees --out $d/pred > logs/predict_$d.log 2>&1
  # Gene trees from the graphs with each repair setting, then scoring.
  for r in none count count+sparse; do
    $PY $S/ldo_builder.py --graph $d/pred/ldo_edges.tsv --genes $d/pred/genes.tsv --species-tree $d/species_tree.nwk \
        --repair $r --trees-out $d/trees_$r.txt > logs/builder_${d}_$r.log 2>&1
    evaluate_ldo $d/trees_$r.txt ldo_${tag}_$r $d
  done
  # count+sparse with the duplication-age test of the count repair (developed on the simphy simulations); it reads
  # the alignments in <d>/trees/.
  $PY $S/ldo_builder.py --graph $d/pred/ldo_edges.tsv --genes $d/pred/genes.tsv --species-tree $d/species_tree.nwk \
      --repair count+sparse --count-age 1.0 --seqs-trees $d/trees --trees-out $d/trees_count+sparse+age.txt \
      > logs/builder_${d}_count+sparse+age.log 2>&1
  evaluate_ldo $d/trees_count+sparse+age.txt ldo_${tag}_count+sparse+age $d
done

# ---------------------------------------------------------------- 4. OrthoFinder
# How OrthoFinder is used, and why the comparison with the LDO pipeline is as fair as we could make it:
#   * Same genes: OrthoFinder gets exactly the genes of the LDO baseline (conversion st), raw sequences, one FASTA
#     per species.
#   * Same species tree: OrthoFinder is given the true species tree (-s), as the LDO reconstruction is, instead of
#     inferring its own.
#   * Its own best output: we score OrthoFinder's ortholog pairs (its Orthologues/ folder, from its gene-tree
#     reconciliation), not its orthogroups, which would mix in paralogs.
#   * Two settings: pooled (all 9 families in one run: OrthoFinder's normal use, it must find the families itself)
#     and, with PERFAM=1, per family (one run per family: the same family knowledge as the LDO pipeline).
#   * Fixed version: OrthoFinder 3 at a pinned commit (setup.sh); DIAMOND search with default settings.
#
# Score normalisation, and the patch. OrthoFinder corrects every hit score for sequence length with a line fitted
# on the hits of each pair of species. For a pair with fewer than 2 hits no line can be fitted, and OrthoFinder
# discards the pair's hits. With small families (as here, and in the simulations) many species pairs share only
# one or two genes, so true orthologs are lost before the analysis starts: a limitation of the input size, not of
# OrthoFinder's method. We therefore also run OrthoFinder with an opt-in patch (orthofinder_waterfall_norm.patch,
# inactive unless the OF_NORM_* variables are set):
#   1. a run with OF_NORM_LOG writes the line fitted for every species pair (normfit.tsv);
#   2. fit_norm_params.py takes the median line over the pairs with >= 30 hits (the reliable fits): 0.636, -0.807;
#   3. the patched run (OF_NORM_MIN=30 OF_NORM_PARAMS=...) uses that median line for the pairs with < 30 hits;
#      pairs with >= 30 hits are treated exactly as in standard OrthoFinder.
# The parameters come from OrthoFinder's own fits on these sequences; no truth is used. The per-family runs reuse
# the line fitted on the pooled run (a single family has too few hits to fit one).
# Effect on the curated truth (8 families, precision / recall / F1):
#   pooled       standard 0.991 / 0.819 / 0.897   ->  patched 0.996 / 0.962 / 0.979   (recall +0.143)
#   per family   standard 0.990 / 0.753 / 0.855   ->  patched 0.996 / 0.962 / 0.979   (recall +0.209)
# Both versions are reported; the patched one is the OrthoFinder baseline we compare with.
#
# Runs: the first run of an input folder (-f) does the DIAMOND all-vs-all search (the slow part) and the standard
# analysis; the other runs reuse its search (-b WorkingDirectory) and redo only the analysis.

# OrthoFinder inputs: one FASTA per species, for all families together (of_in/pooled) and per family
# (of_in/fam/STxxx), and the species tree pruned to the species of each input folder (<folder>.sp.nwk, for -s).
step "OrthoFinder: proteomes and species trees"
$PY $S/make_proteomes.py --genes st/ldo_cur/genes.tsv --seqs st/seqs --out of_in > logs/proteomes.log
$PY $S/prune_species_tree.py st/species_tree.nwk of_in/pooled of_in/fam/ST* > logs/prune.log

# Score OrthoFinder results folders:
#   $1 key of the eval files (of_<variant>)   $2... Results_* folders (their predicted pairs are merged)
# writes eval/<key>_curated.tsv (curated-tree pairs) and eval/<key>_st006.tsv (ST006 against the curators' list).
of_eval() {
  local key=$1; shift
  rm -f eval/${key}_curated.tsv eval/${key}_st006.tsv
  $PY $S/evaluate_orthofinder.py --truth st/ldo_cur/orthologs.tsv --genes st/ldo_cur/genes.tsv --run "$@" \
      --label "$key/curated" --per-family eval/${key}_curated.tsv >> logs/eval_of.log
  $PY $S/evaluate_orthofinder.py --truth st/curator_orthologs.tsv --genes st/ldo_cur/genes.tsv --run "$@" \
      --families ST006 --label "$key/st006" --per-family eval/${key}_st006.tsv >> logs/eval_of.log
}

# Pooled, standard: DIAMOND search of all 9 families together and the analysis with the true species tree.
step "OrthoFinder: pooled run (DIAMOND all-vs-all + standard analysis, ~7 min)"
mkdir -p of_out          # OrthoFinder needs the parent of -o to exist (and exits 0 when it does not)
$OF -f of_in/pooled -o of_out/pooled -s of_in/pooled.sp.nwk -X -t $T -a $T > logs/of_pooled.log 2>&1
R=$(ls -d of_out/pooled/Results_* | head -1)
W=$R/WorkingDirectory

# Pooled, patched (steps 1-3 above), from the saved search: the fit log (-og: stop after the orthogroups), the
# median parameters, then the patched analysis.
step "OrthoFinder: fit the normalisation parameters, then the patched run"
OF_NORM_LOG=$OUT/of_out/normfit.tsv $OF -b $W -s of_in/pooled.sp.nwk -n fitlog -X -t $T -a $T -og > logs/of_pooled_fitlog.log 2>&1
NORM=$($PY $S/fit_norm_params.py of_out/normfit.tsv --min-hits 30)
echo "median fitted parameters (b,a) over species pairs with >= 30 hits: $NORM" | tee logs/norm_params.txt
OF_NORM_MIN=30 OF_NORM_PARAMS=$NORM $OF -b $W -s of_in/pooled.sp.nwk -n fixnorm30 -X -t $T -a $T > logs/of_pooled_fixnorm30.log 2>&1

# Score the two pooled analyses.
of_eval of_pooled         $R
of_eval of_pooled_fixnorm $W/OrthoFinder/Results_fixnorm30

# Optional (PERFAM=1, about 15 extra minutes): one OrthoFinder run per family instead of the 9 families pooled.
# Purpose: a check on the pooled OrthoFinder results, not needed for the main comparison.
#   - In a pooled run, genes of different families can be grouped together; one run per family rules this out.
#   - Each family alone has few genes per species pair, which is where OrthoFinder's standard score normalisation
#     loses the most hits (see above); the patched run shows that the fixed parameters repair this.
# For each family: (1) the standard run, DIAMOND search + analysis with the family's pruned true species tree;
# (2) the patched analysis from that search, with the parameters $NORM fitted on the pooled run.
# Output: the rows "OrthoFinder, per family, ..." of summary.md (swisstree/expected_summary.md was made with PERFAM=1;
# without it, those two rows are missing and everything else is the same).
if [ "${PERFAM:-0}" = 1 ]; then
  step "OrthoFinder: one run per family"
  for d in of_in/fam/ST*/; do
    f=$(basename $d)
    $OF -f $d -o of_out/fam_$f -s of_in/fam/$f.sp.nwk -X -t $T -a $T > logs/of_fam_$f.log 2>&1
    Wf=$(ls -d of_out/fam_$f/Results_*/WorkingDirectory | head -1)
    OF_NORM_MIN=30 OF_NORM_PARAMS=$NORM $OF -b $Wf -s of_in/fam/$f.sp.nwk -n fixnorm30 -X -t $T -a $T \
        > logs/of_fam_${f}_fixnorm30.log 2>&1
  done
  # Score the per-family runs, all families together.
  of_eval of_perfam         of_out/fam_ST*/Results_*
  of_eval of_perfam_fixnorm of_out/fam_ST*/Results_*/WorkingDirectory/OrthoFinder/Results_fixnorm30
fi

# ---------------------------------------------------------------- 5. report
# Remove the temporary scoring files and write summary.md from the eval/ files.
step "summary"
rm -rf eval/tmp eval/tmp_st006.txt
$PY $S/summarize.py --eval eval --out summary.md > /dev/null
echo "done: $OUT/summary.md"
