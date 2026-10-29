#!/bin/bash
# LDO reconstruction on Ensembl Compara gene trees: prepare the downloaded trees, compute the true orthologs (with the
# dubious nodes read as duplications, then as speciations), predict LDO graphs from the protein alignments, run
# ldo_builder.py, run OrthoFinder, and score both on ortholog pairs (and the reconstructed trees on clusters).
#
#   source env/env.sh                                   # PY - see README
#   python scripts/ensembl_download.py --out DATA ...   # needs network access to Ensembl (see README)
#   bash ensembl/ensembl_run.sh DATA OUT
#
# DATA is the folder written by ensembl_download.py: DATA/trees/*.json and DATA/compara_species_tree.nh.
# Results go to OUT; the report to OUT/summary.md.
#
# Environment variables:
#   SPECIES_NEWICK  Compara species tree (default DATA/compara_species_tree.nh)
#   MAX_GENES       families with more genes are skipped (default 1500)
#   REPAIRS         ldo_builder.py --repair settings (default "none count count+sparse")
#   COUNT_AGE       also run count+sparse with --count-age COUNT_AGE (default 1.0; "" = skip), reported as
#                   count+sparse+age
#   ORTHOFINDER     1 (default) = also run OrthoFinder (pooled: all families in one run, standard and patched), as
#                   for SwissTree; 0 = skip. OF and THREADS (default 2) as in the other pipelines.
#
# Dubious nodes: Ensembl Compara marks some duplication nodes as "dubious" (low duplication-consistency score: the
# two subtrees share few species, so the node may well be a speciation in a tree with errors). The predictions do not
# depend on them (they are made from the sequences only); only the truth does. Every evaluation is therefore done
# twice, with the dubious nodes read as Dup (tag "dup") and as Spec (tag "spec"), and the report also gives the
# scores on the families without dubious nodes, where both readings agree.
set -euo pipefail

# Arguments, programs and defaults of the environment variables above (S = the scripts folder).
DATA=$(realpath "$1")
OUT=$(realpath -m "${2:-ensembl_out}")
S=$(cd "$(dirname "$0")/.." && pwd)/scripts
PY=${PY:-python3}; OF=${OF:-orthofinder}; T=${THREADS:-2}
SPECIES_NEWICK=${SPECIES_NEWICK:-$DATA/compara_species_tree.nh}
REPAIRS=${REPAIRS:-"none count count+sparse"}
COUNT_AGE=${COUNT_AGE-1.0}
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$OUT/eval" "$OUT/logs"
cd "$OUT"

# Print a step title with the time.
step() { echo "== $(date +%H:%M:%S) $*"; }

# 1. Prepared trees (OUT/trees: leaves with aligned sequence and a taxon of the species tree), species tree with
#    taxon-id leaves (OUT/species_tree.nwk), per-family report (OUT/prepare_report.tsv).
step "prepare"
$PY $S/ensembl_prepare.py --trees "$DATA/trees" --species-newick "$SPECIES_NEWICK" --out . \
    --max-genes "${MAX_GENES:-1500}" | tee logs/prepare.log

# 2. Truth: the ortholog pairs of the Ensembl trees, dubious nodes read as Dup (truth_dup/) and as Spec (truth_spec/).
step "true orthologs"
for dub in dup spec; do
  $PY $S/true_orthologs.py --trees trees --out truth_$dub --dubious $dub > logs/true_orthologs_$dub.log 2>&1
done

# 3. Predicted LDO graph from the protein alignments (Poisson distance).
step "predict LDO graphs"
$PY $S/predict_ldo.py --trees trees --out pred > logs/predict.log 2>&1

# Score one reconstruction under both readings of the dubious nodes:
#   $1 trees file of ldo_builder.py   $2 name of the setting (none, count, ...)
# writes eval/ldo_<setting>_<dub>.tsv (per family, evaluate_orthologs.py) and eval/clusters_<setting>_<dub>.tsv / .txt
# (cluster_analysis.py: per family / printed totals).
evaluate() {
  for dub in dup spec; do
    $PY $S/evaluate_orthologs.py --trees $1 --truth truth_$dub/orthologs.tsv --out eval/tmp > /dev/null
    mv eval/tmp/per_family.tsv eval/ldo_$2_$dub.tsv
    $PY $S/cluster_analysis.py --pred $1 --true-trees trees --dubious $dub --out eval/clusters_$2_$dub.tsv \
        > eval/clusters_$2_$dub.txt
  done
}

# 4. Reconstruction with each repair setting, then scoring.
for r in $REPAIRS; do
  step "ldo_builder.py --repair $r"
  $PY $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
      --repair $r --trees-out trees_$r.txt > logs/builder_$r.log 2>&1
  evaluate trees_$r.txt $r
done
# count+sparse with the duplication-age test of the count repair (alignment distances from OUT/trees).
if [ -n "$COUNT_AGE" ]; then
  step "ldo_builder.py --repair count+sparse --count-age $COUNT_AGE"
  $PY $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
      --repair count+sparse --count-age $COUNT_AGE --seqs-trees trees --trees-out trees_count+sparse+age.txt \
      > logs/builder_count+sparse+age.log 2>&1
  evaluate trees_count+sparse+age.txt count+sparse+age
fi

# 5. OrthoFinder, as in swisstree/reproduce.sh (see the explanation there): the same genes (sequences without gaps),
#    the true species tree (-s), scored on OrthoFinder's ortholog pairs. Pooled: all families in one run (DIAMOND
#    search + standard analysis), then from the saved search the normalisation fit log, the median parameters of the
#    species pairs with >= 30 hits, and the patched analysis (pairs with < 30 hits use those parameters).
if [ "${ORTHOFINDER:-1}" = 1 ]; then
  step "OrthoFinder: proteomes, pooled run (DIAMOND all-vs-all + standard analysis)"
  $PY $S/make_proteomes.py --genes truth_dup/genes.tsv --seqs seqs --out of_in > logs/proteomes.log
  $PY $S/prune_species_tree.py species_tree.nwk of_in/pooled > logs/prune.log
  mkdir -p of_out
  $OF -f of_in/pooled -o of_out/pooled -s of_in/pooled.sp.nwk -X -t $T -a $T > logs/of_pooled.log 2>&1
  R=$(ls -d of_out/pooled/Results_* | head -1)
  W=$R/WorkingDirectory
  step "OrthoFinder: fit the normalisation parameters, then the patched run"
  OF_NORM_LOG=$OUT/of_out/normfit.tsv $OF -b $W -s of_in/pooled.sp.nwk -n fitlog -X -t $T -a $T -og \
      > logs/of_pooled_fitlog.log 2>&1
  NORM=$($PY $S/fit_norm_params.py of_out/normfit.tsv --min-hits 30)
  echo "median fitted parameters (b,a) over species pairs with >= 30 hits: $NORM" | tee logs/norm_params.txt
  OF_NORM_MIN=30 OF_NORM_PARAMS=$NORM $OF -b $W -s of_in/pooled.sp.nwk -n fixnorm30 -X -t $T -a $T \
      > logs/of_pooled_fixnorm30.log 2>&1
  # Score both analyses under both readings of the dubious nodes: eval/of_<variant>_<dub>.tsv.
  for dub in dup spec; do
    $PY $S/evaluate_orthofinder.py --truth truth_$dub/orthologs.tsv --genes truth_$dub/genes.tsv --run $R \
        --label pooled/$dub --per-family eval/of_pooled_$dub.tsv >> logs/eval_of.log
    $PY $S/evaluate_orthofinder.py --truth truth_$dub/orthologs.tsv --genes truth_$dub/genes.tsv \
        --run $W/OrthoFinder/Results_fixnorm30 --label pooled_fixnorm/$dub --per-family eval/of_pooled_fixnorm_$dub.tsv \
        >> logs/eval_of.log
  done
fi

# 6. Report.
step "summary"
rm -rf eval/tmp
$PY $S/ensembl_summarize.py --out . > /dev/null
echo "done: $OUT/summary.md"
