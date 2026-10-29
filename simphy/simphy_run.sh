#!/bin/bash
# LDO reconstruction and OrthoFinder on simphy simulations: label the gene trees with simphy's events, predict LDO graphs
# from the sequences, run ldo_builder.py, run OrthoFinder, and score both on ortholog pairs (and LDO on clusters).
#
#   source env/env.sh                         # PY, OF (patched OrthoFinder) - see README
#   bash simphy_run.sh SIMPHY_ROOT OUT
#
# SIMPHY_ROOT is searched for replicate folders (those containing s_tree.trees, e.g.
# simpphy_p10_dl1_t0_s300/ssim_..._seed3001/1). Results go to OUT/<setting>/<seedNNNN>/, the report to OUT/summary.md.
#
# Environment variables:
#   OF_MODES   OrthoFinder runs to do, any of: pooled perfam perfam_fixnorm   (default: "pooled"; "" = no OrthoFinder)
#              pooled         one run per replicate with all its families (gene ids are made unique per family)
#              perfam         one run per family (unpatched normalisation)
#              perfam_fixnorm one run per family, species pairs with < 30 hits normalised with the median parameters
#                             fitted in the pooled run of the same replicate (needs "pooled" too)
#   PERFAM_FAMILIES  families for the per-family runs, e.g. 1-10 (default: all; ~15 s per family per variant)
#   JOBS       parallel per-family OrthoFinder runs (default 2), THREADS threads for the pooled run (default 2)
#   REPAIRS    ldo_builder.py --repair settings (default "none count count+sparse")
#   COUNT_AGE  also run count+sparse with --count-age COUNT_AGE (default 1.0; "" = skip), reported as count+sparse+age
#   ONLY       regular expression: only replicates whose path matches it
set -euo pipefail

# Arguments, programs and defaults of the environment variables above (S = the scripts folder).
ROOT=$(realpath "$1")
OUT=$(realpath -m "${2:-simphy_out}")
S=$(cd "$(dirname "$0")/.." && pwd)/scripts
PY=${PY:-python3}; OF=${OF:-orthofinder}
OF_MODES=${OF_MODES-pooled}; JOBS=${JOBS:-2}; T=${THREADS:-2}
REPAIRS=${REPAIRS:-"none count count+sparse"}
COUNT_AGE=${COUNT_AGE-1.0}
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$OUT"

# The replicates: every folder under ROOT that contains s_tree.trees (filtered by ONLY), in sorted order.
mapfile -t RUNS < <(find "$ROOT" -name s_tree.trees -printf '%h\n' | sort | grep -E "${ONLY:-.}")
echo "${#RUNS[@]} replicates; OrthoFinder: ${OF_MODES:-none}"

for RUN in "${RUNS[@]}"; do
  # Output folder of the replicate: OUT/<setting folder>/<seedNNNN>.
  rel=${RUN#$ROOT/}
  setting=$(echo "$rel" | cut -d/ -f1)
  seed=$(echo "$rel" | grep -o 'seed[0-9]*' | head -1)
  W="$OUT/$setting/$seed"
  # A replicate finished in an earlier run (it has a 'done' file) is skipped; an unfinished one is redone from scratch.
  if [ -e "$W/done" ]; then echo "skip $setting/$seed (done)"; continue; fi
  echo "== $(date +%H:%M:%S) $setting/$seed"
  rm -rf "$W"; mkdir -p "$W/eval" "$W/logs"
  # one replicate, in a background subshell so that set -e applies inside it (bash ignores set -e in a
  # subshell that is part of a || list); a failure is reported and the loop goes on
  (
  set -e
  cd "$W"

  # 1. labelling (simphy's own events) + inputs: labelled JSON trees with the alignments (trees/), species tree, OrthoFinder inputs (of_in/)
  $PY $S/simphy_prepare.py --run "$RUN" --out . > logs/prepare.log
  # 2. truth: ortholog pairs of the labelled trees (ldo/orthologs.tsv) and genes.tsv
  $PY $S/true_orthologs.py --trees trees --out ldo > logs/true_orthologs.log 2>&1
  # 3. predicted LDO graph (Jukes-Cantor distance, DNA)
  $PY $S/predict_ldo.py --trees trees --out pred --model jc > logs/predict.log 2>&1
  # 4. reconstruction + scoring, for each repair setting: trees (ldo_trees_<r>.txt), ortholog pairs (eval/ldo_<r>.tsv)
  #    and clusters (eval/clusters_<r>.tsv)
  for r in $REPAIRS; do
    $PY $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
        --repair $r --trees-out ldo_trees_$r.txt > logs/builder_$r.log 2>&1
    $PY $S/evaluate_orthologs.py --trees ldo_trees_$r.txt --truth ldo/orthologs.tsv --out eval/tmp > logs/eval_ldo_$r.log
    cp eval/tmp/per_family.tsv eval/ldo_$r.tsv
    $PY $S/cluster_analysis.py --pred ldo_trees_$r.txt --true-trees trees --out eval/clusters_$r.tsv > logs/clusters_$r.log
  done
  # count+sparse with the duplication-age test of the count repair (--count-age, uses the alignments in trees/),
  # scored the same way
  if [ -n "$COUNT_AGE" ]; then
    r=count+sparse+age
    $PY $S/ldo_builder.py --graph pred/ldo_edges.tsv --genes pred/genes.tsv --species-tree species_tree.nwk \
        --repair count+sparse --count-age $COUNT_AGE --seqs-trees trees --model jc --trees-out ldo_trees_$r.txt > logs/builder_$r.log 2>&1
    $PY $S/evaluate_orthologs.py --trees ldo_trees_$r.txt --truth ldo/orthologs.tsv --out eval/tmp > logs/eval_ldo_$r.log
    cp eval/tmp/per_family.tsv eval/ldo_$r.tsv
    $PY $S/cluster_analysis.py --pred ldo_trees_$r.txt --true-trees trees --out eval/clusters_$r.tsv > logs/clusters_$r.log
  fi
  rm -rf eval/tmp

  # 5. OrthoFinder (DNA: -d -S blastn; true species tree -s)
  # How OrthoFinder is used, and why the comparison is as fair as we could make it:
  #   * Same genes: the sequences of the LDO pipeline (the simphy alignments without gaps), one FASTA per species.
  #   * Same species tree: OrthoFinder is given the true species tree (-s), as the LDO reconstruction is.
  #   * DNA search: the sequences are DNA, so OrthoFinder searches with blastn (-d -S blastn; DIAMOND is protein-only).
  #   * Its own best output: we score OrthoFinder's ortholog pairs (its Orthologues/ folder), not its orthogroups.
  #   * pooled: one run per replicate with all its 100 families (gene ids made unique per family). This is standard
  #     OrthoFinder (OF_NORM_LOG only logs the normalisation fits, it changes nothing) and the baseline we
  #     recommend: with 100 families every species pair has many hits.
  #   * perfam / perfam_fixnorm: one run per family (the same family knowledge as the LDO pipeline). Here OrthoFinder's
  #     score normalisation hurts it: a species pair with fewer than 2 hits cannot be fitted and its hits are
  #     discarded, and in a small family most pairs have only a few hits (a single-copy family loses all of them).
  #     perfam_fixnorm runs the opt-in patch (orthofinder_waterfall_norm.patch): species pairs with < 30 hits use the
  #     median line fitted on the pairs with >= 30 hits of the replicate's pooled run (fit_norm_params.py); pairs with
  #     >= 30 hits are unchanged. The parameters come from OrthoFinder's own fits on these sequences; no truth is used.
  # Effect, on the 249 families with per-family runs (families 1-10 of each replicate; precision / recall / F1):
  #   pooled           0.573 / 0.711 / 0.634
  #   perfam           0.470 / 0.557 / 0.510
  #   perfam_fixnorm   0.514 / 0.614 / 0.560   (recall +0.057; single-copy dl0 families: recall 0.000 -> 0.643)
  # All three are reported; the pooled run is the OrthoFinder baseline we compare with.
  mkdir -p of_out
  NORM=""
  # Pooled run: all the families of the replicate in one run. OF_NORM_LOG makes the (patched) OrthoFinder log the
  # normalisation fit of every species pair; score it, then keep the median parameters of the pairs with >= 30 hits
  # (NORM, also saved in of_out/norm_params.txt) for the perfam_fixnorm runs.
  if [[ " $OF_MODES " == *" pooled "* ]]; then
    rm -f of_out/normfit.tsv
    OF_NORM_LOG=$W/of_out/normfit.tsv $OF -f of_in/pooled -o of_out/pooled -d -S blastn -s of_in/pooled.sp.nwk \
        -X -t $T -a $T > logs/of_pooled.log 2>&1
    $PY $S/evaluate_orthofinder.py --truth ldo/orthologs.tsv --genes ldo/genes.tsv --run of_out/pooled/Results_* \
        --label pooled --per-family eval/of_pooled.tsv > logs/eval_of_pooled.log
    NORM=$($PY $S/fit_norm_params.py of_out/normfit.tsv --min-hits 30 2>/dev/null || true)
    echo "$NORM" > of_out/norm_params.txt
  fi
  # Per-family runs (perfam, perfam_fixnorm): one OrthoFinder run per family.
  if [[ " $OF_MODES " == *" perfam"* ]]; then
    # The families to run: those of PERFAM_FAMILIES (e.g. "1-10,15"; empty = all) that have >= 2 species.
    fams=$($PY - "$PERFAM_FAMILIES" <<'EOF'
import os, sys
spec = sys.argv[1]
fams = sorted(d for d in os.listdir("of_in/fam") if os.path.isdir(os.path.join("of_in/fam", d)))
# OrthoFinder needs >= 2 species; a one-species family has no ortholog pairs to score anyway
fams = [d for d in fams if len([x for x in os.listdir(os.path.join("of_in/fam", d)) if x.endswith(".fa")]) >= 2]
if spec:
    want = set()
    for part in spec.split(","):
        a, _, b = part.partition("-")
        want |= set(range(int(a), int(b or a) + 1))
    fams = [f for f in fams if int(f) in want]
print(" ".join(fams))
EOF
)
    export OF W
    for variant in perfam perfam_fixnorm; do
      [[ " $OF_MODES " == *" $variant "* ]] || continue
      if [ $variant = perfam_fixnorm ] && [ -z "$NORM" ]; then
        echo "  perfam_fixnorm skipped: needs the pooled run (fitted parameters)"; continue
      fi
      export VARIANT=$variant NORM
      # JOBS families at a time (xargs -P), one thread each, into of_out/<variant>_<family>; with perfam_fixnorm
      # the pairs with < 30 hits use the pooled parameters NORM. A failed family run does not stop the others (it
      # is counted as no prediction by the scoring).
      echo $fams | tr ' ' '\n' | xargs -P $JOBS -I{} bash -c '
        f={}; o=$W/of_out/${VARIANT}_$f
        if [ $VARIANT = perfam_fixnorm ]; then export OF_NORM_MIN=30 OF_NORM_PARAMS=$NORM; fi
        $OF -f $W/of_in/fam/$f -o $o -d -S blastn -s $W/of_in/fam/$f.sp.nwk -X -t 1 -a 1 > $o.log 2>&1 || true'
      # Score the variant on the families run, their pairs pooled.
      runs=$(ls -d of_out/${variant}_*/Results_* 2>/dev/null || true)
      if [ -n "$runs" ]; then
        $PY $S/evaluate_orthofinder.py --truth ldo/orthologs.tsv --genes ldo/genes.tsv --run $runs \
            --families $fams --label $variant --per-family eval/of_$variant.tsv > logs/eval_of_$variant.log
      fi
    done
  fi
  # The replicate is complete.
  touch done
  ) &
  # Wait for the replicate; on failure note it in OUT/failed.txt and go on with the next one.
  if ! wait $!; then echo "  FAILED $setting/$seed (see $W/logs)"; echo "$setting/$seed" >> "$OUT/failed.txt"; fi
done

# Report over all the replicates (summary.md, all_families.tsv), and the list of failed replicates if any.
$PY $S/simphy_summarize.py --out "$OUT" > /dev/null
[ -s "$OUT/failed.txt" ] && echo "failed replicates (not in the summary): $(tr '\n' ' ' < "$OUT/failed.txt")"
echo "done: $OUT/summary.md"
