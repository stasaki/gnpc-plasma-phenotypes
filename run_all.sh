#!/usr/bin/env bash
# Run every analysis step in order, from the repository root:
#   bash run_all.sh            # all steps, including CatBoost training (needs a GPU)
#   TRAIN=0 bash run_all.sh    # skip training; use runs already in ../output/04_classifier_runs
#   CLASSIFIER_RUNS=/path/to/runs TRAIN=0 bash run_all.sh   # skip training; read runs from another folder
#
# Logs: ../output/logs/<step>.log; a summary table is printed at the end.
# A failing step is reported and the remaining steps still run.
#
# Fig. 6c/6d use SpeakEasy clusters, which have no fixed seed. To reproduce the published
# clusters, put the cached SE_res_*.rds files in ../data/risk_profile_clusters/; they are
# copied into ../output/04_risk_profiles/ before that step (without them, new clusters are
# computed).
set -u
cd "$(dirname "$0")"
LOG=../output/logs
mkdir -p "$LOG"
TRAIN=${TRAIN:-1}
if [ "$TRAIN" = "1" ] && [ -n "${CLASSIFIER_RUNS:-}" ]; then
  echo "CLASSIFIER_RUNS is set: use TRAIN=0 with it (training would write into that folder)"; exit 1
fi
SUMMARY=()

run() {   # run <name> <command...>
  local name=$1; shift
  local t0=$(date +%s)
  echo ">>> $name"
  if "$@" > "$LOG/$name.log" 2>&1; then status=OK; else status=FAILED; fi
  SUMMARY+=("$(printf '%-36s %-7s %6ds' "$name" "$status" $(( $(date +%s) - t0 )))")
  echo "    $status (log: $LOG/$name.log)"
}
rmd() { Rscript -e "rmarkdown::render('$1', knit_root_dir = getwd(), output_dir = '../output/$2')"; }

run 01_validation                 rmd 01_validation/01_validation.Rmd 01_validation
run 02_01_ad_assoc_by_site        rmd 02_ad_association/01_ad_assoc_by_site.Rmd 02_ad_assoc_by_site
run 02_02_ad_assoc_site_adjusted  rmd 02_ad_association/02_ad_assoc_site_adjusted.Rmd 02_ad_assoc_site_adjusted
run 02_03_apoe_adjusted_assoc     Rscript 02_ad_association/03_apoe_adjusted_assoc.R
run 02_04_ad_assoc_by_age         Rscript 02_ad_association/04_ad_assoc_by_age.R
run 03_01_age_trajectory          Rscript 03_trajectories/01_age_trajectory.R
run 03_02_age_trajectory_protein  Rscript 03_trajectories/02_age_trajectory_protein.R
run 03_03_table_s8                Rscript 03_trajectories/03_table_s8.R
if [ "$TRAIN" = "1" ]; then
  run 04_01_train_classifiers     Rscript 04_ad_classifier/01_train_classifiers.R
else
  echo ">>> 04_01_train_classifiers skipped (TRAIN=0)"
fi
run 04_02_collect_results         Rscript 04_ad_classifier/02_collect_results.R
run 04_03_compare_model_perf      Rscript 04_ad_classifier/03_compare_model_perf.R
run 04_04_variable_contributions  Rscript 04_ad_classifier/04_variable_contributions.R
if ls ../data/risk_profile_clusters/SE_res_*.rds > /dev/null 2>&1; then
  mkdir -p ../output/04_risk_profiles
  cp -n ../data/risk_profile_clusters/SE_res_*.rds ../output/04_risk_profiles/
  chmod u+w ../output/04_risk_profiles/SE_res_*.rds   # the step re-saves SE_res_ad_nci.rds
  echo "    using cached clusters from ../data/risk_profile_clusters/"
fi
run 04_05_risk_profiles           Rscript 04_ad_classifier/05_risk_profiles.R
run 04_06_fig6d_individuals       Rscript 04_ad_classifier/06_fig6d_individuals.R
run 04_07_predict_cdr_progression Rscript 04_ad_classifier/07_predict_cdr_progression.R
run 05_01_pdftd_specificity       Rscript 05_ad_specificity/01_pdftd_specificity.R
run 05_02_ad_classifier_oof       Rscript 05_ad_specificity/02_ad_classifier_oof_specificity.R
run 05_03_demographic_persite     Rscript 05_ad_specificity/03_demographic_accuracy_persite.R

# record the environment, with the analysis packages loaded so their versions are listed
Rscript -e 'for (p in c("mgcv", "nlme", "Matrix", "catboost", "speakeasyR", "data.table", "dplyr", "ggplot2", "future", "furrr", "lme4", "rmarkdown")) suppressMessages(requireNamespace(p, quietly = TRUE)); writeLines(capture.output(sessionInfo()), "../output/logs/sessionInfo.txt")'

echo; echo "Summary"; printf '%s\n' "${SUMMARY[@]}" | tee "$LOG/_summary.txt"
