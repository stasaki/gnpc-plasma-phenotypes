# Settings shared by every script.
#
# Expected layout (paths are relative to this repository, so nothing needs editing):
#   <project>/
#     code/     this repository
#     data/     ProteomicsV1_3ms.rds, predicted_phenotypes.rds, excluded_contributors.txt
#     output/   written by the scripts, one subfolder per step (see step_dir())
#
# Run every script from the repository root, for example
#   Rscript 03_trajectories/01_age_trajectory.R
#   Rscript -e 'rmarkdown::render("01_validation/01_validation.Rmd", knit_root_dir = getwd())'

# ---- inputs (all stay on the AD Workbench) ----------------------------------
# Harmonized GNPC V1.3 SomaScan data (SummarizedExperiment)
PROTEOMICS_RDS <- normalizePath("../data/ProteomicsV1_3ms.rds", mustWork = FALSE)
# Predicted phenotypes + annotation tables (see R/load_inputs.R)
PREDICTED_RDS  <- normalizePath("../data/predicted_phenotypes.rds", mustWork = FALSE)
# GNPC contributor codes left out of every analysis (see Methods), one per line
EXCLUDED_CONTRIBUTORS <- trimws(readLines("../data/excluded_contributors.txt"))
EXCLUDED_CONTRIBUTORS <- EXCLUDED_CONTRIBUTORS[nzchar(EXCLUDED_CONTRIBUTORS)]
# Tag in the classifier run folder names: "ex" + the excluded codes
EX_TAG <- paste0("ex", paste(EXCLUDED_CONTRIBUTORS, collapse = ""))

# ---- outputs -------------------------------------------------------------------
# The classifier summaries select runs by matching words in their file paths, so this
# path must not contain EX_TAG, "p11", "clean" or "cross_validation".
dir.create("../output", showWarnings = FALSE)
OUTPUT_DIR <- normalizePath("../output")

# Parallel workers for the trajectory fits: the cores this R session may actually use
# (availableCores() respects container/cgroup limits; detectCores() counts the whole
# machine), capped at 8 because R allows at most 128 connections. The number of
# workers does not change the results.
N_CORES <- min(8L, as.integer(parallelly::availableCores()))

source("R/load_inputs.R")

# Optional: read already trained classifiers from another folder instead of
# OUTPUT_DIR/04_classifier_runs, e.g. to rerun the later steps without retraining:
#   CLASSIFIER_RUNS=/path/to/runs TRAIN=0 bash run_all.sh
# That folder is only read; 04_ad_classifier/01_train_classifiers.R refuses to run while it is set.
CLASSIFIER_RUNS <- Sys.getenv("CLASSIFIER_RUNS")

# Output folder of one analysis step, created on first use.
step_dir <- function(step) {
  if (step == "04_classifier_runs" && nzchar(CLASSIFIER_RUNS)) return(normalizePath(CLASSIFIER_RUNS, mustWork = TRUE))
  d <- file.path(OUTPUT_DIR, step)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
