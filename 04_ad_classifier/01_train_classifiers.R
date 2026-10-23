# Train the GNPC AD classifiers (CatBoost, leave-one-contributor-out) used for Fig. 6,
# Extended Data Fig. 8 and 10, and Supplementary Tables 10-12.
# Runs on the AD Workbench (needs the GNPC data).
#   Inputs : PROTEOMICS_RDS, PREDICTED_RDS (read by each config in 04_ad_classifier/param_config/)
#   Configs: run_param8/9/10/12 (predicted phenotypes / all proteins / 120 PCs / SHAP-selected proteins)
#            and run_param14-17 (same, with disease-free controls)
#   Outputs: step_dir("04_classifier_runs")/p<N>_<label>/... (one folder per config)
# Run from the repository root:
#   Rscript 04_ad_classifier/01_train_classifiers.R

source("config.R")
if (nzchar(CLASSIFIER_RUNS)) stop("CLASSIFIER_RUNS is set: unset it before training (training would write into that folder)")
options(stringsAsFactors = FALSE)
source("R/vendor/run_pipeline.R")



library(data.table)
library(tidyverse)
library(caret)
library(parallel)
source("R/vendor/ml_func1.R")
source("R/vendor/catboost_pipeline.R")
source("R/ml_session.R")

mainfunc = function(param_config) {
  save(param_config, file = paste0(param_config$OUTLOC,
                                   "/param_config.RData"))
  expdata = param_config$prep_data(param_config)
  ml_session(expdata, param_config)
}


  

  param_config = run_pipeline(mainfunc, param_location = normalizePath("04_ad_classifier/param_config/"),
                              param_suffix = "run_param", OUTDIR = step_dir("04_classifier_runs"),
                              PARAM = NULL, return_first = F, SGE = F, MCLAPPLY = F,
                              MC.cores = 4, main_script = normalizePath("04_ad_classifier/01_train_classifiers.R"),
                              sge_options = paste0("-wd ", step_dir("04_classifier_runs"),
                                                   " -V -S /bin/bash -pe smp 4"), 
                              run_pnum = c(8,9,10,12,14,15,16,17))

