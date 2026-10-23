# Model setting for ml_session() (R/ml_session.R): binary CatBoost classifier of "trait",
# one model per held-out fold ("fold" = contributor), with class weights.
# Read by run_pipeline() via param_loc in 04_ad_classifier/param_config/run_param*.R.

PARAM$combn = c("setting")
PARAM <- list(
  label="trait",
  target = "trait",
  model_type = "binary", #multiclass,binary
  test_fold = NULL,
  valid_ratio = 0.2,
  fold_col = "fold",
  setting = c("default_with_weights"), # other options: "tuning_with_weights", "default_no_weights", "tuning_no_weights"
  results_dir = NULL
)

PARAM$select_samples <- function(expdata,PARAM){
  
  return(expdata)
}
