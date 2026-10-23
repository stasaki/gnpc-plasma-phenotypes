# Loaders for the two inputs of the GNPC analysis code
# ============================================================================
#   ProteomicsV1_3ms.rds       harmonized GNPC V1.3 SomaScan data (SummarizedExperiment)
#   predicted_phenotypes.rds   a list, built on the AD Workbench:
#     $prediction  one row per sample: predicted scores, class probabilities and
#                  predicted class labels from the ROSMAP-trained TabPFN models
#     $tables      the small annotation tables the analysis scripts read
#
# Both files stay on the AD Workbench (GNPC data use agreement).
# ============================================================================

read_predicted_phenotypes <- function(path) readRDS(path)

# The prediction table in the two forms the analysis scripts use.
#   "full"            every column.
#   "no_class_labels" drops each classification trait's predicted class label
#                     and its first class-probability column (the reference
#                     class), leaving one column per continuous trait and the
#                     remaining class probabilities.
# Column names: continuous traits and class labels are bare IDs ("T0001");
# probability columns are "<ID>_<class>", in the class order of the model.
get_prediction <- function(pp, version = c("full", "no_class_labels")) {
  version <- match.arg(version)
  pred <- pp$prediction
  if (version == "full") return(pred)
  prob_cols   <- setdiff(grep("_", colnames(pred), value = TRUE), "sample_id")
  prob_trait  <- sub("_.+", "", prob_cols)
  label_cols  <- intersect(colnames(pred), unique(prob_trait))
  first_class <- prob_cols[!duplicated(prob_trait)]
  pred[, !colnames(pred) %in% c(label_cols, first_class)]
}

# Predicted class labels of one trait (e.g. "T0002" sex, "T0607" APOE genotype)
# as sample_id + prediction, the layout of TabPFN's *_predictions.csv.gz files.
get_class_labels <- function(pp, trait) {
  out <- data.table::as.data.table(pp$prediction[, c("sample_id", trait)])
  data.table::setnames(out, trait, "prediction")
  out
}

get_table <- function(pp, name) {
  if (!name %in% names(pp$tables)) stop("no table '", name, "' in predicted_phenotypes.rds")
  pp$tables[[name]]
}
