# Classifier config run_param10 (output folder p10_one_contributoer_out_pca_<EX_TAG>), read by
# 04_ad_classifier/01_train_classifiers.R via run_pipeline(). Runs on the AD Workbench.
# Predictors: first 120 principal components of all proteins.
# Samples: GNPC plasma, excluded contributors (config.R) removed, AD vs all non-AD participants;
# folds = contributor (leave-one-contributor-out).
# Inputs: PROTEOMICS_RDS, PREDICTED_RDS (config.R).

# ---------------------------------------------------------------------------------
# Script Configuration
# ---------------------------------------------------------------------------------

# General Configurations
PARAM = list()
PARAM$label = paste0("one_contributoer_out_pca_", EX_TAG)

PARAM$combn = c("trait")
PARAM$OVERWRITE = FALSE

PARAM$trait = list()
PARAM$trait[["ad"]] = list(trait = "ad",
                           param_loc = "04_ad_classifier/param_analysis/ml_categorical/")

# ---------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------


PARAM$prep_data <- function(PARAM) {
  
  library(SummarizedExperiment)
  pp <- read_predicted_phenotypes(PREDICTED_RDS)
  prediction_stacking <- get_prediction(pp, "no_class_labels")
  prediction <- readRDS(PROTEOMICS_RDS)
  metadata <- colData(prediction) %>% as.data.frame()
  
  metadata = metadata %>%
    filter(sample_id%in%prediction_stacking$sample_id)%>%
    filter(!sample_matrix %in% c("CSF", "Serum"))%>%
    filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)%>%
    dplyr::select(IID=sample_id,trait = ad,fold=contributor_code)%>%
    filter(trait %in%c(0, 1))
  
  
  prediction=prediction[,prediction$sample_id%in%metadata$IID]
  
  
  indx <- rowSums(assays(prediction)[["data"]] < 0) < 100
  prediction <- prediction[indx, ]
  
  # Replace negative values with NA and log-transform
  assays(prediction)[["data"]][assays(prediction)[["data"]] < 0] <- NA
  assays(prediction)[["data"]] <- log(assays(prediction)[["data"]])
  
  prediction = assay(prediction,"data")
  prediction = as.matrix(prediction)%>%t()
  
  # Impute missing values with the column (protein) mean; PCA cannot handle missing values
  for(i in 1:ncol(prediction)) {
    prediction[is.na(prediction[,i]), i] <- mean(prediction[,i], na.rm = TRUE)
  }
  
  # Zero-variance columns (a check only; nothing is removed)
  zero_var_cols <- which(apply(prediction, 2, function(x) sd(x) == 0))
  length(zero_var_cols)        # Number of zero-variance columns
  colnames(prediction)[zero_var_cols]  # Names of zero-variance columns
  
  
  # Run PCA
  # Each protein is centered and scaled before PCA
  library(irlba)
  
  # Apply PCA using irlba (much faster for large matrices when you only need some PCs)
  pca_result <- prcomp_irlba(prediction, n = 120, center = TRUE, scale = TRUE)
  
  # Extract the first 120 principal components
  # Check if we have at least 120 components available
  n_components <- min(120, ncol(pca_result$x))
  pc_scores <- pca_result$x[, 1:n_components]
  rownames(pc_scores) = rownames(prediction)
  
  
  
  
  if(!identical(metadata$IID,rownames(pc_scores))){
    stop()
  }
  
  data <- metadata %>%
    cbind(.,pc_scores)
  
  
  data = data%>%arrange(IID)
  
  # Check the balance
  table(data$fold, data$trait)
  data=as.data.frame(data)
  return(data)
}


# ---------------------------------------------------------------------------------
# Pipeline flag
# ---------------------------------------------------------------------------------
PARAM$ml = list()
PARAM$ml$flag = TRUE
