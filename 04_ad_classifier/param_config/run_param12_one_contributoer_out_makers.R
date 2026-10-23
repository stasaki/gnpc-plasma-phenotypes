# Classifier config run_param12 (output folder p12_one_contributoer_out_markers_<EX_TAG>), read by
# 04_ad_classifier/01_train_classifiers.R via run_pipeline(). Runs on the AD Workbench.
# Predictors: SHAP-selected proteins (top 5-100).
# Samples: GNPC plasma, excluded contributors (config.R) removed, AD vs all non-AD participants;
# folds = contributor (leave-one-contributor-out).
# Inputs: PROTEOMICS_RDS, PREDICTED_RDS (config.R).

# ---------------------------------------------------------------------------------
# Script Configuration
# ---------------------------------------------------------------------------------

# General Configurations
PARAM = list()
PARAM$label = paste0("one_contributoer_out_markers_", EX_TAG)

PARAM$combn = c("trait","top_n")
PARAM$OVERWRITE = FALSE

PARAM$trait = list()
PARAM$trait[["ad"]] = list(trait = "ad",
                              param_loc = "04_ad_classifier/param_analysis/ml_categorical/")
PARAM$top_n = c(100,50,25,10,5)

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
  
  # Impute missing values with the column (protein) mean
  for(i in 1:ncol(prediction)) {
    prediction[is.na(prediction[,i]), i] <- mean(prediction[,i], na.rm = TRUE)
  }
  

  if(!identical(metadata$IID,rownames(prediction))){
    stop()
  }
  
  top_markers = get_table(pp, "top_markers")[[paste0("top_",PARAM$top_n)]]
  indx = colnames(prediction)%in%top_markers
  data <- metadata %>%
    cbind(.,prediction[,indx])


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
