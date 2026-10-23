# Classifier config run_param8 (output folder p8_one_contributoer_out_nesting_<EX_TAG>), read by
# 04_ad_classifier/01_train_classifiers.R via run_pipeline(). Runs on the AD Workbench.
# Predictors: predicted phenotypes (all, or one domain), filtered by ROSMAP performance (metric_th).
# Samples: GNPC plasma, excluded contributors (config.R) removed, AD vs all non-AD participants;
# folds = contributor (leave-one-contributor-out).
# Inputs: PROTEOMICS_RDS, PREDICTED_RDS (config.R).

# ---------------------------------------------------------------------------------
# Script Configuration
# ---------------------------------------------------------------------------------

# General Configurations
PARAM = list()
PARAM$label = paste0("one_contributoer_out_nesting_", EX_TAG)

PARAM$combn = c("trait","feature","metric_th")
PARAM$OVERWRITE = FALSE

PARAM$trait = list()
PARAM$trait[["ad"]] = list(trait = "ad",
                           param_loc = "04_ad_classifier/param_analysis/ml_categorical/")
PARAM$feature = c("all_1st_tier","Clinical.Diagnosis","Cognition",
  "Demographics", "Genetics", "Medical.Conditions",
  "Medications","Motor.and.Gait","Pathology","ex.Clinical.Diagnosis")
PARAM$metric_th = c(0,0.1,0.01,0.2,0.3)
# ---------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------


PARAM$prep_data <- function(PARAM) {
  library(SummarizedExperiment)
  metadata <- colData(readRDS(PROTEOMICS_RDS)) %>% as.data.frame()
  pp <- read_predicted_phenotypes(PREDICTED_RDS)
  prediction <- get_prediction(pp, "no_class_labels")
  
  
  feature_annot = get_table(pp, "feature_groups")
  feature_annot = feature_annot%>%group_by(Category,group)%>%summarise(target_id = list(target_id))%>%filter(group!="exclude")
  feature_annot = bind_rows(tibble(Category="ex.Clinical.Diagnosis",
                                   group = "",
                                   target_id = list(feature_annot%>%filter(group=="")%>%
                                                      filter(Category!="Clinical Diagnosis")%>%pull(target_id)%>%unlist()%>%
                                                      unique())),
                            feature_annot)
  feature_annot = bind_rows(tibble(Category="all_1st_tier",
                                   group = "",
                                   target_id = list(feature_annot%>%filter(group=="")%>%pull(target_id)%>%unlist()%>%
                                                      unique())),
                            feature_annot)
  
  
  names(feature_annot$target_id) = paste0(feature_annot$Category,"__",feature_annot$group)%>%make.names()%>%gsub("__$","",.)
  
  feature_use = feature_annot$target_id[names(feature_annot$target_id)==PARAM$feature]%>%unlist()
  
  # Keep predicted phenotypes whose median ROSMAP performance
  # (weighted_n_auprc for classification, Pearson r for regression) is above metric_th
  get_table(pp, "rosmap_classification_performance")%>%
    group_by(`Target ID`)%>%
    summarise(metric=median(weighted_n_auprc))%>%
    bind_rows(.,
              get_table(pp, "rosmap_regression_performance")%>%
                group_by(`Target ID`)%>%
                summarise(metric=median(Pearson)))%>%
    filter(metric> PARAM$metric_th)%>%pull(`Target ID`) -> feature_passed
  
  indx = feature_use%>%gsub("_.+","",.) %in%  feature_passed
  feature_use = feature_use[indx]
  
  prediction = prediction[,c("sample_id",feature_use)]
  
  data <- metadata %>%
    filter(!sample_matrix %in% c("CSF", "Serum"))%>%
    filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)%>%
    dplyr::select(IID=sample_id,trait = ad,fold=contributor_code)%>%
    filter(trait %in%c(0, 1))%>%
    inner_join(prediction, by = c("IID"="sample_id"))
  
  
  # Check the balance
  table(data$fold, data$trait)
  
  data=as.data.frame(data)
  class_labels <- c("0" = "NCI", "1" = "MCI", "2" = "dementia")
  return(data)
}


# ---------------------------------------------------------------------------------
# Pipeline flag
# ---------------------------------------------------------------------------------
PARAM$ml = list()
PARAM$ml$flag = TRUE
