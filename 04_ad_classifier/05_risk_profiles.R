# ==============================================================================
# Fig. 6c: AD risk profiles. Clusters individuals by their SHAP vectors from the
# AD classifier (speakeasyR) and tests how each profile differs.
# Inputs : predicted_phenotypes.rds (target_annotation_updated);
#          classifier run p8_one_contributoer_out_nesting_<EX_TAG> in
#          step_dir("04_classifier_runs") (01_train_classifiers.R).
# Outputs: step_dir("04_risk_profiles"). The per-person cluster files SE_res_*.rds
#          are cached there; clustering runs only if they are missing.
# Runs on the AD Workbench. From the repo root:
#   Rscript 04_ad_classifier/05_risk_profiles.R
# ==============================================================================

# ==============================================================================
# SHAP values of the AD classifier (predictors = predicted phenotypes)
# Goal: cluster participants by their SHAP vectors and describe each cluster
# ==============================================================================

# Setup -----------------------------------------------------------------------
source("config.R")
pp <- read_predicted_phenotypes(PREDICTED_RDS)

# Load required libraries
library(tidyverse)
library(data.table)
library(caret)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(broom)

# Source custom functions
source("R/vendor/ml_func1.R")
source("R/vendor/catboost_pipeline.R")


target_annotations = get_table(pp, "target_annotation_updated")%>%
  dplyr::select(-label)%>%
  dplyr::rename(label=plot_label)

gen_lookup = function(code_str){
  parts <- strsplit(code_str%>%gsub("\\.$","",.), ";\\s*")[[1]]
  lookup <- setNames(
    sub(".*=\\s*", "", parts),
    sub("\\s*=.*", "", parts)
  )
  return(lookup)
}
target_annotations%>%
  rowwise()%>%
  mutate(class_lookup = gen_lookup(Class_label)%>%list()) -> target_annotations
indx = which(target_annotations$Target=="apoe_genotype")
names(target_annotations$class_lookup[[indx]]) = as.character(1:6)


# Define paths ----------------------------------------------------------------
base_path <- file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01/ml/p1_trait"))
submodels_path <- file.path(base_path, "submodels")

# Load and process data -------------------------------------------------------

# 1. Load sample IDs
sample_id <- readRDS(file.path(base_path, "sample_id.rds"))

# 2. Load and process SHAP values
shap_files <- list.files(
  path = submodels_path,
  pattern = "shap_values.rds",
  recursive = TRUE,
  full.names = TRUE
)

shap_values <- shap_files %>%
  map(readRDS) %>%
  bind_rows()

# Process SHAP values matrix
rownames(shap_values) <- sample_id[as.numeric(shap_values$sample_id)]
shap_values <- shap_values %>%
  select(-sample_id, -actual) %>%
  as.matrix() %>%
  t()

data.frame(feature=rownames(shap_values))%>%
  mutate(target_id_no_class = gsub("_class.+","",feature))%>%
  inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
  filter(grepl("class",feature))%>%
  mutate(class_code = gsub(".+class_","",feature)%>%as.numeric())%>%
  group_by(target_id_no_class)%>%
  mutate(N_class = dplyr::n())%>%
  rowwise()%>%
  mutate(class_code = ifelse(N_class>2,class_code+1,class_code),
         class_code=as.character(class_code))%>%
  dplyr::select(-N_class)%>%
  rowwise()%>%
  mutate(class_label = class_lookup%>% unlist() %>% .[class_code])%>%
  bind_rows(.,
            data.frame(feature=rownames(shap_values))%>%
              mutate(target_id_no_class = gsub("_class.+","",feature))%>%
              inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
              filter(!grepl("class",feature)))%>%
  mutate(class_label=replace(class_label,class_label%in%c("Yes","Limbic-type","Present"),NA))%>%
  mutate(label = replace(label,!is.na(class_label),paste0(label," - ", class_label)[!is.na(class_label)]))%>%
  dplyr::select(Category,label,feature) -> feature_annotation



# 3. Load and process predictions
prediction_files <- list.files(
  path = submodels_path,
  pattern = "predictions.rds",
  recursive = TRUE,
  full.names = TRUE
)

predictions <- prediction_files %>%
  map(readRDS) %>%
  map(~ filter(.x, split == "test")) %>%
  bind_rows()

predictions$IID <- sample_id[as.numeric(predictions$sample_id)]

# 4. Load raw input data and create metadata
load(file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01/param_config.RData")))
# Use the data-preparation function of this repository's config file rather than the copy
# saved inside the run (runs trained elsewhere hold their own paths). Same code, same data.
cfg <- new.env(); cfg$PARAM <- NULL
sys.source("04_ad_classifier/param_config/run_param8_one_contributoer_out_nesting.R", envir = cfg)
param_config$prep_data <- cfg$PARAM$prep_data
rawinput <- param_config$prep_data(param_config)

metadata <- rawinput %>%
  select(IID:fold) %>%
  inner_join(predictions, by = "IID") %>%
  select(IID, trait, fold, predicted_class, predicted_prob)

# Align data ------------------------------------------------------------------
# Match sample order between metadata and SHAP values
sample_index <- match(metadata$IID, colnames(shap_values))
shap_values <- shap_values[, sample_index]


# clustering ----------
library(speakeasyR)
library(proxy)
if(!file.exists(file.path(step_dir("04_risk_profiles"), "SE_res_all.rds"))){
  # Cosine similarity matrix (1 - cosine distance)
  dist_matrix <- dist(shap_values%>%t(), method = "cosine")
  dist_matrix = as.matrix(dist_matrix)
  dist_matrix = 1-dist_matrix
  res_all = speakeasyR::cluster(dist_matrix,independent_runs = 10,min_clust = 100,target_partitions = 20,subcluster = 3,max_threads = 10)
  res_all = t(res_all)%>%as.data.frame()
  colnames(res_all) = c("cluster1","cluster2","cluster3")
  res_all$IID = colnames(shap_values)
  
  saveRDS(res_all,file=file.path(step_dir("04_risk_profiles"), "SE_res_all.rds"))
  
}

if(!file.exists(file.path(step_dir("04_risk_profiles"), "SE_res_ad.rds"))){
  
  metadata%>%
    filter(trait==1,predicted_class==1)%>%pull(IID) -> keep
  dist_matrix <- dist(shap_values[,keep]%>%t(), method = "cosine")
  dist_matrix = as.matrix(dist_matrix)
  dist_matrix = 1-dist_matrix
  res_ad = speakeasyR::cluster(dist_matrix,independent_runs = 10,min_clust = 100,target_partitions = 20,subcluster = 3,max_threads = 10)
  
  res_ad = t(res_ad)%>%as.data.frame()
  colnames(res_ad) = c("cluster1","cluster2","cluster3")
  res_ad$IID = keep
  
  saveRDS(res_ad,file=file.path(step_dir("04_risk_profiles"), "SE_res_ad.rds"))
}

if(!file.exists(file.path(step_dir("04_risk_profiles"), "SE_res_nci.rds"))){
  metadata%>%
    filter(trait==0,predicted_class==0)%>%pull(IID) -> keep
  dist_matrix <- dist(shap_values[,keep]%>%t(), method = "cosine")
  dist_matrix = as.matrix(dist_matrix)
  dist_matrix = 1-dist_matrix
  res_nci = speakeasyR::cluster(dist_matrix,independent_runs = 10,min_clust = 100,target_partitions = 20,subcluster = 3,max_threads = 10)
  
  res_nci = t(res_nci)%>%as.data.frame()
  colnames(res_nci) = c("cluster1","cluster2","cluster3")
  res_nci$IID = keep
  
  saveRDS(res_nci,file=file.path(step_dir("04_risk_profiles"), "SE_res_nci.rds"))
  
}
SE_res_ad <- readRDS(file.path(step_dir("04_risk_profiles"), "SE_res_ad.rds"))
SE_res_nci <- readRDS(file.path(step_dir("04_risk_profiles"), "SE_res_nci.rds"))

SE_res_ad%>%
  mutate(cluster1=paste0("ad",cluster1),
         cluster2=paste0("ad",cluster2),
         cluster3=paste0("ad",cluster3))%>%
  bind_rows(.,SE_res_nci%>%
              mutate(cluster1=paste0("nci",cluster1),
                     cluster2=paste0("nci",cluster2),
                     cluster3=paste0("nci",cluster3)))%>%
  saveRDS(.,file=file.path(step_dir("04_risk_profiles"), "SE_res_ad_nci.rds"))


one_sample_t_test <- function(data, test_cluster) {
  
  # Calculate grand mean for this feature
  grand_mean <- mean(data$value, na.rm = TRUE)
  
  # Get data for target cluster
  cluster_data <- data %>%
    filter(cluster1 == test_cluster)
  
  # Perform one-sample t-test against grand mean
  test_result <- t.test(cluster_data$value, mu = grand_mean, alternative = "two.sided")
  
  # Determine direction
  cluster_mean <- mean(cluster_data$value, na.rm = TRUE)
  direction <- ifelse(cluster_mean > grand_mean, "above average", "below average")
  
  # Return results
  tibble(
    cluster = test_cluster,
    p_value = test_result$p.value,
    direction = direction,
    cluster_mean = cluster_mean,
    grand_mean = grand_mean,
    t_statistic = test_result$statistic,
    n_cluster = nrow(cluster_data)
  )
}

# cluster characterization 
# Function to process each clustering result file
process_cluster_file <- function(file_path, file_label) {
  
  cat("Processing:", file_label, "\n")
  
  # Read clustering results
  res <- readRDS(file_path)
  
  # Process metadata with clustering results
  metadata_processed <- metadata %>% 
    inner_join(res, by = "IID") %>%
    mutate(
      cluster1 = paste0("g", cluster1),
      cluster2 = paste0("g", cluster2),
      cluster1 = ifelse(cluster1 %in% names(table(cluster1)[table(cluster1) < 150]), 
                        NA, cluster1),
      cluster2 = ifelse(cluster2 %in% names(table(cluster2)[table(cluster2) < 150]), 
                        NA, cluster2)
    )%>%
    ungroup()
  
  # Cluster summary statistics
  cluster_summary <- metadata_processed %>%
    drop_na(cluster1) %>%
    ungroup() %>%
    dplyr::select(trait, fold, predicted_class, predicted_prob, cluster1) %>%
    unique() %>%
    group_by(cluster=cluster1) %>%
    summarise(
      N_sample = dplyr::n(),
      trait_one_ratio = sum(trait)/dplyr::n(),
      predicted_class_one_ratio = sum(as.numeric(predicted_class))/dplyr::n(),
      predicted_prob_mean = mean(predicted_prob),
      fold = list(table(fold) %>% as.data.frame()),
      .groups = "drop"
    ) %>%
    unnest(fold) %>%
    pivot_wider(names_from = fold, values_from = Freq, names_prefix = "fold_")%>%
    ungroup()
  
  # SHAP values analysis
  df_cluster_shap <- shap_values %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "key", values_to = "value") %>%
    inner_join(metadata_processed, by = c("key" = "IID")) %>%
    inner_join(feature_annotation, by = c("rowname" = "feature"))%>%
    ungroup()
  
  df_cluster_shap%>%
    group_by(Category,key,trait, fold,  predicted_class, predicted_prob, cluster1, cluster2, cluster3)%>%
    summarise(value=sum(value))%>%
    ungroup() -> df_cluster_shap_cat
  
  # Raw input analysis
  df_cluster_raw <- rawinput %>%
    as.data.frame() %>%
    dplyr::select(-trait, -fold) %>%
    pivot_longer(-IID, names_to = "key", values_to = "value") %>%
    inner_join(metadata_processed, by = "IID") %>%
    inner_join(feature_annotation, by = c("key" = "feature"))%>%
    ungroup()
  
  # Statistical tests for SHAP values - Cluster1
  result_shap_c1 <- df_cluster_shap %>%
    drop_na(cluster1) %>%
    group_by(Category, label, rowname) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "SHAP",
      file_source = file_label,
      feature = rowname,
      cluster_type = "cluster1"
    )
  
  result_shap_c1_cat <- df_cluster_shap_cat %>%
    ungroup()%>%
    drop_na(cluster1) %>%
    group_by(Category) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "SHAP_category",
      file_source = file_label,
      cluster_type = "cluster1"
    )
  
  # Statistical tests for SHAP values - Cluster2
  result_shap_c2 <- df_cluster_shap %>%
    drop_na(cluster2) %>%
    dplyr::select(-cluster1)%>%
    dplyr::rename(cluster1 = cluster2) %>%  # Rename for consistency with test function
    group_by(Category, label, rowname) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "SHAP",
      file_source = file_label,
      feature = rowname,
      cluster_type = "cluster2"
    )
  
  result_shap_c2_cat <- df_cluster_shap_cat %>%
    drop_na(cluster2) %>%
    dplyr::select(-cluster1)%>%
    dplyr::rename(cluster1 = cluster2) %>%  # Rename for consistency with test function
    group_by(Category) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "SHAP_category",
      file_source = file_label,
      cluster_type = "cluster2"
    )
  
  
  # Statistical tests for raw values - Cluster1
  result_raw_c1 <- df_cluster_raw %>%
    drop_na(cluster1) %>%
    group_by(Category, label, key) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "Raw",
      file_source = file_label,
      feature = key,
      cluster_type = "cluster1"
    )
  
  # Statistical tests for raw values - Cluster2
  result_raw_c2 <- df_cluster_raw %>%
    drop_na(cluster2) %>%
    dplyr::select(-cluster1)%>%
    dplyr::rename(cluster1 = cluster2) %>%  # Rename for consistency with test function
    group_by(Category, label, key) %>%
    summarise(
      test_results = list(
        unique(cluster1) %>%
          map_dfr(~one_sample_t_test(cur_data(), .x))
      ),
      .groups = "drop"
    ) %>%
    unnest(test_results) %>%
    mutate(
      analysis_type = "Raw",
      file_source = file_label,
      feature = key,
      cluster_type = "cluster2"
    )
  
  # Combine all statistical results
  result_shap <- bind_rows(result_shap_c1, result_shap_c2)
  result_raw <- bind_rows(result_raw_c1, result_raw_c2)
  result_shap_cat <- bind_rows(result_shap_c1_cat, result_shap_c2_cat)
  
  # Cluster characterization for cluster2 as well
  cluster2_summary <- metadata_processed %>%
    drop_na(cluster2) %>%
    ungroup() %>%
    dplyr::select(trait, fold, predicted_class, predicted_prob, cluster2) %>%
    unique() %>%
    group_by(cluster = cluster2) %>%
    summarise(
      N_sample = dplyr::n(),
      trait_one_ratio = sum(trait)/dplyr::n(),
      predicted_class_one_ratio = sum(as.numeric(predicted_class))/dplyr::n(),
      predicted_prob_mean = mean(predicted_prob),
      fold = list(table(fold) %>% as.data.frame()),
      .groups = "drop"
    ) %>%
    unnest(fold) %>%
    pivot_wider(names_from = fold, values_from = Freq, names_prefix = "fold_")
  
  return(list(
    cluster1_summary = cluster_summary %>% mutate(file_source = file_label, cluster_type = "cluster1"),
    cluster2_summary = cluster2_summary %>% mutate(file_source = file_label, cluster_type = "cluster2"),
    shap_results = result_shap,
    raw_results = result_raw,
    result_shap_cat=result_shap_cat,
    metadata = metadata_processed
  ))
}

# Process all four files
files_to_process <- list(
  "SE_res_all" = file.path(step_dir("04_risk_profiles"), "SE_res_all.rds"),
  "SE_res_ad" = file.path(step_dir("04_risk_profiles"), "SE_res_ad.rds"),
  "SE_res_nci" = file.path(step_dir("04_risk_profiles"), "SE_res_nci.rds"),
  "SE_res_ad_nci" = file.path(step_dir("04_risk_profiles"), "SE_res_ad_nci.rds")
)

# Process each file
all_results <- map2(files_to_process, names(files_to_process), process_cluster_file)


# Combine results across all files
combined_cluster1_summary <- map_dfr(all_results, ~.x$cluster1_summary)
combined_cluster2_summary <- map_dfr(all_results, ~.x$cluster2_summary)
combined_shap_results <- map_dfr(all_results, ~.x$shap_results)
combined_raw_results <- map_dfr(all_results, ~.x$raw_results)
combined_shap_results_cat <- map_dfr(all_results, ~.x$result_shap_cat)
# Comprehensive cluster characterization
cluster_characterization <- bind_rows(
  combined_cluster1_summary,
  combined_cluster2_summary
)
fwrite(cluster_characterization,file=file.path(step_dir("04_risk_profiles"), "cluster_characterization.txt"),sep = "\t")
fwrite(combined_shap_results,file=file.path(step_dir("04_risk_profiles"), "combined_shap_results.txt"),sep = "\t")
fwrite(combined_raw_results,file=file.path(step_dir("04_risk_profiles"), "combined_raw_results.txt"),sep = "\t")
fwrite(combined_shap_results_cat,file=file.path(step_dir("04_risk_profiles"), "combined_shap_results_cat.txt"),sep = "\t")


file_path=file.path(step_dir("04_risk_profiles"), "SE_res_ad_nci.rds")

# Read clustering results
res <- readRDS(file_path)

# Process metadata with clustering results
metadata_processed <- metadata %>% 
  inner_join(res, by = "IID") %>%
  mutate(
    cluster1 = paste0("g", cluster1),
    cluster2 = paste0("g", cluster2),
    cluster1 = ifelse(cluster1 %in% names(table(cluster1)[table(cluster1) < 150]), 
                      NA, cluster1),
    cluster2 = ifelse(cluster2 %in% names(table(cluster2)[table(cluster2) < 150]), 
                      NA, cluster2)
  )%>%
  ungroup()

# SHAP values analysis
df_cluster_shap <- shap_values %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(-rowname, names_to = "key", values_to = "value") %>%
  inner_join(metadata_processed, by = c("key" = "IID")) %>%
  inner_join(feature_annotation, by = c("rowname" = "feature"))%>%
  ungroup()

df_cluster_shap%>%
  group_by(Category,key,trait, fold,  predicted_class, predicted_prob, cluster1, cluster2, cluster3)%>%
  summarise(value=sum(value))%>%
  ungroup() -> df_cluster_shap_cat

df_cluster_shap_cat%>%
  dplyr::select(key,Category,value)%>%
  spread(key,value)%>%
  as.data.frame()-> shap_mat_cat

rownames(shap_mat_cat) = shap_mat_cat$Category
shap_mat_cat$Category = NULL
shap_mat_cat = as.matrix(shap_mat_cat)

metadata_processed%>%
  drop_na(cluster1)%>%
  arrange(cluster1)%>%pull(IID) -> IID_ordered

metadata_processed = metadata_processed%>%
  mutate(cluster1 = c( "AD1",  "AD2",  "AD3",  "Non-AD1", "Non-AD2")[match(cluster1,c( "gad1",  "gad2",  "gad3",  "gnci2", "gnci4"))])


shap_mat_cat = shap_mat_cat[,IID_ordered]

CATEGORY_LEVELS <- c("Pathology: Alzheimer's disease", "Pathology: Tangles" , "Pathology: Amyloid", "Pathology: TDP-43",  "Pathology: Vascular", "Pathology: Others","Cognition", 
                     "Motor and Gait","Disabilities", "Medical Conditions", "Lifestyle", "Medications", 
                     "Genetics", "Demographics")
rownames(shap_mat_cat)[rownames(shap_mat_cat)=="Pathology"] = "Pathology: Others"
shap_mat_cat = shap_mat_cat[CATEGORY_LEVELS,]
cluster1_anno <- metadata_processed$cluster1[match(colnames(shap_mat_cat), metadata_processed$IID)]
fold_anno <- metadata_processed$fold[match(colnames(shap_mat_cat), metadata_processed$IID)]
# Colors for the 5 subtypes
subtype_colors = c(
  "AD1" = "#D73027",      # Dark red
  "AD2" = "#FC8D59",      # Orange-red  
  "AD3" = "#FEE08B",      # Light orange/yellow
  "Non-AD1" = "#4575B4",  # Dark blue
  "Non-AD2" = "#74ADD1"   # Light blue
)

# Create annotation with predefined colors
ha = HeatmapAnnotation(
  Subtypes = cluster1_anno,
  col = list(Subtypes = subtype_colors),
  show_legend = FALSE
)
set.seed(12364)
p = Heatmap(
  shap_mat_cat,
  name = "SHAP Value",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  top_annotation = ha,
  column_split = cluster1_anno, # split columns by subtype
  row_names_gp = gpar(fontfamily = "Arial"),        # Row names font
  column_names_gp = gpar(fontfamily = "Arial"),     # Column names font
  heatmap_legend_param = list(
    labels_gp = gpar(fontfamily = "Arial"),         # Legend labels
    title_gp = gpar(fontfamily = "Arial"),           # Legend title
    legend_width = unit(10, "cm")         # More space for legend labels
  ),
  col = colorRamp2(c(-0.5,0,0.5), c("#4078BD", "white", "#B83A3F")),
  row_names_max_width = unit(10, "cm"),  # Increase space for row labels
)
png(file.path(step_dir("04_risk_profiles"), "shap_heatmap_cat_all.png"), 
    width = 12, height = 4, units = "in", res = 300)
draw(p,padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()
