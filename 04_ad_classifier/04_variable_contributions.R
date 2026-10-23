# ==============================================================================
# Fig. 6b: SHAP feature importance of the AD classifier and UMAP of the
# predicted phenotypes by SHAP-value correlation.
# Inputs : predicted_phenotypes.rds (target_annotation, target_annotation_updated);
#          classifier run p8_one_contributoer_out_nesting_<EX_TAG> in
#          step_dir("04_classifier_runs") (01_train_classifiers.R).
# Outputs: step_dir("04_variable_contributions")
# Runs on the AD Workbench. From the repo root:
#   Rscript 04_ad_classifier/04_variable_contributions.R
# ==============================================================================

# ==============================================================================
# SHAP values of the AD classifier (predictors = predicted phenotypes)
# Goal: rank predicted phenotypes by mean |SHAP| and map them by SHAP-value correlation
# ==============================================================================

# Setup -----------------------------------------------------------------------
source("config.R")
pp <- read_predicted_phenotypes(PREDICTED_RDS)

invisible(step_dir("04_variable_contributions"))

# Load required libraries
library(tidyverse)
library(data.table)
library(caret)
library(parallel)
library(ComplexHeatmap)
library(circlize)

# Source custom functions
source("R/vendor/ml_func1.R")
source("R/vendor/catboost_pipeline.R")

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


target_annotations = get_table(pp, "target_annotation")

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


shap_values%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  gather(IID,shap_value,-feature) -> ds

ds%>%
  group_by(feature)%>%
  summarise(shap_value=mean(abs(shap_value)))%>%
  arrange(desc(shap_value)) -> shap_rank

fwrite(shap_rank,file=file.path(step_dir("04_variable_contributions"), "shap_rank.txt"))


taget_annotation_updated = get_table(pp, "target_annotation_updated")

shap_values%>%
  as.data.frame()%>%
  rownames_to_column("feature")%>%
  mutate(target_id_no_class = gsub("_class.+","",feature))%>%
  inner_join(.,taget_annotation_updated%>%dplyr::select(`Target ID`,Category),by=c("target_id_no_class"="Target ID"))%>%
  gather(IID,shap_value,-feature,-Category,-target_id_no_class)%>%
  group_by(Category,IID)%>%
  summarise(shap_value = sum(shap_value))%>%
  group_by(Category)%>%
  summarise(shap_value=mean(abs(shap_value)))%>%
  arrange(desc(shap_value)) -> shap_rank_category
fwrite(shap_rank_category,file=file.path(step_dir("04_variable_contributions"), "shap_rank_category.txt"))


shap_rank%>%
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
            shap_rank%>%
              mutate(target_id_no_class = gsub("_class.+","",feature))%>%
              inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
              filter(!grepl("class",feature)))%>%
  arrange(desc(shap_value))%>%
  mutate(class_label=replace(class_label,class_label%in%c("Yes","Limbic-type","Present"),NA))%>%
  mutate(label = replace(label,!is.na(class_label),paste0(label," - ", class_label)[!is.na(class_label)]))%>%
  dplyr::select(Category,label,shap_value)%>%
  head(20)->  vis_data


# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For number formatting

# Use the provided color scheme
category_colors <- c(
  "Pathology" = "#E69F00",
  "Clinical Diagnosis" = "#56B4E9",
  "Cognition" = "#009E73", 
  "Motor and Gait" = "#F0E442",
  "Medical Conditions" = "#0072B2", 
  "Medications" = "#D55E00",
  "Genetics" = "#CC79A7", 
  "Demographics" = "#999999"
)

# Create the plot
ggplot(vis_data, aes(y = reorder(label, shap_value), x = shap_value, fill = Category)) +
  # Add the bars
  geom_bar(stat = "identity", width = 0.7) +
  # Add thin separating lines for better readability
  geom_hline(yintercept = seq(0.5, nrow(vis_data) - 0.5, by = 1), 
             color = "white", linewidth = 0.5) +
  # Color according to category
  scale_fill_manual(values = category_colors) +
  # Set x-axis limits with no padding on the left
  scale_x_continuous(
    expand = c(0, 0), 
    limits = c(0, max(vis_data$shap_value) * 1.05),
    labels = number_format(accuracy = 0.001)
  ) +
  # Add labels
  labs(
    title = "Feature Importance Based on SHAP Values",
    y = NULL,
    x = "Mean absolute SHAP value",
    fill = "Feature Category"
  ) +
  # Apply a publication-ready theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.line.x = element_line(color = "gray50"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) -> p

# Save the plot with high resolution
save_plot <- function(plot, filename, width = 10, height = 10) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}
output_dir=step_dir("04_variable_contributions")
save_plot(p, "top10_contributors.png", width = 10, height = 5)
save_plot(p, "top10_contributors.pdf", width = 10, height = 5)


shap_rank%>%
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
            shap_rank%>%
              mutate(target_id_no_class = gsub("_class.+","",feature))%>%
              inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
              filter(!grepl("class",feature)))%>%
  arrange(desc(shap_value))%>%
  mutate(class_label=replace(class_label,class_label%in%c("Yes","Limbic-type","Present"),NA))%>%
  mutate(label = replace(label,!is.na(class_label),paste0(label," - ", class_label)[!is.na(class_label)]))%>%
  dplyr::select(feature,Category,label) -> target_annotations_slim

ds %>%
  dplyr::select(IID, feature, shap_value) %>%
  spread(feature, shap_value) %>%
  dplyr::select(-IID) %>%
  cor() -> shap_cor

fwrite(shap_cor,file=file.path(step_dir("04_variable_contributions"), "shap_cor.txt"),row.names = T,col.names = T)
  
# Keep features with no missing correlation (e.g. constant SHAP values)
indx = is.na(shap_cor) %>% rowSums() == 0
shap_cor = shap_cor[indx, indx]

# Convert correlation matrix to a distance matrix
# Using 1-abs(correlation) as the distance measure
shap_dist <- 1 - abs(shap_cor)

# Apply UMAP to the distance matrix
library(umap)
set.seed(42)
umap_result <- umap(shap_dist, n_neighbors = 15, min_dist = 0.1)

# Create a dataframe for plotting
library(dplyr)
library(tidyr)
umap_df <- data.frame(
  x = umap_result$layout[,1],
  y = umap_result$layout[,2],
  feature = colnames(shap_cor)
) %>%
  inner_join(., target_annotations_slim, by = "feature") %>%
  mutate(
    is_top20 = label %in% vis_data$label,
    point_size = ifelse(is_top20, 4, 2.5),
    point_alpha = ifelse(is_top20, 1, 0.7),
    label = ifelse(is_top20, gsub(".+:", "", label), "")
  )

# Category colors
library(ggrepel)
library(ggplot2)
library(viridis)
library(showtext)
library(ggforce)

# Enhanced color palette
category_colors <- c(
  "Pathology" = "#E69F00",
  "Clinical Diagnosis" = "#56B4E9",
  "Cognition" = "#009E73", 
  "Motor and Gait" = "#F0E442",
  "Medical Conditions" = "#0072B2", 
  "Medications" = "#D55E00",
  "Genetics" = "#CC79A7", 
  "Demographics" = "#999999"
)

# Plot the results
p <- ggplot(umap_df, aes(-x, -y, color = Category)) +
  # Add a light grey background
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "grey97", alpha = 0.4) +
  
  # Add points with dynamic size and alpha
  geom_point(aes(size = point_size, alpha = point_alpha), shape = 16) +
  
  # Add labels
  geom_text_repel(
    aes(label = label, size = ifelse(is_top20, 4, 0)), 
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.size = 0.2,
    segment.color = "grey50",
    max.overlaps = 50,
    min.segment.length = 0.1,
    force = 5
  ) +
  
  # Set color scales
  scale_color_manual(values = category_colors, name = "Feature Category") +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_size_identity() +
  scale_alpha_identity() +
  
  # Add coordinate equal aspect ratio
  coord_equal() +
  
  # Enhanced theme
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 20)),
    plot.caption = element_text(size = 9, color = "grey30", margin = margin(t = 15)),
    axis.title = element_text(size = 13, face = "bold", color = "grey20"),
    axis.text = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = "grey90"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_line(color = "grey92", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.7),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
  ) +
  
  # Labels
  labs(
    title = "Feature Relationships Based on SHAP Value Correlation",
    subtitle = "UMAP of features, distance = 1 - |SHAP-value correlation|",
    caption = paste0("Analysis based on SHAP values from ", nrow(umap_df), " features. Top 20 features highlighted."),
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )
p

output_dir=step_dir("04_variable_contributions")
save_plot(p, "umap_shap.png", width = 10, height = 8)
save_plot(p, "umap_shap.pdf", width = 10, height = 8)
# Fig 6b source data: one row per phenotype feature (no per-participant data).
# The plot draws (-x, -y).
fwrite(umap_df, file.path(output_dir, "umap_df.csv"))

