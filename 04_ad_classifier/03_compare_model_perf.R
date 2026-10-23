# Compare AD classifier performance across predictor sets: Fig. 6a, Extended Data Fig. 8
# (panel a = figs1, b = figs3, c = figs2) and Supplementary Table 10, plus the paired
# Wilcoxon tests of Fig. 6a. Runs on the AD Workbench.
#   Input : step_dir("04_collect_result")/fold_metrics.rds (made by 04_ad_classifier/02_collect_results.R)
#   Output: step_dir("04_compare_model_perf"): fig6a.{pdf,txt}, fig6a_paired_test.txt,
#           figs1.{pdf,txt}, figs2.{pdf,txt}, figs3.{pdf,txt}
# Run from the repository root:
#   Rscript 04_ad_classifier/03_compare_model_perf.R

source("config.R")
library(tidyverse)
library(data.table)
outdir = step_dir("04_compare_model_perf")

fold_metrics = readRDS(file.path(step_dir("04_collect_result"),"fold_metrics.rds"))
fold_metrics = fold_metrics%>%
  filter(!fold %in% paste0("fold_", EXCLUDED_CONTRIBUTORS))%>%
  filter(fold!="fold_N")

indx = fold_metrics$file_name%>%grepl("cross_validation",.)
fold_metrics$cross_valiation_type = "one-contributor-out"
fold_metrics$cross_valiation_type[indx] = "10-fold"
fold_metrics$control = "Non-AD Controls"
indx = fold_metrics$file_name%>%grepl("clean",.)
fold_metrics$control[indx] = "Disease-Free Controls"

fold_metrics = fold_metrics%>%
  filter(grepl(EX_TAG,file_name))

fold_metrics = fold_metrics%>%
  filter(!grepl(paste0(EX_TAG, "_$"),input_type))

fold_metrics$feature = ""
fold_metrics$feature  = fold_metrics$file_name%>%gsub("/ml/p1_trait/cv_summary/evaluation_summary.rds","",.)%>%
  gsub("./job_output/","",.)%>%gsub("/trait..ad__feature..|/trait..ad","",.)%>%
  gsub(".+one_contributor_out_","",.)%>%
  gsub(".+_one_contributoer_out_","",.)%>%
  gsub(paste0("nesting_", EX_TAG),"",.)%>%
  gsub(paste0("_", EX_TAG),"",.)%>%gsub("_clean_control","",.)

fold_metrics = fold_metrics%>%
  dplyr::select(fold,control,feature,normalized_pr_auc,somers_d,auc,n_var_no_class)%>%
  drop_na()%>%
  group_by(control,feature,n_var_no_class)%>%nest()

fold_metrics$model_th=NA
fold_metrics$model_th[grepl("metric",fold_metrics$feature)] = fold_metrics$feature%>%gsub(".+\\.\\.","",.)%>%.[grepl("metric",fold_metrics$feature)]
fold_metrics$feature = fold_metrics$feature%>%
  gsub("markers__top_n..(.*)$","SHAP-selected (top=\\1)",.)%>%
  gsub("__metric.+","",.)%>%
  gsub("ex.Clinical.Diagnosis","All phenotypes",.)%>%
  gsub("raw","All proteins",.)%>%
  gsub("pca","Principal components",.)%>%
  gsub("modules","Co-abundance modules",.)%>%
  gsub("\\."," ",.)
fold_metrics = fold_metrics%>%
  filter(!feature%in%c( "Clinical Diagnosis","all_1st_tier"))

fold_metrics = fold_metrics%>%
  filter(!((model_th%in%c("0","0.1","0.2","0.3"))&(!feature%in%c("All phenotypes")) ))%>%
  mutate(feature=replace(feature,feature=="All phenotypes",paste0(feature," (r>",model_th,")")))%>%
  mutate(feature=paste0(feature," (",n_var_no_class,")"))%>%
  ungroup()%>%
  dplyr::select(-n_var_no_class)
fold_metrics$category ="Domain-specific"
fold_metrics$category[grepl("All pheno",fold_metrics$feature)] = "All domains"
fold_metrics$category[grepl("Principal|SHAP|Co-ab|All pro",fold_metrics$feature)] = "Protein-based"
fold_metrics$category = factor(fold_metrics$category,levels = c("All domains", "Protein-based", "Domain-specific"))
category = rep(factor(
  c(
    "All domains",
    rep("Protein-based", 3),
    rep("Domain-specific", 7)
  ),
  levels = c("All domains", "Protein-based", "Domain-specific")
), each = 9)

##### Fig. 6a ####
results_by_fold = fold_metrics%>%
  filter(control=="Non-AD Controls")%>%
  filter(!((!grepl("r>0.01",feature))&grepl("All domain",category)))%>%
  filter(!((!grepl("top=5)|All proteins|Principal",feature))&grepl("Protein-based",category)))%>%
  mutate(feature = gsub("\\(top.+\\) ","",feature))%>%
  mutate(feature = gsub("\\(r>.+\\) ","",feature))%>%
  dplyr::select(-control,-model_th)%>%
  unnest()%>%
  dplyr::select(model=feature,fold,nAUPRC=normalized_pr_auc,somers_d,auc,category)

library(ggforce)
library(patchwork)
# Calculate summary statistics
results_summary <- results_by_fold %>%
  group_by(model, category) %>%
  summarise(
    auc_mean = mean(auc),
    nAUPRC_mean = mean(nAUPRC),
    nAUPRC_sd = sd(nAUPRC),
    nAUPRC_ci = 1.96 * sd(nAUPRC) / sqrt(n()),
    somers_d_mean = mean(somers_d),
    somers_d_sd = sd(somers_d),
    somers_d_ci = 1.96 * sd(somers_d) / sqrt(n()),
    .groups = "drop"
  )

# Calculate ranking by mean nAUPRC within category
results_summary <- results_summary %>%
  group_by(category) %>%
  mutate(rank_within = rank(-nAUPRC_mean)) %>%
  ungroup() %>%
  arrange(category, desc(rank_within)) %>%
  mutate(model = fct_inorder(model))

# Pivot to long format
results_long <- results_summary %>%
  pivot_longer(
    cols = c(nAUPRC_mean, somers_d_mean),
    names_to = "metric",
    values_to = "mean"
  ) %>%
  mutate(
    sd = if_else(metric == "nAUPRC_mean", nAUPRC_sd, somers_d_sd),
    ci = if_else(metric == "nAUPRC_mean", nAUPRC_ci, somers_d_ci),
    metric = factor(metric, 
                    levels = c("nAUPRC_mean", "somers_d_mean"),
                    labels = c("Prevalence-normalized AUCPR", "Somers' D (normalized AUCROC)"))
  ) %>%
  select(model, category, metric, mean, sd,ci)

# =============================================================================
# PLOT with pointrange (mean ± 95% CI)
# =============================================================================

p_nAUPRC <- results_long %>%
  filter(metric == "Prevalence-normalized AUCPR") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  
  labs(x = "Prevalence-normalized AUCPR")

p_somers <- results_long %>%
  filter(metric == "Somers' D (normalized AUCROC)") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  labs(x = "Somers' D (normalized AUCROC)")

# =============================================================================
# COMBINE
# =============================================================================

p_combined <- p_nAUPRC + p_somers +
  plot_layout(ncol = 2)
p_combined
# =============================================================================
# SAVE
# =============================================================================

ggsave(file.path(outdir,"fig6a.pdf"), p_combined, width = 9, height = 4)
fwrite(results_summary,file.path(outdir,"fig6a.txt"),sep = "\t")


# Paired Wilcoxon tests (two-sided, exact): P, n_pairs and rank-biserial r
results_by_fold%>%filter(grepl("All pheno|SHAP",model))%>%arrange(fold)%>%
  dplyr::select(category,fold,nAUPRC)%>%
  spread(category,nAUPRC) -> input_for_paried_test
wilcox.test(input_for_paried_test$`All domains`,input_for_paried_test$`Protein-based`,paired = TRUE,alternative = "two.sided",exact = TRUE)%>%
  broom::tidy()%>%
  mutate(metric = "nAUPRC",contrast = "All domain vs SHAP-selected",
         n_pairs = nrow(input_for_paried_test),
         r_effect = 1 - (2 * statistic) / (n_pairs * (n_pairs + 1) / 2)) -> res_paired

results_by_fold%>%filter(grepl("All",model))%>%arrange(fold)%>%
  dplyr::select(category,fold,nAUPRC)%>%
  spread(category,nAUPRC) -> input_for_paried_test
wilcox.test(input_for_paried_test$`All domains`,input_for_paried_test$`Protein-based`,paired = TRUE,alternative = "two.sided",exact = TRUE)%>%
  broom::tidy()%>%
  mutate(metric = "nAUPRC",contrast = "All domain vs All proteins",
         n_pairs = nrow(input_for_paried_test),
         r_effect = 1 - (2 * statistic) / (n_pairs * (n_pairs + 1) / 2))%>%
  bind_rows(res_paired,.)-> res_paired

results_by_fold%>%filter(grepl("All pheno|SHAP",model))%>%arrange(fold)%>%
  dplyr::select(category,fold,somers_d)%>%
  spread(category,somers_d) -> input_for_paried_test
wilcox.test(input_for_paried_test$`All domains`,input_for_paried_test$`Protein-based`,paired = TRUE,alternative = "two.sided",exact = TRUE)%>%
  broom::tidy()%>%
  mutate(metric = "somers_d",contrast = "All domain vs SHAP-selected",
         n_pairs = nrow(input_for_paried_test),
         r_effect = 1 - (2 * statistic) / (n_pairs * (n_pairs + 1) / 2))%>%
  bind_rows(res_paired,.)-> res_paired

results_by_fold%>%filter(grepl("All",model))%>%arrange(fold)%>%
  dplyr::select(category,fold,somers_d)%>%
  spread(category,somers_d) -> input_for_paried_test
wilcox.test(input_for_paried_test$`All domains`,input_for_paried_test$`Protein-based`,paired = TRUE,alternative = "two.sided",exact = TRUE)%>%
  broom::tidy()%>%
  mutate(metric = "somers_d",contrast = "All domain vs All proteins",
         n_pairs = nrow(input_for_paried_test),
         r_effect = 1 - (2 * statistic) / (n_pairs * (n_pairs + 1) / 2))%>%
  bind_rows(res_paired,.)-> res_paired
fwrite(res_paired,file.path(outdir,"fig6a_paired_test.txt"),sep = "\t")


####### Feature selection gradient (Extended Data Fig. 8a) ######
results_by_fold = fold_metrics%>%
  filter(control=="Non-AD Controls")%>%
  filter(!((!grepl("r>0.01",feature))&grepl("All domain",category)))%>%
  filter(((!grepl("Principal|abundance",feature))&grepl("Protein-based",category)))%>%
  mutate(feature = gsub("top=","Top-",feature))%>%
  dplyr::select(-control,-model_th)%>%
  arrange(gsub(".+\\(|\\)","",feature)%>%as.numeric())%>%
  mutate(feature=factor(feature,levels = feature))%>%
  unnest()%>%
  dplyr::select(model=feature,fold,nAUPRC=normalized_pr_auc,somers_d,category)

library(ggforce)
library(patchwork)
# Calculate summary statistics
results_summary <- results_by_fold %>%
  group_by(model, category) %>%
  summarise(
    nAUPRC_mean = mean(nAUPRC),
    nAUPRC_sd = sd(nAUPRC),
    nAUPRC_ci = 1.96 * sd(nAUPRC) / sqrt(n()),
    somers_d_mean = mean(somers_d),
    somers_d_sd = sd(somers_d),
    somers_d_ci = 1.96 * sd(somers_d) / sqrt(n()),
    .groups = "drop"
  )

# Pivot to long format
results_long <- results_summary %>%
  pivot_longer(
    cols = c(nAUPRC_mean, somers_d_mean),
    names_to = "metric",
    values_to = "mean"
  ) %>%
  mutate(
    sd = if_else(metric == "nAUPRC_mean", nAUPRC_sd, somers_d_sd),
    ci = if_else(metric == "nAUPRC_mean", nAUPRC_ci, somers_d_ci),
    metric = factor(metric, 
                    levels = c("nAUPRC_mean", "somers_d_mean"),
                    labels = c("Prevalence-normalized AUCPR", "Somers' D (normalized AUCROC)"))
  ) %>%
  select(model, category, metric, mean, sd,ci)

# =============================================================================
# PLOT with pointrange (mean ± 95% CI)
# =============================================================================

p_nAUPRC <- results_long %>%
  filter(metric == "Prevalence-normalized AUCPR") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  
  labs(x = "Prevalence-normalized AUCPR")

p_somers <- results_long %>%
  filter(metric == "Somers' D (normalized AUCROC)") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  labs(x = "Somers' D (normalized AUCROC)")

# =============================================================================
# COMBINE
# =============================================================================

p_combined <- p_nAUPRC + p_somers +
  plot_layout(ncol = 2)
p_combined
# =============================================================================
# SAVE
# =============================================================================

ggsave(file.path(outdir,"figs1.pdf"), p_combined, width = 9, height = 4)
fwrite(results_summary,file.path(outdir,"figs1.txt"),sep = "\t")


####### Sensitivity to the phenotype threshold (Extended Data Fig. 8c) ######
results_by_fold = fold_metrics%>%
  filter(control=="Non-AD Controls")%>%
  filter(((!grepl("r>0\\)",feature))&grepl("All domain",category)))%>%
  dplyr::select(-control,-model_th)%>%
  mutate(feature=factor(feature,levels = feature))%>%
  unnest()%>%
  dplyr::select(model=feature,fold,nAUPRC=normalized_pr_auc,somers_d,category)

library(ggforce)
library(patchwork)
# Calculate summary statistics
results_summary <- results_by_fold %>%
  group_by(model, category) %>%
  summarise(
    nAUPRC_mean = mean(nAUPRC),
    nAUPRC_sd = sd(nAUPRC),
    nAUPRC_ci = 1.96 * sd(nAUPRC) / sqrt(n()),
    somers_d_mean = mean(somers_d),
    somers_d_sd = sd(somers_d),
    somers_d_ci = 1.96 * sd(somers_d) / sqrt(n()),
    .groups = "drop"
  )

# Pivot to long format
results_long <- results_summary %>%
  pivot_longer(
    cols = c(nAUPRC_mean, somers_d_mean),
    names_to = "metric",
    values_to = "mean"
  ) %>%
  mutate(
    sd = if_else(metric == "nAUPRC_mean", nAUPRC_sd, somers_d_sd),
    ci = if_else(metric == "nAUPRC_mean", nAUPRC_ci, somers_d_ci),
    metric = factor(metric, 
                    levels = c("nAUPRC_mean", "somers_d_mean"),
                    labels = c("Prevalence-normalized AUCPR", "Somers' D (normalized AUCROC)"))
  ) %>%
  select(model, category, metric, mean, sd,ci)

# =============================================================================
# PLOT with pointrange (mean ± 95% CI)
# =============================================================================

p_nAUPRC <- results_long %>%
  filter(metric == "Prevalence-normalized AUCPR") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  
  labs(x = "Prevalence-normalized AUCPR")

p_somers <- results_long %>%
  filter(metric == "Somers' D (normalized AUCROC)") %>%
  ggplot(aes(x = mean, y = model, color = category)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  scale_color_manual(
    values = c(
      "All domains" = "#2E86AB",
      "Protein-based" = "#F18F01",
      "Domain-specific" = "#A23B72"
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  labs(x = "Somers' D (normalized AUCROC)")

# =============================================================================
# COMBINE
# =============================================================================

p_combined <- p_nAUPRC + p_somers +
  plot_layout(ncol = 2)
p_combined
# =============================================================================
# SAVE
# =============================================================================


ggsave(file.path(outdir,"figs2.pdf"), p_combined, width = 9, height = 3)
fwrite(results_summary,file.path(outdir,"figs2.txt"),sep = "\t")


####### Sensitivity to the control definition (Extended Data Fig. 8b) ######
results_by_fold = fold_metrics%>%
  filter(!((!grepl("r>0.01",feature))&grepl("All domain",category)))%>%
  filter(!((!grepl("top=5)|All proteins|Principal",feature))&grepl("Protein-based",category)))%>%
  mutate(feature = gsub("\\(top.+\\) ","",feature))%>%
  mutate(feature = gsub("\\(r>.+\\) ","",feature))%>%
  dplyr::select(-model_th)%>%
  unnest()%>%
  dplyr::select(control,model=feature,fold,nAUPRC=normalized_pr_auc,somers_d,category)

library(ggforce)
library(patchwork)
# Calculate summary statistics
results_summary <- results_by_fold %>%
  group_by(control,model, category) %>%
  summarise(
    nAUPRC_mean = mean(nAUPRC),
    nAUPRC_sd = sd(nAUPRC),
    nAUPRC_ci = 1.96 * sd(nAUPRC) / sqrt(n()),
    somers_d_mean = mean(somers_d),
    somers_d_sd = sd(somers_d),
    somers_d_ci = 1.96 * sd(somers_d) / sqrt(n()),
    .groups = "drop"
  )

# Calculate ranking by mean nAUPRC within category
results_summary <- results_summary %>%
  group_by(category) %>%
  mutate(rank_within = rank(-nAUPRC_mean)) %>%
  ungroup() %>%
  arrange(category, desc(rank_within)) %>%
  mutate(model = fct_inorder(model))

# Pivot to long format
results_long <- results_summary %>%
  pivot_longer(
    cols = c(nAUPRC_mean, somers_d_mean),
    names_to = "metric",
    values_to = "mean"
  ) %>%
  mutate(
    sd = if_else(metric == "nAUPRC_mean", nAUPRC_sd, somers_d_sd),
    ci = if_else(metric == "nAUPRC_mean", nAUPRC_ci, somers_d_ci),
    metric = factor(metric, 
                    levels = c("nAUPRC_mean", "somers_d_mean"),
                    labels = c("Prevalence-normalized AUCPR", "Somers' D (normalized AUCROC)"))
  ) %>%
  select(control,model, category, metric, mean, sd,ci)

# =============================================================================
# PLOT with pointrange (mean ± 95% CI)
# =============================================================================

p_nAUPRC <- results_long %>%
  filter(metric == "Prevalence-normalized AUCPR") %>%
  ggplot(aes(x = mean, y = model, color = control)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    position = position_dodge(width=0.4),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  
  labs(x = "Prevalence-normalized AUCPR")

p_somers <- results_long %>%
  filter(metric == "Somers' D (normalized AUCROC)") %>%
  ggplot(aes(x = mean, y = model, color = control)) +
  
  geom_pointrange(
    aes(xmin = mean - ci, xmax = mean + ci),
    position = position_dodge(width=0.4),
    size = 0.4,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  facet_col(
    facets = vars(category),
    scales = "free_y",
    space = "free"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.3, "lines")
  ) +
  labs(x = "Somers' D (normalized AUCROC)")

# =============================================================================
# COMBINE
# =============================================================================

p_combined <- p_nAUPRC + p_somers +
  plot_layout(ncol = 2)
p_combined
# =============================================================================
# SAVE
# =============================================================================

ggsave(file.path(outdir,"figs3.pdf"), p_combined, width = 9, height = 4)
fwrite(results_summary,file.path(outdir,"figs3.txt"),sep = "\t")
