# 03_trajectories/01_age_trajectory.R
# Age trajectories of predicted phenotypes in GNPC plasma (GAMs): Fig. 5a-d,
# Extended Data Fig. 5-6, and gam_stats_all_variants.csv for Supplementary Table 8.
# Inputs:  PROTEOMICS_RDS, PREDICTED_RDS (see config.R)
# Outputs: step_dir("03_age_trajectory") -- GAM statistics, onset ages, fitted
#          trajectories (.rds), Source Data CSVs and figures.
# Runs on the AD Workbench. From the repository root:
#   Rscript 03_trajectories/01_age_trajectory.R
# REPLOT_ONLY (default FALSE) fits every GAM and writes the statistics. The
# published panels were re-plotted from the saved fits (REPLOT_ONLY <- TRUE)
# after a full run; that pass does not touch the statistics tables.

source("config.R")
pp <- read_predicted_phenotypes(PREDICTED_RDS)

if(T){
  
  
  library(dplyr)
  library(data.table)
  library(pROC)
  library(ggplot2)
  library(tidyverse)
  library(PRROC)  # For precision-recall curves
  library(SummarizedExperiment)
  
  # Load proteomics data
  expdata <- readRDS(PROTEOMICS_RDS)
  metadata <- colData(expdata) %>% as.data.frame()
  
  
  # Load prediction data
  prediction <- get_prediction(pp, "no_class_labels")
  
  # Join metadata and prediction data
  data <- metadata %>%
    inner_join(prediction, by = "sample_id")
  
  data = data%>%
    filter(!contributor_code%in%EXCLUDED_CONTRIBUTORS)
  data$contributor_code = paste0(data$contributor_code,"_",data$sample_matrix)
  
  data$control =  !rowAnys(data%>%
                             select(all_of(c("pd", "als", "ad", "ftd", "stroke")))%>%
                             as.matrix() == "1")
  
  data%>%
    filter(grepl("Plasma",sample_matrix)) -> data_plasma
  
  # Standardize APOE genotype order (e.g., "32" becomes "23")
  data_plasma$APOE[data_plasma$APOE %in% "32"] <- "23"
  
  data_plasma <- data_plasma %>%
    mutate(APOE = na_if(as.character(APOE), "-1"),
           APOE4_count = recode(APOE, `22` = 0L, `23` = 0L, `24` = 1L, `33` = 0L, `34` = 1L, `44` = 2L, .default = NA_integer_),
           sex = na_if(as.character(sex), "-1"),
           ad = na_if(as.character(ad), "-1"))
  
  data_plasma = data_plasma%>%
    filter((!!sym("ad") == 0 & control) | !!sym("ad") == 1)
  
  # Predicted sex (class labels of T0002)
  sex_predictions <- get_class_labels(pp, "T0002")
  
  # Predicted APOE genotype (class labels of T0607)
  apoe_predictions <- get_class_labels(pp, "T0607")
  
  data_plasma = data_plasma%>%
    inner_join(.,sex_predictions%>%
                 dplyr::rename(sex_pred = prediction),by="sample_id")%>%
    inner_join(.,apoe_predictions%>%
                 dplyr::rename(APOE_pred = prediction)%>%
                 mutate(APOE4_count_pred = recode(APOE_pred, `22` = 0L, `23` = 0L, `24` = 1L, `33` = 0L, `34` = 1L, `44` = 2L, .default = NA_integer_),),
               by="sample_id")
  data_plasma$sex = 2-data_plasma$sex%>%as.numeric()
  c("T0726","T0727","T0734","T0049","T0019","T0003",
    "T0009","T0017","T0593","T0048","T0561","T0725",
    "T0722","T0730","T0552","T0047","T0546","T0560","T0559","T0675","T0691","T0692","T0693") -> selected_target
  c("T0615","T0720","T0593","T0732","T0616","T0706","T0713","T0614",
    "T0725","T0018","T0545","T0719","T0013","T0715","T0007","T0611","T0722",
    "T0012","T0004","T0010","T0005","T0015","T0675","T0003","T0691","T0011",
    "T0017","T0693","T0019","T0009","T0014",
    "T0692","T0001","T0049") -> selected_target
  flip_vars <- c("T0003", "T0009","T0019","T0593","T0012","T0011","T0010",
                 "T0004","T0005","T0006","T0007","T0008","T0595",
                 "T0614","T0615","T0616","T0617","T0618","T0619")
  flip_vars = intersect(flip_vars,selected_target)
  
  taget_annotation = get_table(pp, "target_annotation_updated")
  taget_annotation%>%
    filter(`Target ID`%in%selected_target)%>%
    mutate(Category = gsub(":.+","",Category))
  
  
  lapply(selected_target, function(x){
    colnames(data_plasma)[grepl(x,data_plasma%>%colnames())]
  })%>%unlist() -> vars
  
  
  data_plasma %>%
    filter(age_at_visit != -1) %>%
    select(SAMPLE.ID, all_of(vars)) %>%
    pivot_longer(cols = -SAMPLE.ID, names_to = "key", values_to = "value")%>%
    mutate(multiple = gsub(".+_","",key)%>%as.numeric())%>%
    mutate(multiple = replace(multiple,is.na(multiple),1))%>%
    mutate(key=gsub("_.+","",key))%>%
    group_by(SAMPLE.ID,key)%>%
    summarise(value=sum(value))%>%
    ungroup()%>%
    spread(key,value)%>%
    inner_join(data_plasma%>%dplyr::select(SAMPLE.ID,age_at_visit,ad,sex,
                                           sex_pred,APOE,APOE4_count,
                                           APOE_pred,APOE4_count_pred,contributor_code),.,by="SAMPLE.ID") -> data_plasma_no_class
  
  vars_no_class = vars%>%gsub("_.+","",.)
  
  # Prepare data: drop invalid ages, flip sign for specified vars
  prepared <- data_plasma_no_class %>%
    filter(age_at_visit != -1) %>%
    mutate(across(all_of(flip_vars), ~ -., .names = "{.col}")) %>%
    dplyr::select(age_at_visit,ad,sex,sex_pred,APOE,APOE4_count,APOE_pred,APOE4_count_pred,contributor_code, all_of(vars_no_class))
  
  # For each variable: mean value in 10-year age bins among non-AD participants;
  # min_bin / max_bin = the lowest and highest bin means (used for 0-1 scaling)
  ref_values <- prepared %>%
    filter(ad=="0")%>%
    dplyr::select(-(ad:APOE4_count_pred))%>%
    pivot_longer(cols = all_of(vars_no_class), names_to = "key", values_to = "value") %>%
    mutate(age_bin = cut(age_at_visit, breaks = c(-Inf,seq(20, 90, 10),Inf),right=FALSE))%>%
    group_by(key, age_bin) %>%
    summarize(ref_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(key)%>%
    filter(ref_mean%in%c(max(ref_mean),min(ref_mean)))%>%
    arrange(ref_mean)%>%
    mutate(group = c("min_bin","max_bin"))%>%
    dplyr::select(-age_bin)%>%
    ungroup()%>%
    pivot_wider(names_from = group, values_from = ref_mean) 
  
  
  # Apply 0-1 scaling
  long_scaled <- prepared %>%
    pivot_longer(cols = all_of(vars_no_class), names_to = "key", values_to = "value") %>%
    left_join(ref_values, by = "key") %>%
    mutate(
      scaled_value = (value - min_bin) / (max_bin - min_bin)
    )
  
  long_scaled = long_scaled%>%
    group_by(key)%>%
    nest()%>%
    inner_join(.,taget_annotation%>%
                 filter(`Target ID`%in%selected_target)%>%
                 mutate(Category = gsub(":.+","",Category))%>%
                 dplyr::select(-Target,-Short.Description,-label,-Class_label),by=c("key"="Target ID"))%>%
    unnest()
  
  effect_sizes <- long_scaled %>%
    filter(age_at_visit<40)%>%
    group_by(key) %>%
    summarize(                                     # max_bin/min_bin constant per key
      range_original = dplyr::first(max_bin) - dplyr::first(min_bin),
      range_sd = (dplyr::first(max_bin) - dplyr::first(min_bin)) / sd(value)
    )%>%
    ungroup()%>%
    unique()%>%
    inner_join(.,taget_annotation%>%
                 filter(`Target ID`%in%selected_target)%>%
                 mutate(Category = gsub(":.+","",Category))%>%
                 dplyr::select(-Target,-Short.Description,-label,-Class_label),by=c("key"="Target ID"))
  
  # Multi-Variable Trajectory Visualization #######
  
  library(tidyverse)
  library(mgcv)
  library(cowplot)
  library(viridis)
  library(RColorBrewer)
  
  
  # Use the provided color scheme
  category_colors <- c(
    "Pathology" = "#E69F00",
    "Pathology: Alzheimer's disease"= "#E69F00",
    "Pathology: Tangles"= "#E69F00" ,
    "Pathology: Amyloid"= "#E69F00", 
    "Pathology: TDP-43"= "#E69F00",  
    "Pathology: Vascular"= "#E69F00",
    "Disabilities" = "#56B4E9",
    "Cognition" = "#009E73", 
    "Motor and Gait" = "#F0E442",
    "Medical Conditions" = "#0072B2", 
    "Lifestyle" = "#0072B2", 
    "Medications" = "#D55E00",
    "Genetics" = "#CC79A7", 
    "Demographics" = "#999999"
  )
}


# 2. FIT GAMs AND GET PREDICTIONS ####

# Age grid for predictions
age_grid <- seq(
  min(long_scaled$age_at_visit, na.rm = TRUE),
  max(long_scaled$age_at_visit, na.rm = TRUE),
  by = 0.5
)

# Fit GAM for each variable and extract predictions
fit_gam_trajectory <- function(data, var_name) {
  
  df_var <- data %>% filter(plot_label == var_name)
  
  if (nrow(df_var) < 50) return(NULL)
  df_var$ad = as.numeric(df_var$ad)
  df_var <- df_var %>%
    mutate(contributor_code = factor(contributor_code))
  
  model <- tryCatch(
    gam(scaled_value ~ s(age_at_visit, k = 10)+contributor_code+sex+ad, data = df_var, method = "REML"),
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NULL)
  
  library(broom)
  gam_summ <- summary(model)
  tidy_param <- tidy(model, parametric = TRUE, conf.int = TRUE, conf.level = 0.95)
  tidy_smooth <- tidy(model, parametric = FALSE)
  tidy_all <- bind_rows(
    tidy_param %>% mutate(component = "parametric"),
    tidy_smooth %>% mutate(component = "smooth")
  )
  
  # Full Spearman correlation statistics
  cor_result <- cor.test(df_var$age_at_visit, df_var$scaled_value, method = "spearman") %>% broom::tidy()
  
  pred_df <- expand.grid(
    age_at_visit = age_grid,
    sex = c(0.5),#mean(df_var$sex_pred, na.rm = TRUE),
    contributor_code = model$xlevels$contributor_code,
    ad = 0
  )
  
  pred <- predict(model, newdata = pred_df, se.fit = TRUE)
  pred_df = pred_df%>%
    mutate(fit = pred$fit,
           se = pred$se.fit,
           lower = pred$fit - 1.96 * pred$se.fit,
           upper = pred$fit + 1.96 * pred$se.fit,
           plot_label = var_name)
  pred_df = pred_df%>%
    group_by(plot_label)%>%
    nest()%>%
    mutate(tidy_all=list(tidy_all))%>%
    mutate(r_squared = gam_summ$r.sq,
           dev_explained = gam_summ$dev.expl,
           n_obs = nrow(df_var),
           spearman_rho = cor_result$estimate,
           spearman_statistic = cor_result$statistic,
           spearman_p_value = cor_result$p.value)%>%
    ungroup()
  return(pred_df)
}



fit_gam_trajectory_apoe <- function(data, var_name,pred_type =c("actual")) {
  
  df_var <- data %>% filter(plot_label == var_name)
  
  if (nrow(df_var) < 50) return(NULL)
  df_var$ad = as.numeric(df_var$ad)
  if(pred_type=="pred"){
    # APOE4_ordered is not used by the model below (it uses the numeric APOE4 count)
    df_var <- df_var %>%
      mutate(APOE4_ordered = ordered(APOE4_count_pred))
    df_var <- df_var %>%
      mutate(contributor_code = factor(contributor_code))
    
    model <- tryCatch(
      gam(scaled_value ~ s(age_at_visit,k = 10) +
            ti(age_at_visit, APOE4_count_pred, k = c(10, 3), bs = c("cr", "cr")) +
            contributor_code+sex+ad+APOE4_count_pred, data = df_var, method = "REML"),
      error = function(e) NULL
    )
    
    
  }else if(pred_type=="actual"){
    model <- tryCatch(
      gam(scaled_value ~ s(age_at_visit, k = 10) + s(age_at_visit, k = 10,by = APOE4_count)+sex+ad+APOE4_count, data = df_var, method = "REML"),
      error = function(e) NULL
    )
  }
  
  if (is.null(model)) return(NULL)
  
  library(broom)
  gam_summ <- summary(model)
  tidy_param <- tidy(model, parametric = TRUE, conf.int = TRUE, conf.level = 0.95)
  tidy_smooth <- tidy(model, parametric = FALSE)
  tidy_all <- bind_rows(
    tidy_param %>% mutate(component = "parametric"),
    tidy_smooth %>% mutate(component = "smooth")
  )
  
  # Full Spearman correlation statistics
  cor_result <- cor.test(df_var$age_at_visit, df_var$scaled_value, method = "spearman") %>% broom::tidy()
  
  if(pred_type=="pred"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0.5),#mean(df_var$sex_pred, na.rm = TRUE),
      APOE4_count_pred = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      contributor_code = model$xlevels$contributor_code,
      ad = 0
    )%>%mutate(APOE4_ordered  = ordered(APOE4_count_pred))
  }else if(pred_type=="actual"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0.5),#mean(df_var$sex_pred, na.rm = TRUE),
      APOE4_count = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      ad = 0
    )
  }
  
  pred <- predict(model, newdata = pred_df, se.fit = TRUE)
  pred_df = pred_df%>%
    mutate(fit = pred$fit,
           se = pred$se.fit,
           lower = pred$fit - 1.96 * pred$se.fit,
           upper = pred$fit + 1.96 * pred$se.fit,
           plot_label = var_name)
  
  pred_df = pred_df%>%
    group_by(plot_label)%>%
    nest()%>%
    mutate(tidy_all=list(tidy_all))%>%
    mutate(r_squared = gam_summ$r.sq,
           dev_explained = gam_summ$dev.expl,
           n_obs = nrow(df_var),
           spearman_rho = cor_result$estimate,
           spearman_statistic = cor_result$statistic,
           spearman_p_value = cor_result$p.value)%>%
    ungroup()%>%
    mutate(pred_type=pred_type)
  return(pred_df)
}

fit_gam_trajectory_sex <- function(data, var_name,pred_type =c("actual")) {
  
  df_var <- data %>% filter(plot_label == var_name)
  
  if (nrow(df_var) < 50) return(NULL)
  df_var$ad = as.numeric(df_var$ad)
  if(pred_type=="pred"){
    # Ordered factor: the by = sex_ordered smooth is the difference from the reference sex
    df_var <- df_var %>%
      mutate(sex_ordered = ordered(sex))
    df_var <- df_var %>%
      mutate(contributor_code = factor(contributor_code))
    
    model <- tryCatch(
      gam(scaled_value ~ s(age_at_visit, k = 10) + 
            s(age_at_visit, k = 10,by = sex_ordered)+contributor_code+APOE4_count_pred+sex_ordered+ad, data = df_var, method = "REML"),
      error = function(e) NULL
    )
    
  }else if(pred_type=="actual"){
    model <- tryCatch(
      gam(scaled_value ~ s(age_at_visit, k = 10) + s(age_at_visit, k = 10,by = APOE4_count)+sex+ad, data = df_var, method = "REML"),
      error = function(e) NULL
    )
  }
  
  if (is.null(model)) return(NULL)
  
  library(broom)
  gam_summ <- summary(model)
  tidy_param <- tidy(model, parametric = TRUE, conf.int = TRUE, conf.level = 0.95)
  tidy_smooth <- tidy(model, parametric = FALSE)
  tidy_all <- bind_rows(
    tidy_param %>% mutate(component = "parametric"),
    tidy_smooth %>% mutate(component = "smooth")
  )
  
  # Full Spearman correlation statistics
  cor_result <- cor.test(df_var$age_at_visit, df_var$scaled_value, method = "spearman") %>% broom::tidy()
  
  if(pred_type=="pred"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0,1),#mean(df_var$sex, na.rm = TRUE),
      APOE4_count_pred = c(0),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      contributor_code = model$xlevels$contributor_code,
      ad = 0
    )%>%mutate(APOE4_ordered  = ordered(APOE4_count_pred),
               sex_ordered = ordered(sex))
  }else if(pred_type=="actual"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0.5),#mean(df_var$sex, na.rm = TRUE),
      APOE4_count = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      ad = 0
    )
  }
  
  pred <- predict(model, newdata = pred_df, se.fit = TRUE)
  pred_df = pred_df%>%
    mutate(fit = pred$fit,
           se = pred$se.fit,
           lower = pred$fit - 1.96 * pred$se.fit,
           upper = pred$fit + 1.96 * pred$se.fit,
           plot_label = var_name)
  pred_df = pred_df%>%
    group_by(plot_label)%>%
    nest()%>%
    mutate(tidy_all=list(tidy_all))%>%
    mutate(r_squared = gam_summ$r.sq,
           dev_explained = gam_summ$dev.expl,
           n_obs = nrow(df_var),
           spearman_rho = cor_result$estimate,
           spearman_statistic = cor_result$statistic,
           spearman_p_value = cor_result$p.value)%>%
    ungroup()%>%
    mutate(pred_type=pred_type)
  return(pred_df)
}


fit_gam_trajectory_apoe_sex <- function(data, var_name,pred_type =c("actual")) {
  
  df_var <- data %>% filter(plot_label == var_name)
  
  if (nrow(df_var) < 50) return(NULL)
  df_var$ad = as.numeric(df_var$ad)
  if(pred_type=="pred"){
    # Ordered factor: the by = sex_ordered smooth is the difference from the reference sex
    df_var <- df_var %>%
      mutate(sex_ordered = ordered(sex))
    df_var <- df_var %>%
      mutate(contributor_code = factor(contributor_code))
    # ANOVA-type tensor decomposition:
    # s(age) main effect
    # s(age, by=sex) age-by-sex smooth
    # ti(age, APOE4) age-by-APOE4 interaction
    # ti(age, APOE4, by=sex) three-way interaction age x APOE4 x sex
    # Include parametric sex and APOE terms to capture main effects
    
    model <- tryCatch(
      mgcv::gam(
        scaled_value ~
          # main effects
          s(age_at_visit, k = 10) +
          contributor_code + sex_ordered + ad + APOE4_count_pred +
          # lower-order interactions
          s(age_at_visit, k = 10, by = sex_ordered) +
          ti(age_at_visit, APOE4_count_pred, k = c(10, 3), bs = c("cr", "cr")) +
          # three-way interaction: age x APOE4 x sex
          ti(age_at_visit, APOE4_count_pred, by = sex_ordered, k = c(10, 3), bs = c("cr", "cr")),
        data = df_var,
        method = "REML"
      ),
      error = function(e) NULL
    )
  }else if(pred_type=="actual"){
    model <- tryCatch(
      gam(scaled_value ~ s(age_at_visit, k = 10) + s(age_at_visit, k = 10,by = APOE4_count)+sex+ad, data = df_var, method = "REML"),
      error = function(e) NULL
    )
  }
  
  if (is.null(model)) return(NULL)
  
  library(broom)
  gam_summ <- summary(model)
  tidy_param <- tidy(model, parametric = TRUE, conf.int = TRUE, conf.level = 0.95)
  tidy_smooth <- tidy(model, parametric = FALSE)
  tidy_all <- bind_rows(
    tidy_param %>% mutate(component = "parametric"),
    tidy_smooth %>% mutate(component = "smooth")
  )
  
  # Full Spearman correlation statistics
  cor_result <- cor.test(df_var$age_at_visit, df_var$scaled_value, method = "spearman") %>% broom::tidy()
  
  if(pred_type=="pred"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0,1),#mean(df_var$sex, na.rm = TRUE),
      APOE4_count_pred = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      contributor_code = model$xlevels$contributor_code,
      ad = 0
    )%>%mutate(APOE4_ordered  = factor(APOE4_count_pred),
               sex_ordered = ordered(sex))
  }else if(pred_type=="actual"){
    pred_df <- expand.grid(
      age_at_visit = age_grid,
      sex = c(0.5),#mean(df_var$sex, na.rm = TRUE),
      APOE4_count = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      ad = 0
    )
  }
  
  pred <- predict(model, newdata = pred_df, se.fit = TRUE)
  pred_df = pred_df%>%
    mutate(fit = pred$fit,
           se = pred$se.fit,
           lower = pred$fit - 1.96 * pred$se.fit,
           upper = pred$fit + 1.96 * pred$se.fit,
           plot_label = var_name)
  pred_df = pred_df%>%
    group_by(plot_label)%>%
    nest()%>%
    mutate(tidy_all=list(tidy_all))%>%
    mutate(r_squared = gam_summ$r.sq,
           dev_explained = gam_summ$dev.expl,
           n_obs = nrow(df_var),
           spearman_rho = cor_result$estimate,
           spearman_statistic = cor_result$statistic,
           spearman_p_value = cor_result$p.value)%>%
    ungroup()%>%
    mutate(pred_type=pred_type)
  return(pred_df)
}

# Get all unique variables
all_vars <- unique(long_scaled$plot_label)

library(furrr)
library(purrr)
library(dplyr)

# 1) Choose core count safely
n_cores <- N_CORES  # respects quotas (SLURM/cgroups)

# 2) Plan a multisession cluster (works on Linux; use multicore if you prefer)
plan(multisession, workers = n_cores)


# Below this age the predicted phenotypes are out-of-training-range extrapolations
# (ROSMAP training age >= ~56); the <EXTRAP_AGE segment is drawn as a FADED solid line
# in the line panels and marked with a boundary line in the heatmap. The same age is
# the cut for the `ge65` variant below, so the figure fades exactly the range the
# sensitivity analysis drops.
EXTRAP_AGE <- 65

# ---- three analyses, ONE refit -------------------------------------------
# full, rank and ge65 are all fitted here, in this run, from the same
# `long_scaled`. That is the point of doing it this way: the primary numbers,
# the two floor-robust re-analyses AND the figure below then all come from a
# single fit, so they can be put side by side in Supplementary Table 8 without
# having to reconcile two runs.
#   full  primary analysis -- min-max scaled outcome, all ages.
#   rank  rank-INT -- the outcome is replaced by its inverse-normal rank scores
#         within phenotype, over the full age range. The test then sees the
#         outcome only through its ORDER, so a monotone compression of the young
#         end is invisible to it. n is unchanged, so this costs no power.
#   ge65  primary outcome, chronological age >= EXTRAP_AGE only. Drops the
#         compressed young participants, but also costs power and shortens the
#         age span over which an age-slope interaction is estimated -- a real
#         effect can stay real and still lose significance here, so read this
#         together with `rank` rather than on its own.
VARIANTS <- c("full", "rank", "ge65")

rank_int <- function(x) {                  # inverse-normal rank (Blom) transform
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(x)))
}

make_variant <- function(data, variant) {
  switch(variant,
    full = data,
    ge65 = data %>% filter(age_at_visit >= EXTRAP_AGE),
    rank = data %>%                        # within phenotype, over the full range
      group_by(plot_label) %>%
      mutate(scaled_value = rank_int(scaled_value)) %>%
      ungroup(),
    stop("unknown variant: ", variant)
  )
}

# The three-way apoe_sex model is not reported (not in Fig 5, not in
# Supplementary Table 8), so fit_all_models() does not fit it.
# fit_gam_trajectory_apoe_sex() is left defined.
fit_all_models <- function(dat) {
  list(
    base = furrr::future_map_dfr(all_vars, ~ fit_gam_trajectory(dat, .x), .progress = TRUE),
    apoe = furrr::future_map_dfr(all_vars, ~ fit_gam_trajectory_apoe(dat, .x, "pred"), .progress = TRUE),
    sex  = furrr::future_map_dfr(all_vars, ~ fit_gam_trajectory_sex(dat, .x, "pred"), .progress = TRUE)
  )
}

# REPLOT_ONLY = TRUE reuses the fits already saved by a previous run instead of
# re-fitting: set it when only the FIGURES need changing (panel styling, colours,
# row order, which phenotypes are shown). The GAMs are the expensive part; data
# preparation above still runs, which is minutes rather than hours.
#   FALSE -> fit everything, write gam_stats_all_variants.csv (the full run)
#   TRUE  -> readRDS the saved trajectories, rebuild the figures only
# The statistics tables are NOT rewritten under REPLOT_ONLY, so a replot can
# never silently change the numbers in Supplementary Table 8.
REPLOT_ONLY <- FALSE

RES_DIR <- step_dir("03_age_trajectory")

if (REPLOT_ONLY) {
  message("REPLOT_ONLY: loading saved fits from ", RES_DIR, " -- no GAMs are fitted")
  need <- file.path(RES_DIR, c("trajectories.rds", "trajectories_apoe.rds", "trajectories_sex.rds"))
  if (!all(file.exists(need))) {
    stop("REPLOT_ONLY needs the saved fits, missing:\n  ",
         paste(need[!file.exists(need)], collapse = "\n  "))
  }
  traj              <- NULL          # guards the stats export in section 8
  trajectories      <- readRDS(need[1])
  trajectories_apoe <- readRDS(need[2])
  trajectories_sex  <- readRDS(need[3])
} else {
  traj <- lapply(setNames(VARIANTS, VARIANTS), function(v) {
    dat <- make_variant(long_scaled, v)
    message(sprintf("fitting variant '%s' (%d rows, %d phenotypes)",
                    v, nrow(dat), dplyr::n_distinct(dat$plot_label)))
    fit_all_models(dat)
  })

  # Everything downstream -- onset ages, all panels, Fig 5 -- is the primary
  # `full` variant. The other two exist only for the statistics table.
  trajectories          <- traj$full$base
  trajectories_apoe     <- traj$full$apoe
  trajectories_sex      <- traj$full$sex
}

# 4) Optional: reset plan to sequential
plan(sequential)


# Add category. The saved .rds already carry Category (they are saved after this
# join), so joining again under REPLOT_ONLY would produce Category.x/Category.y
# and break every panel that facets or colours by Category.
category_map = long_scaled%>%ungroup()%>%dplyr::select(Category,plot_label)%>%unique()
add_category <- function(df) {
  if ("Category" %in% names(df)) df else left_join(df, category_map, by = "plot_label")
}
trajectories      <- add_category(trajectories)
trajectories_apoe <- add_category(trajectories_apoe)
trajectories_sex  <- add_category(trajectories_sex)

# 3. ESTIMATE ONSET AGES (% Abnormal Method) ######

# Reference levels shown as guide lines in panel A and as crossing points in panel C,
# so the two panels correspond (same set + same shades in both).
LEVELS <- c(0.5)
lvcol  <- setNames(c("grey30"), as.character(LEVELS))
# Panel C's onset markers use their own, lighter fill: on the heatmap they sit as
# filled circles with a black outline, and grey30 reads as a dark blob against
# the tiles. Panel A's dashed guide line keeps lvcol -- a thin line needs the
# darker grey to stay visible. The two still mark the same LEVELS.
lvfill_c <- setNames(c("grey75"), as.character(LEVELS))
# EXTRAP_AGE is defined in section 2, where it is also the cut for the `ge65`
# sensitivity variant.

estimate_onset <- function(data, trajectories, var_name, threshold = 0.2, pct_cutoff = 0.10) {

  df_var <- data %>% filter(plot_label == var_name)
  
  result <- df_var %>%
    mutate(age_bin = cut(age_at_visit, breaks = seq(10, 100, 10), right=FALSE,
                         labels = seq(20, 100, 10))) %>%
    filter(!is.na(age_bin)) %>%
    mutate(abnormal = scaled_value > threshold) %>%
    group_by(age_bin) %>%
    summarize(
      pct_abnormal = mean(abnormal, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(age_bin = as.numeric(as.character(age_bin)))
  
  # Find the first age bin from which the fraction abnormal stays above pct_cutoff
  above_threshold <- result$pct_abnormal > pct_cutoff
  onset <- NA_real_
  for (i in seq_along(result$age_bin)) {
    if (all(above_threshold[i:length(above_threshold)])) {
      onset <- result$age_bin[i]
      break
    }
  }
  
  
  # Also get GAM-based onset for comparison
  traj <- trajectories %>% filter(plot_label == var_name)%>%
    dplyr::select(-tidy_all)%>%
    unnest()%>%
    group_by(age_at_visit,spearman_rho)%>%summarise(fit=mean(fit))%>%
    ungroup()%>%
    mutate(fit = (fit - min(fit)) / (max(fit) - min(fit)))
  
  above_threshold <- traj$fit > threshold
  
  onset_gam <- NA_real_
  for (i in seq_along(traj$age_at_visit)) {
    if (all(above_threshold[i:length(above_threshold)])) {
      onset_gam <- traj$age_at_visit[i]
      break
    }
  }

  # crossing age for each display LEVEL (first age from which fit stays above the level)
  cross_lv <- sapply(LEVELS, function(L) {
    ab <- traj$fit > L; a <- NA_real_
    for (i in seq_along(traj$age_at_visit))
      if (all(ab[i:length(ab)])) { a <- traj$age_at_visit[i]; break }
    a
  })

  tibble(
    plot_label = var_name,
    onset_pct = onset,
    onset_gam = onset_gam,
    max_pct_abnormal = max(result$pct_abnormal, na.rm = TRUE),
    spearman_rho = traj$spearman_rho[1],
    !!!setNames(as.list(cross_lv), paste0("onset_L", seq_along(LEVELS)))
  )
}

# Estimate onset for all variables
onset_results <- map_df(all_vars, ~ estimate_onset(long_scaled,trajectories, .x, threshold = 0.5,pct_cutoff = 0.50))

# Add category and sort by onset
onset_results <- onset_results %>%
  left_join(category_map, by = "plot_label") %>%
  arrange(onset_pct)

# 4. PANEL A (Fig. 5a): Key Representative Trajectories ####

# Select key variables (one or two per category)
key_vars <- c( "Body mass index","Tangle density", "ADL score","Number of medications" ,
               "Global cognitive function","Motor function composite")[-1]
category_colors_vars = category_colors
indx = match(names(category_colors_vars),effect_sizes$Category[match(key_vars,effect_sizes$plot_label)])
names(category_colors_vars)[!is.na(indx)] = key_vars[indx[!is.na(indx)]]

# Filter trajectories
traj_key <- trajectories %>%
  filter(plot_label %in% key_vars) %>%
  mutate(plot_label = factor(plot_label, levels = key_vars))%>%
  dplyr::select(-tidy_all)%>%
  unnest()%>%
  group_by(Category,plot_label,age_at_visit)%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))


panel_a <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = plot_label, fill = plot_label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(data = function(d) d[d$age_at_visit >= EXTRAP_AGE, ], linewidth = 1.2) +
  geom_line(data = function(d) d[d$age_at_visit <= EXTRAP_AGE, ], linewidth = 1.2, alpha = 0.30) +
  geom_hline(yintercept = LEVELS, linetype = "dashed", color = lvcol, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
  scale_color_manual(values = category_colors_vars) +
  scale_fill_manual(values = category_colors_vars) +
  scale_x_continuous(breaks = seq(20, 100, 10)) +
  scale_y_continuous(limits = c(-0.3, 1.2), breaks = seq(0, 1., 0.2)) +
  labs(
    x = "Age (years)",
    y = "Scaled value",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none",#c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 10, 5, 5)
  )

panel_a

# 6. PANEL C (Fig. 5b): Heatmap of Trajectories (selected phenotypes) #####
# Create heatmap matrix 
age_bins <- seq(20, 90, 1)
key_vars = c("Macroscopic Cerebral Infarcts","Lewy Body disease","Hippocampal sclerosis",
             "Body mass index","ADL score","MMSE","Global cognitive function",
             "Motor function composite","IADL score","Parkinsonian signs (UPDRS)",
             "Mobility disability","Number of medications","Global burden of AD pathology",
             "Diastolic blood pressure","Cerebral Atherosclerosis","Amyloid",
             "Tangle density","Nigral Neuronal Loss","TDP-43","Systolic blood pressure",
             "Gross Cerebral Infarcts","Cerebral amyloid angiopathy","Arteriolosclerosis")
heatmap_df <- trajectories %>%
  filter(!is.na(Category)) %>%
  filter(spearman_rho>0.1)%>%
  dplyr::select(-tidy_all)%>%
  filter(plot_label%in%key_vars)%>%
  unnest()%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))%>%
  group_by(plot_label, Category) %>%
  reframe(
    age_bins = list(age_bins),
    fit_interp = list(approx(age_at_visit, fit, xout = unlist(age_bins))$y)
  ) %>%
  unnest(cols = c(age_bins, fit_interp)) %>%
  dplyr::rename(age = age_bins, fit = fit_interp)

# Order by onset age
var_order <- onset_results %>%
  filter(spearman_rho>0.1)%>%
  arrange(desc(onset_gam)) %>%
  filter(plot_label%in%key_vars)%>%
  pull(plot_label)

heatmap_df <- heatmap_df %>%
  filter(plot_label %in% var_order) %>%
  mutate(plot_label = factor(plot_label, levels = rev(var_order)))


# crossing points at the SAME levels as panel A's guide lines (one point per level)
onset_points_long <- onset_results %>%
  filter(plot_label %in% var_order) %>%
  dplyr::select(plot_label, dplyr::all_of(paste0("onset_L", seq_along(LEVELS)))) %>%
  tidyr::pivot_longer(-plot_label, names_to = "lvid", values_to = "age") %>%
  mutate(level = LEVELS[as.integer(gsub("onset_L", "", lvid))],
         plot_label = factor(plot_label, levels = rev(var_order))) %>%
  filter(!is.na(age))

panel_c <- ggplot(heatmap_df, aes(x = age, y = plot_label, fill = fit)) +
  geom_tile() +
  # extrapolation boundary (heatmap has no curve to fade)
  geom_vline(xintercept = EXTRAP_AGE, color = "gray40", linewidth = 0.4) +
  # one geom_point layer per LEVEL, at the same LEVELS as panel A's dashed lines
  lapply(LEVELS, function(L) geom_point(
    data = dplyr::filter(onset_points_long, level == L),
    aes(x = age, y = plot_label), shape = 21, fill = lvfill_c[as.character(L)],
    color = "black", size = 2.6, stroke = 0.7, inherit.aes = FALSE)) +
  # hjust offsets the label by that many text-widths from the point; -0.9 moves
  # the number clear of the marker circle.
  geom_text(data = dplyr::filter(onset_points_long, level == 0.5),
            aes(x = age, y = plot_label, label = round(age, 0)),
            hjust = -0.9, size = 2.2, color = "black", inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "#F7F7F7", 
    high = "#B2182B",
    midpoint = 0.3,
    limits = c(-0.3, 1.3),
    name = "Scaled value"
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(20, 90, 10)) +
  labs(
    x = "Age (years)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    legend.key.height = unit(.5, "cm"),
    legend.title = element_text(size = 7),  # title font size
    legend.text  = element_text(size = 7),   # item labels font size
    plot.margin = margin(5, 10, 5, 5),
    axis.title = element_text(size = 10)
  )
panel_c

# Add category color strip on left side
category_strip <- onset_results %>%
  filter(plot_label %in% var_order) %>%
  mutate(plot_label = factor(plot_label, levels = rev(var_order))) %>%
  ggplot(aes(x = 1, y = plot_label, fill = Category)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = category_colors) +
  theme_void() +
  theme(legend.position = "none")

# APOE panel (Fig. 5c) ####

key_vars = c("Macroscopic Cerebral Infarcts","Lewy Body disease","Hippocampal sclerosis",
             "Body mass index","ADL score","MMSE","Global cognitive function",
             "Motor function composite","IADL score","Parkinsonian signs (UPDRS)",
             "Mobility disability","Number of medications","Global burden of AD pathology",
             "Diastolic blood pressure","Cerebral Atherosclerosis","Amyloid",
             "Tangle density","Nigral Neuronal Loss","TDP-43","Systolic blood pressure",
             "Gross Cerebral Infarcts","Cerebral amyloid angiopathy","Arteriolosclerosis")
trajectories_apoe%>%
  dplyr::select(-data)%>%
  unnest()%>%
  filter(grepl("\\,",term))%>%
  arrange((p.value))%>%
  filter(plot_label%in%key_vars)%>%
  filter(p.value<=0)%>%
  pull(plot_label)%>%.[1:min(10,length(.))]%>%unique() -> key_vars
key_vars = c("Tangle density","Amyloid", "Global burden of AD pathology",
             "Global cognitive function", "Motor function composite",
             "Mobility disability","Number of medications")

# Map category colours to the key variables
category_colors_vars = category_colors
indx = match(names(category_colors_vars),effect_sizes$Category[match(key_vars,effect_sizes$plot_label)])
names(category_colors_vars)[!is.na(indx)] = key_vars[indx[!is.na(indx)]]

# Filter trajectories
traj_key <- trajectories_apoe %>%
  filter(plot_label %in% key_vars) %>%
  mutate(plot_label = factor(plot_label, levels = key_vars))%>%
  dplyr::select(-tidy_all)%>%
  unnest()%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))


pvals_df <- trajectories_apoe%>%
  dplyr::select(-data)%>%
  unnest()%>%
  filter(grepl("\\,",term))%>%
  arrange(desc(statistic))%>%
  filter(plot_label%in%key_vars)

# Choose annotation positions per facet
# Here we place text near the top-right of each panel.
annot_positions <- traj_key %>%
  group_by(plot_label) %>%
  summarize(x_pos = 50,
            y_pos = max(upper %||% fit, na.rm = TRUE)) %>%
  ungroup()

# Build labels with exact P values
pval_labels <- pvals_df %>%
  left_join(annot_positions, by = "plot_label") %>%
  mutate(
    # p only (no F). mgcv returns 0 (or a tiny negative) when its smooth-term P
    # is below its numerical accuracy -- seen for true P up to ~2.2e-6 -- so a P
    # of 0 is printed as "p < 3e-6" (not "< 2e-16"); other P values as computed
    p_text = ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g"))),
    label_text = p_text
  )%>%
  mutate(plot_label = factor(plot_label, levels = key_vars))



panel_d <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(APOE4_count_pred),
                                fill = factor(APOE4_count_pred))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(data = function(d) d[d$age_at_visit >= EXTRAP_AGE, ], linewidth = 1.2) +
  geom_line(data = function(d) d[d$age_at_visit <= EXTRAP_AGE, ], linewidth = 1.2, alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_color_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  scale_fill_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  labs(
    x = "Age (years)",
    y = "Scaled value",
    color = NULL,
    fill = NULL
  ) +
  facet_wrap(~plot_label,ncol=7,scales="free_y")+
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none",#c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    axis.title = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )+
  geom_text(data = pval_labels,
            aes(x = x_pos, y = y_pos, label = label_text),
            inherit.aes = FALSE,
            hjust = 1, vjust = 1,
            size = 2,
            lineheight = 0.95,
            color = "black")


panel_d

# APOE panel, all phenotypes (Extended Data Fig. 5) ####
trajectories_apoe%>%
  arrange(Category)%>%pull(plot_label)%>%setdiff(.,c("Age","ADPRS ensemble (averaging)" )) -> key_vars
key_vars = key_vars[!grepl("slope",key_vars)]

# Map category colours to the key variables
category_colors_vars = category_colors
indx = match(names(category_colors_vars),effect_sizes$Category[match(key_vars,effect_sizes$plot_label)])
names(category_colors_vars)[!is.na(indx)] = key_vars[indx[!is.na(indx)]]

# Filter trajectories
traj_key <- trajectories_apoe %>%
  filter(plot_label %in% key_vars) %>%
  mutate(plot_label = factor(plot_label, levels = key_vars))%>%
  dplyr::select(-tidy_all)%>%
  unnest()%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))


pvals_df <- trajectories_apoe%>%
  dplyr::select(-data)%>%
  unnest()%>%#View
  filter(grepl("\\,",term))%>%
  arrange(desc(statistic))%>%
  filter(plot_label%in%key_vars)

# Choose annotation positions per facet
# Here we place text near the top-right of each panel.
annot_positions <- traj_key %>%
  group_by(plot_label) %>%
  summarize(x_pos = 50,
            y_pos = max(upper %||% fit, na.rm = TRUE)) %>%
  ungroup()

# Build labels with exact P values
pval_labels <- pvals_df %>%
  left_join(annot_positions, by = "plot_label") %>%
  mutate(
    # p only (no F). mgcv returns 0 (or a tiny negative) when its smooth-term P
    # is below its numerical accuracy -- seen for true P up to ~2.2e-6 -- so a P
    # of 0 is printed as "p < 3e-6" (not "< 2e-16"); other P values as computed
    p_text = ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g"))),
    label_text = p_text
  )%>%
  ungroup()%>%
  arrange(p.value)%>%
  mutate(plot_label = factor(plot_label, levels = plot_label))

traj_key$plot_label = factor(traj_key$plot_label,levels = pval_labels$plot_label%>%levels)

panel_d_sub <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(APOE4_count_pred),
                                    fill = factor(APOE4_count_pred))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(data = function(d) d[d$age_at_visit >= EXTRAP_AGE, ], linewidth = 1.2) +
  geom_line(data = function(d) d[d$age_at_visit <= EXTRAP_AGE, ], linewidth = 1.2, alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_color_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  scale_fill_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  labs(
    x = "Age (years)",
    y = "Scaled value",
    color = NULL,
    fill = NULL
  ) +
  facet_wrap(~plot_label+Category,ncol=6,scales="free_y")+
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none",#c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    axis.title = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )+
  geom_text(data = pval_labels,
            aes(x = x_pos, y = y_pos, label = label_text),
            inherit.aes = FALSE,
            hjust = 1, vjust = 1,
            size = 2,
            lineheight = 0.95,
            color = "black")


panel_d_sub

# Sex x age panel (Fig. 5d) ####
key_vars = c("Macroscopic Cerebral Infarcts","Lewy Body disease","Hippocampal sclerosis",
             "Body mass index","ADL score","MMSE","Global cognitive function",
             "Motor function composite","IADL score","Parkinsonian signs (UPDRS)",
             "Mobility disability","Number of medications","Global burden of AD pathology",
             "Diastolic blood pressure","Cerebral Atherosclerosis","Amyloid",
             "Tangle density","Nigral Neuronal Loss","TDP-43","Systolic blood pressure",
             "Gross Cerebral Infarcts","Cerebral amyloid angiopathy","Arteriolosclerosis")
trajectories_sex%>%
  dplyr::select(-data)%>%
  unnest()%>%
  filter(grepl("\\:",term))%>%
  arrange((p.value))%>%
  filter(plot_label%in%key_vars)%>%
  filter(p.value<=0)%>%
  pull(plot_label)%>%.[1:min(10,length(.))]%>%unique() -> key_vars
key_vars=c(key_vars,"Amyloid","Global cognitive function")
key_vars = c("Tangle density","Amyloid", "Global burden of AD pathology",
             "Global cognitive function", "Motor function composite",
             "Mobility disability","Number of medications")
# Filter trajectories
traj_key <- trajectories_sex %>%
  filter(plot_label %in% key_vars) %>%
  mutate(plot_label = factor(plot_label, levels = key_vars))%>%
  dplyr::select(-tidy_all)%>%
  unnest()%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))



pvals_df <- trajectories_sex%>%
  dplyr::select(-data)%>%
  unnest()%>%
  filter(grepl("\\:",term))%>%
  arrange(desc(statistic))%>%
  filter(plot_label%in%key_vars)

# Choose annotation positions per facet
# Here we place text near the top-right of each panel.
annot_positions <- traj_key %>%
  group_by(plot_label) %>%
  summarize(x_pos = 50,
            y_pos = max(upper %||% fit, na.rm = TRUE)) %>%
  ungroup()

# Build labels with exact P values
pval_labels <- pvals_df %>%
  left_join(annot_positions, by = "plot_label") %>%
  mutate(
    # p only (no F). mgcv returns 0 (or a tiny negative) when its smooth-term P
    # is below its numerical accuracy -- seen for true P up to ~2.2e-6 -- so a P
    # of 0 is printed as "p < 3e-6" (not "< 2e-16"); other P values as computed
    p_text = ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g"))),
    label_text = p_text
  )%>%
  mutate(plot_label = factor(plot_label, levels = key_vars))



panel_e <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(sex),
                                fill = factor(sex))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(data = function(d) d[d$age_at_visit >= EXTRAP_AGE, ], linewidth = 1.2) +
  geom_line(data = function(d) d[d$age_at_visit <= EXTRAP_AGE, ], linewidth = 1.2, alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_color_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  scale_fill_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  labs(
    x = "Age (years)",
    y = "Scaled value",
    color = NULL,
    fill = NULL
  ) +
  facet_wrap(~plot_label,ncol=7,scales="free_y")+
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none",#c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    axis.title = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )+
  geom_text(data = pval_labels,
            aes(x = x_pos, y = y_pos, label = label_text),
            inherit.aes = FALSE,
            hjust = 1, vjust = 1,
            size = 2,
            lineheight = 0.95,
            color = "black")

panel_e


# Sex x age panel, all phenotypes (Extended Data Fig. 6) ####
trajectories_sex%>%
  arrange(Category)%>%pull(plot_label)%>%setdiff(.,c("Age","ADPRS ensemble (averaging)" )) -> key_vars
key_vars = key_vars[!grepl("slope",key_vars)]

# Filter trajectories
traj_key <- trajectories_sex %>%
  filter(plot_label %in% key_vars) %>%
  mutate(plot_label = factor(plot_label, levels = key_vars))%>%
  dplyr::select(-tidy_all)%>%
  unnest()%>%
  filter(contributor_code=="A_EDTA Plasma")%>%
  group_by(plot_label)%>%
  mutate(lower=(lower - min(fit)) / (max(fit) - min(fit)),
         upper=(upper - min(fit)) / (max(fit) - min(fit)),
         fit = (fit - min(fit)) / (max(fit) - min(fit)))



pvals_df <- trajectories_sex%>%
  dplyr::select(-data)%>%
  unnest()%>%
  filter(grepl("\\:",term))%>%
  arrange(desc(statistic))%>%
  filter(plot_label%in%key_vars)

# Choose annotation positions per facet
# Here we place text near the top-right of each panel.
annot_positions <- traj_key %>%
  group_by(plot_label) %>%
  summarize(x_pos = 50,
            y_pos = max(upper %||% fit, na.rm = TRUE)) %>%
  ungroup()

# Build labels with exact P values
pval_labels <- pvals_df %>%
  left_join(annot_positions, by = "plot_label") %>%
  mutate(
    # p only (no F). mgcv returns 0 (or a tiny negative) when its smooth-term P
    # is below its numerical accuracy -- seen for true P up to ~2.2e-6 -- so a P
    # of 0 is printed as "p < 3e-6" (not "< 2e-16"); other P values as computed
    p_text = ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g"))),
    label_text = p_text
  )%>%
  arrange(p.value)%>%
  mutate(plot_label = factor(plot_label, levels = plot_label))

traj_key$plot_label = factor(traj_key$plot_label,levels = pval_labels$plot_label%>%levels)


panel_e_sub <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(sex),
                                    fill = factor(sex))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(data = function(d) d[d$age_at_visit >= EXTRAP_AGE, ], linewidth = 1.2) +
  geom_line(data = function(d) d[d$age_at_visit <= EXTRAP_AGE, ], linewidth = 1.2, alpha = 0.30) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_color_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  scale_fill_manual(
    values = c("0" = "#009E73", "1" = "#E69F00", "2" = "#D55E00"),
    labels = c("0" = "Non-carrier", "1" = "Heterozygote", "2" = "Homozygote"),
    name = "APOE ε4"
  ) +
  labs(
    x = "Age (years)",
    y = "Scaled value",
    color = NULL,
    fill = NULL
  ) +
  facet_wrap(~plot_label+Category,ncol=6,scales="free_y")+
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none",#c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    axis.title = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )+
  geom_text(data = pval_labels,
            aes(x = x_pos, y = y_pos, label = label_text),
            inherit.aes = FALSE,
            hjust = 1, vjust = 1,
            size = 2,
            lineheight = 0.95,
            color = "black")

panel_e_sub

# 7. COMBINE PANELS ####
# Top row, right: category strip + heatmap
top_row_right <- plot_grid(
  category_strip,
  panel_c,
  ncol = 2,
  rel_widths = c(0.03, 1),
  align = "h"
)

top_row <- plot_grid(
  panel_a,
  top_row_right,
  nrow = 1,
  rel_widths = c(0.7,1)
)


final_figure = plot_grid(
  top_row,
  panel_d,
  panel_e,
  nrow = 3,
  rel_heights = c(1.1,1,1)
)

# 8. SAVE #####
print(final_figure)
dir.create(RES_DIR, showWarnings = FALSE)

# ORDER MATTERS: the fits take hours, the figures take seconds, and a broken
# graphics device on the workbench kills the script at the first ggsave. Data
# first, pictures last, so a device problem can never cost a run.
#
# The workbench box throws "Graphics API version mismatch" on PNG: ggplot2 uses
# ragg for PNG when it is installed, and the installed ragg (1.4.0) was built
# against a different R graphics engine than the running R (4.5.1). It cannot be
# rebuilt there -- ragg needs freetype2/libpng/libtiff/libjpeg/libwebp headers,
# which are not available. So on a PNG failure fall back to base R's png device,
# which is part of R itself and therefore always matches. PDF goes through
# grDevices and is unaffected either way.
png_base <- function(filename, width, height, ...) {
  grDevices::png(filename, width = width, height = height,
                 units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib",
                 bg = "white")
}

safe_ggsave <- function(filename, plot, ...) {
  ok <- tryCatch({ ggsave(filename, plot, ...); TRUE },
                 error = function(e) {
                   message("  ggsave failed for ", basename(filename), ": ",
                           conditionMessage(e))
                   FALSE
                 })
  if (ok || !grepl("\\.png$", filename, ignore.case = TRUE)) return(invisible(ok))

  tryCatch({ ggsave(filename, plot, ..., device = png_base)
             message("  -> recovered via grDevices::png")
             invisible(TRUE) },
           error = function(e) {
             warning("both PNG devices failed for ", basename(filename), ": ",
                     conditionMessage(e), call. = FALSE)
             invisible(FALSE)
           })
}

saveRDS(onset_results,file=file.path(RES_DIR, "onset_results.rds"))
# the trajectory .rds are written in the guarded block below, so a REPLOT_ONLY
# pass cannot overwrite the saved fits with the copies it just read back in

# Export GAM statistics as CSV tables
out_dir <- RES_DIR

# Helper: extract tidy stats + model-level summaries from trajectory objects
extract_gam_stats <- function(traj_df) {
  traj_df %>%
    dplyr::select(plot_label, tidy_all, r_squared, dev_explained, n_obs,
                  spearman_rho, spearman_statistic, spearman_p_value) %>%
    unnest(tidy_all)
}

# The three analyses in one tidy table. `variant` is full / rank / ge65 and
# `model` is base / apoe / sex, so Supplementary Table 8 can take its primary,
# rank-INT and age >= 65 columns from this single file -- same run, same data
# preparation, directly comparable. The per-model CSVs written below are the
# `full` variant.
if (is.null(traj)) {
  message("REPLOT_ONLY: statistics tables left untouched")
} else {
  gam_stats_all <- purrr::imap_dfr(traj, function(fits, variant) {
    purrr::imap_dfr(fits[c("base", "apoe", "sex")], function(traj_df, model) {
      extract_gam_stats(traj_df) %>% mutate(model = model, .before = 1)
    }) %>% mutate(variant = variant, .before = 1)
  })
  write.csv(gam_stats_all,
            file.path(out_dir, "gam_stats_all_variants.csv"), row.names = FALSE)

  write.csv(extract_gam_stats(trajectories),
            file.path(out_dir, "gam_stats_combined.csv"), row.names = FALSE)
  write.csv(extract_gam_stats(trajectories_apoe),
            file.path(out_dir, "gam_stats_apoe.csv"), row.names = FALSE)
  write.csv(extract_gam_stats(trajectories_sex),
            file.path(out_dir, "gam_stats_sex.csv"), row.names = FALSE)

  saveRDS(trajectories,file=file.path(RES_DIR, "trajectories.rds"))
  saveRDS(trajectories_apoe,file=file.path(RES_DIR, "trajectories_apoe.rds"))
  saveRDS(trajectories_sex,file=file.path(RES_DIR, "trajectories_sex.rds"))
}
write.csv(onset_results,
          file.path(out_dir, "onset_results.csv"), row.names = FALSE)

# panel_c's plotting input, so the heatmap can be redrawn without the
# trajectory .rds and without re-fitting. Aggregate only:
# per phenotype, the 0-1 normalised fitted curve interpolated onto 20-90 by 1 yr.
# Row order and the onset circles come from onset_results.csv above.
write.csv(heatmap_df,
          file.path(out_dir, "panelC_heatmap_df.csv"), row.names = FALSE)

# Source Data: exactly what each published panel draws -- the plot's own
# data (GAM predictions on the age grid, 0-1 scaled as plotted) plus the data of
# any layer that brings its own table (P-value labels, onset markers). Aggregate
# only: no per-participant rows, so these may leave the workbench.
write_source_data <- function(p, stem) {
  flat <- function(d) d %>% ungroup() %>% dplyr::select(where(~ !is.list(.x)))
  write.csv(flat(p$data), file.path(out_dir, paste0("source_", stem, "_data.csv")),
            row.names = FALSE)
  own <- Filter(function(l) is.data.frame(l$data), p$layers)
  for (i in seq_along(own)) {
    geom <- tolower(sub("^Geom", "", class(own[[i]]$geom)[1]))
    write.csv(flat(own[[i]]$data),
              file.path(out_dir, sprintf("source_%s_%s%d.csv", stem, geom, i)),
              row.names = FALSE)
  }
}
write_source_data(panel_a,     "fig5a")
write_source_data(panel_c,     "fig5b")
write_source_data(panel_d,     "fig5c")
write_source_data(panel_e,     "fig5d")
write_source_data(panel_d_sub, "figS7")   # Extended Data Fig. 5
write_source_data(panel_e_sub, "figS8")   # Extended Data Fig. 6

cat("data written to", out_dir, "-- figures next\n")

# ---- figures last (see safe_ggsave above) ----
safe_ggsave(file.path(RES_DIR, "figure_trajectories_combined.pdf"), final_figure,
            width = 11, height = 6, dpi = 300)
safe_ggsave(file.path(RES_DIR, "figure_trajectories_combined.png"), final_figure,
            width = 11, height = 6, dpi = 300)

safe_ggsave(file.path(RES_DIR, "figure_trajectories_apoe_all.pdf"), panel_d_sub,
            width = 11, height = 6, dpi = 300)
safe_ggsave(file.path(RES_DIR, "figure_trajectories_apoe_all.png"), panel_d_sub,
            width = 11, height = 6, dpi = 300)


safe_ggsave(file.path(RES_DIR, "figure_trajectories_sex_all.pdf"), panel_e_sub,
            width = 11, height = 6, dpi = 300)
safe_ggsave(file.path(RES_DIR, "figure_trajectories_sex_all.png"), panel_e_sub,
            width = 11, height = 6, dpi = 300)
