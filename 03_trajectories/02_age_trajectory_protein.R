# 03_trajectories/02_age_trajectory_protein.R
# Age trajectories (GAMs) of measured plasma proteins (hard-coded SomaScan
# aptamers linked to amyloid, tangles and global AD pathology) in GNPC:
# Extended Data Fig. 7.
# Inputs:  PROTEOMICS_RDS, PREDICTED_RDS (predicted sex and APOE; see config.R)
# Outputs: step_dir("03_age_trajectory_protein") -- GAM statistics, fitted
#          trajectories (.rds), Source Data CSVs and figures.
# Runs on the AD Workbench. From the repository root:
#   Rscript 03_trajectories/02_age_trajectory_protein.R

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
  amyloid = c("seq.15321.8","seq.9599.6","seq.22782.80","seq.10462.14","seq.17685.9",
              "seq.15511.37","seq.14318.1","seq.8997.4","seq.10980.11","seq.25082.3")
  tangle = c("seq.15321.8","seq.8997.4","seq.9599.6",
             "seq.14114.18","seq.10980.11","seq.17671.58",
             "seq.21887.2","seq.19622.7","seq.22782.80","seq.13102.1")
  gpath = c("seq.15321.8","seq.9599.6","seq.14114.18","seq.22782.80","seq.10980.11",
            "seq.11293.14","seq.25086.83","seq.13640.5","seq.7932.23","seq.3066.12")
  al_vars = c(amyloid,tangle,gpath)%>%unique()
  indx = rownames(expdata)%in%al_vars
  expdata=expdata[indx,]
  metadata <- colData(expdata) %>% as.data.frame()
  
  # Protein levels: negative values set to NA, then natural log
  assay(expdata,"data")%>%t() -> prediction
  prediction[prediction<0]=NA
  prediction=log(prediction)
  prediction%>%
    as.data.frame()%>%
    rownames_to_column(var="sample_id") -> prediction
  
  
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
  
  # Prepare data: drop invalid ages
  prepared <- data_plasma %>%
    filter(age_at_visit != -1) %>%
    dplyr::select(age_at_visit,ad,sex,sex_pred,APOE,APOE4_count,APOE_pred,APOE4_count_pred,contributor_code, all_of(al_vars))
  
  # Long format. No 0-1 scaling here: scaled_value (set below) is the log protein level
  long_scaled <- prepared %>%
    pivot_longer(cols = all_of(al_vars), names_to = "key", values_to = "value") 
  
  # Multi-Variable Trajectory Visualization #######
  
  library(tidyverse)
  library(mgcv)
  library(cowplot)
  library(viridis)
  library(RColorBrewer)
  long_scaled$scaled_value=long_scaled$value
  
  long_scaled = long_scaled%>%
    group_by(key)%>%
    nest()%>%
    inner_join(.,    rowData(expdata)%>%as.data.frame()%>%
                 dplyr::select(probeid,symbol)%>%
                 mutate(amyloid = ifelse(probeid%in%amyloid,"A",""),
                        tangles = ifelse(probeid%in%tangle,"T",""),
                        gpath = ifelse(probeid%in%gpath,"G",""))%>%
                 mutate(symbol=paste0(symbol," (",probeid,"),",amyloid,tangles,gpath))%>%
                 dplyr::select(probeid,symbol),by=c("key"="probeid"))%>%
    mutate(key=symbol)%>%
    dplyr::select(-symbol)%>%
    unnest()
}


# 2. FIT GAMs AND GET PREDICTIONS ####

# Age grid for predictions
age_grid <- seq(
  min(long_scaled$age_at_visit, na.rm = TRUE),
  max(long_scaled$age_at_visit, na.rm = TRUE),
  by = 0.5
)

# Fit GAM for each variable and extract predictions

fit_gam_trajectory_apoe <- function(data, var_name,pred_type =c("actual")) {
  
  df_var <- data %>% filter(key == var_name)
  
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
      sex = c(0.5),#mean(df_var$sex, na.rm = TRUE),
      APOE4_count_pred = c(0,1,2),#mean(df_var$APOE4_count_pred, na.rm = TRUE),
      contributor_code = model$xlevels$contributor_code,
      ad = 0
    )%>%mutate(APOE4_ordered  = ordered(APOE4_count_pred))
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

fit_gam_trajectory_sex <- function(data, var_name,pred_type =c("actual")) {
  
  df_var <- data %>% filter(key == var_name)
  
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


# Get all unique variables
all_vars <- unique(long_scaled$key)

library(furrr)
library(purrr)
library(dplyr)

# 1) Choose core count safely
n_cores <- N_CORES  # respects quotas (SLURM/cgroups)

# 2) Plan a multisession cluster (works on Linux; use multicore if you prefer)
plan(multisession, workers = n_cores)


# Fit all trajectories
trajectories_apoe <- furrr::future_map_dfr(all_vars, ~ fit_gam_trajectory_apoe(long_scaled, .x,"pred"),.progress = TRUE)
trajectories_sex <- furrr::future_map_dfr(all_vars, ~ fit_gam_trajectory_sex(long_scaled, .x,"pred"),.progress = TRUE)

# 4) Optional: reset plan to sequential
plan(sequential)



# Filter trajectories
traj_key <- trajectories_apoe %>%
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
  arrange(desc(statistic))

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
    # mgcv returns 0 when its smooth-term P is below its accuracy (true P up to ~2.2e-6): print "p < 3e-6"
    p_text = paste0("F = ", signif(statistic, 3), "\n", ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g")))),
    label_text = p_text
  )



panel_d <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(APOE4_count_pred),
                                fill = factor(APOE4_count_pred))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
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

# Filter trajectories
traj_key <- trajectories_sex %>%
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
  arrange(desc(statistic))

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
    # mgcv returns 0 when its smooth-term P is below its accuracy (true P up to ~2.2e-6): print "p < 3e-6"
    p_text = paste0("F = ", signif(statistic, 3), "\n", ifelse(p.value <= 0, "p < 3e-6", paste0("p = ", formatC(p.value, digits = 2, format = "g")))),
    label_text = p_text
  )



panel_e <- ggplot(traj_key, aes(x = age_at_visit, y = fit, color = factor(sex),
                                fill = factor(sex))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
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

# 7. COMBINE PANELS ####
# Two rows: APOE panel, then sex panel


final_figure = plot_grid(
  panel_d,
  panel_e,
  nrow = 2,
  rel_heights = c(1,1)
)

# 8. SAVE #####
print(final_figure)
RES_DIR <- step_dir("03_age_trajectory_protein")
ggsave(file.path(RES_DIR, "figure_trajectories_combined.pdf"), final_figure,
       width = 11, height = 6, dpi = 300)
ggsave(file.path(RES_DIR, "figure_trajectories_combined.png"), final_figure,
       width = 11, height = 6, dpi = 300)



saveRDS(trajectories_apoe,file=file.path(RES_DIR, "trajectories_apoe.rds"))
saveRDS(trajectories_sex,file=file.path(RES_DIR, "trajectories_sex.rds"))

# Export GAM statistics as CSV tables
out_dir <- RES_DIR

extract_gam_stats <- function(traj_df) {
  traj_df %>%
    dplyr::select(plot_label, tidy_all, r_squared, dev_explained, n_obs,
                  spearman_rho, spearman_statistic, spearman_p_value) %>%
    unnest(tidy_all)
}

write.csv(extract_gam_stats(trajectories_apoe),
          file.path(out_dir, "gam_stats_apoe.csv"), row.names = FALSE)
write.csv(extract_gam_stats(trajectories_sex),
          file.path(out_dir, "gam_stats_sex.csv"), row.names = FALSE)

# Source Data for Extended Data Fig. 7: the plot's own data (GAM predictions on the age
# grid, 0-1 scaled as plotted) plus the P-value label table. Aggregate only.
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
write_source_data(panel_d, "figS9a")
write_source_data(panel_e, "figS9b")

