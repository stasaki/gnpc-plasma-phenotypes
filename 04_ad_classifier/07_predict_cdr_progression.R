# ==============================================================================
# Extended Data Fig. 10 and Supplementary Table 12: tests whether the baseline
# per-domain predicted AD risk (domain models of the AD classifier) is associated
# with the change in CDR over follow-up. Linear mixed models in GNPC participants
# with repeat visits.
# Inputs : ProteomicsV1_3ms.rds; classifier runs p8_one_contributoer_out_nesting_<EX_TAG>
#          in step_dir("04_classifier_runs") (01_train_classifiers.R).
# Outputs: step_dir("04_predict_future")
# Runs on the AD Workbench. From the repo root:
#   Rscript 04_ad_classifier/07_predict_cdr_progression.R
# ==============================================================================

source("config.R")

library(tidyverse)
library(tidyverse)
library(dplyr)
library(data.table)
library(pROC)
library(ggplot2)
library(knitr)
library(PRROC)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)  # For colorblind-friendly palettes
library(hrbrthemes)  # For nice typography
library(ggtext)  # For rich text formatting
library(forcats)   # For factor reordering
library(r2glmm)
library(broom.mixed)
library(patchwork)
library(ggplot2)
library(dplyr)

# metadata ####
ProteomicsV1_3ms <- readRDS(PROTEOMICS_RDS)
metadata = colData(ProteomicsV1_3ms)%>%as.data.frame()

keep=metadata$IID%>%.[duplicated(.)]
metadata = metadata%>%
  filter(IID%in%keep)
metadata%>%
  group_by(IID)%>%
  summarise(N=dplyr::n())%>%pull(N)%>%table

# domain specific prediction values ####
list.files(file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/")))%>%
  .[grepl("0.01",.)]%>%
  .[!grepl(".all_1st_tier|\\.Clinical.Diagnosis",.)]%>%
  lapply(., function(x){
    base_path <- paste0(file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/")),x,"/ml/p1_trait")
    submodels_path <- file.path(base_path, "submodels")
    # 1. Load sample IDs
    sample_id <- readRDS(file.path(base_path, "sample_id.rds"))
    prediction_files <- list.files(
      path = submodels_path,
      pattern = "predictions.rds",
      recursive = TRUE,
      full.names = TRUE
    )
    
    predictions <- prediction_files %>%
      map(readRDS) %>%
      bind_rows()%>%
      filter(split=="test")
    predictions$sample_id <- sample_id[as.numeric(predictions$sample_id)]
    predictions$model =gsub(".+feature..","",x)%>%gsub("__.+","",.)
    predictions%>%return()
  })%>%
  bind_rows() -> predictions_all
metadata = metadata%>%
  left_join(.,predictions_all%>%
              dplyr::select(sample_id,predicted_prob,model)%>%
              spread(model,predicted_prob),by=c("SAMPLE.ID"="sample_id"))

metadata$sex = metadata$sex -1 

metadata = metadata%>%
  filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS,
         sample_matrix%in%c("Citrate Plasma","EDTA Plasma"))%>%
  filter(age_at_visit!= -1,
         cdr != -1)

trait = "cdr"
primary_predictors = c("Cognition","Demographics",
                       "Genetics","Medical.Conditions",
                       "Medications","Motor.and.Gait","Pathology")
covs = c("sex","contributor_code","age_at_visit")


ds = metadata%>%.[,c("IID","SAMPLE.ID","ad",trait,covs,primary_predictors)]%>%
  drop_na()

keep=ds$IID%>%.[duplicated(.)]
ds = ds%>%
  filter(IID%in%keep)

ds = ds%>%
  group_by(IID)%>%
  filter(sum(duplicated(age_at_visit))==0)%>%
  ungroup()

ds%>%
  group_by(IID)%>%
  summarise(N=dplyr::n())%>%pull(N)%>%table


# lmer model ####
library(lme4)
library(lmerTest)
library(MuMIn)

# First, create a baseline time variable (time since first visit)
ds <- ds %>%
  group_by(IID) %>%
  mutate(
    time = age_at_visit - min(age_at_visit),  # years from baseline
  ) %>%
  ungroup()

# Extract baseline values for each predictor
baseline_vars <- ds %>%
  filter(time == 0) %>%
  dplyr::select(IID, cdr_bl = cdr, Cognition_bl = Cognition, Demographics_bl = Demographics,
                Genetics_bl = Genetics, Medical.Conditions_bl = Medical.Conditions,
                Medications_bl = Medications, Motor.and.Gait_bl = Motor.and.Gait,
                Pathology_bl = Pathology, baseline_age=age_at_visit)

# Merge back
df <- inner_join(ds, baseline_vars, by = "IID")


library(tableone)
metadata%>%
  filter(SAMPLE.ID%in%df$SAMPLE.ID)%>%
  group_by(IID) %>%
  mutate(
    time = age_at_visit - min(age_at_visit),  # years from baseline
    n_followup = dplyr::n(),
    year_followup = mean(time)
  ) %>%
  ungroup()%>%
  mutate(baseline = time==0)%>%
  mutate(years_of_education=replace(years_of_education,years_of_education==-1,NA))%>%
  mutate(bmi=replace(bmi,bmi==-1,NA))%>%
  mutate(cdr=replace(cdr,cdr==-1,NA))%>%
  mutate(cognitive_test_score=replace(cognitive_test_score,cognitive_test_battery!="MMSE",NA),
         cognitive_test_score=replace(cognitive_test_score,cognitive_test_score>30,NA),
         MMSE=cognitive_test_score)%>%
  mutate(race=replace(race,race==-1,NA)) -> meta_demo


race_mapping = c("1"="American_Indian/Alaska_Native",
                 "2"="Black_or_African_American",
                 "3"="Native_Hawaiian_or_Other_Pacific_Islander",
                 "4"="Asian",
                 "5"="White",
                 "66"="Other(Group_C)",
                 "77"="Non_Caucasian_or_Black(Group_B)",
                 "88"="Non_Caucasian(Group_A)",
                 "99"="Unknown",
                 "-1"="Not_Stated")
meta_demo$race = race_mapping[meta_demo$race]
meta_demo$ad = meta_demo$ad==1
colnames(meta_demo)[colnames(meta_demo)=="cdr"] = "CDR"

CreateTableOne(data = meta_demo%>%filter(baseline),     
               vars = c("sex","ad","age_at_visit","n_followup","year_followup",
                        "contributor_code","race","APOE","years_of_education","CDR","MMSE"),
               factorVars = c("ad","sex","contributor_code","race","APOE"))-> tb1

# Get the printed table as a matrix
mat <- print(tb1, missing = TRUE, printToggle = FALSE, showAllLevels = TRUE)

# Get missing counts
missing_n <- sapply(meta_demo%>%filter(baseline)%>%.[,c("sex","ad","age_at_visit","n_followup","year_followup",
                                                        "contributor_code","race","APOE","years_of_education","CDR","MMSE")],
                    function(x) sum(is.na(x)))

# Loop through your variables and update matching rows
for (var in names(missing_n)) {
  # Find rows that start with this variable name
  idx <- grep(paste0("^", var), rownames(mat))
  if (length(idx) > 0) {
    # Only update the first match (the main variable row)
    pct <- mat[idx[1], "Missing"]
    mat[idx[1], "Missing"] <- paste0(missing_n[var], " (", pct, ")")
  }
}

write.table(print(mat, quote = FALSE), 
            file = file.path(step_dir("04_predict_future"), "table_one_baseline.txt"), 
            sep = "\t", 
            row.names = TRUE, 
            col.names = TRUE, 
            quote = FALSE)


CreateTableOne(data = meta_demo%>%filter(!baseline),     
               vars = c("sex","ad","age_at_visit","n_followup","year_followup",
                        "contributor_code","race","APOE","years_of_education","CDR","MMSE"),
               factorVars = c("ad","sex","contributor_code","race","APOE"))-> tb1

# Get the printed table as a matrix
mat <- print(tb1, missing = TRUE, printToggle = FALSE, showAllLevels = TRUE)

# Get missing counts
missing_n <- sapply(meta_demo%>%filter(!baseline)%>%.[,c("sex","ad","age_at_visit","n_followup","year_followup",
                                                         "contributor_code","race","APOE","years_of_education","CDR","MMSE")],
                    function(x) sum(is.na(x)))

# Loop through your variables and update matching rows
for (var in names(missing_n)) {
  # Find rows that start with this variable name
  idx <- grep(paste0("^", var), rownames(mat))
  if (length(idx) > 0) {
    # Only update the first match (the main variable row)
    pct <- mat[idx[1], "Missing"]
    mat[idx[1], "Missing"] <- paste0(missing_n[var], " (", pct, ")")
  }
}


# Save as tab-delimited table
write.table(print(mat, quote = FALSE), 
            file = file.path(step_dir("04_predict_future"), "table_one_fu.txt"), 
            sep = "\t", 
            row.names = TRUE, 
            col.names = TRUE, 
            quote = FALSE)

all_terms = c("cdr_bl","baseline_age", "sex", "contributor_code", "Cognition_bl","Demographics_bl",
              "Genetics_bl", "Medical.Conditions_bl","Medications_bl", "Motor.and.Gait_bl", "Pathology_bl")

all_terms%>%
  lapply(., function(term_ex){
    # One model per baseline term: the term and its interaction with time,
    # plus age, sex and contributor (and their interactions with time)
    model <- lmer(
      paste0("cdr_change ~ time +baseline_age+contributor_code+sex+",
             paste0(all_terms%>%intersect(term_ex),collapse = "+"),"+",
             paste0("time:(baseline_age+contributor_code+sex+",paste0(all_terms%>%intersect(term_ex),collapse = "+"),")"),"+(0 + time| IID)")%>%as.formula(),
      data = df %>%
        group_by(IID) %>%
        arrange(time)%>%
        mutate(cdr_change = cdr - dplyr::first(cdr)) %>%  # change from baseline
        ungroup()%>%
        mutate(
          baseline_age = scale(baseline_age)[,1],
          cdr_bl = scale(cdr_bl)[,1],
          Cognition_bl = scale(Cognition_bl)[,1],
          Demographics_bl = scale(Demographics_bl)[,1],
          Genetics_bl = scale(Genetics_bl)[,1],
          Medical.Conditions_bl = scale(Medical.Conditions_bl)[,1],
          Medications_bl = scale(Medications_bl)[,1],
          Motor.and.Gait_bl = scale(Motor.and.Gait_bl)[,1],
          Pathology_bl = scale(Pathology_bl)[,1]
        ),
      REML = FALSE
    )
    # sample sizes per univariate model
    n_obs_uni <- nrow(model@frame)
    n_groups_uni <- length(unique(df$IID))
    broom::tidy(model,conf.int = TRUE)%>%
      as.data.frame()%>%
      mutate(ex_term = term_ex,
             n_obs = n_obs_uni,
             n_groups = n_groups_uni)%>%
      return()
  })%>%
  bind_rows() -> uni_model_res
uni_model_res%>%
  fwrite(.,file=file.path(step_dir("04_predict_future"), "lmer_uni_var_table.txt"),quote = F,sep = "\t")

# Extract coefficients
coef_df <- uni_model_res%>%
  filter(effect == "fixed")%>%
  filter(grepl("time:.+_bl",term))

# Time interactions (effect on rate of change)
time_int_df <- coef_df %>%
  mutate(
    predictor = gsub("time:|_bl_z|_z", "", term),
    predictor = gsub("_|\\.", " ", predictor),
    direction = ifelse(estimate > 0, "Faster decline", "Slower decline"),
    significant = p.value < 0.05
  ) %>%
  filter(!grepl("batch|Batch|contributor", predictor, ignore.case = TRUE)) %>%
  arrange(estimate)%>%
  mutate(predictor=gsub("cdr bl","CDR",predictor),
         predictor=gsub("sex","Sex",predictor),
         predictor=gsub("baseline age","Age",predictor),
         predictor=gsub(" bl"," model",predictor))

# Fit mixed model with baseline predictors affecting slope
model <- lmer(
  cdr_change ~ time + baseline_age + sex + contributor_code +cdr_bl+
    Cognition_bl + 
    Medical.Conditions_bl + Medications_bl + Motor.and.Gait_bl + Pathology_bl +
    time:(cdr_bl+baseline_age + sex + contributor_code +Cognition_bl + Medical.Conditions_bl + 
            Medications_bl + Motor.and.Gait_bl + Pathology_bl) +
    (0 + time| IID),
  data = df %>%
    group_by(IID) %>%
    arrange(time)%>%
    mutate(cdr_change = cdr - dplyr::first(cdr)) %>%  # change from baseline
    ungroup()%>%
    mutate(
      baseline_age = scale(baseline_age)[,1],
      cdr_bl = scale(cdr_bl)[,1],
      Cognition_bl = scale(Cognition_bl)[,1],
      Demographics_bl = scale(Demographics_bl)[,1],
      Genetics_bl = scale(Genetics_bl)[,1],
      Medical.Conditions_bl = scale(Medical.Conditions_bl)[,1],
      Medications_bl = scale(Medications_bl)[,1],
      Motor.and.Gait_bl = scale(Motor.and.Gait_bl)[,1],
      Pathology_bl = scale(Pathology_bl)[,1]
    ),
  REML = FALSE
)
# CIs, R², sample sizes for multivariate model
broom::tidy(model, conf.int = TRUE)%>%
  fwrite(.,file=file.path(step_dir("04_predict_future"), "lmer_multi_var_table.txt"),quote = F,sep = "\t")
summary(model)
r2_multi <- r.squaredGLMM(model)
n_obs_multi <- nrow(model@frame)
n_groups_multi <- length(unique(df$IID))
tibble(R2m = r2_multi[1, "R2m"], R2c = r2_multi[1, "R2c"],
       n_obs = n_obs_multi, n_groups = n_groups_multi) %>%
  fwrite(., file = file.path(step_dir("04_predict_future"), "lmer_multi_var_summary.txt"), quote = FALSE, sep = "\t")




time_int_df = time_int_df%>%
  left_join(.,r2beta(model, method = "nsj", partial = TRUE),by=c("term"="Effect"))
fwrite(time_int_df,file=file.path(step_dir("04_predict_future"), "time_int_df.txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)


# Forest plot
p1 <- ggplot(time_int_df, aes(x = estimate, y = reorder(predictor, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, color = "gray40", orientation = "y") +
  geom_point(aes(color = direction, shape = significant), size = 3) +
  scale_color_manual(values = c("Faster decline" = "firebrick", "Slower decline" = "steelblue")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),guide="none") +  # filled vs open
  labs(
    x = "Effect on CDR Change Rate\n(per 1 SD, per year)",
    y = NULL,
    color = NULL,
    shape = "p < 0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")+
  theme(axis.text.y = element_text(face = "bold"))

p1

# R² bar - use full data, geom_col handles NA automatically
time_int_df$Rsq=time_int_df$Rsq*100
p_r2 <- ggplot(time_int_df, aes(y = reorder(predictor, estimate), x = Rsq)) +
  geom_col(fill = "steelblue", width = 0.5, na.rm = TRUE) +
  geom_text(aes(label = ifelse(is.na(Rsq), "", sprintf("%.1f%%", Rsq))),
            hjust = -0.1, size = 3, na.rm = TRUE)+
  xlim(0, max(time_int_df$Rsq, na.rm = TRUE) * 1.4) +
  labs(x = "Unique R²", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 0))

wrap_plots(p1, p_r2, widths = c(3, 1))

# Get the top 2 significant predictors
top_preds <- time_int_df %>%
  filter(predictor!="CDR")%>%
  filter(significant) %>%
  arrange(p.value) %>%
  slice_head(n = 2) %>%
  pull(predictor)

# Map back to variable names
top_preds_z <- paste0(gsub(" ", ".", top_preds%>%gsub(" model","",.)), "_bl")

library(ggeffects)
# Generate predictions for each
pred_list <- lapply(top_preds_z, function(pred) {
  ggpredict(model, terms = c("time [0:10 by=1]", paste0(pred, " [-1, 0, 1]"))) %>%
    as.data.frame()%>%
    mutate(predictor = gsub("_bl", " model", pred),
           predictor = gsub("\\.", " ", predictor))
})

pred_all <- bind_rows(pred_list)


# Get p-values for top predictors
pval_df <- time_int_df %>%
  mutate(
    pval_label = case_when(
      p.value < 0.001 ~ sprintf("p = %.1e", p.value),  # 0.000002 → p = 2.0e-06
      TRUE ~ sprintf("p = %.3f", p.value)               # 0.001 → p = 0.001
    )
  ) %>%
  dplyr::select(predictor, pval_label, estimate)

# Merge with predictions
pred_all <- pred_all %>%
  left_join(pval_df, by = "predictor")
# Extended Data Fig. 10c source data: model-predicted CDR change at -1/0/+1 SD (aggregate)
fwrite(pred_all, file = file.path(step_dir("04_predict_future"), "figS12c_pred_all.txt"),
       quote = F, sep = "\t")

# Plot with p-values
p4 <- ggplot(pred_all%>%mutate(predictor=factor(predictor,levels=pred_all %>% distinct(predictor, pval_label)%>%pull(predictor))),
             aes(x = x, y = predicted, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  facet_wrap(~predictor, scales = "free_y",ncol=1) +
  scale_color_manual(values = c("-1" = "steelblue", "0" = "gray50", "1" = "firebrick"),
                     labels = c("-1 SD", "Mean", "+1 SD")) +
  scale_fill_manual(values = c("-1" = "steelblue", "0" = "gray50", "1" = "firebrick"),
                    labels = c("-1 SD", "Mean", "+1 SD")) +
  labs(
    x = "Time (years)",
    y = "Predicted CDR Change",
    color = "Baseline risk",
    fill = "Baseline risk"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")+
  theme(strip.text = element_text(face = "bold"))

p_final <- wrap_plots(p1, p_r2, p4, widths = c(3, 1, 2.6))
ggsave(filename = paste0(file.path(step_dir("04_predict_future"), "lmer.png")),plot = p_final,width = 13*0.8,height = 6.5*0.8,bg = "white")
ggsave(filename = paste0(file.path(step_dir("04_predict_future"), "lmer.pdf")),plot = p_final,width = 13*0.8,height = 6.5*0.8,bg = "white")
