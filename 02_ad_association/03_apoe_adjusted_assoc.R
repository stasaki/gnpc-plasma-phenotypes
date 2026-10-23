# ---- GNPC code release ---------------------------------------------------------
# Part of: Supplementary Table 7 (APOE/age/sex-adjusted columns).
# Makes: APOE (+ age, sex)-adjusted association of each predicted phenotype with AD
#   in GNPC plasma (observed and predicted covariates), with summary and figures.
# Inputs:  PROTEOMICS_RDS and PREDICTED_RDS (set in config.R).
# Outputs: step_dir("02_apoe_adjusted").
# Runs on the AD Workbench (the GNPC data stay there). From the repository root:
#   Rscript 02_ad_association/03_apoe_adjusted_assoc.R
# ---------------------------------------------------------------------------------
# APOE (+age,sex)-adjusted predicted-phenotype -> AD associations
# ============================================================================
# Extends the predicted-phenotype -> AD association (02_ad_assoc_site_adjusted.Rmd,
# `ad ~ prediction + site`) with TWO covariate sets:
#   Model A (site only): ad ~ scale(prediction) + site
#   Model B (OBSERVED):  ad ~ scale(prediction) + APOE4_count      + age_z      + sex_obs  + site
#   Model C (PREDICTED): ad ~ scale(prediction) + APOE4_count_pred + age_pred_z + sex_pred + site
# Tests whether the prediction's OR changes after adjustment and whether the
# prediction adds signal beyond the covariates (Nagelkerke dR^2 of full vs
# covariate-only; LR test). Two-sided. AD cases vs
# clean controls, plasma, excluded contributors removed.
#
# Model B (OBSERVED) = the literal "control for APOE genotype, age, sex" reading;
# tests for confounding by observed demographics. Model C (PREDICTED) is the sharper "is the
# phenotype model just a proteomic-demographic proxy?" test: the dominant covariates
# in a proteomic model are the proteome's OWN age/sex/APOE readouts (predicted age
# T0001, predicted sex T0002, predicted e4-dose from T0607) -- same prediction.rds
# namespace as the tested phenotypes, so it is apples-to-apples. Predicted covariates
# are collinear with the phenotype prediction BY DESIGN; a pure proxy would lose its
# signal -> we also report VIF of the prediction in Model C. Attenuation can be
# larger under C if predicted age/APOE are themselves associated with AD.
#
# Observed covariates live in colData(ProteomicsV1_3ms): APOE genotype, age_at_visit,
# sex (e4-dose / cleaning recipe as in 03_trajectories/01_age_trajectory.R).
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(tidyr)
  library(data.table); library(ggplot2); library(matrixStats)
})

source("config.R")
OUT_DIR  <- step_dir("02_apoe_adjusted")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
.logf <- file(file.path(OUT_DIR, "run_log.txt"), open = "wt"); sink(.logf, split = TRUE)

# AD-relevant phenotypes highlighted in the figure (all phenotypes are tested)
HEADLINE <- c("tangsqrt_est_8reg", "amylsqrt_est_8reg", "braaksc", "ceradsc",
              "niareagansc", "adnc", "cogn_global", "cts_mmse30", "gpath", "dcfdx")

# ---- load ------------------------------------------------------------------
expdata    <- readRDS(PROTEOMICS_RDS)
metadata   <- as.data.frame(colData(expdata))
pp         <- read_predicted_phenotypes(PREDICTED_RDS)
prediction <- get_prediction(pp, "full")
tannot     <- get_table(pp, "target_annotation")
id2name <- setNames(tannot$Target, tannot$`Target ID`)   # T-id -> name

data <- metadata %>% inner_join(prediction, by = "sample_id")
data <- data %>% filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)

cat("\n=== sample_matrix ===\n"); print(table(data$sample_matrix, useNA = "always"))
data <- data %>% filter(grepl("plasma", tolower(trimws(as.character(sample_matrix)))))   # plasma only
data$site <- factor(paste0(data$contributor_code, "_", data$sample_matrix))

# ---- observed covariates (recipe as in 03_trajectories/01_age_trajectory.R) ----
data$APOE <- as.character(data$APOE)
data$APOE[data$APOE == "32"] <- "23"; data$APOE[data$APOE == "42"] <- "24"; data$APOE[data$APOE == "43"] <- "34"
data <- data %>%
  mutate(APOE = na_if(APOE, "-1"),
         APOE4_count = recode(APOE, `22` = 0L, `23` = 0L, `24` = 1L, `33` = 0L,
                              `34` = 1L, `44` = 2L, .default = NA_integer_),
         age_obs = suppressWarnings(as.numeric(as.character(age_at_visit))),
         age_obs = ifelse(age_obs < 0, NA, age_obs),
         sex_obs = suppressWarnings(as.numeric(as.factor(na_if(as.character(sex), "-1")))))
cat("\n[APOE] e4-dose counts:\n"); print(table(data$APOE4_count, useNA = "always"))

# ---- PREDICTED covariates (Model C): the proteome's own age/sex/APOE axes ----
# Same prediction.rds namespace as the tested phenotypes. Classification base
# columns (T0002 sex, T0607 APOE) hold a predicted class LABEL (non-numeric);
# recode -> numeric sex / e4-dose. T0001 = predicted age (numeric).
AGE_PID <- "T0001"; SEX_PID <- "T0002"; APOE_PID <- "T0607"
PRED_COV_IDS <- c(AGE_PID, SEX_PID, APOE_PID)        # excluded from Model C when self-tested
apoe_dose <- function(x) { a <- as.character(x)
a[a == "32"] <- "23"; a[a == "42"] <- "24"; a[a == "43"] <- "34"
recode(a, `22` = 0L, `23` = 0L, `24` = 1L, `33` = 0L, `34` = 1L, `44` = 2L, .default = NA_integer_) }
data$age_pred <- if (AGE_PID %in% names(data)) suppressWarnings(as.numeric(as.character(data[[AGE_PID]]))) else NA_real_
data$sex_pred <- if (SEX_PID %in% names(data)) suppressWarnings(as.numeric(as.factor(as.character(data[[SEX_PID]])))) else NA_real_
data$APOE4_count_pred <- if (APOE_PID %in% names(data)) apoe_dose(data[[APOE_PID]]) else NA_integer_
cat("\n[predicted covariates] non-missing -- age:", sum(is.finite(data$age_pred)),
    "sex:", sum(is.finite(data$sex_pred)), "e4-dose:", sum(!is.na(data$APOE4_count_pred)), "\n")
have_predcov <- all(PRED_COV_IDS %in% names(data)) &&
  all(c(sum(is.finite(data$age_pred)), sum(is.finite(data$sex_pred)), sum(!is.na(data$APOE4_count_pred))) > 0)
if (!have_predcov) warning("predicted covariates unavailable; Model C will be skipped")

# ---- AD cases vs clean controls (as in 02_ad_assoc_site_adjusted.Rmd) ------
dx <- c("ad", "pd", "ftd", "als", "stroke")
DXm <- sapply(dx, function(d) { v <- suppressWarnings(as.numeric(as.character(data[[d]]))); ifelse(v %in% c(0, 1), v, NA) })
control <- rowSums(DXm == 1, na.rm = TRUE) == 0
data$y <- ifelse(DXm[, "ad"] == 1, 1L, ifelse(DXm[, "ad"] == 0 & control, 0L, NA))

md0 <- data %>%
  filter(!is.na(y), !is.na(APOE4_count), !is.na(age_obs), !is.na(sex_obs)) %>%
  mutate(age_z = scale(age_obs)[, 1],
         age_pred_z = as.numeric(scale(age_pred)))    # predicted-age covariate (Model C)
# keep sites with both cases and controls (so site is estimable)
keep <- md0 %>% group_by(site) %>% summarise(a = any(y == 1), c = any(y == 0), .groups = "drop") %>%
  filter(a & c) %>% pull(site)
md0 <- md0 %>% filter(site %in% keep) %>% mutate(site = droplevels(site))
cat(sprintf("\n=== analysis set: %d (AD=%d, control=%d) across %d sites ===\n",
            nrow(md0), sum(md0$y == 1), sum(md0$y == 0), nlevels(md0$site)))
fwrite(data.frame(n = nrow(md0), n_AD = sum(md0$y == 1), n_control = sum(md0$y == 0),
                  n_sites = nlevels(md0$site)), file.path(OUT_DIR, "analysis_set.csv"))

# ---- predictors: numeric predicted-phenotype columns (T-ids) ---------------
# Classification-target BASE columns hold a non-numeric predicted class label
# (not scalable); keep continuous predictions + class-PROBABILITY columns only.
pred_cols <- grep("^T[0-9]+", colnames(prediction), value = TRUE)
pred_cols <- intersect(pred_cols, colnames(md0))
is_num <- vapply(pred_cols, function(c) is.numeric(md0[[c]]), logical(1))
cat("dropped non-numeric (categorical class-label) columns:", sum(!is_num), "\n")
pred_cols <- pred_cols[is_num]
cat("predicted-phenotype columns tested:", length(pred_cols), "\n")

# safe T-id -> name (single-bracket = NA on miss; class cols map via their base id)
base_ids   <- sub("_class_[0-9]+$", "", pred_cols)
cls_suffix <- ifelse(grepl("_class_", pred_cols), sub(".*(_class_[0-9]+)$", "\\1", pred_cols), "")
nm         <- unname(id2name[base_ids])
pred_names <- ifelse(is.na(nm), pred_cols, paste0(nm, cls_suffix))
names(pred_names) <- pred_cols

# ---- Restrict to the Fig 4g predicted-phenotype panel ----
# = feature_annot[!group %in% c("2nd tier","exclude")]$target_id, MINUS the four
# clinical-diagnosis models that Fig 4g excludes (together with the min-class<60
# filter): dcfdx (T0046; circular — AD is *defined* from dcfdx),
# r_pd (T0648), r_stroke (T0649), r_depres (T0650; the last three rare in ROSMAP).
# Result: 108 predictors (69 continuous + 39 categorical class columns; 98 phenotypes),
# matched by target_id. Counts then reflect the disease-association panel.
FIG4G_PANEL <- c(
  "T0003","T0004","T0005","T0006","T0007","T0008","T0019","T0614","T0615","T0616","T0617","T0618",
  "T0619","T0001","T0002_class_1","T0613","T0691","T0692","T0693","T0606_class_1","T0606_class_2",
  "T0606_class_3","T0606_class_4","T0606_class_5","T0607_class_1","T0607_class_2","T0607_class_3",
  "T0607_class_4","T0607_class_5","T0611","T0049","T0047","T0048","T0653_class_1","T0654_class_1",
  "T0655_class_1","T0656_class_1","T0657_class_1","T0658_class_1","T0659_class_1","T0660_class_1",
  "T0661_class_1","T0662_class_1","T0663_class_1","T0664_class_1","T0665_class_1","T0666_class_1",
  "T0667_class_1","T0668_class_1","T0669_class_1","T0670_class_1","T0671_class_1","T0672_class_1",
  "T0673_class_1","T0675","T0009","T0010","T0011","T0012","T0013","T0014","T0015","T0016","T0017",
  "T0018","T0546_class_1","T0552","T0553","T0554","T0555","T0556","T0557","T0558","T0559","T0560",
  "T0561","T0726_class_1","T0727_class_1","T0727_class_2","T0727_class_3","T0728_class_1","T0730",
  "T0545","T0593","T0594","T0595","T0721","T0724","T0732","T0705","T0706","T0707","T0708","T0709",
  "T0710","T0711","T0712","T0725","T0713","T0714","T0715","T0716","T0717","T0718","T0719","T0720",
  "T0722","T0734_class_1")
keep_panel <- pred_cols %in% FIG4G_PANEL
pred_cols  <- pred_cols[keep_panel]
pred_names <- pred_names[pred_cols]
cat("restricted to Fig 4g predicted-phenotype panel:", length(pred_cols), "of", length(FIG4G_PANEL), "\n")

nagelkerke <- function(full, null) {
  n <- length(full$y); L0 <- as.numeric(logLik(null)); L1 <- as.numeric(logLik(full))
  (1 - exp(2 * (L0 - L1) / n)) / (1 - exp(2 * L0 / n))
}

# covariate-only models (same rows for every predictor -> fit once)
cov_only <- glm(y ~ APOE4_count + age_z + sex_obs + site, binomial, md0)
r2_cov   <- nagelkerke(cov_only, glm(y ~ 1, binomial, md0))
# predicted-covariate-only baseline (Model C); guarded if predicted covs missing
cov_only_pred <- if (have_predcov)
  tryCatch(glm(y ~ APOE4_count_pred + age_pred_z + sex_pred + site, binomial, md0), error = function(e) NULL) else NULL
r2_cov_pred   <- if (!is.null(cov_only_pred)) nagelkerke(cov_only_pred, glm(y ~ 1, binomial, model.frame(cov_only_pred))) else NA_real_

get_term <- function(fit, term = "pred_z") {
  co <- summary(fit)$coefficients
  if (!(term %in% rownames(co))) return(setNames(rep(NA_real_, 5), c("beta", "se", "p", "lo", "hi")))
  ci <- confint.default(fit)[term, ]   # Wald CI (fast; avoids profile-likelihood cost)
  setNames(c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "Pr(>|z|)"], ci[1], ci[2]),
           c("beta", "se", "p", "lo", "hi"))
}

# ---- per-predictor: Model A (site-only) vs Model B (adjusted) ---------------
res <- lapply(pred_cols, function(tid) {
  md <- md0 %>% mutate(pred_z = scale(.data[[tid]])[, 1])
  A <- tryCatch(glm(y ~ pred_z + site, binomial, md), error = function(e) NULL)
  B <- tryCatch(glm(y ~ pred_z + APOE4_count + age_z + sex_obs + site, binomial, md), error = function(e) NULL)
  if (is.null(A) || is.null(B)) return(NULL)
  a <- get_term(A); b <- get_term(B)
  lr <- anova(cov_only, B, test = "LRT")
  
  # ---- Model C: PREDICTED covariates (skip when the predictor IS a covariate) ----
  is_cov <- tid %in% PRED_COV_IDS
  C <- if (have_predcov && !is_cov)
    tryCatch(glm(y ~ pred_z + APOE4_count_pred + age_pred_z + sex_pred + site, binomial, md), error = function(e) NULL) else NULL
  cc_  <- if (is.null(C)) setNames(rep(NA_real_, 5), c("beta","se","p","lo","hi")) else get_term(C)
  vifC <- if (is.null(C)) NA_real_ else
    tryCatch(1 / (1 - summary(lm(pred_z ~ APOE4_count_pred + age_pred_z + sex_pred + site, md))$r.squared), error = function(e) NA_real_)
  lrC  <- if (is.null(C) || is.null(cov_only_pred)) NA_real_ else
    tryCatch(anova(cov_only_pred, C, test = "LRT")$`Pr(>Chi)`[2], error = function(e) NA_real_)
  dR2C <- if (is.null(C)) NA_real_ else nagelkerke(C, glm(y ~ 1, binomial, md)) - r2_cov_pred
  
  data.frame(
    target_id = tid, name = pred_names[[tid]],
    beta_unadj = a["beta"], or_unadj = exp(a["beta"]),
    or_unadj_lo = exp(a["lo"]), or_unadj_hi = exp(a["hi"]), p_unadj = a["p"],
    # Model B (observed covariates)
    beta_adj = b["beta"], or_adj = exp(b["beta"]),
    or_adj_lo = exp(b["lo"]), or_adj_hi = exp(b["hi"]), p_adj = b["p"],
    attenuation = b["beta"] / a["beta"],                      # 1 = no attenuation by adjustment
    delta_nagelkerke_R2 = nagelkerke(B, glm(y ~ 1, binomial, md)) - r2_cov,  # variance added beyond OBS covs
    lr_p_addpred = lr$`Pr(>Chi)`[2],                          # LR test: does the prediction add to the OBS covariate-only model?
    # Model C (predicted covariates)
    beta_adjpred = cc_["beta"], or_adjpred = exp(cc_["beta"]),
    or_adjpred_lo = exp(cc_["lo"]), or_adjpred_hi = exp(cc_["hi"]), p_adjpred = cc_["p"],
    attenuation_pred = cc_["beta"] / a["beta"],
    delta_nagelkerke_pred = dR2C,                             # variance added beyond PRED covs
    lr_p_addpred_vs_predcov = lrC,                            # LR test: does the prediction add to the PRED covariate-only model?
    vif_pred = vifC,                                          # collinearity of pred with predicted covs
    row.names = NULL)
}) %>% bind_rows()

res <- res %>% mutate(
  p_adj_fdr = p.adjust(p_adj, "BH"),
  p_adjpred_fdr = p.adjust(p_adjpred, "BH"),
  headline = name %in% HEADLINE) %>%
  arrange(p_adj)
fwrite(res, file.path(OUT_DIR, "apoe_adjusted_assoc.csv"))
cat("\n=== highlighted phenotypes: unadjusted vs OBSERVED-adj (B) vs PREDICTED-adj (C) ===\n")
print(res %>% filter(headline) %>%
        select(name, or_unadj, p_unadj, or_adj, p_adj, attenuation,
               or_adjpred, p_adjpred, attenuation_pred, vif_pred), digits = 3)

# concordance of effect sizes with Model A -- B (observed) and C (predicted)
cc   <- with(res, cor(beta_unadj, beta_adj,     use = "complete.obs"))
ccP  <- with(res, cor(beta_unadj, beta_adjpred, use = "complete.obs"))
# RETENTION: of the predictors significant in Model A (Fig 4g,
# unadjusted), how many remain significant after adjustment? (not "sig in both" / N).
sigA      <- res$p_unadj < 0.05
n_sigA    <- sum(sigA, na.rm = TRUE)
retain_B  <- sum(sigA & res$p_adj     < 0.05, na.rm = TRUE)   # sig in A AND still sig in B
retain_C  <- sum(sigA & res$p_adjpred < 0.05, na.rm = TRUE)   # sig in A AND still sig in C
sig_both  <- with(res, mean(p_unadj < 0.05 & p_adj     < 0.05, na.rm = TRUE))
sig_bothP <- with(res, mean(p_unadj < 0.05 & p_adjpred < 0.05, na.rm = TRUE))
cat(sprintf("\n[A unadj/Fig4g] significant: %d / %d predictors\n", n_sigA, nrow(res)))
cat(sprintf("[B observed ] retained sig: %d / %d (%.0f%%) | beta concordance r = %.3f | median attenuation = %.3f\n",
            retain_B, n_sigA, 100 * retain_B / n_sigA, cc, median(res$attenuation, na.rm = TRUE)))
cat(sprintf("[C predicted] retained sig: %d / %d (%.0f%%) | beta concordance r = %.3f | median attenuation = %.3f | median VIF = %.2f\n",
            retain_C, n_sigA, 100 * retain_C / n_sigA, ccP, median(res$attenuation_pred, na.rm = TRUE),
            median(res$vif_pred, na.rm = TRUE)))
fwrite(data.frame(
  n_panel = nrow(res), n_sig_A = n_sigA,
  retained_B = retain_B, pct_retained_B = 100 * retain_B / n_sigA,
  beta_concordance = cc, median_attenuation = median(res$attenuation, na.rm = TRUE),
  cov_only_nagelkerke_R2 = r2_cov,
  retained_C = retain_C, pct_retained_C = 100 * retain_C / n_sigA,
  beta_concordance_pred = ccP, median_attenuation_pred = median(res$attenuation_pred, na.rm = TRUE),
  predcov_only_nagelkerke_R2 = r2_cov_pred, median_vif_pred = median(res$vif_pred, na.rm = TRUE),
  pct_sig_both = 100 * sig_both, pct_sig_both_pred = 100 * sig_bothP),
  file.path(OUT_DIR, "summary.csv"))

# ---- figures ---------------------------------------------------------------
lim <- range(c(res$beta_unadj, res$beta_adj), na.rm = TRUE)
p_sc <- ggplot(res, aes(beta_unadj, beta_adj)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = headline), alpha = 0.6, size = 1) +
  geom_point(data = subset(res, headline), colour = "#D55E00", size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D55E00"), guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "log-OR, site-only (Model A)", y = "log-OR, + APOE + age + sex",
       title = sprintf("Predicted-phenotype -> AD associations: site-only vs + APOE/age/sex (r = %.2f)", cc)) +
  theme_bw(base_size = 9)
ggsave(file.path(OUT_DIR, "Fig_apoe_adjusted_scatter.pdf"), p_sc, width = 6, height = 6)

# Model C scatter: unadjusted vs PREDICTED-covariate-adjusted log-OR (if available)
if (any(is.finite(res$beta_adjpred))) {
  limP <- range(c(res$beta_unadj, res$beta_adjpred), na.rm = TRUE)
  p_scP <- ggplot(res, aes(beta_unadj, beta_adjpred)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(colour = "grey60", alpha = 0.6, size = 1) +
    geom_point(data = subset(res, headline), colour = "#0072B2", size = 1.8) +
    coord_equal(xlim = limP, ylim = limP) +
    labs(x = "log-OR, site-only (Model A)", y = "log-OR, + PREDICTED age/sex/APOE",
         title = sprintf("AD associations adjusted for the proteomic demographic axes (r = %.2f)", ccP)) +
    theme_bw(base_size = 9)
  ggsave(file.path(OUT_DIR, "Fig_apoe_adjusted_scatter_predcov.pdf"), p_scP, width = 6, height = 6)
}

# forest plot of the highlighted phenotypes: site-only vs +OBSERVED vs +PREDICTED covariates (tidy from CSV)
hdl <- bind_rows(
  res %>% filter(headline) %>% transmute(name, model = "site-only",        or = or_unadj,   lo = or_unadj_lo,   hi = or_unadj_hi),
  res %>% filter(headline) %>% transmute(name, model = "+observed covs",   or = or_adj,     lo = or_adj_lo,     hi = or_adj_hi),
  res %>% filter(headline) %>% transmute(name, model = "+predicted covs",  or = or_adjpred, lo = or_adjpred_lo, hi = or_adjpred_hi)) %>%
  filter(is.finite(or))
p_for <- ggplot(hdl, aes(or, reorder(name, or), colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
  labs(x = "Odds ratio for AD (per SD of predicted phenotype)", y = NULL, colour = NULL,
       title = "Highlighted predicted phenotypes: AD association vs observed & predicted covariates") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "Fig_apoe_adjusted_forest.pdf"), p_for, width = 7, height = 5)

cat("\n=== exported files ===\n"); print(list.files(OUT_DIR))
cat("\nDONE.\n")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
sink(); close(.logf)
