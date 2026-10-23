# ==============================================================================
# Extended Data Fig. 9b (group values): AD-specificity of the predicted
# phenotypes in pure AD, PD, FTD and ALS groups (site-adjusted contrasts), and
# accuracy of predicted age and sex in each group (positive control).
# Inputs : ProteomicsV1_3ms.rds; predicted_phenotypes.rds (prediction "full",
#          target_annotation).
# Outputs: step_dir("05_pdftd_specificity")
# Runs on the AD Workbench. From the repo root:
#   Rscript 05_ad_specificity/01_pdftd_specificity.R
# ==============================================================================

# AD-specificity of the predicted phenotype layer vs PD / FTD / ALS
# ============================================================================
# Also fits a predicted-phenotype AD classifier and tests its specificity
# (how it scores PD/FTD/ALS compared with AD and controls) — see the "AD-classifier SPECIFICITY" block below.
#
# Question: are the AD-related predicted phenotypes lower in those diseases than
# in AD? Are the age and sex predictions equally accurate in every group?
#
# Compares pure diagnosis groups on:
#   - AD-relevant predictions (pathology: tangles/amyloid/Braak/CERAD/ADNC/predicted
#     AD dx; and AD-type cognition): AD vs PD/FTD/ALS
#   - DEMOGRAPHIC signatures (age, sex): prediction accuracy within each disease group
#     (a positive control; their levels differ between groups with the groups' demographics)
#
# Design: pure neurodegenerative groups (exactly one of AD/PD/FTD/ALS; stroke
# does not affect group membership), plasma only,
# excluded contributors (config.R) removed. Disease contrasts are estimated WITHIN shared sites
# (sites containing both groups) so group != cohort/batch. Two-sided tests.
# NB: ALS cases come from few sites, so few sites are shared with AD. The
# AD-vs-ALS contrast therefore rests on few sites (see ad_shared_sites.csv).
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(tidyr)
  library(data.table); library(ggplot2); library(matrixStats)
})
source("config.R")

# ---- paths -----------------------------------------------------------------
OUT_DIR   <- step_dir("05_pdftd_specificity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Capture ALL console diagnostics to disk (run_log.txt).
.logf <- file(file.path(OUT_DIR, "run_log.txt"), open = "wt"); sink(.logf, split = TRUE)

# ---- phenotype panel (resolved to workbench target_ids by NAME below) ------
# ad_dir: +1 if higher value = more AD pathology; -1 if reverse-coded (CERAD,
# NIA-Reagan: low = more AD). Used to orient effect sizes so positive always
# means "AD has more AD pathology". Demographic axis is not oriented (NA).
PANEL <- tribble(
  ~name,                ~axis,           ~label,                ~ad_dir,
  "tangsqrt_est_8reg",  "AD pathology",  "Tangles (8-region)",   1,
  "amylsqrt_est_8reg",  "AD pathology",  "Amyloid (8-region)",   1,
  "braaksc",            "AD pathology",  "Braak stage",          1,
  "ceradsc",            "AD pathology",  "CERAD (low=more AD)",  -1,
  "niareagansc",        "AD pathology",  "NIA-Reagan (low=more AD)", -1,
  "adnc",               "AD pathology",  "ADNC",                 1,
  "cogn_global",        "Cognition",     "Global cognition",    -1,
  "age_at_visit",       "Demographic",   "Predicted age",        NA,
  "msex",               "Demographic",   "Predicted sex",        NA
)
DX_AD_PREDICTION_NAME <- "dcfdx"   # multiclass clinical dx; AD-class prob added if resolvable

# observed demographics in GNPC metadata (for the accuracy-based control)
OBS_AGE_COL <- "age_at_visit"   # observed age column in colData
OBS_SEX_COL <- "sex"            # observed sex column in colData

# ---- load ------------------------------------------------------------------
pp <- read_predicted_phenotypes(PREDICTED_RDS)
expdata  <- readRDS(PROTEOMICS_RDS)
metadata <- as.data.frame(colData(expdata))
prediction <- get_prediction(pp, "full")
tannot   <- get_table(pp, "target_annotation")   # columns include: Target (name), 'Target ID', Class_label

# name -> workbench target_id
name2id <- setNames(tannot$`Target ID`, tannot$Target)

resolve <- function(nm) {
  tid <- name2id[[nm]]
  if (is.null(tid) || is.na(tid) || !(tid %in% colnames(prediction))) return(NA_character_)
  tid
}
# USE.NAMES=FALSE: keep tids UNNAMED, else select(all_of()) renames the data
# columns to the phenotype names and downstream `.data[[tid]]` lookups fail.
PANEL$tid <- vapply(PANEL$name, resolve, character(1), USE.NAMES = FALSE)

# predicted AD-diagnosis probability from the multiclass dcfdx predictor (best effort)
ad_prob_col <- NA_character_
dcfdx_tid <- name2id[[DX_AD_PREDICTION_NAME]]
if (!is.na(dcfdx_tid)) {
  cls <- grep(paste0("^", dcfdx_tid, "_class_"), colnames(prediction), value = TRUE)
  # Class_label for dcfdx encodes class order; pick the AD class by string match.
  lab <- tannot$Class_label[match(dcfdx_tid, tannot$`Target ID`)]
  cat("\n[dcfdx] classes found:", paste(cls, collapse = ", "), "\n[dcfdx] Class_label:", lab, "\n")
  # heuristic: if labels list AD first/contains 'AD', map by index of 'AD'
  if (length(cls) >= 1 && !is.na(lab)) {
    parts <- trimws(strsplit(lab, "[,;|]")[[1]])
    adi <- which(grepl("AD|Alzheimer", parts, ignore.case = TRUE))
    if (length(adi) == 1 && adi <= length(cls)) ad_prob_col <- cls[adi]
  }
  cat("[dcfdx] -> AD-prob column chosen:", ifelse(is.na(ad_prob_col), "NONE (resolve manually)", ad_prob_col), "\n")
}

# ---- build analysis frame --------------------------------------------------
keep_cols <- na.omit(PANEL$tid)
pred_sel <- prediction %>% select(sample_id, all_of(keep_cols), any_of(na.omit(ad_prob_col)))

data <- metadata %>% inner_join(pred_sel, by = "sample_id")

# drop the excluded contributors (config.R)
data <- data %>% filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)

# PLASMA ONLY. The models are plasma-trained; CSF/serum samples are out-of-
# distribution and must be removed. Filter on sample_matrix directly (the raw
# contributor_code is just "Q"/"U"/..., with no "_CSF"/"_Serum" suffix).
cat("\n=== sample_matrix values (pre-filter) ===\n")
print(table(data$sample_matrix, useNA = "always"))
data <- data %>% filter(grepl("plasma", tolower(trimws(as.character(sample_matrix)))))
cat("plasma-only rows:", nrow(data), "\n")
data$site <- paste0(data$contributor_code, "_", data$sample_matrix)

# diagnosis matrix: only 0/1 are valid; -1 (unknown) and any other code (e.g. 2) -> NA
dx_all <- c("ad", "pd", "ftd", "als", "stroke")
nd     <- c("ad", "pd", "ftd", "als")               # neurodegenerative comparators
DX <- data %>% select(all_of(dx_all)) %>%
  mutate(across(everything(), ~ { v <- suppressWarnings(as.numeric(as.character(.x)))
  ifelse(v %in% c(0, 1), v, NA) }))
npos_all <- rowSums(DX == 1, na.rm = TRUE)            # across all 5 (for control)
npos_nd  <- rowSums(DX[nd] == 1, na.rm = TRUE)        # across neurodegenerative only

# PURE groups: exactly one neurodegenerative dx; STROKE does not gate purity
# (vascular, not a comparator). control = none of the 5 positive.
group <- rep(NA_character_, nrow(data))
group[npos_all == 0] <- "control"
for (d in nd) group[npos_nd == 1 & DX[[d]] %in% 1] <- toupper(d)
data$group <- factor(group, levels = c("control", "AD", "PD", "FTD", "ALS"))
# Multi-neurodegenerative (npos_nd>1) and stroke-only individuals are dropped.

ana <- data %>% filter(!is.na(group))
cat("\n=== group sizes (pure, plasma) ===\n"); print(table(ana$group))
saveRDS(table(ana$group), file.path(OUT_DIR, "group_sizes.rds"))
fwrite(as.data.frame(table(group = ana$group)), file.path(OUT_DIR, "group_sizes.csv"))

# resolved panel (phenotype name -> target_id), written for checking
panel_full <- bind_rows(
  PANEL %>% transmute(name, axis, label, ad_dir, tid),
  if (!is.na(ad_prob_col)) tibble(name = "dcfdx_AD_prob", axis = "AD pathology",
                                  label = "Predicted AD diagnosis", ad_dir = 1, tid = ad_prob_col) else NULL
) %>% filter(!is.na(tid))
fwrite(panel_full, file.path(OUT_DIR, "resolved_panel.csv"))
cat("\n=== resolved panel ===\n"); print(panel_full)

# ---- site structure: a disease comparison is only identifiable within sites --
# that contain BOTH groups (otherwise group = cohort/batch). Diagnose + restrict.
sites_per_group <- ana %>% group_by(group) %>%
  summarise(n = dplyr::n(), n_sites = n_distinct(site), .groups = "drop")
cat("\n=== sites per group (plasma) ===\n"); print(sites_per_group)
fwrite(sites_per_group, file.path(OUT_DIR, "sites_per_group.csv"))

cat("\n=== AD-shared sites per comparator (the identifiable contrast) ===\n")
shared_sites <- lapply(c("PD", "FTD", "ALS"), function(cmp) {
  s_ad <- unique(ana$site[ana$group == "AD"]); s_c <- unique(ana$site[ana$group == cmp])
  shared <- intersect(s_ad, s_c)
  cat(sprintf("  AD vs %-4s: AD sites=%d, %s sites=%d, SHARED=%d, n in shared (AD+%s)=%d\n",
              cmp, length(s_ad), cmp, length(s_c), length(shared), cmp,
              sum(ana$site %in% shared & ana$group %in% c("AD", cmp))))
  data.frame(comparator = cmp, ad_sites = length(s_ad), cmp_sites = length(s_c),
             shared_sites = length(shared),
             n_in_shared = sum(ana$site %in% shared & ana$group %in% c("AD", cmp)))
}) %>% bind_rows()
fwrite(shared_sites, file.path(OUT_DIR, "ad_shared_sites.csv"))

# ---- group comparisons: predicted phenotype ~ group + site (WITHIN shared sites) --
comparators <- c("PD", "FTD", "ALS")
stats <- list()
for (i in seq_len(nrow(panel_full))) {
  tid <- panel_full$tid[i]
  for (cmp in comparators) {
    md <- ana %>%
      filter(group %in% c("AD", cmp)) %>%
      transmute(y = .data[[tid]], group = factor(group, levels = c(cmp, "AD")), site) %>%
      drop_na()
    # keep only sites containing BOTH groups -> within-site (cohort-free) contrast
    sh <- md %>% group_by(site) %>%
      summarise(a = any(group == "AD"), c = any(group == cmp), .groups = "drop") %>%
      filter(a & c) %>% pull(site)
    md <- md %>% filter(site %in% sh)
    n_sites <- length(sh)
    if (nrow(md) < 20 || n_distinct(md$group) < 2 || n_sites < 1) next
    md$y <- scale(md$y)[, 1]
    md$site <- droplevels(factor(md$site))
    fit <- if (n_sites > 1) lm(y ~ group + site, md) else lm(y ~ group, md)
    co <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))
    rn <- "groupAD"
    if (!(rn %in% rownames(co))) next
    ad_dir <- panel_full$ad_dir[i]; if (is.na(ad_dir)) ad_dir <- 1
    est <- co[rn, "Estimate"]; lo <- ci[rn, 1]; hi <- ci[rn, 2]
    stats[[length(stats) + 1]] <- tibble(
      phenotype = panel_full$label[i], axis = panel_full$axis[i],
      contrast = paste0("AD_vs_", cmp),
      n_AD = sum(md$group == "AD"), n_cmp = sum(md$group == cmp), n_shared_sites = n_sites,
      std_mean_diff = est,                    # AD minus comparator, in SD units (raw)
      std_error = co[rn, "Std. Error"],
      # AD-oriented: positive = AD has MORE AD-pathology (handles reverse-coded scores)
      effect_ad_oriented = est * ad_dir,
      ci_low = ifelse(ad_dir > 0, lo, -hi), ci_high = ifelse(ad_dir > 0, hi, -lo),
      t = co[rn, "t value"], df = fit$df.residual,
      p_value = co[rn, "Pr(>|t|)"]
    )
  }
}
stats <- bind_rows(stats) %>%
  group_by(contrast) %>% mutate(p_fdr = p.adjust(p_value, "BH")) %>% ungroup() %>%
  arrange(factor(axis, levels = c("AD pathology", "Cognition", "Demographic")), phenotype, contrast)

fwrite(stats, file.path(OUT_DIR, "specificity_stats.csv"))
saveRDS(stats, file.path(OUT_DIR, "specificity_stats.rds"))
cat("\n=== specificity stats (AD vs comparator, SD units) ===\n"); print(as.data.frame(stats), digits = 3)

# ---- descriptive group means ----------------------------------------------
long <- ana %>%
  select(group, all_of(panel_full$tid)) %>%
  pivot_longer(-group, names_to = "tid", values_to = "value") %>%
  left_join(panel_full, by = "tid")
# Full distribution summary per group x phenotype -> redraw box/violin LOCALLY
# (geom_boxplot(stat="identity") with these quantiles; no raw data needed).
desc <- long %>% group_by(axis, label, group) %>%
  summarise(n = sum(!is.na(value)),
            mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            min = suppressWarnings(min(value, na.rm = TRUE)),
            p05 = quantile(value, .05, na.rm = TRUE), p10 = quantile(value, .10, na.rm = TRUE),
            p25 = quantile(value, .25, na.rm = TRUE), median = median(value, na.rm = TRUE),
            p75 = quantile(value, .75, na.rm = TRUE), p90 = quantile(value, .90, na.rm = TRUE),
            p95 = quantile(value, .95, na.rm = TRUE),
            max = suppressWarnings(max(value, na.rm = TRUE)),
            se = sd(value, na.rm = TRUE) / sqrt(pmax(1, sum(!is.na(value)))),
            # Tukey whisker ends (most extreme value within 1.5*IQR) -> exact boxplot
            whisker_lo = { q1 <- quantile(value, .25, na.rm = TRUE); q3 <- quantile(value, .75, na.rm = TRUE)
            f <- q1 - 1.5 * (q3 - q1); suppressWarnings(min(value[!is.na(value) & value >= f])) },
            whisker_hi = { q1 <- quantile(value, .25, na.rm = TRUE); q3 <- quantile(value, .75, na.rm = TRUE)
            f <- q3 + 1.5 * (q3 - q1); suppressWarnings(max(value[!is.na(value) & value <= f])) },
            .groups = "drop")
fwrite(desc, file.path(OUT_DIR, "distribution_summary.csv"))

# Individual OUTLIER points (beyond 1.5*IQR) per group x phenotype -> outlier dots
outliers <- long %>% filter(!is.na(value)) %>% group_by(axis, label, group) %>%
  filter({ q1 <- quantile(value, .25); q3 <- quantile(value, .75); iqr <- q3 - q1
  value < q1 - 1.5 * iqr | value > q3 + 1.5 * iqr }) %>%
  ungroup() %>% select(axis, label, group, value)
cat("\n[boxplot] outlier points exported:", nrow(outliers), "rows\n")
fwrite(outliers, file.path(OUT_DIR, "boxplot_outliers.csv"))

# ---- demographic ACCURACY control (predicted vs OBSERVED, per group) -------
# Positive control. Predicted age/sex levels differ between groups because the
# groups' demographics differ, so the check is accuracy: predicted vs OBSERVED
# age (r) and sex (AUC) within each group.
pred_age_tid <- PANEL$tid[PANEL$name == "age_at_visit"]
pred_sex_tid <- PANEL$tid[PANEL$name == "msex"]

auc_fn <- function(score, label01) {           # manual AUC, no pROC dependency
  ok <- !is.na(score) & !is.na(label01); score <- score[ok]; label01 <- label01[ok]
  np <- sum(label01 == 1); nn <- sum(label01 == 0)
  if (np == 0 || nn == 0) return(NA_real_)
  r <- rank(score); a <- (sum(r[label01 == 1]) - np * (np + 1) / 2) / (np * nn)
  max(a, 1 - a)                                 # orientation-free
}
fisher_ci <- function(r, n) {                   # 95% CI for a correlation
  if (is.na(r) || is.na(n) || n < 4) return(c(NA, NA))
  z <- atanh(r); se <- 1 / sqrt(n - 3); tanh(c(z - 1.96 * se, z + 1.96 * se))
}
auc_ci <- function(auc, np, nn) {               # Hanley-McNeil 95% CI for AUC
  if (is.na(auc) || np == 0 || nn == 0) return(c(NA, NA))
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (np - 1) * (q1 - auc^2) + (nn - 1) * (q2 - auc^2)) / (np * nn))
  c(max(0, auc - 1.96 * se), min(1, auc + 1.96 * se))
}

sx <- tolower(trimws(as.character(ana[[OBS_SEX_COL]])))
cat("\n[sex] observed coding counts (top):\n"); print(utils::head(sort(table(sx, useNA = "always"), decreasing = TRUE), 8))
# Semantic mapping first (handles m/f, male/female, 0/1/2, man/woman). Strictly 0/1/NA.
male_set   <- c("m", "male", "1", "man",   "males")
female_set <- c("f", "female", "0", "2", "woman", "females")
sex01 <- ifelse(sx %in% male_set, 1, ifelse(sx %in% female_set, 0, NA_real_))
# fallback: if semantic map labels <50%, use the two most-common levels
if (mean(!is.na(sex01)) < 0.5) {
  lev <- names(sort(table(sx[!sx %in% c("", "na")]), decreasing = TRUE))
  sex01 <- rep(NA_real_, length(sx))
  if (length(lev) >= 2) { sex01[sx == lev[1]] <- 0; sex01[sx == lev[2]] <- 1 }
  cat("[sex] semantic map weak; fell back to 2 most-common levels:", paste(lev[1:2], collapse = " , "), "\n")
}
# orient so higher predicted P(male) -> 1
if (!is.na(pred_sex_tid) && !all(is.na(sex01)) &&
    suppressWarnings(cor(ana[[pred_sex_tid]], sex01, use = "complete.obs")) < 0) sex01 <- 1 - sex01
# QC: a site whose observed sex is negatively correlated with predicted sex
# (r < -0.3, n >= 10) is treated as having flipped sex coding. Its labels are
# excluded from the sex metric (excluded, not corrected).
if (!is.na(pred_sex_tid)) {
  flip <- ana %>% transmute(site, .p = .data[[pred_sex_tid]], .s = sex01) %>%
    group_by(site) %>%
    summarise(r = suppressWarnings(cor(.p, .s, use = "complete.obs")),
              n = sum(!is.na(.s)), .groups = "drop") %>%
    filter(!is.na(r) & r < -0.3 & n >= 10) %>% pull(site)
  if (length(flip)) {
    cat("[sex] sites with inverted coding (excluded from sex metric):", paste(flip, collapse = ", "), "\n")
    sex01[ana$site %in% flip] <- NA_real_
  }
}
cat("[sex] labeled fraction:", round(mean(!is.na(sex01)), 3), "\n")
cat("[sex] labeled n per group:\n"); print(tapply(!is.na(sex01), ana$group, sum))

acc <- ana %>%
  mutate(.obs_age = suppressWarnings(as.numeric(.data[[OBS_AGE_COL]])),
         .obs_age = ifelse(.obs_age < 0, NA, .obs_age),
         .obs_sex = sex01,
         .pred_age = if (!is.na(pred_age_tid)) .data[[pred_age_tid]] else NA_real_,
         .pred_sex = if (!is.na(pred_sex_tid)) .data[[pred_sex_tid]] else NA_real_) %>%
  group_by(group) %>%
  summarise(n = dplyr::n(),
            age_n   = sum(!is.na(.pred_age) & !is.na(.obs_age)),
            age_r   = suppressWarnings(cor(.pred_age, .obs_age, use = "complete.obs")),
            age_mae = mean(abs(.pred_age - .obs_age), na.rm = TRUE),
            sex_np  = sum(.obs_sex == 1, na.rm = TRUE), sex_nn = sum(.obs_sex == 0, na.rm = TRUE),
            sex_auc = auc_fn(.pred_sex, .obs_sex), .groups = "drop") %>%
  rowwise() %>%
  mutate(age_r_lo = fisher_ci(age_r, age_n)[1], age_r_hi = fisher_ci(age_r, age_n)[2],
         sex_auc_lo = auc_ci(sex_auc, sex_np, sex_nn)[1],
         sex_auc_hi = auc_ci(sex_auc, sex_np, sex_nn)[2]) %>%
  ungroup()
fwrite(acc, file.path(OUT_DIR, "demographic_accuracy_by_group.csv"))
cat("\n=== demographic accuracy by group (positive control) ===\n"); print(as.data.frame(acc), digits = 3)

# ====================================================================================
# Predicted-phenotype AD-classifier SPECIFICITY
# Train an AD-vs-control classifier on the predicted-phenotype panel; compare its
# scores in PD/FTD/ALS with those in AD and in controls.
# ====================================================================================
clf_feats <- panel_full %>% filter(axis %in% c("AD pathology", "Cognition")) %>% pull(tid)
clf_feats <- intersect(clf_feats, colnames(ana))
Xall <- scale(as.matrix(ana[, clf_feats]))                       # standardized features (all groups)
Xdf <- as.data.frame(Xall); colnames(Xdf) <- paste0("f", seq_along(clf_feats))
ff <- as.formula(paste("y ~", paste(colnames(Xdf), collapse = " + ")))

tridx <- which(ana$group %in% c("AD", "control"))
dtr <- Xdf[tridx, , drop = FALSE]; dtr$y <- as.integer(ana$group[tridx] == "AD")
set.seed(1); fold <- sample(rep_len(1:5, nrow(dtr)))             # 5-fold CV out-of-fold AD prob
oof <- rep(NA_real_, nrow(dtr))
for (k in 1:5) {
  m <- glm(ff, binomial, dtr[fold != k, ]); oof[fold == k] <- predict(m, dtr[fold == k, ], type = "response")
}
cv_auc <- auc_fn(oof, dtr$y)
mfull <- glm(ff, binomial, dtr)                                  # refit on AD+control
ana$ad_clf_prob <- as.numeric(predict(mfull, Xdf, type = "response"))   # AD prob for ALL groups
cat(sprintf("\n[AD classifier] AD-vs-control classifier CV AUC = %.3f (features: %d predicted phenotypes)\n",
            cv_auc, length(clf_feats)))

# AD-classifier probability by group (distribution summary for local replot)
clf_grp <- ana %>% group_by(group) %>%
  summarise(n = dplyr::n(), mean = mean(ad_clf_prob), sd = sd(ad_clf_prob),
            p25 = quantile(ad_clf_prob, .25), median = median(ad_clf_prob),
            p75 = quantile(ad_clf_prob, .75),
            whisker_lo = { q1 <- quantile(ad_clf_prob, .25); q3 <- quantile(ad_clf_prob, .75)
            suppressWarnings(min(ad_clf_prob[ad_clf_prob >= q1 - 1.5 * (q3 - q1)])) },
            whisker_hi = { q1 <- quantile(ad_clf_prob, .25); q3 <- quantile(ad_clf_prob, .75)
            suppressWarnings(max(ad_clf_prob[ad_clf_prob <= q3 + 1.5 * (q3 - q1)])) },
            .groups = "drop")
fwrite(clf_grp, file.path(OUT_DIR, "ad_classifier_prob_by_group.csv"))

# specificity AUCs: AD vs each group (does the AD classifier separate AD from other
# diseases?) and non-AD diseases vs control (~0.5 if the classifier does not score that disease above controls)
clf_spec <- bind_rows(
  lapply(c("control", "PD", "FTD", "ALS"), function(g) {
    idx <- ana$group %in% c("AD", g)
    data.frame(contrast = paste0("AD_vs_", g), kind = "AD-discrimination",
               auc = auc_fn(ana$ad_clf_prob[idx], as.integer(ana$group[idx] == "AD")),
               n_AD = sum(ana$group == "AD"), n_other = sum(ana$group == g))
  }),
  lapply(c("PD", "FTD", "ALS"), function(g) {        # AD classifier firing for non-AD disease vs control
    idx <- ana$group %in% c("control", g)
    data.frame(contrast = paste0(g, "_vs_control"), kind = "non-AD firing",
               auc = auc_fn(ana$ad_clf_prob[idx], as.integer(ana$group[idx] == g)),
               n_AD = NA, n_other = sum(ana$group == g))
  }))
clf_spec$cv_auc_AD_vs_control <- cv_auc
fwrite(clf_spec, file.path(OUT_DIR, "ad_classifier_specificity.csv"))
cat("\n=== AD-classifier specificity ===\n"); print(clf_spec, digits = 3)

# ---- figure: predicted phenotype distributions by group --------------------
long$axis  <- factor(long$axis, levels = c("AD pathology", "Cognition", "Demographic"))
long$group <- factor(long$group, levels = c("control", "AD", "PD", "FTD", "ALS"))
PAL <- c(control = "#999999", AD = "#D55E00", PD = "#0072B2", FTD = "#009E73", ALS = "#CC79A7")

p <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.7) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = "Predicted value",
       title = "Predicted phenotypes across pure diagnosis groups (GNPC plasma)",
       subtitle = "CERAD/NIA-Reagan reverse-coded: low = more AD") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7.5), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "FigSx_pdftd_specificity.pdf"), p, width = 11, height = 7)

# ---- figure: AD-oriented effect sizes ---------------------------------------
fp <- stats %>% filter(axis %in% c("AD pathology", "Cognition"))
fp$phenotype <- factor(fp$phenotype, levels = rev(unique(fp$phenotype)))
p_forest <- ggplot(fp, aes(effect_ad_oriented, phenotype, colour = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  position = position_dodge(width = 0.55), size = 0.3) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y") +
  labs(x = "AD-oriented effect (SD): AD − comparator  (positive = more AD pathology in AD)",
       y = NULL, colour = NULL,
       title = "AD-specificity of predicted phenotypes (site-adjusted)") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "FigSx_pdftd_specificity_forest.pdf"), p_forest, width = 8, height = 6)

# ---- figure: demographic accuracy by group (positive control) -------------
accl <- acc %>% select(group, age_r, sex_auc) %>%
  pivot_longer(-group, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, age_r = "Age: predicted-vs-observed r", sex_auc = "Sex: AUC"))
p_acc <- ggplot(accl, aes(group, value, fill = group)) +
  geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Accuracy of predicted age and sex by diagnosis group (positive control)") +
  theme_bw(base_size = 9)
ggsave(file.path(OUT_DIR, "FigSx_demographic_accuracy.pdf"), p_acc, width = 7, height = 4)

# ---- figure: AD-classifier probability by group ----------------------------
clf_long <- ana %>% transmute(group = factor(group, levels = c("control", "AD", "PD", "FTD", "ALS")),
                              ad_clf_prob)
p_clf <- ggplot(clf_long, aes(group, ad_clf_prob, fill = group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.7) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = "Predicted-phenotype AD-classifier probability",
       title = sprintf("AD classifier probability by diagnosis group (CV AUC %.2f, AD vs control)", cv_auc),
       subtitle = "Pure diagnosis groups, plasma") +
  theme_bw(base_size = 9)
ggsave(file.path(OUT_DIR, "FigSx_ad_classifier_specificity.pdf"), p_clf, width = 6, height = 4)

cat("\n=== exported files ===\n")
print(list.files(OUT_DIR))
cat("\nDONE. Outputs in", OUT_DIR, "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
sink(); close(.logf)   # close the run_log.txt capture
