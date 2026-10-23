# ==============================================================================
# Extended Data Fig. 9b (per-site points): accuracy of predicted age (r) and
# sex (AUC) per contributor site within each pure diagnosis group.
# Inputs : ProteomicsV1_3ms.rds; predicted_phenotypes.rds (prediction "full",
#          target_annotation); demographic_accuracy_by_group.csv from
#          step_dir("05_pdftd_specificity") (01_pdftd_specificity.R, run it first).
# Outputs: step_dir("05_demographic_persite")
# Runs on the AD Workbench. From the repo root:
#   Rscript 05_ad_specificity/03_demographic_accuracy_persite.R
# ==============================================================================

# Extended Data Fig. 9b: demographic accuracy PER CONTRIBUTOR SITE
# ============================================================================
# Individual data points for ED Fig. 9b (bars of predicted-vs-observed age r and
# sex AUC per disease group). A correlation or an AUC has no per-person value,
# so the "individual points" are contributor sites: one r and one AUC per site
# within each disease group, as in panel 9a.
#
# Sample set = 01_pdftd_specificity.R exactly (plasma only, excluded contributors out,
# pure diagnosis groups, sites with inverted sex coding excluded from the sex
# metric). The by-group rows written here must reproduce
# demographic_accuracy_by_group.csv; the script stops if they do not.
#
# Only aggregate values leave the workbench: sites below MIN_N are not written
# (their count is), so no row describes a handful of people.
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(data.table)
})
source("config.R")

OUT_DIR  <- step_dir("05_demographic_persite")

MIN_N <- 10   # age: >= 10 people with both ages; sex: >= 10 of EACH sex (as panel 9a)

.logf <- file(file.path(OUT_DIR, "run_log_06_persite.txt"), open = "wt"); sink(.logf, split = TRUE)

# ---- load (as 01_pdftd_specificity.R) --------------------------------------
pp         <- read_predicted_phenotypes(PREDICTED_RDS)
expdata    <- readRDS(PROTEOMICS_RDS)
metadata   <- as.data.frame(colData(expdata))
prediction <- get_prediction(pp, "full")
tannot     <- get_table(pp, "target_annotation")
name2id    <- setNames(tannot$`Target ID`, tannot$Target)
pred_age_tid <- name2id[["age_at_visit"]]
pred_sex_tid <- name2id[["msex"]]
stopifnot(pred_age_tid %in% colnames(prediction), pred_sex_tid %in% colnames(prediction))

data <- metadata %>%
  inner_join(prediction %>% dplyr::select(sample_id, all_of(c(pred_age_tid, pred_sex_tid))),
             by = "sample_id") %>%
  filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS) %>%
  filter(grepl("plasma", tolower(trimws(as.character(sample_matrix)))))
data$site <- paste0(data$contributor_code, "_", data$sample_matrix)

dx_all <- c("ad", "pd", "ftd", "als", "stroke")
nd     <- c("ad", "pd", "ftd", "als")
DX <- data %>% dplyr::select(all_of(dx_all)) %>%
  mutate(across(everything(), ~ { v <- suppressWarnings(as.numeric(as.character(.x)))
  ifelse(v %in% c(0, 1), v, NA) }))
npos_all <- rowSums(DX == 1, na.rm = TRUE)
npos_nd  <- rowSums(DX[nd] == 1, na.rm = TRUE)
group <- rep(NA_character_, nrow(data))
group[npos_all == 0] <- "control"
for (d in nd) group[npos_nd == 1 & DX[[d]] %in% 1] <- toupper(d)
data$group <- factor(group, levels = c("control", "AD", "PD", "FTD", "ALS"))
ana <- data %>% filter(!is.na(group))
cat("=== group sizes (pure, plasma) ===\n"); print(table(ana$group))

# ---- observed sex, oriented and QC'd (as 01) -------------------------------
sx <- tolower(trimws(as.character(ana$sex)))
male_set   <- c("m", "male", "1", "man",   "males")
female_set <- c("f", "female", "0", "2", "woman", "females")
sex01 <- ifelse(sx %in% male_set, 1, ifelse(sx %in% female_set, 0, NA_real_))
stopifnot(mean(!is.na(sex01)) >= 0.5)   # no fallback coding here (unlike 01); stop rather than guess
if (suppressWarnings(cor(ana[[pred_sex_tid]], sex01, use = "complete.obs")) < 0) sex01 <- 1 - sex01
flip <- ana %>% transmute(site, .p = .data[[pred_sex_tid]], .s = sex01) %>%
  group_by(site) %>%
  summarise(r = suppressWarnings(cor(.p, .s, use = "complete.obs")),
            n = sum(!is.na(.s)), .groups = "drop") %>%
  filter(!is.na(r) & r < -0.3 & n >= 10) %>% pull(site)
cat("[sex] sites with inverted coding (excluded from sex metric):", paste(flip, collapse = ", "), "\n")
sex01[ana$site %in% flip] <- NA_real_

d <- ana %>%
  transmute(group, site,
            obs_age  = { a <- suppressWarnings(as.numeric(age_at_visit)); ifelse(a < 0, NA, a) },
            pred_age = .data[[pred_age_tid]],
            obs_sex  = sex01,
            pred_sex = .data[[pred_sex_tid]])

# Directed AUC (P(male) higher in males). Per site a folded max(a, 1-a) would hide
# a site the model gets backwards; at the group level the check below requires
# the directed AUC to match 01's folded AUC.
auc_dir <- function(score, y) {
  ok <- !is.na(score) & !is.na(y); score <- score[ok]; y <- y[ok]
  np <- sum(y == 1); nn <- sum(y == 0)
  if (np == 0 || nn == 0) return(NA_real_)
  (sum(rank(score)[y == 1]) - np * (np + 1) / 2) / (np * nn)
}

summ <- function(x) x %>% summarise(
  n       = dplyr::n(),
  age_n   = sum(!is.na(pred_age) & !is.na(obs_age)),
  age_r   = if (age_n >= 3) cor(pred_age, obs_age, use = "complete.obs") else NA_real_,
  sex_np  = sum(obs_sex == 1 & !is.na(pred_sex), na.rm = TRUE),
  sex_nn  = sum(obs_sex == 0 & !is.na(pred_sex), na.rm = TRUE),
  sex_auc = auc_dir(pred_sex, obs_sex),
  .groups = "drop")

# ---- check: group level reproduces demographic_accuracy_by_group.csv -------
by_group <- d %>% group_by(group) %>% summ()
old <- fread(file.path(step_dir("05_pdftd_specificity"), "demographic_accuracy_by_group.csv"))
chk <- merge(as.data.table(by_group)[, group := as.character(group)], old,
             by = "group", suffixes = c("", "_old"))
cat("\n=== check vs demographic_accuracy_by_group.csv ===\n")
print(chk[, .(group, age_n, age_n_old, age_r, age_r_old, sex_auc, sex_auc_old)], digits = 6)
stopifnot(nrow(chk) == nrow(old),
          all(chk$age_n == chk$age_n_old),
          all(abs(chk$age_r - chk$age_r_old) < 1e-10),
          all(abs(chk$sex_auc - chk$sex_auc_old) < 1e-10))

# ---- per site within each group --------------------------------------------
per <- d %>% group_by(group, site) %>% summ() %>%
  mutate(age_ok = age_n >= MIN_N,
         sex_ok = sex_np >= MIN_N & sex_nn >= MIN_N & !site %in% flip)

excluded <- per %>% group_by(group) %>%
  summarise(sites_total = dplyr::n(),
            age_sites_used = sum(age_ok), age_sites_excluded = sum(!age_ok),
            sex_sites_used = sum(sex_ok), sex_sites_excluded = sum(!sex_ok),
            sex_sites_inverted = sum(site %in% flip), .groups = "drop")

# Only sites that pass MIN_N leave the workbench; the other metric is blanked.
out <- per %>% filter(age_ok | sex_ok) %>%
  mutate(age_r   = ifelse(age_ok, age_r, NA_real_),
         age_n   = ifelse(age_ok, age_n, NA_integer_),
         sex_auc = ifelse(sex_ok, sex_auc, NA_real_),
         sex_np  = ifelse(sex_ok, sex_np, NA_integer_),
         sex_nn  = ifelse(sex_ok, sex_nn, NA_integer_)) %>%
  dplyr::select(group, site, age_n, age_r, sex_np, sex_nn, sex_auc) %>%
  arrange(group, site)

fwrite(out,      file.path(OUT_DIR, "demographic_accuracy_persite.csv"))
fwrite(excluded, file.path(OUT_DIR, "demographic_accuracy_persite_sites.csv"))
cat("\n=== sites used / excluded (MIN_N =", MIN_N, ") ===\n"); print(as.data.frame(excluded))
cat("\n=== per-site values ===\n"); print(as.data.frame(out), digits = 3)

sink(); close(.logf)
