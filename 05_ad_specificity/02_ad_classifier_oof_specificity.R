# ==============================================================================
# Extended Data Fig. 9a: AD-specificity of the AD classifier of Fig. 6a. Uses its
# saved out-of-fold AD probability (one-contributor-out) in pure AD, PD, FTD,
# ALS and control groups; pooled and within-site (stratified) AUCs.
# Inputs : ProteomicsV1_3ms.rds; classifier run p8_one_contributoer_out_nesting_<EX_TAG>
#          in step_dir("04_classifier_runs") (04_ad_classifier/01_train_classifiers.R).
# Outputs: step_dir("05_ad_classifier_oof")
# Runs on the AD Workbench. From the repo root:
#   Rscript 05_ad_specificity/02_ad_classifier_oof_specificity.R
# ==============================================================================

# AD-classifier specificity using the existing CatBoost AD classifier of Fig 6,
# NOT a fresh logistic regression.
# ============================================================================
# Reuses the trained CatBoost AD model and its saved out-of-fold (OOF) AD
# probability (01_pdftd_specificity.R fits its own logistic-regression classifier).
#
# The AD classifier of Fig 6 (one-contributor-out):
#   p8_one_contributoer_out_nesting_<EX_TAG>/
#       trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01/ml/p1_trait/
#   -> per-sample OOF prob in */predictions.rds ($predicted_prob) + sample_id.rds
#
# Fig 6a configuration:
#   control == "Non-AD Controls"  (NOT the "...clean..." Disease-Free variant)
#   feature  == "All phenotypes"  == ex.Clinical.Diagnosis, metric_th 0.01
#   split    == one-contributor-out (excluded contributors and fold N not scored)
# => primary analysis below; CLEAN_CONTROLS = FALSE matches Fig 6.
#
# DESIGN:
#   Fig 6a used "Non-AD Controls", so pure PD/FTD/ALS were in that non-AD control
#   pool and ALSO have a saved OOF predicted_prob. We simply READ the saved OOF for
#   every group (one-contributor-out => each sample's prob is held out from its own
#   contributor). No model is re-applied.
#
#   Coverage is partial by design: contributors with too few AD cases/controls to
#   form a one-contributor-out fold contribute NO OOF, so some samples in every
#   group (AD, control, PD, FTD, ALS) have no prob. We analyse ONLY the genuine
#   held-out OOF subset. The coverage per group is printed below.
#
#   CAVEAT: PD/FTD/ALS were non-AD *training negatives*, so a low AD-prob is the
#   deployed classifier's held-out behaviour, not a fully naive-to-disease
#   specificity.
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(tidyr)
  library(data.table); library(purrr)
})
source("config.R")

# ---- config ----------------------------------------------------------------
CLEAN_CONTROLS <- FALSE   # FALSE = Fig-6a "Non-AD Controls" (primary); TRUE = "clean" run

JOBROOT <- file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG))
FEAT_DIR_NAME <- if (CLEAN_CONTROLS)
  "trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01_clean_control" else
    "trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01"
P1     <- file.path(JOBROOT, FEAT_DIR_NAME, "ml/p1_trait")

OUT_DIR  <- step_dir("05_ad_classifier_oof")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
.logf <- file(file.path(OUT_DIR, "run_log.txt"), open = "wt"); sink(.logf, split = TRUE)

# ---- metadata + pure diagnosis groups (same definition as 01_pdftd_specificity) ----
expdata  <- readRDS(PROTEOMICS_RDS)
metadata <- as.data.frame(colData(expdata))
metadata <- metadata %>%
  filter(grepl("plasma", tolower(trimws(as.character(sample_matrix))))) %>%   # plasma only
  filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)

dx_all <- c("ad", "pd", "ftd", "als", "stroke"); nd <- c("ad", "pd", "ftd", "als")
DX <- metadata %>% dplyr::select(all_of(dx_all)) %>%
  mutate(across(everything(), ~ ifelse(.x %in% c(0, 1), .x, NA_integer_)))
npos_all <- rowSums(DX == 1, na.rm = TRUE)
npos_nd  <- rowSums(DX[nd] == 1, na.rm = TRUE)
metadata$group <- NA_character_
metadata$group[npos_all == 0] <- "control"
for (d in nd) metadata$group[npos_nd == 1 & DX[[d]] %in% 1] <- toupper(d)

# ---- (A) SAVED OOF predicted_prob (all groups) -----------------------------
# Per-fold held-out
# predictions live in submodels/*/predictions.rds; the OOF prediction is the
# split == "test" row (each sample held out from its own contributor's fold).
sample_id   <- readRDS(file.path(P1, "sample_id.rds"))
pred_files  <- list.files(file.path(P1, "submodels"), pattern = "predictions.rds",
                          recursive = TRUE, full.names = TRUE)
oof <- pred_files %>% map(readRDS) %>% map(~ dplyr::filter(.x, split == "test")) %>% bind_rows()
oof$sample_id <- sample_id[as.numeric(oof$sample_id)]
oof <- oof %>% dplyr::select(sample_id, ad_clf_prob = predicted_prob) %>%
  distinct(sample_id, .keep_all = TRUE)
cat("\n[OOF] samples with saved AD-classifier prob:", nrow(oof), "\n")

ana <- metadata %>% filter(!is.na(group)) %>%
  left_join(oof, by = c("SAMPLE.ID" = "sample_id"))
cat("[OOF] coverage by group (has saved OOF prob):\n")
print(ana %>% group_by(group) %>% summarise(n = dplyr::n(), n_oof = sum(!is.na(ad_clf_prob)), .groups = "drop"))

# ---- (B) restrict to the OOF subset ----------------------------------------
# All groups (incl. PD/FTD/ALS) carry a saved OOF prob,
# because Fig 6a used "Non-AD Controls" (pure PD/FTD/ALS were in that pool). Missing
# samples come from contributors dropped by the one-contributor-out design (too few
# AD cases/controls to form a fold). We analyse only the
# genuine held-out OOF predictions; no model application needed.
ana <- ana %>% filter(!is.na(ad_clf_prob))
cat("\n[final] AD-classifier OOF analysis set (held-out predictions only):\n")
print(ana %>% group_by(group) %>% summarise(n_oof = dplyr::n(), .groups = "drop"))

# ---- helpers: directed AUC + Hanley-McNeil 95% CI ---------------------------
# NB: unlike auc_fn() in 01_pdftd_specificity.R, the AUC here is NOT folded with
# max(a, 1-a) and its CI is NOT clamped at 0.5:
#   * folding is one-directional, so it pushes every estimate to >= 0.5 and hides
#     a classifier that scores a group BELOW controls (directed AUC < 0.5).
#   * clamping lo at 0.5 makes the "non-AD firing" rows one-sided by construction, so
#     their CIs could never cross below 0.5.
# The direction is part of what is tested here, so it is kept.
auc_directed <- function(score, label01) {     # NO max(a, 1-a): direction is kept
  ok <- !is.na(score) & !is.na(label01); score <- score[ok]; label01 <- label01[ok]
  np <- sum(label01 == 1); nn <- sum(label01 == 0)
  if (np == 0 || nn == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label01 == 1]) - np * (np + 1) / 2) / (np * nn)
}
hm_se <- function(auc, np, nn) {               # Hanley-McNeil SE (approx)
  if (is.na(auc) || np == 0 || nn == 0) return(NA_real_)
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  sqrt((auc * (1 - auc) + (np - 1) * (q1 - auc^2) + (nn - 1) * (q2 - auc^2)) / (np * nn))
}
auc_fn <- function(score, label01) {
  ok <- !is.na(score) & !is.na(label01); score <- score[ok]; label01 <- label01[ok]
  np <- sum(label01 == 1); nn <- sum(label01 == 0)
  if (np == 0 || nn == 0) return(c(auc = NA_real_, lo = NA_real_, hi = NA_real_))
  a  <- auc_directed(score, label01)
  # Hanley-McNeil is derived for AUC >= 0.5 and is not symmetric under a -> 1-a,
  # so evaluate the SE at the folded value (same number either way, hence
  # invariant to which group is called "positive") and apply it around the
  # DIRECTED estimate. Clamp to [0,1], not [0.5,1].
  se <- hm_se(max(a, 1 - a), np, nn)
  c(auc = a, lo = max(0, a - 1.96 * se), hi = min(1, a + 1.96 * se))
}

# ---- specificity AUCs WITH CIs ---------------------------------------------
mk <- function(g, kind) {
  idx <- ana$group %in% (if (kind == "AD-discrimination") c("AD", g) else c("control", g))
  pos <- if (kind == "AD-discrimination") "AD" else g
  a <- auc_fn(ana$ad_clf_prob[idx], as.integer(ana$group[idx] == pos))
  data.frame(contrast = if (kind == "AD-discrimination") paste0("AD_vs_", g) else paste0(g, "_vs_control"),
             kind = kind, auc = a[1], auc_lo = a[2], auc_hi = a[3],
             n_pos = sum(ana$group[idx] == pos), n_neg = sum(idx) - sum(ana$group[idx] == pos))
}
clf_spec <- bind_rows(
  lapply(c("control", "PD", "FTD", "ALS"), mk, kind = "AD-discrimination"),
  lapply(c("PD", "FTD", "ALS"), mk, kind = "non-AD firing")
)
fwrite(clf_spec, file.path(OUT_DIR, "ad_classifier_specificity_oof.csv"))
cat("\n=== AD-classifier specificity (OOF-based, with 95% CI) ===\n"); print(clf_spec, digits = 3)

# ============================================================================
# WITHIN-SITE (STRATIFIED) AUC — site-confound-free version of the block above
# ============================================================================
# WHY: the pooled AUC above compares an AD sample at site X against a PD sample
# at site Y, so contributor and group are tangled. 01_pdftd_specificity.R
# already restricts its standardized-mean-difference contrasts to sites holding
# BOTH groups; this does the
# same for the classifier AUC. Only pairs from the SAME site are counted, so the
# site contrast is removed by construction (not merely adjusted for).
#
#   AUC_strat = sum_s w_s * AUC_s ,  w_s = (np_s * nn_s) / sum_t (np_t * nn_t)
#
# i.e. the fraction of concordant pairs among all WITHIN-site pairs.
#
# Both this section and the pooled block above use the shared directed AUC
# (auc_directed / hm_se in the helpers section) -- see the folding note there.

# site = contributor x matrix, matching 01_pdftd_specificity.R (plasma-only
# here, so this is effectively contributor, but keep the definition identical).
ana$site <- paste0(ana$contributor_code, "_", ana$sample_matrix)

# A "SOLID" site has >= MIN_PER_GROUP_SITE samples in BOTH groups. The primary
# stratified AUC, its bootstrap CI, and the heterogeneity test are ALL computed on
# solid sites only, so the three agree with each other.
#
# WHY: pair-count weighting gives a site weight proportional to np*nn, so a site
# with a huge majority class and a handful of minority-class samples can still
# draw real weight. And the unweighted per-site range is set by very small sites.
# Filtering before the estimator removes both problems.
MIN_N_SHARED       <- 20   # same estimability guard as 01_pdftd_specificity.R
MIN_PER_GROUP_SITE <- 10   # "solid" site: >= this many in BOTH groups
B_BOOT             <- 2000

# auc_directed() and hm_se() are defined once in the helpers section above and
# shared with the pooled block, so both blocks use the same directed AUC.

# per-site AUC table for one contrast. Sites lacking either group are dropped
# (not shared). Every remaining shared site is kept in the table for transparency
# and flagged `solid`; `w` is its weight in the all-sites estimate and `w_solid`
# its weight in the primary (solid-only) estimate.
per_site_auc <- function(d, pos_group) {
  out <- lapply(split(seq_len(nrow(d)), d$site), function(ix) {
    y <- as.integer(d$group[ix] == pos_group)
    np <- sum(y == 1); nn <- sum(y == 0)
    if (np == 0 || nn == 0) return(NULL)       # not a shared site -> excluded
    a <- auc_directed(d$ad_clf_prob[ix], y)
    data.frame(site = d$site[ix][1], np = np, nn = nn, auc = a,
               se = hm_se(a, np, nn), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(out)) return(NULL)
  out$pairs   <- out$np * out$nn
  out$w       <- out$pairs / sum(out$pairs)
  out$solid   <- out$np >= MIN_PER_GROUP_SITE & out$nn >= MIN_PER_GROUP_SITE
  out$w_solid <- ifelse(out$solid, out$pairs / sum(out$pairs[out$solid]), 0)
  out[order(-out$pairs), ]
}

# two-stage cluster bootstrap over the SOLID sites: resample sites with
# replacement, then subjects within each drawn site. Stage 1 carries the
# between-site variance the pooled Hanley-McNeil CI ignores; stage 2 keeps the
# within-site sampling noise a site-only bootstrap would drop.
#
# The site set here MUST match the one used for the point estimate, or the CI
# would describe a different quantity than the value it is drawn around.
boot_strat_auc <- function(d, pos_group, keep_sites, B = B_BOOT, seed = 1) {
  set.seed(seed)
  cl <- split(seq_len(nrow(d)), d$site)
  cl <- cl[names(cl) %in% keep_sites]
  if (length(cl) < 1) return(c(lo = NA, hi = NA, valid = 0))
  reps <- replicate(B, {
    num <- 0; den <- 0
    for (s in sample(names(cl), length(cl), replace = TRUE)) {
      ix <- cl[[s]]; ib <- sample(ix, length(ix), replace = TRUE)
      y  <- as.integer(d$group[ib] == pos_group)
      np <- sum(y == 1); nn <- sum(y == 0)
      if (np == 0 || nn == 0) next             # degenerate redraw -> drop this site
      a <- auc_directed(d$ad_clf_prob[ib], y)
      num <- num + a * np * nn; den <- den + np * nn
    }
    if (den == 0) NA_real_ else num / den
  })
  c(lo = unname(quantile(reps, .025, na.rm = TRUE)),
    hi = unname(quantile(reps, .975, na.rm = TRUE)),
    valid = sum(!is.na(reps)))
}

mk_strat <- function(g, kind) {
  cmp      <- if (kind == "AD-discrimination") "AD" else "control"
  pos      <- if (kind == "AD-discrimination") "AD" else g
  contrast <- if (kind == "AD-discrimination") paste0("AD_vs_", g) else paste0(g, "_vs_control")
  d <- ana[ana$group %in% c(cmp, g), c("site", "group", "ad_clf_prob")]
  
  ps <- per_site_auc(d, pos)
  base <- data.frame(
    contrast = contrast, kind = kind,
    n_sites_pos = length(unique(d$site[d$group == pos])),
    n_sites_neg = length(unique(d$site[d$group != pos])),
    n_shared_sites = 0L, n_solid_sites = 0L,
    n_pos_solid = 0L, n_neg_solid = 0L,
    estimable = FALSE, reason = "no shared site",
    auc_pooled_directed = auc_directed(d$ad_clf_prob, as.integer(d$group == pos)),
    auc_stratified = NA_real_,                 # primary estimate: solid sites only
    auc_strat_lo = NA_real_, auc_strat_hi = NA_real_, boot_valid = NA_integer_,
    auc_strat_allsites = NA_real_,             # sensitivity: every shared site
    auc_solid_min = NA_real_, auc_solid_max = NA_real_,
    eff_sites = NA_real_, max_w = NA_real_, top_site = NA_character_,
    single_site = NA, Q = NA_real_, Q_df = NA_integer_, Q_p = NA_real_, I2 = NA_real_,
    stringsAsFactors = FALSE)
  if (is.null(ps)) return(base)
  
  base$n_shared_sites     <- nrow(ps)
  base$auc_strat_allsites <- sum(ps$w * ps$auc)
  
  sd_ps <- ps[ps$solid, , drop = FALSE]        # solid sites drive everything below
  base$n_solid_sites <- nrow(sd_ps)
  base$n_pos_solid   <- sum(sd_ps$np)
  base$n_neg_solid   <- sum(sd_ps$nn)
  n_solid <- base$n_pos_solid + base$n_neg_solid
  
  if (nrow(sd_ps) < 1) {
    base$reason <- sprintf("no site has >=%d in both groups (%d shared site(s), %d samples)",
                           MIN_PER_GROUP_SITE, nrow(ps), sum(ps$np) + sum(ps$nn))
    attr(base, "per_site") <- transform(ps, contrast = contrast)
    return(base)
  }
  if (n_solid < MIN_N_SHARED) {
    base$reason <- sprintf("only %d samples in %d solid site(s) (<%d)",
                           n_solid, nrow(sd_ps), MIN_N_SHARED)
    attr(base, "per_site") <- transform(ps, contrast = contrast)
    return(base)
  }
  
  base$estimable      <- TRUE
  base$reason         <- ""
  base$auc_stratified <- sum(sd_ps$w_solid * sd_ps$auc)
  base$auc_solid_min  <- min(sd_ps$auc); base$auc_solid_max <- max(sd_ps$auc)
  # inverse-Simpson effective number of sites: 10 sites means little if one
  # contributor carries most of the weight.
  base$eff_sites   <- 1 / sum(sd_ps$w_solid^2)
  base$max_w       <- max(sd_ps$w_solid)
  base$top_site    <- sd_ps$site[which.max(sd_ps$w_solid)]
  base$single_site <- nrow(sd_ps) < 2          # CI then reflects within-site noise only
  
  bt <- boot_strat_auc(d, pos, keep_sites = sd_ps$site)
  base$auc_strat_lo <- bt["lo"]; base$auc_strat_hi <- bt["hi"]
  base$boot_valid   <- as.integer(bt["valid"])
  
  # Cochran Q across the same solid sites. Diagnostic only -- the point estimate
  # stays pair-count weighted, not inverse-variance weighted.
  h <- sd_ps[!is.na(sd_ps$se) & sd_ps$se > 0, ]
  if (nrow(h) >= 2) {
    wi <- 1 / h$se^2; fixed <- sum(wi * h$auc) / sum(wi)
    Q  <- sum(wi * (h$auc - fixed)^2); dfQ <- nrow(h) - 1
    base$Q <- Q; base$Q_df <- dfQ
    base$Q_p <- pchisq(Q, dfQ, lower.tail = FALSE)
    base$I2  <- max(0, (Q - dfQ) / Q)
  }
  attr(base, "per_site") <- transform(ps, contrast = contrast)
  base
}

strat_list <- c(lapply(c("control", "PD", "FTD", "ALS"), mk_strat, kind = "AD-discrimination"),
                lapply(c("PD", "FTD", "ALS"), mk_strat, kind = "non-AD firing"))
clf_strat  <- bind_rows(strat_list)
per_site   <- bind_rows(lapply(strat_list, attr, "per_site"))

fwrite(clf_strat, file.path(OUT_DIR, "ad_classifier_specificity_oof_stratified.csv"))
fwrite(per_site,  file.path(OUT_DIR, "ad_classifier_specificity_oof_persite.csv"))

cat("\n=== AD-classifier specificity: WITHIN-SITE (stratified, solid sites only) ===\n")
cat("    solid site = >=", MIN_PER_GROUP_SITE, "in BOTH groups; the point estimate,\n",
    "   the bootstrap CI and the Q/I2 test all use the same site set.\n")
print(clf_strat[, c("contrast", "n_shared_sites", "n_solid_sites", "n_pos_solid",
                    "n_neg_solid", "estimable", "auc_pooled_directed",
                    "auc_stratified", "auc_strat_lo", "auc_strat_hi",
                    "auc_strat_allsites", "auc_solid_min", "auc_solid_max",
                    "eff_sites", "max_w", "I2", "reason")], digits = 3)
cat("\n--- per-site AUCs (all shared sites; `solid` flags those used) ---\n")
print(as.data.frame(per_site), digits = 3)

ne <- clf_strat$contrast[!clf_strat$estimable]
if (length(ne)) cat("\n[NOT ESTIMABLE within site]:", paste(ne, collapse = ", "),
                    "-- the pooled AUC for these is not site-adjusted.\n")
ss <- clf_strat$contrast[clf_strat$estimable & clf_strat$single_site %in% TRUE]
if (length(ss)) cat("[EFFECTIVELY SINGLE-SITE]:", paste(ss, collapse = ", "),
                    "-- one solid contributor; the bootstrap redraws that same site\n",
                    "   every replicate, so the CI covers within-site noise ONLY and\n",
                    "   understates uncertainty (single-contributor estimate).\n")
wc <- clf_strat$contrast[clf_strat$estimable & !(clf_strat$single_site %in% TRUE) &
                           clf_strat$max_w > 0.5]
if (length(wc)) cat("[WEIGHT CONCENTRATED >50% IN ONE SITE]:", paste(wc, collapse = ", "),
                    "-- see eff_sites alongside n_solid_sites.\n")

# ---- AD-classifier prob by group (Tukey box from quantiles) ---------------
clf_grp <- ana %>% filter(!is.na(ad_clf_prob)) %>% group_by(group) %>%
  summarise(n = dplyr::n(), mean = mean(ad_clf_prob), sd = sd(ad_clf_prob),
            p25 = quantile(ad_clf_prob, .25), median = median(ad_clf_prob), p75 = quantile(ad_clf_prob, .75),
            whisker_lo = { q1 <- quantile(ad_clf_prob, .25); q3 <- quantile(ad_clf_prob, .75)
            suppressWarnings(min(ad_clf_prob[ad_clf_prob >= q1 - 1.5 * (q3 - q1)])) },
            whisker_hi = { q1 <- quantile(ad_clf_prob, .25); q3 <- quantile(ad_clf_prob, .75)
            suppressWarnings(max(ad_clf_prob[ad_clf_prob <= q3 + 1.5 * (q3 - q1)])) },
            .groups = "drop")
fwrite(clf_grp, file.path(OUT_DIR, "ad_classifier_prob_by_group_oof.csv"))
cat("\n=== AD-classifier prob by group (OOF) ===\n"); print(as.data.frame(clf_grp), digits = 3)

cat("\nDONE. Outputs in", OUT_DIR, "\n")
cat("Within-site outputs: ad_classifier_specificity_oof_stratified.csv (within-site AUC)\n",
    "                     ad_classifier_specificity_oof_persite.csv    (per-site detail)\n")
sink(); close(.logf)
