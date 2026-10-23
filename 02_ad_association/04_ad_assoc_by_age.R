# ---- GNPC code release ---------------------------------------------------------
# Part of: Supplementary Table 7 (age-split columns, <65 / >=65).
# Makes: GNPC predicted-phenotype -> AD association (Fig. 4g model), recomputed by
#   age stratum (all, <65, >=65).
# Inputs:  PROTEOMICS_RDS and PREDICTED_RDS (set in config.R).
# Outputs: step_dir("02_ad_assoc_by_age").
# Runs on the AD Workbench (the GNPC data stay there). From the repository root:
#   Rscript 02_ad_association/04_ad_assoc_by_age.R
# ---------------------------------------------------------------------------------
# Fig 4g cross-cohort concordance, STRATIFIED BY AGE (GNPC side)
# ============================================================================
# Fig 4g compares, per predicted phenotype, the ROSMAP MEASURED-phenotype->AD beta
# (x, computed in ROSMAP, an older cohort) vs the GNPC PREDICTED-phenotype->AD beta
# (y, this script). Question: do the associations differ in younger people, where
# predicted age has a floor near 65? Here we recompute the GNPC y-axis SEPARATELY
# for age<65 and age>=65, so each stratum can be compared with the same ROSMAP
# reference (e.g. Spearman rho; betas that shrink toward 0 in the young). That
# comparison is not part of this script.
#
# GNPC model = EXACT Fig 4g spec (02_ad_assoc_site_adjusted.Rmd):
#   site-adjusted logistic  glm(ad ~ scale(predictor) + contributor_code, binomial)
#   plasma pooled ("all_plasma"), excluded contributors (config.R) removed, AD cases vs clean controls.
# Only the sample is split by observed age; the model is otherwise identical.
# Two-sided Wald z. n and n_case reported per stratum (fewer AD cases give less
# precise betas).
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(tidyr); library(data.table)
})

source("config.R")
OUT_DIR  <- step_dir("02_ad_assoc_by_age")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
.logf <- file(file.path(OUT_DIR, "run_log.txt"), open = "wt"); sink(.logf, split = TRUE)

AGE_THRESH <- 65
TARGET <- "ad"
CSF_SERUM <- c("Q_CSF", "U_Serum", "N_CSF", "T_CSF")   # the all_plasma exclusion (Fig 4g)

# ---- load + build the Fig-4g GNPC data frame (as in 02_ad_assoc_site_adjusted.Rmd) ---
expdata    <- readRDS(PROTEOMICS_RDS)
metadata   <- as.data.frame(colData(expdata))
pp         <- read_predicted_phenotypes(PREDICTED_RDS)
prediction <- get_prediction(pp, "full")
# drop the base predicted column when a *_class* version exists (as in Fig 4g)
drop <- colnames(prediction) %in% (grep("class", colnames(prediction), value = TRUE) %>% gsub("_class.+", "", .))
prediction <- prediction[, !drop]

data <- metadata %>% inner_join(prediction, by = "sample_id") %>% filter(!contributor_code %in% EXCLUDED_CONTRIBUTORS)
data$contributor_code <- paste0(data$contributor_code, "_", data$sample_matrix)
# clean control = free of every neuro dx (== "1" string match, as Fig 4g)
dis <- sapply(c("pd", "als", "ad", "ftd", "stroke"), function(d) as.character(data[[d]]) == "1")
data$control <- rowSums(dis, na.rm = TRUE) == 0
data$age_obs <- { a <- suppressWarnings(as.numeric(as.character(data$age_at_visit))); ifelse(a < 0, NA, a) }

# plasma pooled ("all_plasma"): exclude the CSF/serum contributor-matrix combos
plasma <- data %>% filter(!contributor_code %in% CSF_SERUM)
cat(sprintf("all_plasma N=%d ; with observed age=%d ; AD(0/1)=%s\n",
            nrow(plasma), sum(!is.na(plasma$age_obs)),
            paste(table(suppressWarnings(as.numeric(as.character(plasma$ad)))), collapse = "/")))

# predictor columns: every predicted-phenotype column, in file order, dropping base
# cols that have a _class
predictor_names <- intersect(names(plasma), setdiff(names(prediction), "sample_id"))
predictor_names <- setdiff(predictor_names,
                           grep("class", predictor_names, value = TRUE) %>% gsub("_class.+", "", .) %>% unique())
predictor_names <- predictor_names[vapply(predictor_names, function(c) is.numeric(plasma[[c]]) || all(plasma[[c]] %in% c(0, 1, NA)), logical(1))]
# NB: this is the FULL predictor set of 02_ad_assoc_site_adjusted.Rmd (~94 phenotypes,
# but discrete ones expand into multiple _class_k columns -> ~379 columns). Fig 4g
# narrows this to <100 at plot time: keep the SHAP features, collapse each discrete
# phenotype's classes to its most-significant one, drop Clinical Diagnosis.
cat(sprintf("predictors (full set; Fig 4g subsets to <100 at plot time): %d\n", length(predictor_names)))

# ---- site-adjusted logistic beta for one predictor on one age subset ---------
assoc_one <- function(sub, predictor) {
  md <- as.data.frame(sub)[, c(TARGET, predictor, "control", "contributor_code")]
  names(md)[1:2] <- c("y", "x")
  md$y <- suppressWarnings(as.numeric(as.character(md$y)))
  md$y[md$y == -1] <- NA
  md <- md[stats::complete.cases(md), ]
  md <- md[(md$y == 0 & md$control) | md$y == 1, ]        # AD cases vs clean controls
  na_row <- data.table(predictor = predictor, beta = NA_real_, std_error = NA_real_,
                       z_value = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                       odds_ratio = NA_real_, n = nrow(md), n_case = sum(md$y == 1, na.rm = TRUE))
  if (nrow(md) < 10 || length(unique(md$y)) < 2) return(na_row)
  if (sd(md$x, na.rm = TRUE) == 0) return(na_row)
  md$x <- as.numeric(scale(md$x))
  md$contributor_code <- droplevels(factor(md$contributor_code))
  fml <- if (nlevels(md$contributor_code) >= 2) y ~ x + contributor_code else y ~ x
  tryCatch({
    m  <- glm(fml, family = binomial(link = "logit"), data = md)
    cc <- summary(m)$coefficients
    b <- cc["x", "Estimate"]; se <- cc["x", "Std. Error"]
    data.table(predictor = predictor, beta = b, std_error = se,
               z_value = cc["x", "z value"], p_value = cc["x", "Pr(>|z|)"],
               ci_low = b - 1.96 * se, ci_high = b + 1.96 * se, odds_ratio = exp(b),
               n = nrow(md), n_case = sum(md$y == 1))
  }, error = function(e) { na_row$note <- conditionMessage(e); na_row })
}

# ---- run the three strata ---------------------------------------------------
STRATA <- list(
  full = rep(TRUE, nrow(plasma)),
  lt65 = !is.na(plasma$age_obs) & plasma$age_obs <  AGE_THRESH,
  ge65 = !is.na(plasma$age_obs) & plasma$age_obs >= AGE_THRESH)

res <- rbindlist(lapply(names(STRATA), function(sn) {
  sub <- plasma[STRATA[[sn]], ]
  ncase <- sum(suppressWarnings(as.numeric(as.character(sub$ad))) == 1 & (sub$control | TRUE), na.rm = TRUE)
  cat(sprintf("\n[stratum %s] N=%d, AD cases=%d -> fitting %d predictors...\n",
              sn, nrow(sub), sum(suppressWarnings(as.numeric(as.character(sub$ad))) == 1, na.rm = TRUE),
              length(predictor_names)))
  out <- rbindlist(lapply(predictor_names, function(p) assoc_one(sub, p)), fill = TRUE)
  out[, stratum := sn]
  out
}), fill = TRUE)
setcolorder(res, c("stratum", "predictor", "beta", "std_error", "z_value", "p_value",
                   "ci_low", "ci_high", "odds_ratio", "n", "n_case"))
fwrite(res, file.path(OUT_DIR, "gnpc_ad_assoc_by_age.csv"))

cat("\n=== GNPC predicted-phenotype -> AD associations by age stratum ===\n")
cat("(full = all ages, the Fig 4g y-axis model; lt65/ge65 = age <65 / >=65)\n")
print(res[, .(n_pheno = .N, n_sig = sum(p_value < 0.05, na.rm = TRUE),
              median_n = as.integer(median(n)), median_ncase = as.integer(median(n_case))), by = stratum])
cat("\n=== exported files ===\n"); print(list.files(OUT_DIR))
cat("\nDONE.\n"); writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
sink(); close(.logf)
