# 03_trajectories/03_table_s8.R
# Supplementary Table 8 (GAM statistics: primary, rank-INT and age >= 65 fits).
# Input:   step_dir("03_age_trajectory")/gam_stats_all_variants.csv, written by
#          03_trajectories/01_age_trajectory.R (REPLOT_ONLY <- FALSE)
# Outputs: step_dir("03_table_s8") -- Table_S8_GAM_stats.csv and .xlsx
#          (or a _caption.txt if openxlsx is not installed).
# Runs on the AD Workbench. From the repository root:
#   Rscript 03_trajectories/03_table_s8.R [path/to/gam_stats_all_variants.csv]

# ---------------------------------------------------------------
# Assembles Supplementary Table S8 (GAM statistics) from ONE input:
#
#   gam_stats_all_variants.csv   (03_trajectories/01_age_trajectory.R)
#
# That file carries every number the table needs, from a single model-fitting
# run: variant = full / rank / ge65 crossed with model = base / apoe / sex.
# Taking all three from one run is the point -- the primary statistics and the
# two re-analyses of the scale floor (compressed young end) are then directly comparable, and a test
# statistic never sits next to a P value from a different fit.
#
#   full   primary analysis: min-max scaled outcome, all ages
#   rank   rank-INT: outcome replaced by inverse-normal rank scores over the
#          full age range, so the test sees the outcome only through its order
#          and a monotone compression of the young end is invisible (n unchanged)
#   ge65   primary outcome, chronological age >= 65 only
# ---------------------------------------------------------------

source("config.R")

library(dplyr)
library(readr)
library(tidyr)

# --- Paths ---
args <- commandArgs(trailingOnly = TRUE)
default_csv <- file.path(step_dir("03_age_trajectory"),
                         "gam_stats_all_variants.csv")
variants_csv <- if (length(args) >= 1) args[1] else default_csv
excel_path   <- file.path(step_dir("03_table_s8"), "Table_S8_GAM_stats.xlsx")

if (!file.exists(variants_csv)) {
  stop("gam_stats_all_variants.csv not found:\n  ", variants_csv, "\n",
       "Produce it by running 03_trajectories/01_age_trajectory.R with\n",
       "REPLOT_ONLY <- FALSE, or pass the path as the first argument.")
}

# --- Which term carries the effect of interest, per model ---
#   base: s(age_at_visit)                   — age smooth (Fig 5a-b)
#   apoe: ti(age_at_visit,APOE4_count_pred) — age x APOE interaction (Fig 5c, Extended Data Fig. 5)
#   sex:  s(age_at_visit):sex_ordered1      — age x sex interaction (Fig 5d, Extended Data Fig. 6)
MAIN_TERM <- c(
  base = "s(age_at_visit)",
  apoe = "ti(age_at_visit,APOE4_count_pred)",
  sex  = "s(age_at_visit):sex_ordered1"
)
MODEL_NAME <- c(base = "Age trajectory",
                apoe = "APOE × Age",
                sex  = "Sex × Age")

message("Reading ", variants_csv)
v <- read_csv(variants_csv, show_col_types = FALSE)

need <- c("variant", "model", "plot_label", "term",
          "n_obs", "r_squared", "dev_explained", "edf", "ref.df", "statistic", "p.value")
missing <- setdiff(need, names(v))
if (length(missing)) {
  stop("gam_stats_all_variants.csv is missing columns: ", paste(missing, collapse = ", "))
}

v <- v %>%
  # phenotypes shown in Extended Data Fig. 5-6: all except Age, the AD PRS and the cognitive slopes
  filter(!plot_label %in% c("Age", "ADPRS ensemble (averaging)"),
         !grepl("slope", plot_label)) %>%
  filter(term == MAIN_TERM[model]) %>%
  mutate(Model = factor(MODEL_NAME[model], levels = MODEL_NAME)) %>%
  dplyr::select(variant, Model, Phenotype = plot_label,
                n_obs, r_squared, dev_explained, edf, ref.df, statistic, p.value)

have <- sort(unique(v$variant))
message(sprintf("  variants present: %s | models: %d | phenotypes: %d",
                paste(have, collapse = ", "),
                dplyr::n_distinct(v$Model), dplyr::n_distinct(v$Phenotype)))
if (!all(c("full", "rank", "ge65") %in% have)) {
  stop("expected variants full/rank/ge65, found: ", paste(have, collapse = ", "))
}

# --- One row per Model x Phenotype, one column block per variant --------------
# Every reported P keeps its own EDF, reference df and test statistic: a P value
# without the statistic it came from is not reportable.
# R2 and deviance explained are given for the primary fit only. They are
# model-level, and under rank-INT they describe a different outcome scale, so
# repeating them per variant would invite a comparison that means nothing.
wide <- v %>%
  pivot_wider(
    id_cols     = c(Model, Phenotype),
    names_from  = variant,
    values_from = c(n_obs, r_squared, dev_explained, edf, ref.df, statistic, p.value)
  )

table_s8 <- wide %>%
  dplyr::select(
    Model, Phenotype,
    n                            = n_obs_full,
    R2                           = r_squared_full,
    `Deviance Explained (%)`     = dev_explained_full,
    EDF                          = edf_full,
    `Ref. DF`                    = ref.df_full,
    `Test Statistic`             = statistic_full,
    `P-value`                    = p.value_full,
    `EDF (rank-INT)`             = edf_rank,
    `Ref. DF (rank-INT)`         = ref.df_rank,
    `Test Statistic (rank-INT)`  = statistic_rank,
    `P-value (rank-INT)`         = p.value_rank,
    `n (age >= 65)`              = n_obs_ge65,
    `EDF (age >= 65)`            = edf_ge65,
    `Ref. DF (age >= 65)`        = ref.df_ge65,
    `Test Statistic (age >= 65)` = statistic_ge65,
    `P-value (age >= 65)`        = p.value_ge65
  ) %>%
  arrange(Model) %>%
  mutate(Model = as.character(Model))

stopifnot(!anyNA(table_s8$`P-value`),
          !anyNA(table_s8$`P-value (rank-INT)`),
          !anyNA(table_s8$`P-value (age >= 65)`))

# mgcv can return a very small NEGATIVE P for a highly significant smooth -- a
# floating-point artefact of the test statistic, not a value. Clamp to 0; these
# are all P values that underflowed, and the table reports that floor as 0
# elsewhere.
p_cols <- grep("^P-value", names(table_s8), value = TRUE)
n_clamped <- sum(sapply(table_s8[p_cols], function(x) sum(x < 0)))
if (n_clamped > 0) {
  message(sprintf("  clamped %d negative P value(s) to 0 (mgcv underflow)", n_clamped))
  table_s8[p_cols] <- lapply(table_s8[p_cols], function(x) pmax(x, 0))
}

# summary(gam)$dev.expl is a PROPORTION, but the column is headed "(%)".
# Convert, so the number matches its own header.
table_s8$`Deviance Explained (%)` <- table_s8$`Deviance Explained (%)` * 100

message(sprintf("Table S8: %d rows x %d cols", nrow(table_s8), ncol(table_s8)))

# ">= 65" is kept ASCII in the code; the delivered header uses the proper symbol.
names(table_s8) <- gsub(">=", "≥", names(table_s8), fixed = TRUE)
p_cols <- grep("^P-value", names(table_s8), value = TRUE)   # names changed above

title <- paste0(
  "Supplementary Table 8. Generalized additive model (GAM) statistics ",
  "for pseudotemporal trajectory analyses of predicted phenotypes in GNPC. ",
  "Models fitted with mgcv::gam(method=\"REML\"). ",
  "Site-specific contributor effects are omitted for brevity. ",
  "All tests are two-sided. ",
  "Alongside the primary fit, the table reports two re-analyses that test ",
  "whether a result depends on the compressed young end of the predicted-value ",
  "scale: rank-INT, in which the outcome is replaced by its inverse-normal rank ",
  "scores over the full age range, so the test sees the outcome only through its ",
  "order and any monotone compression is invisible (n unchanged); and age ≥ 65, ",
  "in which the primary outcome is re-fitted on chronological age ≥ 65 only, ",
  "which drops the compressed young participants but also costs power and ",
  "shortens the age span over which an age-slope interaction is estimated. ",
  "All three use the same model specification and come from the same fitting ",
  "run, so they are directly comparable. ",
  "Values are given to three significant figures. P values for smooth terms are ",
  "as computed by mgcv; where mgcv returned 0 (its test is not resolved below ",
  "about 2e-6), P is reported as < 3e-6."
)

# CSV twin: the easier form to diff between runs. It stays fully numeric -- it
# is the machine-readable form.
csv_path <- sub("\\.xlsx$", ".csv", excel_path)
write_csv(table_s8, csv_path)
message("Saved to ", csv_path)

# The Excel copy is the one a reader opens, so P values are formatted there:
# 3 significant figures, and a P of 0 reported as "< 3e-6". mgcv computes
# smooth-term P with psum.chisq() (Davies' algorithm, absolute tol = 2e-05) and
# returns 0 (occasionally a tiny negative, clamped above) once it cannot resolve
# the tail -- in testing, for true P as large as ~2.2e-6. "0", or "< 2.2e-16",
# would be a claim the fit cannot support. This turns the P columns into text,
# which is why it is not done to the CSV.
P_ZERO_BOUND <- "< 3e-6"
fmt_p <- function(p) ifelse(p <= 0, P_ZERO_BOUND, formatC(p, format = "g", digits = 3))
n_floored <- sum(sapply(table_s8[p_cols], function(x) sum(x <= 0)))
message(sprintf("  %d P value(s) of 0 -> \"%s\"", n_floored, P_ZERO_BOUND))

table_x <- table_s8
table_x[p_cols] <- lapply(table_x[p_cols], fmt_p)

# Everything else to 3 significant figures, matching the P values. These stay
# NUMERIC so the sheet is still sortable -- only the P columns have to be text,
# because of the "< 3e-6" floor. The sample counts are left alone: they are
# exact integers, and signif() would turn 13793 into 13800.
num_cols <- setdiff(names(table_x)[sapply(table_x, is.numeric)],
                    grep("^n( |$)", names(table_x), value = TRUE))
table_x[num_cols] <- lapply(table_x[num_cols], signif, digits = 3)

# openxlsx is not guaranteed on the workbench; the CSV above is the fallback.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Table S8")
  openxlsx::writeData(wb, "Table S8", x = title, startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::writeData(wb, "Table S8", x = table_x, startRow = 2, colNames = TRUE)
  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
  message("Saved to ", excel_path)
} else {
  writeLines(title, sub("\\.xlsx$", "_caption.txt", excel_path))
  message("openxlsx not installed; wrote the CSV and the caption only.")
}
