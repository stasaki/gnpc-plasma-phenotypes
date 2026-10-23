# ==============================================================================
# Fig. 6d: SHAP and per-domain predicted AD probability of one typical
# individual per risk profile, drawn over the profile's distribution (violins).
# Inputs : ProteomicsV1_3ms.rds; predicted_phenotypes.rds (target_annotation);
#          classifier runs p8_one_contributoer_out_nesting_<EX_TAG> in
#          step_dir("04_classifier_runs") (01_train_classifiers.R);
#          SE_res_ad_nci.rds from step_dir("04_risk_profiles") (05_risk_profiles.R).
# Outputs: step_dir("04_fig6d")
# Runs on the AD Workbench. From the repo root:
#   Rscript 04_ad_classifier/06_fig6d_individuals.R
# ==============================================================================

# ==============================================================================
# Fig 6d — violin + individual point
# ==============================================================================
# For each feature category the figure shows the SHAP DISTRIBUTION across every
# member of a risk profile as a violin, with the representative individual's
# phi_i laid over it as a single large point.
#
#   left  panel : SHAP (log-odds)            violin = profile, point = the case
#                                            point red/blue split at 0
#   right panel : predicted AD probability   violin = profile, point = the case
#                                            point red/blue split at 0.5
#
# The representative individual of each profile is the member with the highest
# total cosine similarity to the other members (see REP_METHOD below).
#
# Figure size: width = 2.36*4, height = 3.12*4  ->  9.44 x 12.48 inches
#
# Outputs:
#   clean_composite_violin.pdf / .png     Fig 6d
#   fig6d_violin_density.csv              violin density curves (AGGREGATED)
#   fig6d_violin_summary.csv              n / mean / sd / median / IQR / range
#   fig6d_case_values.csv                 the 5 individuals' phi_i and probs
#   fig6d_cases.csv                       profile labels + n_profile
#   fig6d_rep_selection_comparison.csv    the three selection rules compared
#
# NOTE ON DISCLOSURE: nothing member-level is written. The violins are exported
# as density curves and summary statistics only, so the figure can be restyled
# without per-participant rows leaving the workbench.
#
# No ggdist / ggbeeswarm dependency — violins are hand-built polygons.
# ==============================================================================

source("config.R")
pp <- read_predicted_phenotypes(PREDICTED_RDS)

library(tidyverse)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(proxy)
library(patchwork)
library(SummarizedExperiment)

OUT_DIR <- step_dir("04_fig6d")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FIG_W <- 2.36 * 4     # 9.44 in
FIG_H <- 3.12 * 4     # 12.48 in

# ==============================================================================
# ---- SHAP LOADING ------------------------------------------------------------
# ==============================================================================
target_annotations = get_table(pp, "target_annotation")

gen_lookup = function(code_str){
  parts <- strsplit(code_str%>%gsub("\\.$","",.), ";\\s*")[[1]]
  lookup <- setNames(sub(".*=\\s*", "", parts), sub("\\s*=.*", "", parts))
  return(lookup)
}
target_annotations%>%
  rowwise()%>%
  mutate(class_lookup = gen_lookup(Class_label)%>%list()) -> target_annotations
indx = which(target_annotations$Target=="apoe_genotype")
names(target_annotations$class_lookup[[indx]]) = as.character(1:6)

# per-category predicted AD probability (the right-hand panel)
list.files(file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/")))%>%
  .[grepl("0.01",.)]%>%
  .[!grepl(".all_1st_tier|\\.Clinical.Diagnosis",.)]%>%
  lapply(., function(x){
    base_path <- paste0(file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/")),x,"/ml/p1_trait")
    submodels_path <- file.path(base_path, "submodels")
    sample_id <- readRDS(file.path(base_path, "sample_id.rds"))
    prediction_files <- list.files(submodels_path, pattern = "predictions.rds", recursive = TRUE, full.names = TRUE)
    predictions <- prediction_files %>% map(readRDS) %>% bind_rows()%>% filter(split=="test")
    predictions$sample_id <- sample_id[as.numeric(predictions$sample_id)]
    predictions$model = gsub(".+feature..","",x)%>%gsub("__.+","",.)
    predictions%>%return()
  })%>%
  bind_rows() -> predictions_all

# the AD classifier (ex Clinical Diagnosis) — source of the SHAP shown in Fig 6d
base_path <- file.path(step_dir("04_classifier_runs"), paste0("p8_one_contributoer_out_nesting_", EX_TAG, "/trait..ad__feature..ex.Clinical.Diagnosis__metric_th..0.01/ml/p1_trait"))
submodels_path <- file.path(base_path, "submodels")

sample_id <- readRDS(file.path(base_path, "sample_id.rds"))

shap_files <- list.files(submodels_path, pattern = "shap_values.rds", recursive = TRUE, full.names = TRUE)
shap_values <- shap_files %>% map(readRDS) %>% bind_rows()

prediction_files <- list.files(submodels_path, pattern = "predictions.rds", recursive = TRUE, full.names = TRUE)
predictions <- prediction_files %>% map(readRDS) %>% bind_rows()%>% filter(split=="test")
rownames(predictions) <- sample_id[as.numeric(predictions$sample_id)]

shap_values = shap_values%>%
  inner_join(.,predictions%>%dplyr::select(predicted_prob,sample_id),by="sample_id")
rownames(shap_values) <- sample_id[as.numeric(shap_values$sample_id)]

shap_values = shap_values%>%
  dplyr::select(-sample_id)%>%
  mutate(sample_id=rownames(.))%>%
  gather(feature,shap_value,-sample_id,-actual,-predicted_prob)

shap_values = shap_values%>%
  group_by(feature)%>%
  nest()%>%
  mutate(target_id_no_class = gsub("_class.+","",feature))%>%
  inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
  filter(grepl("class",feature))%>%
  mutate(class_code = gsub(".+class_","",feature)%>%as.numeric())%>%
  group_by(target_id_no_class)%>%
  mutate(N_class = dplyr::n())%>%
  rowwise()%>%
  mutate(class_code = ifelse(N_class>2,class_code+1,class_code),
         class_code=as.character(class_code))%>%
  dplyr::select(-N_class)%>%
  rowwise()%>%
  mutate(class_label = class_lookup%>% unlist() %>% .[class_code])%>%
  bind_rows(.,
            shap_values%>%
              group_by(feature)%>%
              nest()%>%
              mutate(target_id_no_class = gsub("_class.+","",feature))%>%
              inner_join(.,target_annotations,by=c("target_id_no_class"="Target ID"))%>%
              filter(!grepl("class",feature)))%>%
  mutate(class_label=replace(class_label,class_label%in%c("Yes","Limbic-type","Present"),NA))%>%
  mutate(label = replace(label,!is.na(class_label),paste0(label," - ", class_label)[!is.na(class_label)]))%>%
  dplyr::select(-class_lookup)%>%
  unnest(cols = c(data))

ProteomicsV1_3ms <- readRDS(PROTEOMICS_RDS)

# ==============================================================================
# ---- PROFILES + representative individuals ----------------------------------
# ==============================================================================
SE_res_ad_nci <- readRDS(file.path(step_dir("04_risk_profiles"), "SE_res_ad_nci.rds"))

profile_levels <- c("ad1","ad2","ad3","nci2","nci4")
profile_pretty <- c(ad1="AD1", ad2="AD2", ad3="AD3", nci2="Non-AD1", nci4="Non-AD2")

profile_members <- lapply(profile_levels, function(p){
  data.frame(profile = p,
             sample_id = SE_res_ad_nci %>% filter(cluster1 == p) %>% pull(IID),
             stringsAsFactors = FALSE)
}) %>% bind_rows()

# ---- HOW THE REPRESENTATIVE IS CHOSEN ---------------------------------------
# "cosine"    = the rule used for Fig 6d: max total cosine similarity to the
#               other profile members. Cosine is SCALE-INVARIANT, so it scores
#               the DIRECTION of the SHAP vector (the relative mix across
#               categories) and is blind to its MAGNITUDE. It picks the individual
#               whose contribution pattern is most typical of the profile; that
#               individual can still have extreme values.
# "euclidean" = medoid under Euclidean distance: central in direction AND
#               magnitude. May pick different individuals.
# "median"    = member closest to the per-category median (robust variant of
#               "euclidean"; may also pick different individuals).
# The run-time diagnostic below reports where each plotted individual sits
# inside its own profile.
REP_METHOD <- "cosine"      # <- change to "euclidean" / "median" to switch

pick_representative <- function(M, method) {
  switch(method,
    cosine    = names(sort(colSums(1 - as.matrix(proxy::dist(M, method = "cosine")))) %>% tail(1)),
    euclidean = names(sort(colSums(as.matrix(stats::dist(M))))[1]),
    median    = { med <- apply(M, 2, median)
                  rownames(M)[which.min(apply(M, 1, function(r) sum(abs(r - med))))] },
    stop("unknown REP_METHOD: ", method))
}

cat_matrix <- function(ids) {
  shap_values %>%
    filter(sample_id %in% ids) %>%
    group_by(sample_id, feature = Category) %>%
    summarise(shap_value = sum(shap_value), .groups = "drop") %>%
    spread(feature, shap_value) %>%
    as.data.frame() %>%
    column_to_rownames("sample_id")
}

reps <- lapply(profile_levels, function(subtype_label){
  keep <- profile_members %>% filter(profile == subtype_label) %>% pull(sample_id)
  M <- cat_matrix(keep)
  data.frame(profile = subtype_label, profile_label = profile_pretty[[subtype_label]],
             sample_id = pick_representative(M, REP_METHOD),
             n_profile = length(keep), stringsAsFactors = FALSE)
}) %>% bind_rows()

rep_ids <- reps$sample_id   # internal use only — never written out

# ---- the two quantities being plotted, for EVERY profile member --------------
# left panel: category-level SHAP
shap_cat_members <- shap_values %>%
  inner_join(profile_members, by = "sample_id") %>%
  group_by(profile, sample_id, Category) %>%
  summarise(value = sum(shap_value), .groups = "drop") %>%
  mutate(panel = "shap")

# right panel: per-category predicted AD probability
prob_members <- predictions_all %>%
  inner_join(profile_members, by = "sample_id") %>%
  mutate(Category = gsub("\\.", " ", model)) %>%
  dplyr::select(profile, sample_id, Category, value = predicted_prob) %>%
  mutate(panel = "prob")

members_long <- bind_rows(shap_cat_members, prob_members)

case_values <- members_long %>%
  inner_join(reps %>% dplyr::select(profile, profile_label, sample_id), by = c("profile","sample_id"))

# ==============================================================================
# ---- PLOT CONSTANTS ---------------------------------------------------------
# ==============================================================================
category_colors <- c(
  "Pathology" = "#E69F00", "Cognition" = "#009E73", "Motor and Gait" = "#F0E442",
  "Medical Conditions" = "#0072B2", "Medications" = "#D55E00",
  "Genetics" = "#CC79A7", "Demographics" = "#999999"
)
CAT_LEVELS  <- rev(names(category_colors))   # Pathology ends up on top
POS_COL     <- "#B83A3F"
NEG_COL     <- "#4078BD"
VIOLIN_FILL <- "grey86"
THRESH      <- 0.5                           # AD decision threshold

stopifnot(all(unique(members_long$Category) %in% CAT_LEVELS))

# ==============================================================================
# ---- HELPERS — hand-built violins (no ggdist / ggbeeswarm on the workbench) --
# ==============================================================================
# violin_density() returns the density curve per group. This table IS the export:
# it is fully aggregated, and violin_poly() rebuilds the identical polygon from
# it, so the figure can be restyled locally without member-level rows.
violin_density <- function(df, levels = CAT_LEVELS, n_grid = 512) {
  df %>%
    mutate(.g = factor(Category, levels = levels)) %>%
    group_by(.g) %>%
    group_modify(~{
      x <- .x$value
      if (length(unique(x)) < 3) return(tibble(x = numeric(), dens = numeric()))
      d <- stats::density(x, bw = "nrd0", n = n_grid, from = min(x), to = max(x))
      tibble(x = d$x, dens = d$y)
    }) %>% ungroup() %>%
    group_by(.g) %>% mutate(w = dens / max(dens)) %>% ungroup() %>%   # scale = "width"
    mutate(ypos = as.numeric(.g))
}

violin_poly <- function(vd, half_width = 0.40) {
  vd %>% group_by(.g) %>% group_modify(~bind_rows(
    tibble(x = .x$x,      y = .x$ypos + .x$w * half_width),
    tibble(x = rev(.x$x), y = rev(.x$ypos - .x$w * half_width)))) %>% ungroup()
}

sign_split <- function(v, at = 0) ifelse(v > at, "positive", "negative")

# A thin tick at the profile MEDIAN. Not a box — one line per category. It
# separates "where this profile sits" from "where this individual sits", so
# profile-level patterns are read off the distribution rather than off the
# single plotted point. Set FALSE to drop it.
SHOW_MEDIAN_TICK <- TRUE

median_tick <- function(memb) {
  if (!SHOW_MEDIAN_TICK) return(NULL)
  m <- memb %>%
    mutate(.g = factor(Category, levels = CAT_LEVELS)) %>%
    group_by(.g) %>% summarise(m = median(value), .groups = "drop") %>%
    mutate(ypos = as.numeric(.g))
  geom_segment(data = m, aes(x = m, xend = m, y = ypos - 0.28, yend = ypos + 0.28),
               color = "grey38", linewidth = 0.7, inherit.aes = FALSE)
}

# Panels have no titles, subtitles or caption (labels are added outside R) and
# an x-axis on every row.

# Fixed x-range. coord_cartesian only clips the display; the run-time check below
# reports whether any violin tail actually falls outside.
SHAP_LIM    <- c(-2, 2)
SHAP_BREAKS <- seq(-1.5, 1.5, 0.5)

shap_panel <- function(p) {
  memb <- members_long %>% filter(profile == p, panel == "shap")
  cse  <- case_values  %>% filter(profile == p, panel == "shap") %>%
    mutate(.g = factor(Category, levels = CAT_LEVELS), ypos = as.numeric(.g),
           sign = sign_split(value, 0))
  poly <- violin_poly(violin_density(memb))
  ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey35") +
    geom_polygon(data = poly, aes(x, y, group = .g), fill = VIOLIN_FILL, color = NA) +
    median_tick(memb) +
    geom_point(data = cse, aes(value, ypos, fill = sign),
               shape = 21, size = 4.8, stroke = 1.1, color = "white") +
    scale_fill_manual(values = c(positive = POS_COL, negative = NEG_COL)) +
    scale_y_continuous(breaks = seq_along(CAT_LEVELS), labels = CAT_LEVELS,
                       limits = c(0.4, length(CAT_LEVELS) + 0.6), expand = c(0, 0)) +
    scale_x_continuous(breaks = SHAP_BREAKS) +
    coord_cartesian(xlim = SHAP_LIM) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title   = element_blank(),
          axis.text.y  = element_text(size = 14, color = "black"),
          axis.text.x  = element_text(size = 12, color = "black"),
          axis.line.y  = element_blank(),
          axis.ticks.y = element_blank())
}

risk_panel <- function(p) {
  memb <- members_long %>% filter(profile == p, panel == "prob")
  cse  <- case_values  %>% filter(profile == p, panel == "prob") %>%
    mutate(.g = factor(Category, levels = CAT_LEVELS), ypos = as.numeric(.g),
           sign = sign_split(value, THRESH))
  poly <- violin_poly(violin_density(memb))
  ggplot() +
    geom_polygon(data = poly, aes(x, y, group = .g), fill = VIOLIN_FILL, color = NA) +
    geom_vline(xintercept = THRESH, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    median_tick(memb) +
    geom_point(data = cse, aes(value, ypos, fill = sign),
               shape = 21, size = 4.8, stroke = 1.1, color = "white") +
    scale_fill_manual(values = c(positive = POS_COL, negative = NEG_COL)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq_along(CAT_LEVELS), labels = NULL,
                       limits = c(0.4, length(CAT_LEVELS) + 0.6), expand = c(0, 0)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          axis.title         = element_blank(),
          panel.grid.major.y = element_line(color = "grey"),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y        = element_blank(),
          axis.text.x        = element_text(size = 12, color = "black"),
          axis.ticks.y       = element_blank(),
          plot.margin        = margin(10, 15, 10, 10))
}

# ==============================================================================
# ---- ASSEMBLE ---------------------------------------------------------------
# ==============================================================================
demog_cols <- c("ad","age_at_visit","APOE","sex","race","years_of_education","bmi",
                "systolic_blood_pressure_sitting","diastolic_blood_pressure_sitting")
demog <- colData(ProteomicsV1_3ms)[colData(ProteomicsV1_3ms)$sample_id %in% rep_ids, ] %>%
  as.data.frame() %>% dplyr::select(sample_id, any_of(demog_cols))

# overall AD probability for each case. NOTE: take it from shap_values, whose
# sample_id is the UUID — `predictions$sample_id` is still the numeric index at
# this point (the UUIDs were assigned to rownames only), so joining on it fails.
overall_prob <- shap_values %>% ungroup() %>% distinct(sample_id, predicted_prob)

case_meta <- reps %>%
  left_join(demog, by = "sample_id") %>%
  left_join(overall_prob, by = "sample_id") %>%
  mutate(subtitle = sprintf("%s yr, %s, APOE %s, predicted %.2f",
                            age_at_visit, ifelse(sex == 1, "Male", ifelse(sex == 2, "Female", "sex unknown")),
                            APOE, predicted_prob))

# one row per profile: (SHAP panel + probability panel), widths c(1.5, 1),
# stacked with wrap_plots(ncol = 1). No annotation layer.
plot_list <- lapply(profile_levels, function(p) {
  (shap_panel(p) + risk_panel(p)) + plot_layout(widths = c(1.5, 1))
})

g <- wrap_plots(plot_list, ncol = 1)

ggsave(plot = g, filename = file.path(OUT_DIR, "clean_composite_violin.png"),
       width = FIG_W, height = FIG_H, bg = "white")
ggsave(plot = g, filename = file.path(OUT_DIR, "clean_composite_violin.pdf"),
       width = FIG_W, height = FIG_H, bg = "white")

# ---- panel label text (printed, never drawn) --------------------------------
cat("\n---- panel labels, not drawn in the figure (top row -> bottom row) ----\n")
case_meta %>%
  dplyr::select(profile_label, n_profile, subtitle) %>%
  arrange(match(profile_label, profile_pretty[profile_levels])) %>%
  as.data.frame() %>% print(row.names = FALSE)

# ==============================================================================
# ---- EXPORTS (aggregated only — no member-level rows) -----------------------
# ==============================================================================
# (1) violin density curves — enough to rebuild the identical violin locally
density_export <- members_long %>%
  group_by(profile, panel) %>%
  group_modify(~violin_density(.x)) %>%
  ungroup() %>%
  left_join(reps %>% dplyr::select(profile, profile_label), by = "profile") %>%
  dplyr::select(profile, profile_label, panel, Category = .g, x, dens, w, ypos)
fwrite(density_export, file.path(OUT_DIR, "fig6d_violin_density.csv"))

# (2) distribution summary of each violin
summary_export <- members_long %>%
  group_by(profile, panel, Category) %>%
  summarise(n = dplyr::n(), mean = mean(value), sd = sd(value),
            min = min(value), q25 = quantile(value, .25), median = median(value),
            q75 = quantile(value, .75), max = max(value), .groups = "drop") %>%
  left_join(reps %>% dplyr::select(profile, profile_label), by = "profile") %>%
  relocate(profile_label, .after = profile)
fwrite(summary_export, file.path(OUT_DIR, "fig6d_violin_summary.csv"))

# (3) the plotted individuals — phi_i and per-category probability, no sample_id
case_export <- case_values %>%
  dplyr::select(profile, profile_label, panel, Category, value) %>%
  arrange(profile, panel, Category)
fwrite(case_export, file.path(OUT_DIR, "fig6d_case_values.csv"))

fwrite(reps %>% dplyr::select(profile, profile_label, n_profile),
       file.path(OUT_DIR, "fig6d_cases.csv"))

# safety check: no sample_id column in anything written
stopifnot(!"sample_id" %in% names(density_export),
          !"sample_id" %in% names(summary_export),
          !"sample_id" %in% names(case_export))

cat("\nProfile sizes:\n"); print(reps %>% dplyr::select(profile_label, n_profile))

# ---- DIAGNOSTIC: is the "representative" actually central? -------------------
# Percentile of the plotted individual within its own profile, per category.
# 50 = sits mid-violin; 95+ = sits in the tail. Also compares the three
# selection rules side by side.
pctl <- case_values %>%
  filter(panel == "shap") %>%
  dplyr::select(profile, profile_label, Category, case_value = value) %>%
  left_join(shap_cat_members %>% dplyr::select(profile, Category, value),
            by = c("profile","Category"), relationship = "many-to-many") %>%
  group_by(profile, profile_label, Category) %>%
  # NB: case_value[1], not first() — SummarizedExperiment is loaded after
  # tidyverse, so BiocGenerics::first masks dplyr::first and has no method for
  # a plain numeric vector.
  summarise(percentile = round(100 * mean(value <= case_value[1])), .groups = "drop") %>%
  spread(Category, percentile)
cat(sprintf("\nPercentile of the plotted individual within its profile (REP_METHOD = '%s'):\n",
            REP_METHOD))
print(as.data.frame(pctl %>% dplyr::select(-profile)), row.names = FALSE)

cat("\nSelection-rule comparison. median_pctl = median percentile across the 7\n",
    "categories (50 = middle of the profile); magnitude_pctl = percentile of the SHAP vector\n",
    "length (how extreme the individual is overall).\n\n")
rule_cmp <- lapply(c("cosine","euclidean","median"), function(mth) {
  lapply(profile_levels, function(p) {
    keep <- profile_members %>% filter(profile == p) %>% pull(sample_id)
    M    <- cat_matrix(keep)
    id   <- pick_representative(M, mth)
    data.frame(method = mth, profile_label = profile_pretty[[p]],
               median_pctl = round(median(sapply(colnames(M), function(k)
                 100 * mean(M[, k] <= M[id, k])))),
               magnitude_pctl = round(100 * mean(
                 sqrt(rowSums(M^2)) <= sqrt(sum(M[id, ]^2)))),
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
}) %>% bind_rows()
print(as.data.frame(rule_cmp), row.names = FALSE)
fwrite(rule_cmp, file.path(OUT_DIR, "fig6d_rep_selection_comparison.csv"))

# Report how much of the SHAP distribution the fixed x-range c(-2, 2) clips.
clipped <- shap_cat_members %>%
  summarise(n = dplyr::n(),
            n_outside = sum(value < SHAP_LIM[1] | value > SHAP_LIM[2]),
            observed_min = min(value), observed_max = max(value))
cat(sprintf("\nSHAP x-range fixed at [%.1f, %.1f].\n",
            SHAP_LIM[1], SHAP_LIM[2]))
cat(sprintf("  observed range: %.2f to %.2f\n", clipped$observed_min, clipped$observed_max))
cat(sprintf("  %d of %d category-level values (%.2f%%) fall outside and are clipped.\n",
            clipped$n_outside, clipped$n, 100 * clipped$n_outside / clipped$n))
if (clipped$n_outside > 0)
  cat("  -> widen SHAP_LIM if the violin tails look cut off.\n")
cat("\nWrote to", OUT_DIR, ":\n",
    " clean_composite_violin.pdf / .png\n",
    " fig6d_violin_density.csv\n  fig6d_violin_summary.csv\n",
    " fig6d_case_values.csv\n  fig6d_cases.csv\n")
