# Plasma protein-based prediction of Alzheimer's disease phenotypes: analysis code

Analysis code for:

> Tasaki S, de Paiva Lopes K, Vialle RA, Ng B, Gaiteri C, The Global Neurodegeneration Proteomics Consortium (GNPC), Bennett DA. *Fusing depth and breadth in plasma proteomics to scale deep Alzheimer's phenotyping to global cohorts.*

## Data

GNPC data: <https://discover.alzheimersdata.org/catalogue/datasets/e2f3536b-d97b-4303-89bd-6864200807a4>

The data are not included here. The analyses were run on the AD Workbench of the AD Data Initiative.

## How to run

Put this repository in `<project>/code/`, the inputs in `<project>/data/`, and run from the repository root:

``` bash
bash run_all.sh            # all steps, incl. classifier training (needs a GPU)
TRAIN=0 bash run_all.sh    # skip training; use runs already in ../output/04_classifier_runs
```

Results go to `<project>/output/`, one subfolder per step, with logs in `output/logs/`. Settings are in `config.R`. Each script starts with a header that lists its inputs and outputs.

## Analysis steps

| Script | Paper items |
|------------------------------------|------------------------------------|
| `01_validation/01_validation.Rmd` | Fig. 4a–e; Extended Data Fig. 2–3 |
| `02_ad_association/01_ad_assoc_by_site.Rmd` | Fig. 4f–g (GNPC side); Extended Data Fig. 4; Supplementary Table 7 |
| `02_ad_association/02_ad_assoc_site_adjusted.Rmd` | Fig. 4f–g (GNPC side); Extended Data Fig. 4; Supplementary Table 7 |
| `02_ad_association/03_apoe_adjusted_assoc.R` | Supplementary Table 7 (APOE-, age- and sex-adjusted columns) |
| `02_ad_association/04_ad_assoc_by_age.R` | Supplementary Table 7 (age \< 65 / ≥ 65 columns) |
| `03_trajectories/01_age_trajectory.R` | Fig. 5a–d; Extended Data Fig. 5–6 |
| `03_trajectories/02_age_trajectory_protein.R` | Extended Data Fig. 7 |
| `03_trajectories/03_table_s8.R` | Supplementary Table 8 |
| `04_ad_classifier/01_train_classifiers.R` | Trains the AD classifiers (CatBoost, leave-one-contributor-out) |
| `04_ad_classifier/02_collect_results.R` | Collects the cross-validation results |
| `04_ad_classifier/03_compare_model_perf.R` | Fig. 6a; Extended Data Fig. 8; Supplementary Table 10 |
| `04_ad_classifier/04_variable_contributions.R` | Fig. 6b |
| `04_ad_classifier/05_risk_profiles.R` | Fig. 6c |
| `04_ad_classifier/06_fig6d_individuals.R` | Fig. 6d |
| `04_ad_classifier/07_predict_cdr_progression.R` | Extended Data Fig. 10; Supplementary Table 12 |
| `05_ad_specificity/01_pdftd_specificity.R` | Extended Data Fig. 9b (group values) |
| `05_ad_specificity/02_ad_classifier_oof_specificity.R` | Extended Data Fig. 9a |
| `05_ad_specificity/03_demographic_accuracy_persite.R` | Extended Data Fig. 9b (per-site points) |

`R/vendor/` holds helper code copied from our internal packages (CatBoost cross-validation, SHAP, and `run_pipeline()`), so no private package is needed.

## Software

- R 4.4.1
- CatBoost (R) 0.25.1: AD classification (GPU)
- mgcv 1.9.1: generalized additive models for the age trajectories
- lme4 1.1.36 / lmerTest 3.1.3: linear mixed-effects models for CDR progression
- umap 0.2.10.0: UMAP of SHAP values
- speakeasyR 0.1.7: clustering of risk attribution profiles

The phenotype prediction models (Python, TabPFN) are at <https://github.com/stasaki/tabpfn-pipeline>.

## License

MIT; see [LICENSE](LICENSE).

## Contact

Shinya Tasaki, Rush Alzheimer's Disease Center, Rush University Medical Center, Chicago, IL, USA.
