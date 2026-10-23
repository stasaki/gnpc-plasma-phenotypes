# Collect the cross-validation results of all AD classifier runs (input to Fig. 6a,
# Extended Data Fig. 8 and Supplementary Table 10; see 03_compare_model_perf.R).
# Runs on the AD Workbench.
#   Input : step_dir("04_classifier_runs") (made by 04_ad_classifier/01_train_classifiers.R)
#   Output: step_dir("04_collect_result"): summary_stats.{rds,txt}, fold_metrics.{rds,txt}
# Run from the repository root:
#   Rscript 04_ad_classifier/02_collect_results.R

source("config.R")
library(tidyverse)
library(data.table)
outdir = step_dir("04_collect_result")

evaluation_summary_files = list.files(path=step_dir("04_classifier_runs"),pattern = "evaluation_summary.rds",recursive = T,full.names = T)
evaluation_summary_files = evaluation_summary_files[grepl(EX_TAG,evaluation_summary_files)]
evaluation_summary_files = evaluation_summary_files[!grepl("p11",evaluation_summary_files)]
metadata = data.frame(file_name = evaluation_summary_files,
                      file_name_shap = evaluation_summary_files%>%gsub("evaluation_summary","shap_summary",.))%>%
  mutate(input_type = file_name%>%dirname()%>%dirname()%>%dirname()%>%dirname()%>%dirname()%>%basename()%>%gsub(".+validation_|.+_out_","",.))%>%
  rowwise()%>%
  mutate(split_type = ifelse(grepl("cross_validation",file_name%>%dirname()%>%dirname()%>%dirname()%>%dirname()%>%dirname()%>%basename()),"cross_validation","contributor_out"))%>%
  mutate(predictor_group = file_name%>%dirname()%>%dirname()%>%dirname()%>%dirname()%>%basename()%>%gsub(".+feature..","",.)%>%gsub("trait..ad","all",.))

# summary_stats = overall performance of each model
lapply(1:nrow(metadata), function(i){
  summary_stats = metadata$file_name[i]%>%readRDS()%>%.$summary_stats%>%
    as.data.frame()%>%
    mutate(file_name = metadata$file_name[i])
  return(summary_stats)
})%>%
  bind_rows() -> summary_stats 

metadata%>%
  ungroup()%>%
  inner_join(.,summary_stats,by="file_name") -> summary_stats
saveRDS(summary_stats,file=file.path(outdir,"summary_stats.rds"))
fwrite(summary_stats,file=file.path(outdir,"summary_stats.txt"),sep="\t",append = F,row.names = F,col.names = T)

# fold_metrics = performance metrics for each model and each fold
lapply(1:nrow(metadata), function(i){
  summary_stats = metadata$file_name[i]%>%readRDS()%>%.$fold_metrics%>%
    as.data.frame()%>%
    mutate(file_name = metadata$file_name[i])
  return(summary_stats)
})%>%
  bind_rows() -> fold_metrics 

# var_numbers = number of predictors in each model (from the SHAP summary),
# counted per column and per trait (class-probability columns merged)
lapply(1:nrow(metadata), function(i){
  summary_stats = metadata$file_name_shap[i]%>%readRDS()%>%.$summary_stats%>%
    as.data.frame()
  n_var = length(summary_stats$feature)
  n_var_no_class = length(summary_stats$feature%>%gsub("_class.+","",.)%>%unique())
  
  data.frame(file_name = metadata$file_name[i],
             n_var = n_var,
             n_var_no_class=n_var_no_class)%>%
    return()
})%>%
  bind_rows() -> var_numbers


metadata%>%
  ungroup()%>%
  inner_join(.,fold_metrics,by="file_name")%>%
  inner_join(.,var_numbers,by="file_name") -> fold_metrics
saveRDS(fold_metrics,file=file.path(outdir,"fold_metrics.rds"))
fwrite(fold_metrics,file=file.path(outdir,"fold_metrics.txt"),sep="\t",append = F,row.names = F,col.names = T)


