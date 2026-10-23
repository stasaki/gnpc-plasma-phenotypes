# ml_session(): runs the CatBoost cross-validation for one classifier config.
# Writes the results under <run folder>/ml/ in step_dir("04_classifier_runs").
# Runs on the AD Workbench.

ml_session <- function(expdata,param_config){
  
  if (is.null(param_config$ml$flag) || !param_config$ml$flag) {
    return()
  }
  
  print("start ml_session")
  param_config$OUTLOC = paste0(param_config$OUTLOC,"/ml/")    
  dir.create(param_config$OUTLOC,showWarnings = F,recursive = T)
  
  
  if (file.exists(paste0(param_config$OUTLOC,"/param_config.RData"))&!param_config$OVERWRITE){
  }
  
  subfunc<-function(param_config,dat_for_ml){
    
    save(param_config,file=paste0(param_config$OUTLOC,"/param_config.RData"))
    saveRDS(dat_for_ml$IID,file=paste0(param_config$OUTLOC,"/sample_id.rds"))
    if(!"fold_col"%in%names(param_config)){
      # make split
      flag=TRUE
      model_total =param_config$N_BSPLIT*param_config$N_BSPLIT_TIMES
      for(n in 1:model_total){
        outloc=paste0(param_config$OUTLOC,"submodels/model",n,"/")
        flag=flag&file.exists(paste0(outloc,"inTrain.rds"))
      }
      
      if (!flag){
        
        set.seed(param_config$split)
        param_config$split_sample(expdata,model_total) -> data_split
        for (n in 1:model_total){
          outloc=paste0(param_config$OUTLOC,"submodels/model",n,"/")
          dir.create(outloc,recursive = T)
          which(dat_for_ml$IID%in%c(data_split%>%
                                      filter(group!=n)%>%.$IID))%>%
            saveRDS(.,file=paste0(outloc,"inTrain.rds"))
        }
      }
      dat_for_ml=dat_for_ml%>%
        dplyr::select(-IID)
      
      mclapply(1:param_config$N_BSPLIT*param_config$N_BSPLIT_TIMES,function(n){
        bsplit_v2(param_config,dat_for_ml,n)
      },mc.cores = 1)
      
      summary_v2(param_config)
    }else{
      dat_for_ml=dat_for_ml%>%
        dplyr::select(-IID)
      for(fold_indx in unique(dat_for_ml[[param_config$fold_col]])){
        param_config$test_fold = fold_indx
        param_config$results_dir = paste0(param_config$OUTLOC,"submodels/model_",fold_indx,"/")
        train_model(dat_for_ml, param_config)
      }
      cv_results <- collect_cv_results(param_config, dat_for_ml)
    }
  }
  
  run_pipeline(subfunc,
               param_location=param_config$trait$param_loc,
               param_suffix ="ml_param",
               OUTDIR= param_config$OUTLOC,
               PARAM=param_config,return_first=F,
               SGE=FALSE,MCLAPPLY=FALSE,MC.cores=8,
               main_script=NULL,sge_options=NULL,
               run_pnum=NULL,expdata)
  
  save(param_config,file=paste0(param_config$OUTLOC,"/param_config.RData"))
  
  print("end")
  return(TRUE)
}
