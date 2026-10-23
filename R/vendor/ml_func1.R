aggr_result<-function(PARAM){
  facsdata = readRDS(PARAM$efile)
  facsdata = PARAM$select_samples(facsdata)
  
  ml_res = PARAM$ml_res
  res_dirs = dirname(list.files(paste0(ml_res),
                                pattern = "preProc.RData",
                                recursive = T,full.names = T))
  dat_for_ml = t(facsdata$data)
  ds = Reduce(function(x, y) merge(x, y, all=T, 
                                   by=c("SAMPLE.ID")), 
              lapply(res_dirs, function(x){
                load(paste0(x,"/preProc.RData"))
                load(paste0(x,"/model_fit/model_fit.RData"))
                indx = match(colnames(preProc$data),colnames(dat_for_ml))
                dat_for_ml_sub = dat_for_ml[,indx]
                indx = (rowSums(is.na(dat_for_ml_sub)))/ncol(dat_for_ml_sub)<.1
                sum(indx)
                dat_for_ml_sub = dat_for_ml_sub[indx,]
                dat_for_ml_norm = predict(preProc, dat_for_ml_sub)
                res  = predict(model_fit,newdata= dat_for_ml_norm,type="prob")
                res = res[,c("RA"),drop=FALSE]
                x = gsub(ml_res,"",x)
                x = gsub("/submodels/",".",x)
                x = gsub("/","",x)
                colnames(res) = x
                res$SAMPLE.ID = rownames(dat_for_ml_norm)
                return(res)
              }))
  rownames(ds)=ds$SAMPLE.ID
  ds = ds[,-1]
  ds = t(ds)
  
  
  facsdata$data = ds
  indx = match(colnames(ds),facsdata$sample.ann$SAMPLE.ID)
  facsdata$sample.ann = facsdata$sample.ann[indx,]
  facsdata$gene.ann = data.frame(Name = rownames(ds),row.names = rownames(ds))
  facsdata$gene.ann$Model = 
    sapply(strsplit(as.character(facsdata$gene.ann$Name),"\\."),function(x)x[1])
  facsdata$gene.ann$Run = 
    sapply(strsplit(as.character(facsdata$gene.ann$Name),"\\."),function(x)x[2])
  
  
  res_dirs = (list.files(paste0(ml_res),
                         pattern = "model.score.rds",
                         recursive = T,full.names = T))
  
  ds = Reduce(function(x, y) rbind(x, y), 
              lapply(res_dirs, function(x){
                res = readRDS(x)
                x = gsub(ml_res,"",dirname(x))
                x = gsub("/summary","",x)
                x = gsub("/","",x)
                rownames(res) = 
                  gsub("/b_validation/evalmodel.RData","",rownames(res))
                rownames(res) = 
                  gsub("submodels/","",rownames(res))
                rownames(res) = paste0(x,".",rownames(res))
                return(res)
              }))
  
  indx = match(rownames(facsdata$gene.ann),rownames(ds))
  facsdata$gene.ann = cbind(facsdata$gene.ann,ds)
  
  return(facsdata)
}


summary_v2 <-function(PARAM){
  
  # collect model performance
  if (file.exists(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))){
    model.score=readRDS(paste0(PARAM$OUTLOC,"/summary/model.score.rds"))
  }else{
    
    flist = list.files(PARAM$OUTLOC,pattern = "evalmodel.rds",recursive = TRUE) 
    model.score = do.call(rbind,
                          lapply(flist, function(x){
                            evalmodel=readRDS(paste0(PARAM$OUTLOC,x))
                            return(evalmodel$model_score)}))
    rownames(model.score)=flist
    dir.create(paste0(PARAM$OUTLOC,"/summary"),showWarnings = F)
    saveRDS(model.score,file=paste0(PARAM$OUTLOC,"/summary/model.score.rds"))
    
    flist = list.files(PARAM$OUTLOC,pattern = "evalmodel_train.rds",recursive = TRUE) 
    model.score = do.call(rbind,
                          lapply(flist, function(x){
                            evalmodel=readRDS(paste0(PARAM$OUTLOC,x))
                            return(evalmodel$model_score)}))
    rownames(model.score)=flist
    dir.create(paste0(PARAM$OUTLOC,"/summary"),showWarnings = F)
    saveRDS(model.score,file=paste0(PARAM$OUTLOC,"/summary/model.score_train.rds"))
    
  }
  
  # if (sum(!is.na(model.score[,2]))>0){
  #   pdf(paste0(PARAM$OUTLOC,"/summary/ACC_KAPPA.pdf"))
  #   plot(model.score)
  #   dev.off()
  # }
  
  
  # collect variable importance
  if (file.exists(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))){
    model.score=readRDS(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))
  }else{
    flist = list.files(PARAM$OUTLOC,pattern = "varImp.rds",recursive = TRUE) 
    if (length(flist)==0){
      return()
    }
    library(caret)
    
    model.score = do.call(rbind,lapply(flist[-2], function(x){
      a=readRDS(paste0(PARAM$OUTLOC,x))
      a$importance$Var1 = rownames(a$importance)
      #a$importance$Var1
      a$importance$Model.name = x
      return(a$importance)}))
    
    model.score$Var1 = gsub("`","",model.score$Var1)
    saveRDS(model.score,file=paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))
  }
  
  # Step 1: Identify the class columns dynamically
  # Assuming Var1 and Model.name are fixed columns, and all other columns are class labels
  model_columns <- c("Var1", "Model.name")
  class_columns <- setdiff(colnames(model.score), model_columns)
  
  # Step 2: Reshape data from wide to long format
  model.score_long <- model.score %>%
    pivot_longer(cols = all_of(class_columns), 
                 names_to = "Class", 
                 values_to = "Importance")
  
  # Step 3: Calculate mean importance by variable across all classes
  model.score_mean <- model.score_long %>%
    group_by(Var1) %>%
    summarise(Mean = mean(Importance, na.rm = TRUE))
  
  # Step 4: Select top 20 variables based on overall mean importance
  top20 <- model.score_mean$Var1[order(model.score_mean$Mean, decreasing = TRUE)[1:20]]
  
  # Step 5: Filter data to include only top 20 variables
  model.score_filtered <- model.score_long %>%
    filter(Var1 %in% top20)
  
  # Step 6: Convert to factor and set levels for proper ordering
  model.score_filtered$Var1 <- factor(model.score_filtered$Var1, levels = top20)
  
  # Step 7: Create visualizations
  # Boxplot with class-specific coloring
  g1 <- ggplot(model.score_filtered, aes(x = Var1, y = Importance, fill = Class)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, angle = 60, hjust = 1, colour = "grey10"),
          legend.position = "right") +
    labs(title = "Variable Importance by Class",
         x = "Variables",
         y = "Importance Score")
  
  # Faceted plot by class
  g2 <- ggplot(model.score_filtered, aes(x = Var1, y = Importance)) +
    geom_boxplot() +
    facet_wrap(~ Class) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, angle = 60, hjust = 1, colour = "grey10")) +
    labs(title = "Variable Importance by Class",
         x = "Variables",
         y = "Importance Score")
  
  # Save plots
  pdf(paste0(PARAM$OUTLOC, "/summary/VarImp.pdf"), width = 12, height = 8)
  print(g1)
  print(g2)
  dev.off()
  
  # 
  # # varimp similarity
  # model.score%>%
  #   mutate(Model.name=gsub("submodels/|/varImp.rds","",Model.name))%>%
  #   spread(Model.name,Overall)%>%
  #   dplyr::select(-Var1)%>%as.matrix()%>%cor()%>%
  #   as.data.frame()%>%
  #   mutate(modelA=rownames(.))%>%
  #   gather(modelB,cor,-modelA)%>%
  #   saveRDS(.,file=paste0(PARAM$OUTLOC,"/summary/varimp_sim.rds"))
  
  return(TRUE)
}

summary1 <-function(PARAM){
  
  # collect model performance
  if (file.exists(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))){
    model.score=readRDS(paste0(PARAM$OUTLOC,"/summary/model.score.rds"))
  }else{
    
    flist = list.files(PARAM$OUTLOC,pattern = "evalmodel.rds",recursive = TRUE) 
    model.score = do.call(rbind,
                          lapply(flist, function(x){
                            evalmodel=readRDS(paste0(PARAM$OUTLOC,x))
                            return(evalmodel$model_score)}))
    rownames(model.score)=flist
    dir.create(paste0(PARAM$OUTLOC,"/summary"),showWarnings = F)
    saveRDS(model.score,file=paste0(PARAM$OUTLOC,"/summary/model.score.rds"))
  }
  
  if (sum(!is.na(model.score[,2]))>0){
    pdf(paste0(PARAM$OUTLOC,"/summary/ACC_KAPPA.pdf"))
    plot(model.score)
    dev.off()
  }
  
  
  # collect variable importance
  if (file.exists(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))){
    model.score=readRDS(paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))
  }else{
    flist = list.files(PARAM$OUTLOC,pattern = "model_fit.rds",recursive = TRUE) 
    library(caret)
    
    model.score = do.call(rbind,lapply(flist, function(x){
      model_fit=readRDS(paste0(PARAM$OUTLOC,x))
      a =varImp(model_fit)
      a$importance$Var1 = rownames(a$importance)
      a$importance$Var1
      a$importance$Model.name = x
      return(a$importance)}))
    
    model.score$Var1 = gsub("`","",model.score$Var1)
    saveRDS(model.score,file=paste0(PARAM$OUTLOC,"/summary/VarImp.rds"))
  }
  
  model.score_mean = group_by(model.score,Var1)%>%summarise(Mean = mean(Overall))
  top20 = model.score_mean$Var1[order(model.score_mean$Mean,decreasing = TRUE)[1:20]]
  indx = model.score$Var1%in%top20
  model.score = model.score[indx,]
  model.score$Var1 = factor(model.score$Var1,levels = top20)
  library(ggthemes)
  g = ggplot(model.score,aes(Var1,Overall))+geom_boxplot()+
    theme(axis.text.x = element_text(size=15,angle = 60, hjust = 1, colour = "grey10"))
  
  pdf(paste0(PARAM$OUTLOC,"/summary/VarImp.pdf"))
  plot(g)
  dev.off()
  
  # train sample similarity
  SUBMODEL_LOC = paste0(PARAM$OUTLOC,"/submodels/")
  if(length( list.dirs(SUBMODEL_LOC,recursive = F,full.names = T))==1){return(TRUE)}
  
  list.dirs(SUBMODEL_LOC,recursive = F,full.names = T)%>%
    lapply(., function(x){
      load(paste0(x,"/inTrain.RData"))
      out=list()
      out[[basename(x)]]=inTrain[,1]
      return(out)
    })->res
  
  lapply(1:length(res), function(i){
    train_sim = data.frame(modelA=names(res[[i]]),
                           modelB=res%>%lapply(.,function(x)names(x))%>%unlist(),
                           n_overlap = NA,
                           n_total=length(res[[i]][[1]]))
    for(j in 1:length(res)){
      train_sim$n_overlap[j]=sum(res[[i]][[1]]%in%res[[j]][[1]])
    }
    return(train_sim)
  })%>%do.call(rbind,.)%>%
    mutate(ratio_overlap=n_overlap/n_total)%>%
    saveRDS(.,file=paste0(PARAM$OUTLOC,"/summary/train_sim.rds"))
  
  
  # varimp similarity
  model.score%>%
    mutate(Model.name=gsub("/model_fit/model_fit.RData|submodels/","",Model.name))%>%
    spread(Model.name,Overall)%>%
    dplyr::select(-Var1)%>%as.matrix()%>%cor()%>%
    as.data.frame()%>%
    mutate(modelA=rownames(.))%>%
    gather(modelB,cor,-modelA)%>%
    saveRDS(.,file=paste0(PARAM$OUTLOC,"/summary/varimp_sim.rds"))
  
  return(TRUE)
}


#Balanced split function
bsplit <-function(PARAM,dat_for_ml,n){
  library(caret)
  
  outloc=paste0(PARAM$OUTLOC,"submodels/model",n,"/")
  if(file.exists(paste0(outloc,"model_fit/model_fit.RData"))){
    return()
  }
  dir.create(outloc,recursive=TRUE,showWarnings = F)
  
  dir.create(paste0(outloc,"b_validation/"),showWarnings = F)
  dir.create(paste0(outloc,"ub_validation/"),showWarnings = F)
  dir.create(paste0(outloc,"model_fit/"),showWarnings = F)
  
  set.seed(n)
  
  if (file.exists(paste0(outloc,"inTrain.RData"))){
    load(paste0(outloc,"inTrain.RData"))
  }else{
    
    # inTrain = createDataPartition(y = dat_for_ml[[PARAM$OUTCOME]],
    #                               p = 1-PARAM$B_RATIO, list=FALSE) 
    # save(inTrain,file=paste0(outloc,"inTrain.RData"))
  }
  
  
  preProc  <- preProcess(dat_for_ml[ inTrain,-1],verbose=FALSE,
                         method=c("knnImpute",PARAM$EXP_CENTER,PARAM$EXP_SCALE))
  save(preProc,file=paste0(outloc,"preProc.RData"))
  
  dat_for_train <- predict(preProc, dat_for_ml[ inTrain,-1])
  dat_for_test     <- predict(preProc, dat_for_ml[ -inTrain,-1])
  ## training data 
  dat_for_train[[PARAM$OUTCOME]] = dat_for_ml[[PARAM$OUTCOME]][inTrain]
  ## test data
  dat_for_test[[PARAM$OUTCOME]] = dat_for_ml[[PARAM$OUTCOME]][-inTrain]
  
  
  dat_for_test_val = dat_for_test[,!colnames(dat_for_test)%in%PARAM$OUTCOME]
  dat_for_test_class= dat_for_test[,PARAM$OUTCOME]
  
  
  
  if (PARAM$N_UBSPLIT==0){
    ## build model
    if (PARAM$MODEL=="xgbTree"){
      model_fit = xgbTree_build(PARAM,dat_for_train)
    }else if(PARAM$MODEL=="glmnet"){
      model_fit = enet_build(PARAM,dat_for_train)
    }else if(PARAM$MODEL=="pls"){
      model_fit = pls_build(PARAM,dat_for_train)
    }
    
    
    # evaluation 
    evalmodel = eval_model(model_fit,
                           dat_for_test_val,
                           dat_for_test_class,
                           metric = PARAM$mmetric)
    
    
    saveRDS(evalmodel,file=paste0(outloc,"b_validation/evalmodel.rds"))
    saveRDS(model_fit,file=paste0(outloc,"model_fit/model_fit.rds"))
  }
  
  
  
  
  #   for(k in 1:PARAM$N_UBSPLIT){
  #     
  #     model_fit=ubsplit(PARAM,dat_for_train,outloc,k)
  #     
  #     # evaluation 
  #     evalmodel = eval_model(model_fit,dat_for_test_val,dat_for_test_class)
  #     
  #     save(evalmodel,file=paste0(outloc,"b_validation/",k,"_evalmodel.RData"))
  #     save(model_fit,file=paste0(outloc,"model_fit/",k,"_modle_fit.RData"))
  #   }
  
}

bsplit_v2 <-function(PARAM,dat_for_ml,n){
  library(caret)
  
  outloc=paste0(PARAM$OUTLOC,"submodels/model",n,"/")
  if(file.exists(paste0(outloc,"evalmodel_train.rds"))){
    return()
  }
  dir.create(outloc,recursive=TRUE,showWarnings = F)
  
  #  dir.create(paste0(outloc,"b_validation/"),showWarnings = F)
  #  dir.create(paste0(outloc,"ub_validation/"),showWarnings = F)
  #  dir.create(paste0(outloc,"model_fit/"),showWarnings = F)
  
  set.seed(n)
  
  if (file.exists(paste0(outloc,"inTrain.rds"))){
    inTrain = readRDS(paste0(outloc,"inTrain.rds"))
  }else{
    
    # inTrain = createDataPartition(y = dat_for_ml[[PARAM$OUTCOME]],
    #                               p = 1-PARAM$B_RATIO, list=FALSE) 
    # save(inTrain,file=paste0(outloc,"inTrain.RData"))
  }
  
  dat_for_train <-dat_for_ml[ inTrain,-1,drop=F]
  dat_for_test  <- dat_for_ml[ -inTrain,-1,drop=F]
  ## training data 
  dat_for_train[[PARAM$OUTCOME]] = dat_for_ml[[PARAM$OUTCOME]][inTrain]
  ## test data
  dat_for_test[[PARAM$OUTCOME]] = dat_for_ml[[PARAM$OUTCOME]][-inTrain]
  
  
  dat_for_test_val = dat_for_test[,!colnames(dat_for_test)%in%PARAM$OUTCOME,drop=F]
  dat_for_test_class= dat_for_test[,PARAM$OUTCOME,drop=F]
  
  
  
  if (PARAM$N_UBSPLIT==0){
    ## build model
    if (PARAM$MODEL=="xgbTree"){
      model_fit = xgbTree_build(PARAM,dat_for_train)
    }else if(PARAM$MODEL=="glmnet"){
      model_fit = enet_build(PARAM,dat_for_train)
    }else if(PARAM$MODEL=="pls"){
      indx=PARAM$pls_ncomp<=(ncol(dat_for_train)-1)
      PARAM$pls_ncomp=PARAM$pls_ncomp[indx]
      model_fit = pls_build(PARAM,dat_for_train)
    }else if(PARAM$MODEL=="xgbTree_mlr"){
      model_fit = xgbTree_mlr_build(PARAM,dat_for_train,outloc)
    }else if(PARAM$MODEL=="xgbTree_mlr3"){
      model_fit = xgbTree_mlr3_build(PARAM,dat_for_train,outloc)
    }
    
    if(PARAM$MODEL=="xgbTree_mlr"){
      #test
      evalmodel = eval_model_mlr(model_fit,
                                 dat_for_test_val,
                                 dat_for_test_class[[PARAM$OUTCOME]],
                                 metric=PARAM$mmetric,
                                 PARAM)
      saveRDS(evalmodel,file=paste0(outloc,"evalmodel.rds"))
      
      # train
      evalmodel = eval_model_mlr(model_fit,
                                 dat_for_train,
                                 dat_for_train[[PARAM$OUTCOME]],
                                 metric = PARAM$mmetric,
                                 PARAM)
      saveRDS(evalmodel,file=paste0(outloc,"evalmodel_train.rds"))
      
      
      # shap
      data=as.matrix(dat_for_ml[,-1,drop=F])
      colnames(data)=make.names(colnames(data))
      data=data[,model_fit$learner.model$feature_names,drop=F]
      shap_contrib <- predict(model_fit$learner.model,
                              data,
                              predcontrib = TRUE, approxcontrib = FALSE)
      
      
      calc_shap = function(shap_contrib){
        rownames(shap_contrib)=rownames(dat_for_ml)
        shap_contrib=shap_contrib[,-which(colnames(shap_contrib)=="BIAS"),drop=F]
        
        shap=data.frame(Var=colnames(shap_contrib),
                        shap=colMeans(abs(shap_contrib)))
        shap$cor=NA
        
        for (i in 1:ncol(shap_contrib)){
          shap$cor[i]=cor(data[,colnames(shap_contrib)[i]],
                          shap_contrib[,colnames(shap_contrib)[i]],
                          method = "spearman")
        }
        return(shap)
      }
      
      if(is.list(shap_contrib)){
        shap = lapply(shap_contrib, function(x){
          calc_shap(x)
        })
      }else{
        shap = calc_shap(shap_contrib)
      }


      
      saveRDS(shap,file=paste0(outloc,"/shap.rds"))
    }else{
      # evaluation 
      evalmodel = eval_model(model_fit,
                             dat_for_test_val,
                             dat_for_test_class,
                             metric = PARAM$mmetric)
      
      saveRDS(evalmodel,file=paste0(outloc,"evalmodel.rds"))
      
      tryCatch({
        saveRDS(varImp(model_fit),file=paste0(outloc,"varImp.rds"))
      }, warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })
      
      evalmodel = eval_model(model_fit,
                             dat_for_train,
                             dat_for_train[,PARAM$OUTCOME,drop=F],
                             metric = PARAM$mmetric)
      
      saveRDS(evalmodel,file=paste0(outloc,"evalmodel_train.rds"))
      
      
      #model_fit$finalModel=NULL
      model_fit$trainingData=NULL
      #model_fit$terms=NULL 
    }
    if(is.null(PARAM$save_model) || !PARAM$save_model){
      return()
    }else{
      saveRDS(model_fit,file=paste0(outloc,"model_fit.rds"))
    }
    
  }
  
}

#Unbalanced split function
ubsplit <- function(PARAM,dat_for_train,outloc,k){
  library(caret)
  nsample = dim(dat_for_train)[1]
  inRandTrain = sample(nsample,round(nsample*(1-PARAM$UB_RATIO)))
  
  ## training data 
  dat_for_randtrain = dat_for_train[inRandTrain,]
  ## test data
  dat_for_randtest= dat_for_train[-inRandTrain,]
  
  dat_for_randtest_val = dat_for_randtest[,!colnames(dat_for_randtest)%in%PARAM$OUTCOME]
  dat_for_randtest_class = dat_for_randtest[,PARAM$OUTCOME]
  
  ## build model
  model_fit = pls_build(PARAM,dat_for_randtrain)
  
  ## evaluation
  evalmodel = eval_model(model_fit,
                         dat_for_randtest_val,
                         dat_for_randtest_class,
                         metric=PARAM$mmetric)
  
  save(evalmodel,file=paste0(outloc,"ub_validation/",k,"_evalmodel.RData"))
  
  return(model_fit)
}

#Model buld for each split
enet_build<-function(PARAM,dat_for_randtrain){
  library(caret)
  # PLSDA

  if(!is.null(PARAM$summaryFunction) && PARAM$summaryFunction=="multiClassSummary"){
    tr = caret::trainControl(method="repeatedcv", 
                             number = PARAM$CV_FOLD,
                             repeats=PARAM$CV_REPEAT,
                             classProbs=PARAM$classProbs,
                             summaryFunction = multiClassSummary)
    #print(tr$summaryFunction)
  }else{
    tr = caret::trainControl(method="repeatedcv", 
                             number = PARAM$CV_FOLD,
                             repeats=PARAM$CV_REPEAT,
                             classProbs=PARAM$classProbs,
                             summaryFunction = defaultSummary)
    #print(tr$summaryFunction)
  }
  
  
  train_grid = expand.grid(alpha=PARAM$glmnet_alpha,
                           lambda=PARAM$glmnet_lambda)
  
  if(!is.null(PARAM$family)){
    model_fit = caret::train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                             data = dat_for_randtrain,
                             method = "glmnet",
                             tuneGrid = train_grid,
                             trControl=tr,
                             family = PARAM$family,
                             metric = PARAM$train_metric,
                             lower.limits=PARAM$glmnet_lower.limits,
                             upper.limits=PARAM$glmnet_upper.limits)
  }else{
    model_fit = caret::train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                             data = dat_for_randtrain,
                             method = "glmnet",
                             tuneGrid = train_grid,
                             trControl=tr,
                             metric = PARAM$train_metric,
                             lower.limits=PARAM$glmnet_lower.limits,
                             upper.limits=PARAM$glmnet_upper.limits)
  }
  return(model_fit)
  
}

pls_build<-function(PARAM,dat_for_randtrain){
  library(caret)
  # PLSDA
  tr = caret::trainControl(method="repeatedcv", 
                    number = PARAM$CV_FOLD,
                    repeats= PARAM$CV_REPEAT,
                    classProbs=PARAM$classProbs,
                    summaryFunction = defaultSummary)
  train_grid = expand.grid(ncomp = PARAM$pls_ncomp)
  model_fit = caret::train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                    data = dat_for_randtrain,
                    method = "pls",
                    tuneGrid = train_grid,
                    trControl=tr,
                    metric = PARAM$train_metric )
  return(model_fit)
  
}

pls_xgbTree_build<-function(PARAM,dat_for_randtrain){
  library(caret)
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = PARAM$CV_FOLD,
    ## repeated ten times
    repeats = 5,
    classProbs=TRUE)
  
  gbmGrid <-  expand.grid(nrounds=(1:20),
                          eta=seq(0.1,0.5,0.1),
                          max_depth=2:5)
  
  model_fit <- train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                     data = dat_for_randtrain,
                     method = "xgbTree",
                     trControl = fitControl,
                     tuneGrid = gbmGrid,
                     ## This last option is actually one
                     ## for gbm() that passes through
                     objective="binary:logistic",
                     nthread=1,
                     metric = "Kappa")
  
  predleafs = predict(model_fit$finalModel,
                      newdata = xgb.DMatrix(as.matrix(dat_for_randtrain[,-ncol(dat_for_randtrain)])),
                      predleaf=TRUE)
  colnames(predleafs) = paste0("xgb",1:ncol(predleafs))
  
  dat_for_randtrain =cbind(dat_for_randtrain,predleafs)
  # PLSDA
  tr = trainControl(method="repeatedcv", 
                    number = PARAM$CV_FOLD,
                    repeats=5,
                    classProbs=TRUE)
  train_grid = expand.grid(ncomp = 1:30)
  model_fit = train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                    data = dat_for_randtrain,
                    method = "pls",
                    tuneGrid = train_grid,
                    trControl=tr,
                    metric = "Kappa")
  return(model_fit)
  
}
library(rminer)
xgbTree_build<-function(PARAM,dat_for_randtrain){
  library(caret)
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = PARAM$CV_FOLD,
    ## repeated ten times
    repeats = PARAM$CV_REPEAT,
    classProbs=PARAM$classProbs)
  
  
  gbmGrid <-  expand.grid(nrounds=PARAM$xgbTree_nrounds,
                          eta=PARAM$xgbTree_eta,
                          max_depth=PARAM$xgbTree_max_depth,
                          gamma =  PARAM$xgbTree_gamma,
                          colsample_bytree = PARAM$xgbTree_colsample_bytree,
                          min_child_weight = PARAM$xgbTree_min_child_weight,
                          subsample=PARAM$xgbTree_subsample)
  
  model_fit <- train(as.formula(paste0(PARAM$OUTCOME," ~ .")),
                     data = dat_for_randtrain,
                     method = "xgbTree",
                     trControl = fitControl,
                     tuneGrid = gbmGrid,
                     ## This last option is actually one
                     ## for gbm() that passes through
                     objective=PARAM$xgbTree_objective,
                     nthread=32,
                     metric = PARAM$train_metric )
  return(model_fit)
}


xgbTree_mlr_build<-function(PARAM,dat_for_randtrain,outloc){
  library(mlr)
  library(xgboost)
  colnames(dat_for_randtrain)=make.names(colnames(dat_for_randtrain))
  


  if (is.factor(dat_for_randtrain[[PARAM$OUTCOME]]) & length(levels(dat_for_randtrain[[PARAM$OUTCOME]]))>2){
    traintask=makeClassifTask( data=dat_for_randtrain,
                               target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("classif.xgboost",predict.type = "prob")
    PARAM$xgbTree_par.vals$objective="multi:softprob"
  }else if (is.factor(dat_for_randtrain[[PARAM$OUTCOME]]) & length(levels(dat_for_randtrain[[PARAM$OUTCOME]]))==2){
    traintask=makeClassifTask( data=dat_for_randtrain,
                            target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("classif.xgboost",predict.type = "prob")
  }else{
    traintask=makeRegrTask( data=dat_for_randtrain,
                            target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("regr.xgboost",predict.type = "response")
  }
  

  lrn$par.vals <- PARAM$xgbTree_par.vals
  
  
  #set resampling strategy
  rdesc <- makeResampleDesc("CV",
                            iters=PARAM$CV_FOLD,
                            predict = "both")
  #search strategy
  ctrl <- makeTuneControlRandom(maxit = PARAM$xgbTree_nsamples)
  
  #set parallel backend
  library(parallel)
  library(parallelMap) 
  parallelStartSocket(cpus = PARAM$xgbTree_ncpus)
  #parameter tuning
  mytune <- mlr::tuneParams(learner = lrn, task = traintask, resampling = rdesc,
                       measures = PARAM$xgbTree_measures,
                       par.set = PARAM$xgbTree_params,
                       control = ctrl, show.info = T)
  parallelStop()
  
  saveRDS(mytune,file=paste0(outloc,"tune_obj.rds"))
  
  res_tune = mlr::generateHyperParsEffectData(mytune,partial.dep = T)
  saveRDS(res_tune,file=paste0(outloc,"tune_res.rds"))
  
  # res_tune$data%>%
  #   ggplot(.,aes(expvar.test.mean,expvar.train.mean))+
  #   geom_point()+geom_abline(slope=1,intercept=0)->p
  # 
  # cowplot::save_plot(paste0(outloc,"tune_res.pdf"),plot = p,base_height = 5,base_width = 5)
  
  lrn = mlr::setHyperPars(lrn, par.vals = mytune$x)
  model_fit = mlr::train(lrn, traintask)
  
  varImp.importance= try(getFeatureImportanceLearner(lrn,model_fit))
  if(inherits(varImp.importance, 'try-error')){
    
  }else{
    varImp=list()
    varImp$importance=getFeatureImportanceLearner(lrn,model_fit)%>%data.frame()
    colnames(varImp$importance)="Overall"
    varImp$model="xgbTree_mlr"
    varImp$calledFrom="getFeatureImportanceLearner"
    saveRDS(varImp,file=paste0(outloc,"varImp.rds"))
  }
  return(model_fit)
}
#Evaluate function
eval_model<-function(model_fit,
                     dat_for_randtest_val,
                     dat_for_randtest_class,
                     metric=c("ACC","KAPPA")){
  library(caret)
  library(rminer)
  evalmodel = list()
  model_predict_res = predict(model_fit, dat_for_randtest_val)
  model_prediction = data.frame(model_predict_res = model_predict_res,
                                dat_for_randtest_class = dat_for_randtest_class)  
  rownames(model_prediction)=rownames(dat_for_randtest_val)
  
  ## prediction_accuracy
  model_score = mmetric(model_predict_res,
                        dat_for_randtest_class[,1],
                        metric=metric)
  
  if(model_fit$modelType=="Classification"){
    evalmodel$model_prob = predict(model_fit, dat_for_randtest_val,type="prob")
  }
  evalmodel$model_prediction = model_prediction
  evalmodel$model_score = model_score
  return(evalmodel)
}


eval_model_mlr<-function(model_fit,
                         dat_for_test_val,
                         dat_for_test_class,
                         metric=c("ACC","KAPPA"),
                         PARAM){
  
  dat_for_test_val[[PARAM$OUTCOME]]=dat_for_test_class
  colnames(dat_for_test_val)=make.names(colnames(dat_for_test_val))
  
  if (is.factor(dat_for_test_val[[PARAM$OUTCOME]])){
    model_fit$learner$predict.type="prob"
    traintask=makeClassifTask( data=dat_for_test_val, target=PARAM$OUTCOME)
  }else{
    traintask=makeRegrTask( data=dat_for_test_val, target=PARAM$OUTCOME)
  }
  
  
  library(rminer)
  evalmodel = list()
  model_predict_res = predict(model_fit, traintask)
  model_prediction = data.frame(model_predict_res = model_predict_res$data$response,
                                dat_for_randtest_class = model_predict_res$data$truth)  
  rownames(model_prediction)=model_predict_res$data%>%rownames()
  
  ## prediction_accuracy
  model_score = mmetric(model_predict_res$data$response,
                        model_predict_res$data$truth,
                        metric=metric)
  
  if (is.factor(dat_for_test_val[[PARAM$OUTCOME]])){
    evalmodel$model_prob = model_predict_res$data
    evalmodel$model_prob$response = NULL
    evalmodel$model_prob$truth = NULL
  }
  
  evalmodel$model_prediction = model_prediction
  evalmodel$model_score = model_score
  return(evalmodel)
}

#dat_for_randtest_val=dat_for_test_val
#dat_for_randtest_class=dat_for_test_class
#mmetric(model_predict_res,dat_for_randtest_class,c("ACC","KAPPA"))


xgbTree_mlr3_build<-function(PARAM,dat_for_randtrain,outloc){
  library(mlr3)
  library(xgboost)
  colnames(dat_for_randtrain)=make.names(colnames(dat_for_randtrain))
  
  
  
  if (is.factor(dat_for_randtrain[[PARAM$OUTCOME]]) & length(levels(dat_for_randtrain[[PARAM$OUTCOME]]))>2){
    traintask=makeClassifTask( data=dat_for_randtrain,
                               target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("classif.xgboost",predict.type = "prob")
    PARAM$xgbTree_par.vals$objective="multi:softprob"
  }else if (is.factor(dat_for_randtrain[[PARAM$OUTCOME]]) & length(levels(dat_for_randtrain[[PARAM$OUTCOME]]))==2){
    traintask=makeClassifTask( data=dat_for_randtrain,
                               target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("classif.xgboost",predict.type = "prob")
  }else{
    traintask=makeRegrTask( data=dat_for_randtrain,
                            target=PARAM$OUTCOME)
    #create learner
    lrn <- makeLearner("regr.xgboost",predict.type = "response")
  }
  
  
  lrn$par.vals <- PARAM$xgbTree_par.vals
  
  
  #set resampling strategy
  rdesc <- makeResampleDesc("CV",
                            iters=PARAM$CV_FOLD,
                            predict = "both")
  #search strategy
  ctrl <- makeTuneControlRandom(maxit = PARAM$xgbTree_nsamples)
  
  #set parallel backend
  library(parallel)
  library(parallelMap) 
  parallelStartSocket(cpus = PARAM$xgbTree_ncpus)
  #parameter tuning
  mytune <- mlr::tuneParams(learner = lrn, task = traintask, resampling = rdesc,
                            measures = PARAM$xgbTree_measures,
                            par.set = PARAM$xgbTree_params,
                            control = ctrl, show.info = T)
  parallelStop()
  
  saveRDS(mytune,file=paste0(outloc,"tune_obj.rds"))
  
  res_tune = mlr::generateHyperParsEffectData(mytune,partial.dep = T)
  saveRDS(res_tune,file=paste0(outloc,"tune_res.rds"))
  
  # res_tune$data%>%
  #   ggplot(.,aes(expvar.test.mean,expvar.train.mean))+
  #   geom_point()+geom_abline(slope=1,intercept=0)->p
  # 
  # cowplot::save_plot(paste0(outloc,"tune_res.pdf"),plot = p,base_height = 5,base_width = 5)
  
  lrn = mlr::setHyperPars(lrn, par.vals = mytune$x)
  model_fit = mlr::train(lrn, traintask)
  
  varImp.importance= try(getFeatureImportanceLearner(lrn,model_fit))
  if(inherits(varImp.importance, 'try-error')){
    
  }else{
    varImp=list()
    varImp$importance=getFeatureImportanceLearner(lrn,model_fit)%>%data.frame()
    colnames(varImp$importance)="Overall"
    varImp$model="xgbTree_mlr"
    varImp$calledFrom="getFeatureImportanceLearner"
    saveRDS(varImp,file=paste0(outloc,"varImp.rds"))
  }
  return(model_fit)
}
