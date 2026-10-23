# run_pipeline() and its helpers from the pipelineUtil package (version 0.1.0),
# copied as plain R code. Identical to the version installed on the AD Workbench

run_pipeline <- 
function (mainfunc, param_location, param_suffix, OUTDIR, PARAM = NULL, 
    return_first = FALSE, SGE = FALSE, MCLAPPLY = FALSE, MC.cores = 8, 
    main_script = NULL, sge_options = NULL, run_pnum = NULL, 
    ...) 
{
    params = list.files(path = param_location, pattern = paste0("^", 
        param_suffix, "*"), recursive = F, full.names = T)
    if (length(params) == 0) {
        return()
    }
    PARAM_master = PARAM
    PARAM_LOCs = c()
    for (i in 1:length(params)) {
        PARAM = PARAM_master
        pname = gsub(paste0(param_suffix, "|\\.R"), "", basename(params[i]))
        pname = gsub("_.+", "", pname)
        if (!is.null(run_pnum)) {
            if (sum(as.numeric(pname) %in% run_pnum) != 1) {
                (next)()
            }
        }
        source(params[i], local = T)
        PARAM = clear_environment(PARAM)
        if (is.null(PARAM$combn)) {
            n_comb = 1
        }
        else {
            lcomb = list_cmb(PARAM, PARAM$combn)
            n_comb = nrow(lcomb)
        }
        for (n in 1:n_comb) {
            PARAM_temp = PARAM
            p_desc = c()
            for (ln in PARAM$combn) {
                if (lcomb[n, ln] == "NULL") {
                  PARAM_temp[[ln]] = NULL
                }
                else {
                  if (is.list(lcomb[n, ln])) {
                    PARAM_temp[[ln]] = lcomb[n, ln][[1]]
                  }
                  else {
                    PARAM_temp[[ln]] = lcomb[n, ln]
                  }
                }
                if (is.list(lcomb[n, ln])) {
                  p_desc = paste(c(p_desc, paste0(ln, "..", names(lcomb[n, 
                    ln]))), collapse = "__")
                }
                else {
                  p_desc = paste(c(p_desc, paste0(ln, "..", lcomb[n, 
                    ln])), collapse = "__")
                }
            }
            print(paste0(n, " p", pname, "_", paste(PARAM_temp[c("label")], 
                collapse = "_"), "_", p_desc))
            PARAM_temp$OUTLOC = paste0(OUTDIR, "/p", pname, "_", 
                PARAM_temp["label"])
            PARAM_temp$OUTLOC = paste0(PARAM_temp$OUTLOC, "/", 
                p_desc, "/")
            dir.create(PARAM_temp$OUTLOC, recursive = T, showWarnings = F)
            PARAM_temp$OUTLOC = paste0(normalizePath(PARAM_temp$OUTLOC), 
                "/")
            PARAM_temp = clear_environment(PARAM_temp)
            if (return_first & !SGE) {
                return(PARAM_temp)
            }
            if (SGE) {
                PARAM_LOC = paste0(PARAM_temp$OUTLOC, "/param_config.rds")
                saveRDS(PARAM_temp, file = PARAM_LOC)
                job_name = paste0("p", pname, "_", PARAM_temp$label, 
                  "_", p_desc)
                job_name = make.names(job_name)
                dir.create("job_scripts", showWarnings = F, recursive = T)
                cmd = paste0(c(". ~/.bash_profile\nuse GCC-5.2\nR --vanilla --slave --args", 
                  PARAM_LOC, getwd(), paste0("< ", main_script)), 
                  collapse = " ")
                write(cmd, file = paste0("job_scripts/", job_name, 
                  n, ".sh"))
                cmd = paste0(c("qsub ", sge_options, paste0("job_scripts/", 
                  job_name, n, ".sh")), collapse = " ")
                if (return_first) {
                  return(cmd)
                }
                system(cmd)
                cmd = paste0("qstat | wc -l")
                NJOBS = as.numeric(system(cmd, intern = TRUE))
                while (NJOBS > 990) {
                  cmd = paste0("qstat | wc -l")
                  NJOBS = as.numeric(system(cmd, intern = TRUE))
                  Sys.sleep(1)
                }
            }
            else if (MCLAPPLY) {
                PARAM_LOC = paste0(PARAM_temp$OUTLOC, "/param_config.rds")
                saveRDS(PARAM_temp, file = PARAM_LOC)
                PARAM_LOCs = c(PARAM_LOCs, PARAM_LOC)
            }
            else {
                mainfunc(PARAM_temp, ...)
            }
        }
    }
    if (MCLAPPLY) {
        library(parallel)
        mclapply(PARAM_LOCs, function(PARAM_LOC) {
            PARAM_temp = readRDS(PARAM_LOC)
            mainfunc(PARAM_temp, ...)
        }, mc.cores = MC.cores)
    }
}

clear_environment <- 
function (param_config) 
{
    indx = which(sapply(param_config, is.function))
    for (i in indx) {
        environment(param_config[[i]]) = environment()
    }
    return(param_config)
}

list_cmb <- 
function (data, list.names) 
{
    nclass = prod((sapply(data[list.names], length)))
    cmb = matrix(NA, nrow = nclass, ncol = length(list.names))
    colnames(cmb) = list.names
    cmb = data.frame(cmb)
    for (ln in list.names) {
        ln.num = length(data[[ln]])
        cmb[, ln] = rep(1:ln.num, nclass/ln.num)
        cmb = cmb[order(cmb[, ln]), , drop = FALSE]
    }
    for (ln in list.names) {
        cmb[[ln]] = data[[ln]][cmb[, ln]]
    }
    return(cmb)
}

