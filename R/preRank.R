preRank <- function(y.vect, id.vect, M, COV=NULL, method = c("LMM","GEE"), parallel=TRUE, ncore = detectCores(), write=FALSE)
{
  if (length(method)>1) method = "LMM"  # By default the model is LMM
  else if (!method %in% c("LMM","GEE")) stop("The method should be 'LMM' or 'GEE'")
   
  M.names = colnames(M) 
  L.M = dim(M)[2]
  
  if (is.null(COV))
  {
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA)
    modelstatement_LMM = y.vect ~ Mone + (1|id.vect)
    modelstatement_GEE = y.vect ~ Mone
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)          
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA, COV = COV)
    modelstatement_LMM = eval(parse(text=(paste0("y.vect ~ Mone +", paste0(paste0("COV.",conf.names),collapse = "+"), "+ (1|id.vect)"))))
    modelstatement_GEE = eval(parse(text=(paste0("y.vect ~ Mone +", paste0(paste0("COV.",conf.names),collapse = "+")))))
  }
  
  if ((method == "LMM")[1])
  {
    cat("Model: LMM\n")
    modelstatement = modelstatement_LMM
    doOne <- function(i, datarun, M){
      datarun$Mone <- M[,i]
      model <- try(lmer(modelstatement, data = datarun))
      if("try-error" %in% class(model)){
        b <- rep(NA, 3)
      } else {
        res=summary(model)$coefficients
        b <- as.numeric(res[2,1:3])
      }
      invisible(b)
    }
  }
  
  if ((method == "GEE")[1])
  {
    cat("Model: GEE\n")
    modelstatement = modelstatement_GEE
    doOne <- function(i, datarun, M){
    datarun$Mone <- M[,i]
    model <- try(geeglm(modelstatement, data = datarun, id = id.vect))
    if("try-error" %in% class(model)){
      b <- rep(NA, 3)
      } else {
        res=summary(model)$coefficients
        b <- as.numeric(res[2,1:3])
      }
      invisible(b)
    }
  }

  if (parallel == TRUE)
  {  
    cat(paste0("Running pre-ranking with ", ncore, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != ncore) registerDoParallel(ncore)
  } else {
    cat(paste0("Running pre-ranking with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }
  
  results <- foreach(n = idiv(L.M, chunks = ncore), M_chunk = iblkcol_lag(M, chunks = ncore),.combine = 'rbind', .packages = c('lme4',"geepack")) %dopar% {
    do.call('rbind',lapply( seq_len(n), doOne, datarun, M_chunk) )
  }
  results = data.frame(results)
  
  rownames(results) = M.names
  results = setNames(results, c("Estimate","Std.Err","Statistic"))
  results$p.value = 2*(1-pnorm(abs(results$Stat)))
  results$BH.FDR=p.adjust(results$p.value,"fdr")              # calculate Benjamini and Hochberg FDR
  results$Bonferroni=p.adjust(results$p.value,"bonferroni")   # calculate Bonferroni adjusted p-value
  results = results[order(results$p.value),]

  if (write == TRUE)
  {
    write.csv(results,paste0("Pre-ranking_",method,".csv"))
    cat(paste0("Pre-ranking results using ",method, " has been saved to ",getwd(),"   (",Sys.time(),")\n"))
  }
  cat(paste0("Done!    (",Sys.time(),")\n"))
  return(results)
}
