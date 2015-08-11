## Pre-ranking genomic marks using generalized estimation equation (GEE)

prGEE <- function(y.vect, id.vect, M, COV=NULL, parallel=TRUE, write=FALSE)
{
  if (parallel == TRUE)
  {  
    num_cores <- detectCores()
    cat(paste0("Running pre-ranking with ", num_cores, " cores in parallel.\n"))
    options(mc.cores = num_cores)
  } else {
    num_cores <- 1
    cat(paste0("Running pre-ranking with single core.\n"))
    options(mc.cores = num_cores)
  }
  
  if (is.null(COV))
  {
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA)
    modelstatement = y.vect ~ Mone
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)          
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA, COV = COV)
    modelstatement = eval(parse(text=(paste0("y.vect ~ Mone +", paste0(paste0("COV.",conf.names),collapse = "+")))))
  }
  
  doOne <- function(i, datarun){
    datarun$Mone <- M[,i]
    model <- try(geeglm(modelstatement, data = datarun, id = id.vect))
    if("try-error" %in% class(model)){
      b <- rep(NA, 4)
    } else {
      res=summary(model)$coefficients
      b <- as.numeric(res[2,])
    }
    invisible(b)
  }
  
  mirnames = colnames(M) 
  
  results <- data.frame(do.call(rbind, mclapply( setNames(seq_len(ncol(M)), colnames(M)), doOne, datarun = datarun))) 
  results = setNames(results, c("Estimate","SE","Wald","pvalue"))
  results = results[order(results$pvalue),]
  results$BH.FDR=p.adjust(results$pvalue,"fdr")              # calculate Benjamini and Hochberg FDR
  results$Bonferroni=p.adjust(results$pvalue,"bonferroni")   # calculate Bonferroni adjusted p-value
  
  if (write == TRUE)
  {
    write.csv(ps_lmm,"Pre-ranking_GMM.csv")
    cat(paste0("Pre-ranking results using GMM has been saved to ",getwd(),".\n"))
  }
  return(results)
}