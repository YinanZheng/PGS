#' Site-by-site Independent Pre-ranking for Genomic Marks
#' 
#' \code{preRank} is used to fit mixed-effect models or generalize estimation equations across genomic marks and rank by p-values.
#' 
#' @param y.vect a vector of dependent variable.
#' @param id.vect a vector of subjuect ID.
#' @param M a data frame or matrix of genomic dataset. Rows represent samples, columns represent genomic marks.
#' @param COV a data frame or matrix of covariates dataset.
#' @param method a character string specifying fitting method. Linear mixed-effect model (\code{"LMM"}) and generalized estimation equation (\code{"GEE"}) method are supported. Default = \code{"LMM"}.
#' @param corstr a character string specifying the correlation structure when \code{method = "GEE"}. The following are permitted: '"independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'. Default = "ar1". (see \code{\link{geeglm}})
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel. Effective when paralle = \code{TRUE}. By default, max number of cores will be used.
#' @param write logical. Export ranking results to csv file in the working directory if \code{TRUE}. Defaul = \code{FALSE}.
#' 
#' @return a list object which contains a vector of strings specifying pre-ranked names of genomic marks (\code{$names}), their corresponding estimates (\code{$estimates}), and related scaling information. If covariates are included, the estimates for each covariates are available.
#' 
#' @seealso see \code{\link{pgsfit}} using the results from \code{\link{preRank}} as input to run PGS. For more details on LMM, GEE, and MLR, see \code{\link{lmer}}, \code{\link{geeglm}}, and \code{\link{lm}}.
#'
#' @examples
#' ### Dataset preview
#' BJdata()
#'
#' ### Convert binary variables into factor type 
#' BJlung$gender = factor(BJlung$gender)
#' BJlung$heat = factor(BJlung$heat)
#' BJlung$cigwear = factor(BJlung$cigwear)
#' 
#' ### Merge miRNA and lung function dataset
#' BJdata <- merge(BJmirna, BJlung, by=c("SID","WD"))
#' 
#' ### Data must be sorted by study subject ID and multiple measurements indicator
#' BJdata <- BJdata[with(BJdata, order(SID, WD)), ]
#' 
#' ### Extract dependent variable (lung function)
#' y.vect<-BJdata$FEV1
#' 
#' ### Extract subjuect ID variable indicating repeated measures             
#' id.vect<-BJdata$SID        
#' 
#' ### Extract microRNA data matrix   
#' M<-BJdata[,3:168]   
#' 
#' ### Extract covariate data matrix          
#' COV<-BJdata[,170:179]
#'            
#' ### LMM site-by-site pre-ranking results
#' prerank_LMM_par = preRank(y.vect, id.vect, M, COV, method = "LMM")

#' ### GEE site-by-site pre-ranking results
#' prerank_GEE_par = preRank(y.vect, id.vect, M, COV, method = "GEE")

#' ### Save the full site-by-site testing results into a csv file in current working directory
#' prerank_LMM_par = preRank(y.vect, id.vect, M, COV, method = "LMM", write = T)

preRank <- function(y.vect, id.vect=NULL, M, COV=NULL, method = c("LMM","GEE","MLR"), corstr = "ar1", parallel=TRUE, ncore = detectCores(), write=FALSE)
{
  if (is.null(id.vect) | all(table(id.vect) == 1))
  {
    cat("No repeated measures detected. Applying single/multiple linear regression model (MLR).\n")
    id.vect = 1:length(y.vect) 
    method = "MLR"
  } else if (length(method) > 1) 
  {
    method = "LMM"  # By default the model is LMM if m >=2
  } else if (!method %in% c("LMM","GEE","MLR")) 
    stop("The method can be 'LMM' (linear mixed-effect model), 'GEE' (generalized estimation equation), or 'MLR' (multiple linear regression).")

  y.vect_scale = scaleto(y.vect)
  M_scale = scaleto(M)
  COV_scale = scaleto(NULL)
  L.M = ncol(M)
  res.COV = NULL
  
  if (is.null(COV))
  { 
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect_scale$d, Mone = NA)
    modelstatement_LMM = y.vect ~ Mone + (1|id.vect)
    modelstatement_MLR = modelstatement_GEE = y.vect ~ Mone
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    COV_scale = scaleto(COV)
    conf.names = colnames(COV)          
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect_scale$d, Mone = NA, COV = COV_scale$d)
    modelstatement_LMM = eval(parse(text=(paste0("y.vect ~ Mone +", paste0(paste0("COV.",conf.names),collapse = "+"), "+ (1|id.vect)"))))
    modelstatement_MLR = modelstatement_GEE = eval(parse(text=(paste0("y.vect ~ Mone +", paste0(paste0("COV.",conf.names),collapse = "+")))))
    modelstatement_LMM.COV = eval(parse(text=(paste0("y.vect ~ ", paste0(paste0("COV.",conf.names),collapse = "+"), "+ (1|id.vect)"))))
    modelstatement_MLR.COV = modelstatement_GEE.COV = eval(parse(text=(paste0("y.vect ~ ", paste0(paste0("COV.",conf.names),collapse = "+")))))  
  }
  
  if (method == "MLR")
  {
    cat("Model: MLR\n")
    modelstatement = modelstatement_MLR
    doOne <- function(i, datarun, Mdat){
      datarun$Mone <- Mdat[,i]
      model <- try(lm(modelstatement, data = datarun))
      if("try-error" %in% class(model)){
        b <- rep(NA, 3)
      } else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3])}
      invisible(b)
    }
    if(!is.null(COV)) res.COV = as.numeric(summary(lm(modelstatement_MLR.COV,data = subset(datarun, select = -Mone)))$coefficients[-1,1])
  }
  
  if (method == "LMM")
  {
    cat("Model: LMM\n")
    modelstatement = modelstatement_LMM
    doOne <- function(i, datarun, Mdat){
      datarun$Mone <- Mdat[,i]
      model <- try(lmer(modelstatement, data = datarun))
      if("try-error" %in% class(model)){
        b <- rep(NA, 3)
      } else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3]) }
      invisible(b)
    }
    if(!is.null(COV)) res.COV = as.numeric(summary(lmer(modelstatement_LMM.COV,data = subset(datarun, select = -Mone)))$coefficients[-1,1])
  }
  
  if (method == "GEE")
  {
    cat("Model: GEE\n")
    modelstatement = modelstatement_GEE
    doOne <- function(i, datarun, Mdat){
    datarun$Mone <- Mdat[,i]
    model <- try(geeglm(modelstatement, data = datarun, id = id.vect, corstr = corstr))
    if("try-error" %in% class(model)){
      b <- rep(NA, 3)
      } else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3])}
      invisible(b)
    }
    if(!is.null(COV)) res.COV = as.numeric(summary(geeglm(modelstatement_GEE.COV,data = subset(datarun, select = -Mone), id = id.vect, corstr = corstr))$coefficients[-1,1])
  }

  if (parallel == TRUE)
  {  
    if(ncore > detectCores())
    {
      cat(paste0("You requested ", ncore, " cores. There are only ", detectCores()," in your machine!"),'\n')
      ncore = detectCores()
    }
    cat(paste0("Running pre-ranking with ", ncore, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != ncore) registerDoParallel(ncore)
  } else {
    cat(paste0("Running pre-ranking with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }
  
  results <- foreach(n = idiv(L.M, chunks = ncore), M_chunk = iblkcol_lag(M_scale$d, chunks = ncore),.combine = 'rbind', .packages = c('lme4',"geepack")) %dopar% {
    do.call('rbind',lapply( seq_len(n), doOne, datarun, M_chunk) )
  }
  
  results = data.frame(results)
  rownames(results) = M_scale$dn
  results = setNames(results, c("Estimate.scaled","Std.Err","Statistic"))
  results$p.value = 2*(1-pnorm(abs(results$Statistic)))
  results$BH.FDR=p.adjust(results$p.value,"fdr")              # calculate Benjamini and Hochberg FDR
  results$Bonferroni=p.adjust(results$p.value,"bonferroni")   # calculate Bonferroni adjusted p-value
  results$center = M_scale$dc
  results$scale = M_scale$ds
  results$y.center = y.vect_scale$dc
  results$y.scale = y.vect_scale$ds
  results = data.frame(Estimate.originalScale = with(results, Estimate.scaled * y.scale / scale), results)
  
  results = results[order(results$p.value), ]
  
  if (write == TRUE)
  {
    write.csv(results,paste0("Pre-ranking_",method,".csv"))
    cat(paste0("Pre-ranking results using ",method, ", ranked by p-values, has been saved to ",getwd(),"   (",Sys.time(),")\n"))
  }
  cat(paste0("Done!    (",Sys.time(),")\n"))
  
  return(list(names = rownames(results), est = results$Estimate.scaled, est.original = results$Estimate.originalScale, est_center = results$center, est_scale = results$scale, 
              names.cov = COV_scale$dn, est.cov = res.COV, est.cov_center = COV_scale$dc, est.cov_scale = COV_scale$ds,
              y.center = y.vect_scale$dc, y.scale = y.vect_scale$ds))
}