#' Sure Independent Screening
#' 
#' \code{sis} is used to conduct sure independent screening across high-dimensional variables.
#' 
#' @param y.vect a vector of dependent variable.
#' @param id.vect a vector of subjuect ID.
#' @param M a data frame or matrix of genomic dataset. Rows represent samples, columns represent variables.
#' @param COV a data frame or matrix of covariates dataset.
#' @param method a character string specifying fitting method. For data contains >= 2 repeated measures, linear mixed-effect model (\code{"LMM"}) and generalized estimation equation (\code{"GEE"}) method are supported. For single measurement data, simple/multiple linear regression (\code{"MLR"}) is available. Default = \code{"LMM"}.
#' @param corstr a character string specifying the correlation structure when \code{method = "GEE"}. The following are permitted: '"independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'. Default = "ar1". (see \code{\link{geeglm}})
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel. Effective when paralle = \code{TRUE}. By default, max number of cores will be used.
#' @param write logical. Export screening results to csv file in the working directory if \code{TRUE}. Defaul = \code{FALSE}.
#' 
#' @return sure independent screening results in a \code{\link{sis.obj}} object.
#' 
#' @seealso see \code{\link{pgsfit}} using the results from \code{\link{sis}} as input to run PGS; see \code{\link{sis.obj}} for class methods; see \code{\link{lmer}}, \code{\link{geeglm}}, and \code{\link{lm}} for more details on \code{"LMM"},\code{"GEE"}, and \code{"MLR"} methods, respectively.
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
#' ### LMM sure independent screening results
#' sis_LMM_par = sis(y.vect, id.vect, M, COV, method = "LMM")
#' 
#' ### GEE sure independent screening results
#' sis_GEE_par = sis(y.vect, id.vect, M, COV, method = "GEE")
#' 
#' ### Save the full site-by-site testing results into a csv file in current working directory
#' sis_LMM_par = sis(y.vect, id.vect, M, COV, method = "LMM", write = T)
#'
#' sis_LMM_par        # print summary of sure independent screening results
#' plot(sis_LMM_par)  # plot histogram of raw p-values and Q-Q plot
#' coef(sis_LMM_par)  # return coefficients from sure independent screening results
#' 
#' #For more information, please visit: https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function
#' 
#' @export
sis <- function(y.vect, id.vect=NULL, M, COV=NULL, method = c("LMM","GEE","MLR"), corstr = "ar1", parallel=TRUE, ncore = detectCores(), write=FALSE)
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

  L.M = ncol(M); M.names = colnames(M)
  res.COV = NULL
  
  if (is.null(COV))
  { 
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA)
    modelstatement_LMM = y.vect ~ Mone + (1|id.vect)
    modelstatement_MLR = modelstatement_GEE = y.vect ~ Mone
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~.,COV))[,-1]
    conf.names = colnames(COV)          
    datarun = data.frame(id.vect = id.vect, y.vect=y.vect, Mone = NA, COV = COV)
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
      if("try-error" %in% class(model)) b <- rep(NA, 3) else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3])}
      invisible(b)
    }
  }
  
  if (method == "LMM")
  {
    cat("Model: LMM\n")
    modelstatement = modelstatement_LMM
    doOne <- function(i, datarun, Mdat){
      datarun$Mone <- Mdat[,i]
      model <- try(lmer(modelstatement, data = datarun))
      if("try-error" %in% class(model)) b <- rep(NA, 3) else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3]) }
      invisible(b)
    }
  }
  
  if (method == "GEE")
  {
    cat("Model: GEE\n")
    modelstatement = modelstatement_GEE
    doOne <- function(i, datarun, Mdat){
    datarun$Mone <- Mdat[,i]
    model <- try(geeglm(modelstatement, data = datarun, id = id.vect, corstr = corstr))
    if("try-error" %in% class(model)) b <- rep(NA, 3) else { res=summary(model)$coefficients; b <- as.numeric(res[2,1:3])}
      invisible(b)
    }
  }

  if (parallel == TRUE & ncore > 1)
  {  
    if(ncore > detectCores())
    {
      cat(paste0("You requested ", ncore, " cores. There are only ", detectCores()," in your machine!"),'\n')
      ncore = detectCores()
    }
    cat(paste0("Running sure independent screening with ", ncore, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != ncore) registerDoParallel(ncore)
  } else {
    cat(paste0("Running sure independent screening with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }
  
  results <- foreach(n = idiv(L.M, chunks = ncore), M_chunk = iblkcol_lag(M, chunks = ncore),.combine = 'rbind', .packages = c('lme4',"geepack")) %dopar% {
    do.call('rbind',lapply( seq_len(n), doOne, datarun, M_chunk) )
  }
  
  results = data.frame(results)
  rownames(results) = M.names
  results = setNames(results, c("Estimate","Std.Err","Statistic"))
  results$p.value = 2*(1-pnorm(abs(results$Statistic)))
  results$BH.FDR=p.adjust(results$p.value,"fdr")              # calculate Benjamini and Hochberg FDR
  results$Bonferroni=p.adjust(results$p.value,"bonferroni")   # calculate Bonferroni adjusted p-value
  results = results[order(results$p.value), ]
  
  if (write == TRUE)
  {
    write.csv(results,paste0("Screening_",method,".csv"))
    cat(paste0("Sure independent screening results using ",method, ", ranked by p-values, has been saved to ",getwd(),"   (",Sys.time(),")\n"))
  }
  cat(paste0("Done!    (",Sys.time(),")\n"))
  
  res = sis.obj(rownames(results), results$Estimate, results$Std.Err, results$p.value, results$BH.FDR, results$Bonferroni, method)
  return(res)
}