#' Penalized Generalized Estimation Equation with Grid Search
#' 
#' \code{pgsfit} is used to fit and determine the best results from penalized GEE method across pre-specified tunning parameter grid.
#' 
#' @param y.vect a vector of dependent variable.
#' @param id.vect a vector of subjuect ID.
#' @param M a data frame or matrix of genomic dataset. Rows represent samples, columns represent genomic marks.
#' @param COV a data frame or matrix of covariates dataset.
#' @param preRank.vect a vector of strings specifying pre-ranked names of genomic marks. see \code{preRank}.
#' @param Pm.vect a numeric vector of tunning parameter Pm, which specifies the number of top ranking genomic marks.       
#' @param lam.vect a numeric vector of tunning parameter lambda, which specifies the penalty parameter.
#' @param fold k-fold cross-validation in calculating grid error. Default = 10.
#' @param eps convergence threshold. By default iteration stops when beta estimation error <1e-5
#' @param iter.n maximum iteration number. Iteration will stop anyway even if the \code{eps} is not met. Default = 50.
#' @param corstr a character string specifying the working correlation structure. The following are permitted: independence (\code{"indep"}), exchangeable (\code{"exch"}), autoregressive(1) (\code{"ar1"}), and unstructured (\code{"un"}). Default = \code{"ar1"}.
#' @param parallel logical. enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel computation. Effective when \code{parallel} = \code{TRUE}. By default, max number of cores will be used.
#' @param seed an integer specifying seed for cross-validation. If not specified \code{pgsfit} will generate one.
#' 
#' @return Fitting and variable penalization results; a list object \code{pgsobj} consisting of:
#'   \item{coefficients}{estimates (shrinked) from the best model}
#'   \item{grid.err}{cross-validation error grid}
#'   \item{lam.sel.vect}{vector of selected lambda}
#'   \item{flag.stop.corr}{estimation error when iteration stop}
#'   \item{iter.n.corr}{iteration times}
#'   \item{Pm.vect}{working vector of tunning parameter Pm}
#'   \item{lam.vect}{working vector of tunning parameter lambda}
#' 
#' @seealso see \code{\link{preRank}} to obtain proper ranked genomic marks. see \code{\link{plotGrid}} to visualize and diagnose fitting results.
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
#' ### In the example we use linear mixed-effect model (default) for pre-ranking, ranked by p-values
#' prerank_LMM_par = preRank(y.vect, id.vect, M, COV)
#' preRank.vect = prerank_LMM_par
#' 
#' ### Initiate Pm sequence: Number of top ranking genomic marks
#' (Pm.vect<-c(seq(10,160,10)))
#' 
#' ### Initiate lambda sequence: Penalty parameter
#' (lam.vect<-seq(0.02, 0.3, 0.02))      
#' 
#' ### If your computer have multiple cores, it is recommended to enable parallel option (default)
#' PGSfit = pgsfit(y.vect, id.vect, M, COV, preRank.vect, Pm.vect, lam.vect, seed = 2015)

pgsfit<-function(
                  y.vect,                
                  id.vect,               
                  M,                     
                  COV = NULL,            
                  preRank.vect,	         
                  Pm.vect,               
                  lam.vect,              
                  fold = 10,             
                  eps = 1e-5,           
                  iter.n = 50,
                  corstr = "ar1",        
                  parallel = TRUE,       
                  ncore = detectCores(), 
                  seed = NULL
                )
{
  cat("Missing value check...\n")
  component_list = c("Dependent variable", "Genomic data matrix", "Covaraites")
  missingCheck = c( length( (y.vect_missing = which(is.na(y.vect)))) , sum(!is.na((M_missing = which(is.na(M), arr.ind = T)[1]))) , (COV_missing = 0) )
  if(!is.null(COV)) missingCheck[3] = sum(!is.na((COV_missing = which(is.na(COV), arr.ind = T)[1])))
  if ( any( as.logical(missingCheck) ) )
  {
    missing_index = unique(na.omit(c(y.vect_missing, M_missing, COV_missing)))
    cat(paste0("Missing value(s) exist in ", paste(component_list[which(as.logical(missingCheck) == TRUE)], collapse = " & " ), ". ID = ", paste(id.vect[missing_index], collapse = ","), ". Any observation(s) with missing value(s) will be dropped.\n"))
    confirm = TRUE
    while(confirm)
    {
      cat("Do you want to continue? (y/n)\n")
      stdinLetter = readLines(n=1)
      if(stdinLetter == "n") stop("PGS stopped. Please double check the data.") 
      else if(stdinLetter != "y") confirm = TRUE else {
        y.vect = y.vect[-missing_index]; M = M[-missing_index,]; COV = COV[-missing_index,]; id.vect = id.vect[-missing_index]
        cat("Missing observation(s) dropped\n")
        confirm = FALSE
      }
    }
  } else { cat("No missing value found.\n") }
  
  id.vect = sort(as.numeric(as.factor(as.character(id.vect)))) # Convert any id format to number.
  indGen_res = indGen_cpp(id.vect)
  
  if (!is.null(COV)) {COV <- as.matrix((model.matrix(~.,as.data.frame(COV)))[,-1]); n.COV = ncol(COV)} else {n.COV = 0}
  if (is.null(seed)) seed = round(as.numeric(Sys.time()))
  
  cat(paste0("The complete working dataset contains ", indGen_res$n, " individuals and a total of ", indGen_res$obs_n, " observations. (seed = ", seed,")\n"))

  L.Pm = length(Pm.vect)
  L.lam = length(lam.vect)
  
  n.marks = ncol(M)
  PGS.limit = min( (indGen_res$obs_n - n.COV - 1), n.marks) # The max number of p (exclude COV) that PGS can handles should not exceed number of samples
  
  preRank.vect = preRank.vect[1:PGS.limit]

  M = M[,match(preRank.vect, colnames(M))]
  
  M_scale = scale(M)
  M = data.frame(M_scale)
  M_center = data.frame(center = attributes(M_scale)[["scaled:center"]])
  M_scale = data.frame(scale = attributes(M_scale)[["scaled:scale"]])
  
  if (max(Pm.vect) > PGS.limit)
  {
    Pm.vect = c(Pm.vect[Pm.vect < PGS.limit], PGS.limit)
    warning(paste0("Besides the covariates, the maximum number of genomic marks that PGS can afford is ", indGen_res$obs_n - n.COV - 1,
                   ". The max number of working genomic marks = ", PGS.limit,
                   ". Your Pm.vect has been truncated to ", paste0(Pm.vect,collapse = " ")))
  }
    
  # Initialize rooms to store results
  grid.err = matrix(nrow = L.Pm, ncol = L.lam)   # Initiate error grid, column = lambda, row = Pm
  lam.sel.vect = rep(NA, L.Pm)  # Initiate optimized lambda for each Pm
  beta.shrink.corr.list = vector("list", L.Pm)   # Initiate penalized Beta (coefficient) estimates for each Pm
  var.sand.corr.list = vector("list", L.Pm)  # Initiate Variance 
  flag.stop.corr.vect = rep(NA, L.Pm) # Initiate stop flag
  iter.n.corr.vect = rep(NA, L.Pm) # Initiate iteration time
  hat.R.list = vector("list", L.Pm) # Initiate working correlation matrix

  if (parallel == TRUE)
  { 
    if(ncore > detectCores())
    {
      warning(paste0("You requested ", ncore, " cores. There are only ", detectCores()," in your machine!"))
      ncore = detectCores()
    }
    cat(paste0("Start running PGS with ", ncore, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != ncore) registerDoParallel(ncore)
  } else {
    cat(paste0("Start Running PGS with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }

#   if(!is.null(COV))
#   {
    res_par <- foreach(M_chunk = iblkcol_cum(M,Pm.vect), .packages = c("PGS") ) %dopar% {
      set.seed(seed)
      if(!is.null(COV)) {x.mat<-as.matrix(cbind(M_chunk, COV))} else {x.mat<-as.matrix(cbind(M_chunk))}
      p = ncol(x.mat)
      one_run_grid_cpp(y.vect, x.mat, id.vect, fold, p, lam.vect, eps, iter.n, corstr)
    }
#   } else {
#     res_par <- foreach(M_chunk = iblkcol_cum(M,Pm.vect), .packages = c("PGS") ) %dopar% {
#       set.seed(seed)
#       x.mat<-as.matrix(cbind(M_chunk))
#       p = ncol(x.mat)
#       one_run_grid_cpp(y.vect, x.mat, id.vect, fold, p, lam.vect, eps, iter.n, corstr)
#     }
#   }

  
  ## Summarize the results
  for(i in 1:L.Pm)
  {
    temp.run = res_par[[i]]
    grid.err[i,] <- temp.run$lam.cv.cor
    lam.sel.vect[i] <- temp.run$lam.sel.cor
    beta.shrink.corr.list[[i]] <- temp.run$beta.shrink.cor
    var.sand.corr.list[[i]] <- temp.run$var.sand.cor
    flag.stop.corr.vect[i] <- temp.run$flag.stop.cor
    iter.n.corr.vect[i] <- temp.run$iter.n.cor
    hat.R.list[[i]] <- temp.run$hat.R
  }
  
  rownames(grid.err) <- Pm.vect  
  colnames(grid.err) <- lam.vect
  
  best.ind<-which(grid.err==min(grid.err), arr.ind=T)

  beta<-beta.shrink.corr.list[[ best.ind[1] ]]
  varb<-var.sand.corr.list[[ best.ind[1] ]]
  se<-sqrt(varb)
  
  sel.names<-c(preRank.vect[1:Pm.vect[best.ind[1]]],colnames(COV))
  coefficients<-data.frame(Estimate = beta, Std.err = se, z = beta / se, P = 2*pnorm(abs(beta / se), lower.tail = F),
                           CI95.lower = beta - qnorm(0.975)*se, CI95.upper = beta + qnorm(0.975)*se)
  rownames(coefficients) = sel.names
  
  M_center_sel = M_center[sel.names,]; M_center_sel[is.na(M_center_sel)] = 0
  M_scale_sel = M_scale[sel.names,]; M_scale_sel[is.na(M_scale_sel)] = 1
  
  coefficients$center = M_center_sel
  coefficients$scale = M_scale_sel

  res = list(coefficients = coefficients,
             grid.err = grid.err,
             lam.sel.vect = lam.sel.vect,
             flag.stop.corr = flag.stop.corr.vect[best.ind[1]],
             iter.n.corr = iter.n.corr.vect[best.ind[1]],
             hat.R.corr = hat.R.list[[ best.ind[1] ]],
             Pm.vect = Pm.vect,
             lam.vect = lam.vect)
  
  cat(paste0("Done!    (",Sys.time(),")\n"))
  cat(paste0("The best selection results achieved using top ",Pm.vect[best.ind[1]]," marks with lambda = ",lam.vect[best.ind[2]]),"\n")
  
  return(res)
}
