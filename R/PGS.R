PGS<-function(
              y.vect,            #### Dependent variable vector
              id.vect,           #### Subjuect id variable to indicate repeated measures
              M,                 #### Genomic data (each row represents samples, each column represents genomic marks)
              COV = NULL,        #### Covariates data
              PS,	               #### Prescreening results. Should have been properly ranked. The first column must be the biomarkers name and must be the same as the column names of biomarkers matrix)
              Pm.vect,           #### Pm vector, a vector of the number of top ranking biomarkers 	       
              lam.vect,          #### lambda vector, a vector of PGEE penalty parameter
              n,                 #### Sample size
              m,                 #### Number of repeated measures
              rho = 0.3,         ####
              fold = 10,         #### 10-fold cross-validation for calculating grid error
              sigcut = 10^-3,    #### If PGS estimate of biomarker is greater than 0.001, then the biomarker is selected 
              eps = 10^-10,      ####
              eps.stop = 0,      ####
              max.step = 50,     ####
              corr_str = "ar1",  #### Specify the correlation structure
              parallel = TRUE,   #### parallel computing model
              seed = NULL        #### seed
             )
{
  if (sum(is.na(y.vect)) + sum(is.na(id.vect)) + sum(is.na(M)) + sum(is.na(COV)) > 0)
    stop("Missing value(s) exist. Please make sure your data is complete and each individual have the same number of multiple measures.")
  
  n.marks = dim(M)[2]
  PGS.limit = min(n * m, n.marks)
  
  if (max(Pm.vect) >= n * m)
  {
    Pm.vect = Pm.vect[Pm.vect < PGS.limit]
    Pm.vect = c(Pm.vect, PGS.limit)
    warning(paste0("The maximum Pm that PGS can afford is n * m = ", n * m,
                   ". The total number of available genomic marks = ", dim(PS)[1],
                   ". Your Pm.vect has been truncated to ", paste0(Pm.vect,collapse = " ")))
  }
    
  if (max(Pm.vect) > n.marks)
  {
    Pm.vect = Pm.vect[Pm.vect <= n.marks]
    Pm.vect = c(Pm.vect, n.marks)
    warning(paste0("The total number of available genomic marks = ", dim(PS)[1],
                   ". Your Pm.vect has been truncated to ", paste0(Pm.vect,collapse = " "))) 
  }
  
  if(is.null( (PS.Morder = as.character(rownames(PS))) ))
  {stop("The row names of the pre-ranking results must be the name of genomic marks.")}
  
  if (is.null(seed)) seed = as.integer(Sys.time())  # Use current time as seed if not specified
  
  if (!is.null(COV))
  {COV <- as.matrix((model.matrix(~.,COV))[,-1])}
  
  M <- as.matrix(M)
  
  M.names=colnames(M)     # Extract genomic marks names
  conf.names=colnames(COV)   # Extract confounder names
  
  M.ps<-M[,match(PS.Morder,M.names)]   # Sort the genomic data based on the prescreening order
  
  L.Pm = length(Pm.vect)
  L.lam = length(lam.vect)
  
  ### Initialize rooms to store results
  grid.err = matrix(nrow = L.Pm, ncol = L.lam)   # Initiate error grid, column = lambda, row = Pm
  lam.sel = rep(NA,L.Pm)  # Initiate optimized lambda for each Pm
  beta.shrink.corr.list = vector("list", L.Pm)   # Initiate penalized Beta (coefficient) estimates for each Pm
  var.sand.corr.list = vector("list", L.Pm)  # Initiate Variance 
  flag.stop.corr.vect = rep(NA, L.Pm)    # Initiate stop flag 
  alpha.corr.list <- vector("list", L.Pm)  # Initiate correlation matrix       
  est.sigma2.corr.vect <- rep(NA, L.Pm)    
  
  if (parallel == TRUE)
  {  
    num_cores <- detectCores()
    cat(paste0("Start running PGS with ", num_cores, " cores in parallel...   (",Sys.time(),")\n"))
    if(getDoParWorkers() != num_cores) registerDoParallel(num_cores)
  } else {
    num_cores <- 1
    cat(paste0("Start Running PGS with single core...   (",Sys.time(),")\n"))
    registerDoSEQ()
  }

  res_par <- foreach(M_chunk = iblkcol_cum(M,Pm.vect), .packages = c("PGS") ) %dopar% {
    set.seed(seed)
    x.mat<-as.matrix(cbind(M_chunk, COV))
    one_run_grid_cpp(y.vect, x.mat, id.vect, n, m, ncol(x.mat), fold, lam.vect, rho, eps, eps.stop, max.step, corr_str)
  }
  
  ## Summarize the results
  for(i in 1:L.Pm)
  {
    temp.run = res_par[[i]]
    grid.err[i,] <- temp.run$lam.cv.vect/fold
    lam.sel[i] <- temp.run$lam.sel.corr
    beta.shrink.corr.list[[i]] <- temp.run$beta.shrink.corr
    var.sand.corr.list[[i]] <- temp.run$var.sand.corr
    flag.stop.corr.vect[i] <- temp.run$flag.stop.corr
    est.sigma2.corr.vect[i] <- temp.run$est.sigma2
    alpha.corr.list[[i]] <- temp.run$alpha.corr
  }
  
  rownames(grid.err) <- Pm.vect  
  colnames(grid.err) <- lam.vect
  
  best.ind<-which(grid.err==min(grid.err),arr.ind=T)

  beta<-beta.shrink.corr.list[[ best.ind[1] ]]
  varb<-var.sand.corr.list[[ best.ind[1] ]]
  se<-sqrt(varb)
  
  sel.names<-c(PS.Morder[1:Pm.vect[best.ind[1]]],conf.names)
  coefficients<-data.frame(Estimate = beta, Std.err = se, z = beta / se, P = 2*pnorm(abs(beta / se), lower.tail = F),
                           CI95.lower = beta - qnorm(0.975)*se, CI95.upper = beta + qnorm(0.975)*se)
  rownames(coefficients) = sel.names

  res = list(coefficients = coefficients,
             grid.err = grid.err,
             lam.sel = lam.sel,
             flag.stop.corr = flag.stop.corr.vect[best.ind[1]],
             alpha.corr = alpha.corr.list[[best.ind[1]]],
             est.sigma2.corr = est.sigma2.corr.vect[best.ind[1]],
             Pm.vect = Pm.vect,
             lam.vect = lam.vect)
  cat(paste0("Done!    (",Sys.time(),")\n"))
  
  cat(paste0(">>> The best selection results achieved using top ",Pm.vect[best.ind[1]]," biomarkers with lambda = ",lam.vect[best.ind[2]]),"\n")
  return(res)
}
