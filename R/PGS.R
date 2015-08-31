PGS<-function(
              y.vect,                #### Dependent variable vector
              id.vect,               #### Subjuect id variable
              M,                     #### Genomic data (each row represents samples, each column represents genomic marks)
              COV = NULL,            #### Covariates data
              preRank.vect,	         #### Preranking names of genomic marks. 1st marks is the most important and last marks is the least important.
              Pm.vect,               #### Pm vector, a vector of the number of top ranking genomic makrs	       
              lam.vect,              #### lambda vector, a vector of PGEE penalty parameter
              fold = 10,             #### k-fold cross-validation for calculating grid error. By default k=10
              eps = 10^-5,           #### Convergence threshold. By default iteration stops when estimation changes <10^-5
              iter.n = 50,           #### Maximum iternation number.
              corstr = "ar1",        #### Specify the working correlation structure
              parallel = TRUE,       #### parallel computing model
              ncore = detectCores(), #### number of cores to run parallel when paralle = TRUE. By default max number of cores will be used.
              seed = NULL            #### set seed number. If not specified PGS will generate one.
             )
{
  cat("Missing value check...\n")
  component_list = c("Dependent variable", "Genomic data matrix", "Covaraites")
  missingCheck = c( length( (y.vect_missing = which(is.na(y.vect)))) , sum(!is.na((M_missing = which(is.na(M), arr.ind = T)[1]))) , sum(!is.na((COV_missing = which(is.na(COV), arr.ind = T)[1]))) )
  if ( any( as.logical(missingCheck) ) )
  {
    missing_index = unique(na.omit(c(y.vect_missing, M_missing, COV_missing)))
    warning(paste0(">>> Missing value(s) exist in ", paste(component_list[which(missingCheck == TRUE)], collapse = " & " ), ". ID = ", paste(id.vect[missing_index], collapse = ","), ". Any observations with missing values will be dropped."))
    y.vect = y.vect[-missing_index]; M = M[-missing_index,]; COV = COV[-missing_index,]; id.vect = id.vect[-missing_index]
  } else {cat(">>> No missing value found.\n")}
  
  # Convert any id format to number.
  id.vect = sort(as.numeric(as.factor(as.character(id.vect))))
  indGen_res = indGen_cpp(id.vect)
  
  if (!is.null(COV))
  {
    COV <- as.matrix((model.matrix(~.,COV))[,-1])
  }
  
  cat(paste0("The complete working dataset contains ", indGen_res$n, " individuals and a total of ", indGen_res$obs_n, " observations. \n"))

  L.Pm = length(Pm.vect)
  L.lam = length(lam.vect)
  
  n.marks = ncol(M)
  n.COV = ncol(COV)
  PGS.limit = min( (indGen_res$obs_n - n.COV - 1), n.marks) # The max number of p (exclude COV) that PGS can handles should not exceed number of samples
  
  preRank.vect = preRank.vect[1:PGS.limit]
  M = M[,match(preRank.vect, colnames(M))]
  
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
  flag.stop.corr.vect = rep(NA, L.Pm)
  iter.n.corr.vect = rep(NA, L.Pm)

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

  res_par <- foreach(M_chunk = iblkcol_cum(M,Pm.vect), .packages = c("PGS") ) %dopar% {
    set.seed(seed)
    x.mat<-as.matrix(cbind(M_chunk, COV))
    p = ncol(x.mat)
    one_run_grid_cpp(y.vect, x.mat, id.vect, fold, p, lam.vect, eps, iter.n, corstr)
  }
  
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
  }
  
  rownames(grid.err) <- Pm.vect  
  colnames(grid.err) <- lam.vect
  
  best.ind<-which(grid.err==min(grid.err),arr.ind=T)

  beta<-beta.shrink.corr.list[[ best.ind[1] ]]
  varb<-var.sand.corr.list[[ best.ind[1] ]]
  se<-sqrt(varb)
  
  sel.names<-c(preRank.vect[1:Pm.vect[best.ind[1]]],colnames(COV))
  coefficients<-data.frame(Estimate = beta, Std.err = se, z = beta / se, P = 2*pnorm(abs(beta / se), lower.tail = F),
                           CI95.lower = beta - qnorm(0.975)*se, CI95.upper = beta + qnorm(0.975)*se)
  rownames(coefficients) = sel.names

  res = list(coefficients = coefficients,
             grid.err = grid.err,
             lam.sel.vect = lam.sel.vect,
             flag.stop.corr = flag.stop.corr.vect[best.ind[1]],
             iter.n.corr = iter.n.corr.vect[best.ind[1]],
             Pm.vect = Pm.vect,
             lam.vect = lam.vect)
  
  cat(paste0("Done!    (",Sys.time(),")\n"))
  cat(paste0(">>> The best selection results achieved using top ",Pm.vect[best.ind[1]]," biomarkers with lambda = ",lam.vect[best.ind[2]]),"\n")
  
  return(res)
}
