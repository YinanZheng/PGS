#' PGS S3 Class and Methods
#' 
#' \code{pgsfit.obj} is a S3 class to store PGS fitting results. Internally used by \code{pgsfit}. \code{\link{pgsfit.obj}} contains:
#' \itemize{
#'   \item{\code{grid.err}: cross-validation error grid}
#'   \item{\code{lam.sel.vect}: vector of selected lambda for each Pm}
#'   \item{\code{beta.shrink}: list of beta after shrinkage for each Pm}
#'   \item{\code{var.sand}: list of sandwich variance of beta for each Pm}
#'   \item{\code{hat.R}: list of estimated working correlation matrix for each Pm}
#'   \item{\code{convergenceError}: convergence error when iteration stopped for each Pm}
#'   \item{\code{iterationNumber}: iteration times when converged for each Pm}
#'   \item{\code{which.best}: index of the best Pm and lambda among convergent results}
#'   \item{\code{which.bestGlobal}: index of the best Pm and lambda among all results (including non-convergent results)}
#'   \item{\code{convergenceThreshold}: threshold of convergence error (i.e. \code{eps})}
#'   \item{\code{maxIteration}: maximum iteration number allowed (i.e. \code{iter.n})}
#'   \item{\code{Pm.vect}: working vector of tunning parameter Pm}
#'   \item{\code{lam.vect}: working vector of tunning parameter lambda}
#'   \item{\code{scale.info}: scaling information of y.vect, M, and COV to transfer beta estimate back to original scale}
#'   \item{\code{sis.name}: name (id) of genomic mark from sure independent screening results, truncated at \code{pm.max}}
#'   \item{\code{convergeMessage}: "Converged!", "Not converged!", or "Conditionally converged!}
#' }   
#' @seealso see \code{\link{pgsfit}} to run PGS algorithm and obtain \code{pgsfit.obj} object.
#' 
#' @examples
#' ## Print PGS object:
#' pgsfit.obj
#' 
#' ## Plot PGS object (a heat map visualizes the grid search errors):
#' plot(pgsfit.obj)
#' 
#' #Other parameters:
#' #IQR.times: cross-validation errors greater than IQR.times * IQR of the errors will be replaced by "NA" (represented by grey blocks). Default = 1.5.
#' #text.size: size of the text in grid. Default = 4.
#' #xlab.size: size of the x-axis labels. Default = 14.
#' #ylab.size: size of the y-axis labels. Default = 14.
#' 
#' ## Return coefficients from PGS object:
#' coef(pgsfit.obj)
#' 
#' #Other parameters:
#' #whcih: an interger between 1 to pm.n specifying results at which Pm level to be returned. Defaul = NULL. If NULL, return the best results.
#' #p.threshold: threshold of p-values to filter out non-significant variables. Default = 0.05.
#' #nonzero: logical. If TRUE, only variables with non-zero beta results are returned. Default = TRUE.
#'
#' #For more information, please visit: https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function


pgsfit.obj <- function(grid.err, lam.sel.vect, beta.shrink.corr.list, var.sand.corr.list, hat.R.list, flag.stop.corr.vect, iter.n.corr.vect, best.ind, bestall.ind, eps, iter.n, Pm.vect, lam.vect, scale.info, sis.dn, ConvergeMessage) {
  out <- list(grid.err = grid.err,
              lam.sel = lam.sel.vect,
              beta.shrink = beta.shrink.corr.list,
              var.sand = var.sand.corr.list,
              hat.R = hat.R.list,
              convergenceError = flag.stop.corr.vect,
              iterationNumber = iter.n.corr.vect,
              which.best = best.ind,
              which.bestGlobal = bestall.ind,
              convergenceThreshold = eps,           
              maxIteration = iter.n,
              Pm.vect = Pm.vect,
              lam.vect = lam.vect,
              scale.info = scale.info,
              sis.name = sis.dn,
              convergeMessage = ConvergeMessage)
  class(out) <- "pgsfit.obj"
  invisible(out)
}

#' @method print pgsfit.obj
#' @rdname pgsfit.obj
#' @export
print.pgsfit.obj <- function(pgsfit.obj) {
#   message("Use 'plot' method to plot grid search results; use 'coef' method to return coefficient matrix.")
  if (pgsfit.obj$convergeMessage == "Conditionally converged!") cat(paste0("Warning: smallest error achieved at ", pgsfit.obj$Pm.vect[pgsfit.obj$which.bestGlobal[1]], " and -ln(lambda) = ",round(-log(pgsfit.obj$lam.vect[pgsfit.obj$which.bestGlobal[2]]),2), " but PGS did not converge within ", pgsfit.obj$maxIteration, " iternations (convergence error = ", pgsfit.obj$convergenceError[pgsfit.obj$which.bestGlobal[1]]," > the preset threshold: ", pgsfit.obj$convergenceThreshold), ". PGS will continue searching among the converged results.\n")
  if (pgsfit.obj$convergeMessage == "Not converged!") cat(paste0("Warning: no converged result!"),"\n")
  if (pgsfit.obj$convergeMessage == "Converged!") cat(paste0("PGS converged!"),"\n")
  cat(paste0("The best selection results achieved using top ",pgsfit.obj$Pm.vect[pgsfit.obj$which.best[1]]," marks with -ln(lambda) = ", round(-log(pgsfit.obj$lam.vect[pgsfit.obj$which.best[2]]),2)),"\n")
}

#' @method plot pgsfit.obj
#' @rdname pgsfit.obj
#' @export
plot.pgsfit.obj <- function(pgsfit.obj, IQR.times = 1.5, text.size = 4, xlab.size = 14, ylab.size = 14)
{
   col_low = "#40004B"; col_high = "#FFFFFF"
   L.Pm = length(pgsfit.obj$Pm.vect); n.COV = length(pgsfit.obj$scale.info$COV_scale.dn)
   
   grid.err = pgsfit.obj$grid.err*10
   grid.err.m = melt(grid.err,id.vars = pgsfit.obj$Pm.vect)
   colnames(grid.err.m) = c("Pm","Lambda","value")
   grid.err.m$Pm <- factor(grid.err.m$Pm,levels = unique(grid.err.m$Pm))
   lambda.log = sprintf("%.2f", round(-log(grid.err.m$Lambda),2))
   
   converge.ind = (which(pgsfit.obj$iterationNumber == pgsfit.obj$maxIteration))
   minlambda.ind = apply(grid.err,1,which.min)
   minerrGlobal.text = as.character(round(grid.err[pgsfit.obj$which.bestGlobal],2))
   minerr.text = as.character(round(grid.err[pgsfit.obj$which.best],2))
   minerr.textall = as.character(round(apply(grid.err,1,min),2))
   
   grid.err.m$Lambda <- factor(lambda.log, levels = unique(lambda.log))
   grid.err.m$small = grid.err.m$plus = grid.err.m$converge = grid.err.m$sig = grid.err.m$nonzero<- ""
   
   if(pgsfit.obj$convergeMessage == "Converged!")
   {
     grid.err.m$small[(minlambda.ind-1) * L.Pm + c(1:L.Pm)] <- minerr.textall
     grid.err.m$plus[(pgsfit.obj$which.best[2]-1) * L.Pm + pgsfit.obj$which.best[1]] <- minerrGlobal.text
   } else {
     grid.err.m$small[(minlambda.ind[-converge.ind]-1) * L.Pm + c(1:L.Pm)[-converge.ind]] <- minerr.textall[-converge.ind]
     grid.err.m$plus[(pgsfit.obj$which.best[2]-1) * L.Pm + pgsfit.obj$which.best[1]] <- minerr.text
     grid.err.m$converge[(minlambda.ind[converge.ind]-1) * L.Pm + converge.ind] <- minerr.textall[converge.ind]
   }
   grid.err.m$plusGlobal <- ""
   grid.err.m$plusGlobal[(pgsfit.obj$which.bestGlobal[2]-1) * L.Pm + pgsfit.obj$which.bestGlobal[1]] <- minerrGlobal.text
   if(any(grid.err.m$converge!="")) message("Note: some results are not converged (text in red).")

   outlier.ind = which(grid.err.m$value > min(grid.err.m$value) + IQR.times * IQR(grid.err.m$value, na.rm = T))
   grid.err.m$value[outlier.ind] = NA
   
   nonzero.vect = sig.vect = rep(NULL, L.Pm)
   for(i in 1:L.Pm)
   {
     temp = pgswrapper(pgsfit.obj, i); mark.ind = 1:(pgsfit.obj$Pm.vect[i] - n.COV)
     nonzero.vect[i] = sum(temp$Estimate[mark.ind]!=0); sig.vect[i] = sum(temp$P[mark.ind]<0.05)
   }
   
   sig.mat = data.frame(Pm = pgsfit.obj$Pm.vect, Lambda = "p<0.05", value = NA, small = "", plus = "", converge = "", plusGlobal = "", sig = as.character(sig.vect), nonzero = "")
   grid.err.m = rbind(grid.err.m, sig.mat)
   
   nonzero.mat = data.frame(Pm = pgsfit.obj$Pm.vect, Lambda = "Nonzero", value = NA, small = "", plus = "", converge = "", plusGlobal = "", sig = "", nonzero = as.character(nonzero.vect))
   grid.err.m = rbind(grid.err.m, nonzero.mat)
   

   plot.title = "Cross-validation Error Grid"
   plot.subtitle =  bquote("Minimal error achieved at" ~ P[m] ~ "=" ~ .(rownames(grid.err)[pgsfit.obj$which.best[1]]) ~ "and -ln(" * lambda * ") =" ~ .(as.character(round(-log(as.numeric(colnames(grid.err)[pgsfit.obj$which.best[2]])),2))) )
   qplot(x=Lambda, y=Pm, data=grid.err.m, fill = value, geom="tile", xlab = bquote("-ln(" * lambda * ")"), ylab = bquote(P[m])) + 
      scale_fill_gradient("CV error (x10)",low=col_low, high=col_high) + 
      geom_text(aes(label=small, angle = 90), color="white", size=text.size) + 
      geom_text(aes(label=plusGlobal, angle = 90), color="yellow", size=text.size) + 
      geom_text(aes(label=plus, angle = 90), color="yellow", size=text.size) + 
      geom_text(aes(label=converge, angle = 90), color="red", size=text.size) + 
      geom_text(aes(label=sig), color="white", size=text.size) + 
      geom_text(aes(label=nonzero), color="white", size=text.size) + 
      ggtitle(bquote(atop(.(plot.title), atop(.(plot.subtitle), "")))) +
      theme(plot.title = element_text(size = rel(1.7))) +
      theme(axis.title = element_text(size = rel(1.5))) +
      theme(axis.text = element_text(size = rel(1.2))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = xlab.size)) +
      theme(axis.text.y = element_text(size = ylab.size)) +
      theme(legend.position="bottom") 
}

#' @method coef pgsfit.obj
#' @rdname pgsfit.obj
#' @export
coef.pgsfit.obj <- function(pgsfit.obj, pm.ind = NULL, p.threshold = 0.05, nonzero = TRUE)
{
  Pm.N = length(pgsfit.obj$Pm.vect)
  if(is.null(pm.ind)) {pm.ind = pgsfit.obj$which.best[1]; cat(paste0("Return the best fitting results using top ", pgsfit.obj$Pm.vect[pm.ind], " marks."),'\n'); coef = pgswrapper(pgsfit.obj, pm.ind)} else {
    if (pm.ind > Pm.N | pm.ind < 1) stop(paste0("pm.ind must between 1 to ",length(pgsfit.obj$Pm.vect),"!")) else {cat("Return fitting results using top", pgsfit.obj$Pm.vect[pm.ind], "marks.",'\n'); coef = pgswrapper(pgsfit.obj, pm.ind)} }
  
  if (nonzero) coef = subset(coef,Estimate != 0)
  coef = subset(coef, P < p.threshold)
  coef$Estimate = with(coef, Estimate * y.scale / scale)
  coef$CI95.lower = with(coef, CI95.lower * y.scale / scale)
  coef$CI95.upper = with(coef, CI95.upper * y.scale / scale)
  coef = subset(coef, select = c(Estimate, CI95.lower, CI95.upper, P))
  return(coef)
}