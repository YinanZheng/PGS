#' Visualize the Cross-Validation Error Grid with Heat Map
#' 
#' \code{plotGrid} is used to visualize the cross-validation (CV) error grid using a heat map.
#' 
#' @param pgsobj a \code{pgsfit} object.
#' @param IQR.times cross-validation errors greater than \code{IQR.times} * IQR of the errors will be replaced by "NA" (represented by grey blocks). Default = 1.5.
#' @param text.size size of the text in grid. Default = 4.
#' @param xaxislab.size size of the x-axis labels. Default = 14.
#' @param yaxislab.size size of the y-axis labels. Default = 14.
#' 
#' @seealso see \code{\link{pgsfit}} to run PGS and obtain \code{pgsobj} object.
#' 
#' @examples
#' PGSfit = pgsfit(y.vect, id.vect, M, COV, preRank_LMM_par, Pm.vect, lam.vect, seed = 1)
#' plotGrid(PGSfit)

plotGrid <- function(pgsobj, IQR.times = 1.5, text.size = 4, xaxislab.size = 14, yaxislab.size = 14)
{
  col_low = "#40004B"
  col_high = "#FFFFFF"
  
  grid.err = pgsobj$grid.err*10
  grid.err.m = melt(grid.err,id.vars = pgsobj$Pm.vect)
  colnames(grid.err.m) = c("Pm","Lambda","value")
  grid.err.m$Pm <- factor(grid.err.m$Pm,levels = unique(grid.err.m$Pm))
  lambda.log = sprintf("%.2f", round(-log(grid.err.m$Lambda),2))
  
  L.Pm = length(pgsobj$Pm.vect) 
  converge.ind = (which(pgsobj$iterationNumber == pgsobj$maxIteration))
  minlambda.ind = apply(grid.err,1,which.min)
  minerr.text = as.character(round((apply(grid.err,1,min)),2))
  
  grid.err.m$Lambda <- factor(lambda.log, levels = unique(lambda.log))
  grid.err.m$small <- ""
  grid.err.m$plus <- ""
  grid.err.m$converge <- ""
  
  if(length(converge.ind) == 0)
  {
    grid.err.m$small[(minlambda.ind-1) * L.Pm + c(1:L.Pm)] <- minerr.text
    grid.err.m$plus[(pgsobj$which.bestLambda-1) * L.Pm + pgsobj$which.bestPm] <- as.character(round(min(grid.err),2))
  } else {
    grid.err.m$small[(minlambda.ind[-converge.ind]-1) * L.Pm + c(1:L.Pm)[-converge.ind]] <- minerr.text[-converge.ind]
    grid.err.m$plus[(pgsobj$which.bestLambda-1) * L.Pm + pgsobj$which.bestPm] <- as.character(round(min(grid.err[-converge.ind,]),2))
    grid.err.m$converge[(minlambda.ind[converge.ind]-1) * L.Pm + converge.ind] <- minerr.text[converge.ind]
  }
  grid.err.m$plusGlobal <- ""
  grid.err.m$plusGlobal[(pgsobj$which.bestLambdaGlobal-1) * L.Pm + pgsobj$which.bestPmGlobal] <- as.character(round(min(grid.err),2))
  
  #Detect outlier
  outlier.ind = which(grid.err.m$value > min(grid.err.m$value) + IQR.times * IQR(grid.err.m$value))
  grid.err.m$value[outlier.ind] = NA
  
  plot.title = "Cross-validation Error Grid"
  plot.subtitle =  bquote("Minimal error achieved at" ~ P[m] ~ "=" ~ .(rownames(grid.err)[pgsobj$which.bestPm]) ~
                                      "and -ln(" * lambda * ") =" ~ .(as.character(round(-log(as.numeric(colnames(grid.err)[pgsobj$which.bestLambda])),2))) )
  qplot(x=Lambda, y=Pm, data=grid.err.m, fill = value, geom="tile", xlab = bquote("-ln(" * lambda * ")"), ylab = bquote(P[m])) + 
    scale_fill_gradient("CV error (x10)",low=col_low, high=col_high) + 
    geom_text(aes(label=small, angle = 90), color="white", size=text.size) + 
    geom_text(aes(label=plusGlobal, angle = 90), color="yellow", size=text.size) + 
    geom_text(aes(label=plus, angle = 90), color="yellow", size=text.size) + 
    geom_text(aes(label=converge, angle = 90), color="red", size=text.size) + 
    ggtitle(bquote(atop(.(plot.title), atop(.(plot.subtitle), "")))) +
    theme(plot.title = element_text(size = rel(1.7))) +
    theme(axis.title = element_text(size = rel(1.5))) +
    theme(axis.text = element_text(size = rel(1.2))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = xaxislab.size)) +
    theme(axis.text.y = element_text(size = yaxislab.size)) +
    theme(legend.position="bottom") 
  }