#' Visualize the Cross-Validation Error Grid with Heat Map
#' 
#' \code{plotGrid} is used to visualize the cross-validation (CV) error grid using a heat map.
#' 
#' @param pgsobj a \code{pgsfit} object.
#' @param colrange a vector of color strings specifying the color range for heat map. The first color will be used for the lowest CV error; the second color will be used for the highest CV error. If not specified, \code{plotGrid} will generate dark to light purple color scheme by default.
#'
#' @seealso see \code{\link{pgsfit}} to run PGS and obtain \code{pgsobj} object.
#' 
#' @examples
#' PGSfit = pgsfit(y.vect, id.vect, M, COV, preRank.vect, Pm.vect, lam.vect, seed = 1)
#' plotGrid(PGSfit)

plotGrid <- function(pgsobj,colrange = NULL)
{
  if (is.null(colrange))
  {
    mypalette = brewer.pal(11,"PRGn")
    col_low = mypalette[1]
    col_high = mypalette[5]
  } else {
    col_low = colrange[1]
    col_high = colrange[2]
  }
  
  grid.err = pgsobj$grid.err
  best.ind=which(grid.err==min(grid.err),arr.ind=T)         # Locate the smallest error in grid
  grid.err.m = melt(grid.err,id.vars = pgsobj$Pm.vect)
  colnames(grid.err.m) = c("Pm","Lambda","value")
  grid.err.m$Pm <- factor(grid.err.m$Pm,levels = unique(grid.err.m$Pm))
  grid.err.m$Lambda <- factor(grid.err.m$Lambda, levels = unique(grid.err.m$Lambda))
  grid.err.m$small <- ""
  grid.err.m$small[(apply(grid.err,1,which.min)-1) * length(pgsobj$Pm.vect) + 1:length(pgsobj$Pm.vect) ] <- "+"
  grid.err.m$plus <- ""
  grid.err.m$plus[which.min(grid.err.m$value)] <- "+"
  
  plot.title = "Cross-validation Error Grid"
  plot.subtitle = paste0("Minimal error with Pm=",rownames(grid.err)[best.ind[1]]," Lambda=",colnames(grid.err)[best.ind[2]])
  qplot(x=Lambda, y=Pm, data=grid.err.m, fill = value, geom="tile") + 
    scale_fill_gradient(low=col_low, high=col_high) + 
    geom_text(aes(label=small), color="white", size=7) + 
    geom_text(aes(label=plus), color="yellow", size=12) + 
    ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), "")))) +
    theme(plot.title = element_text(size = rel(1.7))) +
    theme(axis.title = element_text(size = rel(1.5))) +
    theme(axis.text = element_text(size = rel(1.2))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(fill=guide_legend(title="CV error"))
}
