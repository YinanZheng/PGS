#' Sure Independent Screening S3 Class and Methods
#' 
#' \code{sis.obj} is a S3 class to store PGS sure independent screening results. Internally used by \code{sis} and \code{pgsfit}. \code{\link{sis.obj}} contains:
#' \itemize{
#'  \item{\code{name}: name (id) of variables}
#'  \item{\code{estimate}: beta estimate}
#'  \item{\code{stderr}: standard error}
#'  \item{\code{p.value}: raw p-value}
#'  \item{\code{bh.fdr}: Benjamini-Hochberg (BH) false discover rate (FDR)}
#'  \item{\code{bonferroni.p}: Bonferroni adjusted p-value}
#'  \item{\code{method}: independent screening method used in \code{sis}}
#' }
#'
#' @seealso see \code{\link{sis}} to run sure independent screening procedures and obtain \code{sis.obj} object.
#' 
#' @examples
#' ## Print sure independent screening object:
#' sis.obj
#' 
#' ## Plot sure independent screening object (histogram of p-values and Q-Q plot):
#' plot(sis.obj)
#' 
#' ## Return sure independent screening results:
#' coefficients = coef(sis.obj)
#' head(coefficients)
#' 
#' #For more information, please visit: https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function

sis.obj <- function(name, estimate, stderr, p.value, bh.fdr, bonferroni.p, method) {
  out <- list(name = name, 
              estimate = estimate, 
              stderr = stderr,
              p.value = p.value, 
              bh.fdr = bh.fdr, 
              bonferroni.p = bonferroni.p,
              method = method)
  class(out) <- "sis.obj"
  invisible(out)
}

#' @method print sis.obj
#' @rdname sis.obj
#' @export
print.sis.obj <- function(sis.obj) {
#   message("Use 'plot' to view distribution of p-values and Q-Q plot")
  cat(paste0("A total of ", length(sis.obj$name)," marks were tested using ", sis.obj$method,"."),'\n')
  cat(paste0(sum(sis.obj$bh.fdr < 0.05)," marks with BH.FDR < 0.05"),'\n')
}

#' @method plot sis.obj
#' @rdname sis.obj
#' @export
plot.sis.obj <- function(sis.obj) {
  p = sis.obj$p.value
  par(mfrow = c(1,2))
  hist(p, xlab = "p-value", breaks = 20, main = "Distribution of p-values", col = "light blue")
  o = -log10(sort(p,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=20,cex=1, main="Q-Q plot of p-values",
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

#' @method coef sis.obj
#' @rdname sis.obj
#' @export
coef.sis.obj <- function(sis.obj) {
  coef = data.frame(name = sis.obj$name,
                    estimate = sis.obj$estimate,
                    stderr = sis.obj$stderr,
                    p.value = sis.obj$p.value,
                    bh.fdr = sis.obj$bh.fdr,
                    bonferroni.p = sis.obj$bonferroni.p, row.names = sis.obj$name)
  return(coef)
}