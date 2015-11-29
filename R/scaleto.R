## Internal function: scaleto

scaleto <- function(dat)
{
  if (is.null(dat)) return(list(dn = NULL, d = NULL, dc = NULL, ds = NULL))
  dat_scale = scale(dat)
  dat_names = names(dat)
  if (any(class(dat) %in% c("matrix","data.frame","data.table"))) {dat_names = colnames(dat);dat = as.matrix(data.frame(dat_scale))} else {dat_names = names(dat);dat = as.numeric(dat_scale)}
  dat_center = as.numeric(attributes(dat_scale)[["scaled:center"]])
  dat_scale = as.numeric(attributes(dat_scale)[["scaled:scale"]])
  return(list(dn = dat_names, d = dat, dc = dat_center, ds = dat_scale))
}