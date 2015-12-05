# Internal function: pgswrapper

pgswrapper <-function(pgsfit.obj, pm.ind)
{
  beta <- pgsfit.obj$beta.shrink[[ pm.ind ]]
  se <- sqrt(pgsfit.obj$var.sand[[ pm.ind ]])
  sel.names<-c(pgsfit.obj$sis.name[1:pgsfit.obj$Pm.vect[pm.ind]],pgsfit.obj$scale.info$COV_scale.dn)
  coef<-data.frame(Estimate = beta, Std.err = se, z = beta / se, P = 2*pnorm(abs(beta / se), lower.tail = F),
                   CI95.lower = beta - qnorm(0.975)*se, CI95.upper = beta + qnorm(0.975)*se)
  rownames(coef) = sel.names
  M_scale_sel = pgsfit.obj$scale.info$M_scale.ds[match(sel.names,pgsfit.obj$scale.info$M_scale.dn)]; M_scale_sel[is.na(M_scale_sel)] = pgsfit.obj$scale.info$COV_scale.ds
  coef$scale = M_scale_sel
  coef$y.scale = pgsfit.obj$scale.info$y.vect_scale.ds
  return(coef)
}