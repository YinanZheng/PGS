% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sis.R
\name{sis}
\alias{sis}
\title{Sure Independent Screening}
\usage{
sis(y.vect, id.vect = NULL, M, COV = NULL, method = c("LMM", "GEE",
  "MLR"), corstr = "ar1", parallel = TRUE, ncore = detectCores(),
  write = FALSE)
}
\arguments{
\item{y.vect}{a vector of dependent variable.}

\item{id.vect}{a vector of subjuect ID.}

\item{M}{a data frame or matrix of genomic dataset. Rows represent samples, columns represent variables.}

\item{COV}{a data frame or matrix of covariates dataset.}

\item{method}{a character string specifying fitting method. For data contains >= 2 repeated measures, linear mixed-effect model (\code{"LMM"}) and generalized estimation equation (\code{"GEE"}) method are supported. For single measurement data, simple/multiple linear regression (\code{"MLR"}) is available. Default = \code{"LMM"}.}

\item{corstr}{a character string specifying the correlation structure when \code{method = "GEE"}. The following are permitted: '"independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'. Default = "ar1". (see \code{\link{geeglm}})}

\item{parallel}{logical. Enable parallel computing feature? Default = \code{TRUE}.}

\item{ncore}{number of cores to run parallel. Effective when paralle = \code{TRUE}. By default, max number of cores will be used.}

\item{write}{logical. Export screening results to csv file in the working directory if \code{TRUE}. Defaul = \code{FALSE}.}
}
\value{
sure independent screening results in a \code{\link{sis.obj}} object.
}
\description{
\code{sis} is used to conduct sure independent screening across high-dimensional variables.
}
\examples{
### Dataset preview
BJdata()

### Convert binary variables into factor type 
BJlung$gender = factor(BJlung$gender)
BJlung$heat = factor(BJlung$heat)
BJlung$cigwear = factor(BJlung$cigwear)

### Merge miRNA and lung function dataset
BJdata <- merge(BJmirna, BJlung, by=c("SID","WD"))

### Data must be sorted by study subject ID and multiple measurements indicator
BJdata <- BJdata[with(BJdata, order(SID, WD)), ]

### Extract dependent variable (lung function)
y.vect<-BJdata$FEV1

### Extract subjuect ID variable indicating repeated measures             
id.vect<-BJdata$SID        

### Extract microRNA data matrix   
M<-BJdata[,3:168]   

### Extract covariate data matrix          
COV<-BJdata[,170:179]
           
### LMM sure independent screening results
sis_LMM_par = sis(y.vect, id.vect, M, COV, method = "LMM")

### GEE sure independent screening results
sis_GEE_par = sis(y.vect, id.vect, M, COV, method = "GEE")

### Save the full site-by-site testing results into a csv file in current working directory
sis_LMM_par = sis(y.vect, id.vect, M, COV, method = "LMM", write = T)

sis_LMM_par        # print summary of sure independent screening results
plot(sis_LMM_par)  # plot histogram of raw p-values and Q-Q plot
coef(sis_LMM_par)  # return coefficients from sure independent screening results

#For more information, please visit: https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function

}
\seealso{
see \code{\link{pgsfit}} using the results from \code{\link{sis}} as input to run PGS; see \code{\link{sis.obj}} for class methods; see \code{\link{lmer}}, \code{\link{geeglm}}, and \code{\link{lm}} for more details on \code{"LMM"},\code{"GEE"}, and \code{"MLR"} methods, respectively.
}

