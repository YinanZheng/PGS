## PGS: Penalized GEE with Grid Search
[![GitHub release](https://img.shields.io/badge/release-v0.2.0-blue.svg)](https://github.com/YinanZheng/PGS/releases)

*PGS* is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) package for association study of high-dimensional repeatedly-measured genomic data. 


## Installation 

In R console,
```r
## Install REMP
library(devtools)
install_github("YinanZheng/PGS",
               dependencies=TRUE)
               
## If SSL cert verification failure
library(RCurl)
library(httr)
set_config( config( ssl_verifypeer = 0L ) )
install_github("YinanZheng/PGS",
               dependencies=TRUE)
```

## Wiki & Examples:

[Wiki: PGS](https://github.com/YinanZheng/PGS/wiki)

[Example: Micro RNA expression and lung function](https://github.com/YinanZheng/PGS/wiki/Example:-miRNA-expression-and-lung-function)

##Citation:
1.	Zheng Y, Fei Z, Zhang W, Starren JB, Liu L, Baccarelli AA, Li Y, Hou L. PGS: a tool for association study of high-dimensional microRNA expression data with repeated measures. Bioinformatics. 2014;30(19):2802-7. doi: 10.1093/bioinformatics/btu396. PubMed PMID: 24947752; PubMed Central PMCID: PMC4173025. http://www.ncbi.nlm.nih.gov/pubmed/24947752

2.	Wang L, Zhou J, Qu A. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. Biometrics. 2012;68(2):353-60. doi: 10.1111/j.1541-0420.2011.01678.x. PubMed PMID: 21955051. http://www.ncbi.nlm.nih.gov/pubmed/21955051

##Contact Package Creator and Maintainer:
Yinan Zheng 

Email: y-zheng@northwestern.edu
