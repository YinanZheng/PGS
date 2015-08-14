#PGS: Penalized GEE with Grid Search
Penalized GEE model for association study of high-dimensional genomic data with repeated (>=2) measures. 

##Install PGS R package:

> Binary package download:

[PGS for Windows](https://github.com/YinanZheng/PGS/releases/download/PGS_v0.0.3/PGS_0.0.3_Windows.zip)

[PGS for Mac OSX](https://github.com/YinanZheng/PGS/releases/download/PGS_v0.0.3/PGS_0.0.3_OSX.tgz)

[PGS for Linux](https://github.com/YinanZheng/PGS/releases/download/PGS_v0.0.3/PGS_0.0.3_LINUX.tar.gz)

> Run the following code in your R session:

_`# Check missing dependencies and install if necessary:`_

    list.of.dependencies <- c("doParallel", "lme4", "geepack", "ggplot2", "reshape2", "RColorBrewer")
    new.packages <- list.of.dependencies[!(list.of.dependencies %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages) else cat("Dependencies are ready!\n")
    
_`# Install PGS:`_

    install.packages("/YOUR PATH/PGS_Version_Platform", repo = NULL)

    library(PGS)
    
##Wiki & Examples:
https://github.com/YinanZheng/PGS/wiki

##Citation:
1.	Zheng Y, Fei Z, Zhang W, Starren JB, Liu L, Baccarelli AA, Li Y, Hou L. PGS: a tool for association study of high-dimensional microRNA expression data with repeated measures. Bioinformatics. 2014;30(19):2802-7. doi: 10.1093/bioinformatics/btu396. PubMed PMID: 24947752; PubMed Central PMCID: PMC4173025. http://www.ncbi.nlm.nih.gov/pubmed/24947752

2.	Wang L, Zhou J, Qu A. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. Biometrics. 2012;68(2):353-60. doi: 10.1111/j.1541-0420.2011.01678.x. PubMed PMID: 21955051. http://www.ncbi.nlm.nih.gov/pubmed/21955051




