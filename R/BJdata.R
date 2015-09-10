#' Blood micro RNA Expression and Lung Function Data in Beijing Air Pollution Study
#' 
#' Blood micro RNA (miRNA) expression data of 166 miRNAs and lung function (FEV1) data of 120 individuals with two measurements in Beijing Air Pollution Study.
#' 
#' @usage 
#' \code{BJdata()}
#' 
#' @format 
#' \code{BJmirna} is a data frame of size 240 by 168.
#' \code{BJlung} is a data frame of size 240 by 13.
#' 
#' @details 
#' The data set \code{BJmirna} contains study subject ID (\code{SID}), indicator of examination days (\code{WD}), and 166 blood miRNA expression data. Blood samples were collection after work one the examination days.
#' The data set \code{BJlung} contains study subject id (\code{SID}), indicator of examination days (\code{WD}), lung function measured by forced expiratory volume in 1 second (\code{FEV1}), air pollution measured by particulate matter up to 2.5 micrometers in size (\code{pm25}), gender (\code{gender}), age (\code{age}), body mass index (\code{bmi}), use of central heating system (\code{heat}), smoke cigarette during examination day (\code{cigwear}), commute time (\code{waytime}), working hours (\code{workhr}), humidity measured by dew point (\code{dwep}), and temparature (\code{temp}). All variables were measured during the examination days. 
#' 
#' @references 
#' Hou L, Barupal J, Zhang W, et al. (2015) Environmental health perspectives.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/26068961}{PubMed})
#' 
#' Baccarelli AA, Zheng Y, Zhang X, et al. (2014) Particle and fibre toxicology.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/25272992}{PubMed})
#'
#' @examples
#' ### Dataset preview
#' BJdata()
#'
#' ### Convert binary variables into factor type. 
#' BJlung$gender = factor(BJlung$gender)
#' BJlung$heat = factor(BJlung$heat)
#' BJlung$cigwear = factor(BJlung$cigwear)
#' 
#' ### Merge miRNA and lung function dataset.
#' BJdata <- merge(BJmirna, BJlung, by=c("SID","WD"))
#' 
#' ### Sort the data by study subject id and multiple measurement indicator.
#' BJdata <- BJdata[with(BJdata, order(SID, WD)), ]

BJdata<-function()
{
  data(BJmirna)
  cat("Preview of PGS built-in dataset 'BJmirna' (truncated):\n")
  print(head(BJmirna[,1:10]))
  cat(paste0("'BJmirna' has ", nrow(BJmirna), " rows and ", ncol(BJmirna), " columns.\n"))
  
  cat("\n")
  
  data(BJlung)
  cat("Preview of PGS built-in dataset 'BJlung':\n")
  print(head(BJlung))
  cat(paste0("'BJlung' has ", nrow(BJlung), " rows and ", ncol(BJlung), " columns.\n"))
}

