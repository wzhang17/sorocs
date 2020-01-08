#' The example data is meant to represent the dataset supplied by the Physician Reliability Study (PRS), 
#' which is explored in Section 5 of the paper. 
#' The 'sampledata' file contains the following variables:
#'
#' @format A data frame with 129 rows and following variables:
#' \describe{
#'   \item{STUDYID}{Subject id}
#'   \item{logREscoremean1}{log value of mean of four Regional Experts¡¯ scores at setting 1}
#'   \item{logREscoremean2}{log value of mean of four Regional Experts¡¯ scores at setting 2}
#'   \item{TN1}{sum of the IE(International Experts)'s diagnoses for positive disease at setting 1}
#'   \item{TN2}{sum of the IE(International Experts)'s diagnoses for positive disease at setting 2}
#'   \item{TN12}{TN1+TN2}
#'   \item{JN1}{number of non-missing IE's diagnoses for positive disease at setting 1}
#'   \item{JN2}{number of non-missing IE's diagnoses for positive disease at setting 2}
#'   \item{JN12}{JN1+JN2}
#'   \item{TNN1}{sum of the IE's diagnoses for severe disease at setting 1}
#'   \item{TNN2}{sum of the IE's diagnoses for severe disease at setting 2}
#'   \item{TNN12}{TNN1+TNN2}
#'   \item{JNN1}{number of non-missing IE's diagnoses for severe disease at setting 1}
#'   \item{JNN2}{number of non-missing IE's diagnoses for severe disease at setting 2}
#'   \item{JNN12}{JNN1+JNN2}
#'    }
#' @source \url{https://doi.org/10.1111/biom.12997}
"asrm"
