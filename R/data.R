#' Known modified DNA nucleotides.
#'
#' A dataset containing the precursor masses of several known DNA adducts.
#'
#' @format A data frame with 349 rows and 7 variables:
#' \describe{
#'   \item{Name}{common adduct name}
#'   \item{Formula}{Lewis formula}
#'   \item{DBID}{Modified DNA DataBase ID / accession}
#'   \item{PCID}{PubChem ID}
#'   \item{Exact.Mass}{the exact compound mass}
#'   \item{Precursor}{the precursor mz}
#'   \item{Aglycone}{the mz of the compound after removal of the ribose}
#'   ...
#' }
"mdnadb"