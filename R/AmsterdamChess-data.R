#' Amsterdam Chess Test (ACT) data
#'
#' Responses and response time data from the Amsterdam Chess Test (ACT).
#' 
#' Variables:
#' \itemize{
#'   \item ELO: Standardized ELO rating (numeric)
#'   \item Y1-Y40: item correct score (1 or 0) for scored items 1 – 40 (numeric)
#'   \item RT1-RT40: response time (seconds) for scored items 1 – 40 (numeric)
#' }
#' Three components of chess expertise are measured: 
#' \itemize{
#'   \item Tactical skill (20 items): item 1-20
#'   \item Positional skill (10 items): item 21-30
#'   \item End-game skill (10 items): item 31-40
#' }
#'
#' @docType data
#'
#' @usage data(AmsterdamChess)
#'
#' @format A dataframe with 259 rows and 81 variables.
#'
#' @keywords datasets
#'
#' @references van der Maas, H. L., & Wagenmakers, E. J. (2005). A psychometric analysis of chess expertise. The American journal of psychology, 118(1), 29-60.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/15822609}{PubMed})
#'
"AmsterdamChess"