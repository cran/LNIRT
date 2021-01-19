#' Credential Form data
#'
#' Responses and response time data from the credential data set of Cizek and Wollack (2016).
#' 
#' Variables:
#' \itemize{
#'   \item EID: Examinee ID (character)
#'   \item FormID: Test form name (character)
#'   \item Flagged: 1/0 variable to indicate whether the test vendor suspects the examinee may have engaged in inappropriate behavior (numeric)
#'   \item Pretest: Pretest item set assigned to candidate (numeric) 
#'   \item Attempt: Count of the attempt number for the candidate. A score of 1 indicates that candidate is a new, first-time examinee. Any examinee sitting for the exam for the fourth time or more is marked as 4+ (character)
#'   \item Country: Country where candidate was educated (character)
#'   \item StateCode: 2-digit code corresponding to the state in which the Candidate applied for licensure (numeric)
#'   \item School_ID: 4-digit code corresponding to the particular institution in which the Candidate received his/her educational training (numeric)
#'   \item Cent_id:  4-digit code corresponding to the particular testing center in which the Candidate sat for the exam (numeric)
#'   \item Tot_time: The number of seconds testing (numeric)
#'   \item iresp.1-170: item responses (1 to 4 or NA) for scored items 1 – 170 (numeric)
#'   \item iresp.171-180: item responses (1 to 4 or NA) for 10 pilot items for pilot set 6 or 9 (numeric)
#'   \item iresp.181-190: item responses (1 to 4 or NA) for 10 pilot items for pilot set 7 or 10 (numeric)
#'   \item iresp.191-200: item responses (1 to 4 or NA) for 10 pilot items for pilot set 8 or 11 (numeric)
#'   \item iraw.1-170: item correct score (1 or 0) for scored items 1 – 170 (numeric)
#'   \item iraw.171-180: item correct score (1 or 0) for 10 pilot items for pilot set 6 or 9 (numeric)
#'   \item iraw.181-190: item correct score (1 or 0) for 10 pilot items for pilot set 7 or 10 (numeric)
#'   \item iraw.191-200: item correct score (1 or 0) for 10 pilot items for pilot set 8 or 11 (numeric)
#'   \item idur.1-170: response time (in seconds) for scored items 1 – 170 (numeric)
#'   \item idur.171-180: response time (in seconds) for 10 pilot items for pilot set 6 or 9 (numeric)
#'   \item idur.181-190: response time (in seconds) for 10 pilot items for pilot set 7 or 10 (numeric)
#'   \item idur.191-200: response time (in seconds) for 10 pilot items for pilot set 8 or 11 (numeric)
#' }
#'
#' @docType data
#'
#' @usage data(CredentialForm1)
#' 
#' @format A dataframe with 1636 rows and 610 variables.
#'
#' @keywords datasets
#'
#' @references Cizek GJ, Wollack JA (eds.) (2016). Handbook of Quantitative Methods for Detecting Cheating on Tests. Routledge.
#' (\href{https://www.taylorfrancis.com/books/handbook-quantitative-methods-detecting-cheating-tests-gregory-cizek-james-wollack/e/10.4324/9781315743097}{Taylor&Francis})
#'
#' @example inst/examples/example.CredentialForm1.R
"CredentialForm1"