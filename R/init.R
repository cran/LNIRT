## Startup functions ------------------------------------

#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  start_message <- c(" \n"
                     , "   LNIRT v", packageDescription("LNIRT")$Version, " (log-normal RT-IRT Modeling)",  "\n"
                     , "   ", rep('-', 41), "\n"
                     , "   J.-P. Fox & K. Klotzke \n", sep = "")
                     #, "   ", rep('-', 20), "\n", sep = "")
  
  packageStartupMessage(start_message)
  invisible()
}

