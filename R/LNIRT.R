#' Log-normal response time IRT modelling
#'
#' @importFrom MASS mvrnorm
#' @importFrom methods hasArg is
#' @importFrom stats ks.test pchisq pgamma pnorm qnorm rbeta
#' rbinom rchisq rgamma rlnorm rnorm runif var cov2cor sd
#' @importFrom utils flush.console packageDescription setTxtProgressBar txtProgressBar
#'
#' @param RT
#' a Person-x-Item matrix of log-response times (time spent on solving an item).
#' @param Y
#' a Person-x-Item matrix of responses.
#' @param data
#' either a list or a simLNIRT object containing the response time and response matrices
#' and optionally the predictors for the item and person parameters.
#' If a simLNIRT object is provided, in the summary the simulated item and time parameters are shown alongside of the estimates.
#' If the required variables cannot be found in the list, or if no data object is given, then the variables are taken
#' from the environment from which LNIRT is called.
#' @param XG
#' the number of MCMC iterations to perform (default: 1000).
#' @param burnin
#' the percentage of MCMC iterations to discard as burn-in period (default: 10).
#' @param XGresid
#' the number of MCMC iterations to perform before residuals are computed (default: 1000).
#' @param guess
#' include guessing parameters in the IRT model (default: false).
#' @param par1
#' use alternative parameterization (default: false).
#' @param residual
#' compute residuals, >1000 iterations are recommended (default: false).
#' @param td
#' estimate the time-discrimination parameter(default: true).
#' @param WL
#' define the time-discrimination parameter as measurement error variance parameter (default: false).
#' @param ident
#' set identification rule (default: 2).
#' @param alpha
#' an optional vector of pre-defined item-discrimination parameters.
#' @param beta
#' an optional vector of pre-defined item-difficulty parameters.
#' @param phi
#' an optional vector of predefined time discrimination parameters.
#' @param lambda
#' an optional vector of predefined time intensity parameters.
#' @param XPA
#' an optional matrix of predictors for the person ability parameters.
#' @param XPT
#' an optional matrix of predictors for the person speed parameters.
#' @param XIA
#' an optional matrix of predictors for the item-difficulty parameters.
#' @param XIT
#' an optional matrix of predictors for the item-intensity parameters.
#' @param MBDY
#' an optional indicator matrix for response missings due to the test design (0: missing by design, 1: not missing by design).
#' @param MBDT
#' an optional indicator matrix for response time missings due to the test design (0: missing by design, 1: not missing by design).
#'
#' @return
#' an object of class LNIRT.
#'
#' @examples
#' \dontrun{
#' # Log-normal response time IRT modelling
#' data <- simLNIRT(N = 500, K = 20, rho = 0.8, WL = FALSE)
#' out <- LNIRT(RT = RT, Y = Y, data = data, XG = 1500, residual = TRUE, WL = FALSE)
#' summary(out) # Print results
#' out$Post.Means$Item.Difficulty # Extract posterior mean estimates
#'
#' library(coda)
#' mcmc.object <- as.mcmc(out$MCMC.Samples$Item.Difficulty) # Extract MCMC samples for coda
#' summary(mcmc.object)
#' plot(mcmc.object)
#' }
#' @export
LNIRT <-
  function(RT,
           Y,
           data,
           XG = 1000,
           burnin = 10,
           XGresid = 1000,
           guess = FALSE,
           par1 = FALSE,
           residual = FALSE,
           td = TRUE,
           WL = FALSE,
           ident = 2,
           alpha,
           beta,
           phi,
           lambda,
           XPA = NULL,
           XPT = NULL,
           XIA = NULL,
           XIT = NULL,
           MBDY = NULL,
           MBDT = NULL) {
    ## ident = 1: Identification : fix mean item difficulty(intensity) and product item (time) discrimination responses and response times
    ## ident = 2: Identification : fix mean ability and speed and product item discrimination responses and response times
    # ident <- 1
    # ident <- 2 # (to investigate person fit using latent scores)
    
    if (XG <= 0) {
      stop("XG must be > 0.")
    }
    if ((burnin <= 0) || (burnin >= 100)) {
      stop("burn-in period must be between 0% and 100%.")
    }
    if (ident != 1 && ident != 2) {
      stop("ident must be 1 or 2.")
    }
    if (residual && (XGresid >= XG || XGresid <= 0)) {
      warning("XGresid must be < XG and > 0. Residuals will not be computed.")
      residual <- FALSE
    }
    
    if (!missing(data) && !is.null(data)) {
      # Try to find RT and Y in the data set first
      tryCatch(
        RT <- eval(substitute(RT), data),
        error = function(e)
          NULL
      )
      tryCatch(
        Y <- eval(substitute(Y), data),
        error = function(e)
          NULL
      )
      
      # Try to find predictors in the data set first
      tryCatch(
        XPA <- eval(substitute(XPA), data),
        error = function(e)
          NULL
      )
      tryCatch(
        XPT <- eval(substitute(XPT), data),
        error = function(e)
          NULL
      )
      tryCatch(
        XIA <- eval(substitute(XIA), data),
        error = function(e)
          NULL
      )
      tryCatch(
        XIT <- eval(substitute(XIT), data),
        error = function(e)
          NULL
      )
      
      # Try to find MBDY and MBDT in the data set first
      tryCatch(
        MBDY <- eval(substitute(MBDY), data),
        error = function(e)
          NULL
      )
      tryCatch(
        MBDT <- eval(substitute(MBDT), data),
        error = function(e)
          NULL
      )
    } else {
      data <- NULL
    }
    Y <- as.matrix(Y)
    RT <- as.matrix(RT)
    
    if ((nrow(Y) != nrow(RT)) || (ncol(Y) != ncol(RT))) {
      stop("Y and RT must be of equal dimension.")
    }
    
    N <- nrow(Y)
    K <- ncol(Y)
    
    if (!is.null(XPA)) {
      XPA <- as.matrix(XPA)
      if (nrow(XPA) != N) {
        stop("nrow(XPA) must be equal to the number of persons.")
      }
    }
    if (!is.null(XPT)) {
      XPT <- as.matrix(XPT)
      if (nrow(XPT) != N) {
        stop("nrow(XPT) must be equal to the number of persons.")
      }
    }
    if (!is.null(XIA)) {
      XIA <- as.matrix(XIA)
      if (nrow(XIA) != K) {
        stop("nrow(XIA) must be equal to the number of items.")
      }
    }
    if (!is.null(XIT)) {
      XIT <- as.matrix(XIT)
      if (nrow(XIT) != K) {
        stop("nrow(XIT) must be equal to the number of items.")
      }
    }
    
    if (!is.null(MBDY)) {
      MBDY <- as.matrix(MBDY)
      if ((nrow(MBDY) != nrow(Y)) || (ncol(MBDY) != ncol(Y))) {
        stop("MBDY and Y must be of equal dimension.")
      }
    }
    if (!is.null(MBDT)) {
      MBDT <- as.matrix(MBDT)
      if ((nrow(MBDT) != nrow(RT)) || (ncol(MBDT) != ncol(RT))) {
        stop("MBDT and RT must be of equal dimension.")
      }
    }
    
    cat (" \n")
    cat ("   LNIRT v", packageDescription("LNIRT")$Version, "\n", sep = "")
    cat ("   ", rep('-', 20), "\n\n", sep = "")
    cat (
      "   * MCMC sampler initialized (XG:",
      XG,
      ", Burnin:",
      paste(burnin, "%", sep = ""),
      ")\n",
      sep = ""
    )
    cat ("   * Binary response matrix loaded (", N, "x", K, ") \n", sep = "")
    cat ("   * Response time matrix loaded (", N, "x", K, ") \n\n", sep = "")
    
    # Initialize progress bar
    cat ("   MCMC progress: \n")
    pb <-
      txtProgressBar(
        min = 1,
        max = XG,
        initial = 1,
        style = 3,
        width = 45,
        char = "="
      )
    
    PNO <- guess # TRUE: guessing included
    #WL <- 1  #time discrimination = 1/sqrt(error variance)
    #td <- 1  #time discrimination fixed to one
    
    if (missing(alpha)) {
      discr <- TRUE  #discrimination estimated
    } else {
      discr <- FALSE  #discrimination known
    }
    
    if (missing(beta)) {
      diffc <- TRUE  #difficulties estimated
    } else {
      diffc <- FALSE  #difficulties known
    }
    
    if (missing(phi)) {
      discrt <- 0 #time discrimination estimated
    } else{
      discrt <- 1 #time discrimination known
    }
    
    if (missing(lambda)) {
      difft <- 0 #time intensity estimated
    } else{
      difft <- 1 #time intensity known
    }
    
    
    
    
    if (is.null(XPA) && is.null(XPT)) {
      nopredictorp <- TRUE
      XPA <- matrix(1, ncol = 1, nrow = N) #default intercept for ability
      XPT <- matrix(1, ncol = 1, nrow = N) #default intercept for speed
    }
    else {
      nopredictorp <- FALSE
      #if (!is.null(XPA)) {#location XPA set to zero (BUT not Dummy coded variables)
      ##	XPA <- matrix(XPA,ncol=ncol(XPA),nrow=N) - matrix(apply(XPA,2,mean),ncol=ncol(XPA),nrow=N,byrow=T)
      #}
      #if (!is.null(XPT)) {#location XPT set to zero (BUT not Dummy coded variables)
      ##	XPT <- matrix(XPT,ncol=ncol(XPT),nrow=N) - matrix(apply(XPT,2,mean),ncol=ncol(XPT),nrow=N,byrow=T)
      #}
      if (is.null(XPA)) {
        XPA <- matrix(1, ncol = 1, nrow = N) #default intercept for ability
      }
      if (is.null(XPT)) {
        XPT <- matrix(1, ncol = 1, nrow = N) #default intercept for speed
      }
    }
    MmuP <- matrix(0, nrow = XG, ncol = c(ncol(XPA) + ncol(XPT)))
    
    
    if (is.null(XIA) & is.null(XIT)) {
      nopredictori <- TRUE
      MmuI <- matrix(0, nrow = XG, ncol = 4)
      XIA <- matrix(NA, nrow = K, ncol = 0)
      XIT <- matrix(NA, nrow = K, ncol = 0)
    }
    else {
      nopredictori <- FALSE
      if (!is.null(XIA)) {
        #location XIA set to zero (BUT not Dummy coded variables)
        ##	XIA <- matrix(XIA,ncol=ncol(XIA),nrow=K) - matrix(apply(XIA,2,mean),ncol=ncol(XIA),nrow=K,byrow=T)
        if (ident == 2) {
          XIA <- cbind(rep(1, K), XIA)
        }
      }
      if (!is.null(XIT)) {
        #location XIT set to zero (BUT not Dummy coded variables)
        ##	XIT <- matrix(XIT,ncol=ncol(XIT),nrow=K) - matrix(apply(XIT,2,mean),ncol=ncol(XIT),nrow=K,byrow=T)
        if (ident == 2) {
          XIT <- cbind(rep(1, K), XIT)
        }
      }
      if (is.null(XIA)) {
        XIA <- matrix(1, ncol = 1, nrow = K) #default intercept for ability
      }
      if (is.null(XIT)) {
        XIT <- matrix(1, ncol = 1, nrow = K) #default intercept for speed
      }
      MmuI <- matrix(0, nrow = XG, ncol = 2 + c(ncol(XIA) + ncol(XIT)))
    }
    
    ## population theta (ability - speed)
    theta <- matrix(rnorm(N * 2), ncol = 2) # 1: ability, 2: speed
    muP <-
      matrix(0, nrow = N, ncol = 2) # Mean estimates for person parameters
    SigmaP <- diag(2) # Person covariance matrix
    muP0 <- matrix(0, 1, 2)
    SigmaP0 <- diag(2)
    
    ## population item (ability - speed)
    ab <-
      matrix(rnorm(K * 4), ncol = 4) # 1: item discrimination 2: item difficulty 3: time discrimination 4: time intensity
    ab[, c(1, 3)] <- 1
    muI <-
      t(matrix(rep(c(1, 0), 2 * K), ncol = K)) # Mean estimates for item parameters
    muI0 <- muI[1,]
    SigmaI <- diag(4) # Item covariance matrix
    SigmaI0 <- diag(4) / 10
    for (ii in 1:4) {
      SigmaI0[ii,] <- SigmaI0[ii,] * rep(c(0.5, 3), 2)
    }
    ifelse(PNO, guess0 <- rep(0.2, K), guess0 <- rep(0, K))
    
    ## storage
    MT <- MT2 <- array(0, dim = c(N, 2))
    MAB <- array(0, dim = c(XG, K, 4))
    Mguess <- matrix(0, XG, K)
    MSP <- array(0, dim = c(XG, 2, 2))
    MSI <- array(0, dim = c(XG, 4, 4))
    Msigma2 <- matrix(0, ncol = K, nrow = XG)
    sigma2 <- rep(1, K)
    
    lZP <- 0
    lZPT <- 0
    lZI <- 0
    lZPAT <- matrix(0, ncol = 1, nrow = N)
    lZIA <- matrix(0, ncol = 1, nrow = K)
    lZPA <- matrix(0, ncol = 1, nrow = N)
    
    EAPl0 <- matrix(0, ncol = K, nrow = N)
    PFl <- matrix(0, ncol = 1, nrow = N)
    PFlp <- matrix(0, ncol = 1, nrow = N)
    IFl <- matrix(0, ncol = 1, nrow = K)
    IFlp <- matrix(0, ncol = 1, nrow = K)
    
    EAPresid <- matrix(0, ncol = K, nrow = N)
    EAPresidA <- matrix(0, ncol = K, nrow = N)
    EAPKS <- matrix(0, ncol = 1, nrow = K)
    EAPKSA <- matrix(0, ncol = 1, nrow = K)
    
    EAPphi <- matrix(0, ncol = 1, nrow = K)
    EAPlambda <- matrix(0, ncol = 1, nrow = K)
    EAPtheta <- matrix(0, ncol = 2, nrow = N)
    EAPsigma2 <- matrix(0, ncol = 1, nrow = K)
    
    if (XG > XGresid) {
      MPF <- matrix(0, ncol = N, nrow = XG - XGresid)
      MPFb <- matrix(0, ncol = N, nrow = XG - XGresid)
      MPFp <- matrix(0, ncol = N, nrow = XG - XGresid)
      MPFbp <- matrix(0, ncol = N, nrow = XG - XGresid)
      
      MPFTb <- matrix(0, ncol = N, nrow = XG - XGresid)
      MPFTbp <- matrix(0, ncol = N, nrow = XG - XGresid)
    }
    EAPbeta <- rep(0, K)
    EAPalpha <- rep(0, K)
    EAPmub <- 0
    EAPsigmab <- 1
    EAPSigmaP <- 1
    EAPmuI <- 0
    EAPSigmaI <- 1
    iis <- 1
    
    ## Missing By Design
    if (is.null(MBDY)) {
      MBDY <-
        matrix(1, ncol = K, nrow = N) ## Y : no missing by design, 0=missing by design,1=not missing by design
      MBDYI <- TRUE #no design missing Y
    } else{
      MBDYI <- FALSE #design missings Y
      Y[MBDY == 0] <- 0 #recode design missings to 0
    }
    if (is.null(MBDT)) {
      MBDT <-
        matrix(1, ncol = K, nrow = N) ## RT: no missing by design, 0=missing by design,1=not missing by design
      MBDTI <- TRUE #no design missing RT
    } else{
      MBDTI <- FALSE #design missings RT
      RT[MBDT == 0] <- 0 #recode design missings to 0
    }
    
    # Indicate which responses are missing
    D <- matrix(1, ncol = 1, nrow = N * K)
    D[which(is.na(Y))] <-
      0	## missing Y values and missing by design
    D <- matrix(D, nrow = N, ncol = K)
    
    # Indicate which RT's are missing
    DT <- matrix(1, ncol = 1, nrow = N * K)
    DT[which(is.na(RT))] <-
      0	## missing RT values and missing by design
    DT <- matrix(DT, nrow = N, ncol = K)
    
    EAPCP1 <- matrix(0, ncol = 1, nrow = N)
    EAPCP2 <- matrix(0, ncol = 1, nrow = N)
    EAPCP3 <- matrix(0, ncol = 1, nrow = N)
    
    
    # Output
    MCMC.Samples <- list()
    MCMC.Samples$Person.Ability <- matrix(NA, nrow = XG, ncol = N)
    MCMC.Samples$Person.Speed <- matrix(NA, nrow = XG, ncol = N)
    
    
    ## Start MCMC algorithm
    
    for (ii in 1:XG) {
      # For each iteration
      if (sum(DT == 0) > 0) {
        # Simulate missing RT's
        if (MBDTI) {
          if (WL) {
            RT <-
              SimulateRT(
                RT = RT,
                zeta = theta[, 2],
                lambda = ab[, 4],
                phi = rep(1, K),
                sigma2 = sigma2,
                DT = DT
              )
          } else {
            RT <-
              SimulateRT(
                RT = RT,
                zeta = theta[, 2],
                lambda = ab[, 4],
                phi = ab[, 3],
                sigma2 = sigma2,
                DT = DT
              )
          }
        } else{
          if (WL) {
            RT <-
              SimulateRTMBD(
                RT = RT,
                zeta = theta[, 2],
                lambda = ab[, 4],
                phi = rep(1, K),
                sigma2 = sigma2,
                DT = DT,
                MBDT = MBDT
              )
          } else {
            RT <-
              SimulateRTMBD(
                RT = RT,
                zeta = theta[, 2],
                lambda = ab[, 4],
                phi = ab[, 3],
                sigma2 = sigma2,
                DT = DT,
                MBDT = MBDT
              )
          }
        }
      }
      # ability test
      if (MBDYI) {
        if (PNO) {
          # If guessing
          if (sum(D == 0) > 0) {
            # Simulate missing responses
            Y <-
              SimulateY(
                Y = Y,
                theta = theta[, 1],
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                guess0 = guess0,
                D = D
              )
          }
          if (par1) {
            SR <-
              DrawS_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                guess0 = guess0,
                theta0 = theta[, 1],
                Y = Y
              )
            ZR <-
              DrawZ_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                theta0 = theta[, 1],
                S = SR,
                D = D
              )
          } else {
            SR <-
              DrawS_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                guess0 = guess0,
                theta0 = theta[, 1],
                Y = Y
              )
            ZR <-
              DrawZ_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                theta0 = theta[, 1],
                S = SR,
                D = D
              )
          }
        } else {
          if (par1) {
            ZR <-
              DrawZ_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                theta0 = theta[, 1],
                S = Y,
                D = D
              )
          } else {
            ZR <-
              DrawZ_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                theta0 = theta[, 1],
                S = Y,
                D = D
              )
          }
        }
      } else {
        if (PNO) {
          # If guessing
          if (sum(D == 0) > 0) {
            # Simulate missing responses
            Y <-
              SimulateYMBD(
                Y = Y,
                theta = theta[, 1],
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                guess0 = guess0,
                D = D,
                MBDY = MBDY
              )
          }
          if (par1) {
            SR <-
              DrawSMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                guess0 = guess0,
                theta0 = theta[, 1],
                Y = Y,
                MBDY = MBDY
              )
            ZR <-
              DrawZMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                theta0 = theta[, 1],
                S = SR,
                D = D,
                MBDY = MBDY
              )
          } else {
            SR <-
              DrawSMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                guess0 = guess0,
                theta0 = theta[, 1],
                Y = Y,
                MBDY = MBDY
              )
            ZR <-
              DrawZMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                theta0 = theta[, 1],
                S = SR,
                D = D,
                MBDY = MBDY
              )
          }
        } else {
          if (par1) {
            ZR <-
              DrawZMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2] * ab[, 1],
                theta0 = theta[, 1],
                S = Y,
                D = D,
                MBDY = MBDY
              )
          } else {
            ZR <-
              DrawZMBD_LNIRT(
                alpha0 = ab[, 1],
                beta0 = ab[, 2],
                theta0 = theta[, 1],
                S = Y,
                D = D,
                MBDY = MBDY
              )
          }
        }
      }
      
      
      # Draw ability parameter
      dum <- Conditional(1, muP, SigmaP, theta)
      if (MBDYI) {
        if (par1) {
          theta[, 1] <-
            DrawTheta_LNIRT(
              alpha0 = ab[, 1],
              beta0 = ab[, 2] * ab[, 1],
              Z = ZR,
              mu = dum$CMU,
              sigma = dum$CVAR
            )
        } else {
          theta[, 1] <-
            DrawTheta_LNIRT(
              alpha0 = ab[, 1],
              beta0 = ab[, 2],
              Z = ZR,
              mu = dum$CMU,
              sigma = dum$CVAR
            )
        }
      } else{
        if (par1) {
          theta[, 1] <-
            DrawThetaMBD_LNIRT(
              alpha0 = ab[, 1],
              beta0 = ab[, 2] * ab[, 1],
              Z = ZR,
              mu = dum$CMU,
              sigma = dum$CVAR,
              MBDY = MBDY
            )
        } else {
          theta[, 1] <-
            DrawThetaMBD_LNIRT(
              alpha0 = ab[, 1],
              beta0 = ab[, 2],
              Z = ZR,
              mu = dum$CMU,
              sigma = dum$CVAR,
              MBDY = MBDY
            )
        }
      }
      
      if (ident == 2) {
        # rescale for identification
        theta[, 1] <- theta[, 1] - mean(theta[, 1])
      }
      
      # Draw speed parameter
      dum <- Conditional(2, muP, SigmaP, theta)
      if (MBDTI) {
        if (par1) {
          theta[, 2] <-
            DrawZeta(
              RT = RT,
              phi = ab[, 3],
              lambda = ab[, 3] * ab[, 4],
              sigma2 = sigma2,
              mu = as.vector(dum$CMU),
              sigmaz = as.vector(dum$CVAR)
            )  ## speed
        } else {
          theta[, 2] <-
            DrawZeta(
              RT = RT,
              phi = ab[, 3],
              lambda = ab[, 4],
              sigma2 = sigma2,
              mu = as.vector(dum$CMU),
              sigmaz = as.vector(dum$CVAR)
            )  ## speed
        }
      } else{
        if (par1) {
          theta[, 2] <-
            DrawZetaMBD(
              RT = RT,
              phi = ab[, 3],
              lambda = ab[, 3] * ab[, 4],
              sigma2 = sigma2,
              mu = as.vector(dum$CMU),
              sigmaz = as.vector(dum$CVAR),
              MBDT = MBDT
            )  ## speed
        } else {
          theta[, 2] <-
            DrawZetaMBD(
              RT = RT,
              phi = ab[, 3],
              lambda = ab[, 4],
              sigma2 = sigma2,
              mu = as.vector(dum$CMU),
              sigmaz = as.vector(dum$CVAR),
              MBDT = MBDT
            )  ## speed
        }
      }
      
      # Rescale for identification
      if (ident == 2) {
        theta[, 2] <- theta[, 2] - mean(theta[, 2])
      }
      
      MCMC.Samples$Person.Ability[ii,] <- theta[, 1]
      MCMC.Samples$Person.Speed[ii,] <- theta[, 2]
      
      MT[1:N, 1:2] <- MT[1:N, 1:2] + theta # for theta
      MT2[1:N, 1:2] <- MT2[1:N, 1:2] + theta ^ 2 # for sd(theta)
      
      # Draw guessing parameter from beta distribution
      if (PNO) {
        if (MBDYI) {
          guess0 <- DrawC_LNIRT(S = SR, Y = Y)
        } else{
          guess0 <- DrawCMBD_LNIRT(S = SR,
                                   Y = Y,
                                   MBDY = MBDY)
        }
      }
      Mguess[ii,] <- guess0
      
      # ab - 1: item discrimination 2: item difficulty 3: time discrimination 4: time intensity
      if (WL) {
        ab1 <- cbind(ab[, 1], ab[, 2], 1 / sqrt(sigma2), ab[, 4])
      } else {
        ab1 <- ab
      }
      
      # Draw item discrimination parameter
      if (discr) {
        dum <-
          Conditional(
            kk = 1,
            Mu = muI,
            Sigma = SigmaI,
            Z = ab1
          )  #discrimination
        if (MBDYI) {
          if (par1 == 1) {
            ab[, 1] <-
              abs(
                DrawAlpha_LNIRT(
                  theta = theta[, 1],
                  beta = ab[, 2],
                  Z = ZR,
                  mu = dum$CMU,
                  sigma = dum$CVAR
                )
              )
          } else {
            ab[, 1] <-
              abs(
                DrawPhi_LNIRT(
                  RT = ZR,
                  lambda = -ab[, 2],
                  zeta = -theta[, 1],
                  sigma2 = rep(1, K),
                  mu = dum$CMU,
                  sigmal = dum$CVAR
                )
              )
          }
        } else{
          if (par1 == 1) {
            ab[, 1] <-
              abs(
                DrawAlphaMBD_LNIRT(
                  theta = theta[, 1],
                  beta = ab[, 2],
                  Z = ZR,
                  mu = dum$CMU,
                  sigma = dum$CVAR,
                  MBDY = MBDY
                )
              )
          } else {
            ab[, 1] <-
              abs(
                DrawAlphaMBD2_LNIRT(
                  theta = theta[, 1],
                  beta = ab[, 2],
                  Z = ZR,
                  mu = dum$CMU,
                  sigma = dum$CVAR,
                  MBDY = MBDY
                )
              )
          }
        }
      } else {
        ab[, 1] <- alpha
      }
      ab[, 1] <- ab[, 1] / (prod(ab[, 1]) ^ (1 / K))
      
      if (WL) {
        ab1 <- cbind(ab[, 1], ab[, 2], 1 / sqrt(sigma2), ab[, 4])
      } else {
        ab1 <- ab
      }
      
      # Draw item difficulty parameter
      if (diffc) {
        dum <-
          Conditional(
            kk = 2,
            Mu = muI,
            Sigma = SigmaI,
            Z = ab1
          )  #difficulty
        if (MBDYI) {
          if (par1 == 1) {
            ab[, 2] <-
              DrawBeta_LNIRT(
                theta = theta[, 1],
                alpha = ab[, 1],
                Z = ZR,
                mu = dum$CMU,
                sigma = dum$CVAR
              )
          } else {
            ab[, 2] <-
              -DrawLambda_LNIRT(ZR,-ab[, 1], theta[, 1], rep(1, K), dum$CMU, dum$CVAR)$lambda
          }
        } else{
          if (par1 == 1) {
            ab[, 2] <-
              DrawBetaMBD_LNIRT(
                theta = theta[, 1],
                alpha = ab[, 1],
                Z = ZR,
                mu = dum$CMU,
                sigma = dum$CVAR,
                MBDY = MBDY
              )
          } else {
            ab[, 2] <-
              DrawBetaMBD2_LNIRT(
                theta = theta[, 1],
                alpha = ab[, 1],
                Z = ZR,
                mu = dum$CMU,
                sigma = dum$CVAR,
                MBDY = MBDY
              )
          }
        }
      } else {
        ab[, 2] <- beta
      }
      if (par1) {
        # ab[,2] <- ab[,2]/ab[,1]
      }
      if (ident == 1) {
        # rescale for identification
        ab[, 2] <- ab[, 2] - mean(ab[, 2])
      }
      
      # Draw time discrimination parameter
      if ((WL) | (!td)) {
        # no time discrimination, 1/(sqrt(error variance)) = discrimination on MVN prior
        ab[, 3] <- 1
      } else {
        if (discrt == 0) {
          dum <- Conditional(3, muI, SigmaI, ab)  #time discrimination
          if (MBDTI) {
            ab[, 3] <-
              abs(DrawPhi_LNIRT(RT, ab[, 4], theta[, 2], sigma2, dum$CMU, dum$CVAR))
          } else{
            ab[, 3] <-
              abs(DrawPhiMBD_LNIRT(RT, ab[, 4], theta[, 2], sigma2, dum$CMU, dum$CVAR, MBDT =
                                     MBDT))
          }
        } else{
          ab[, 3] <- phi
        }
        ab[, 3] <- ab[, 3] / (prod(ab[, 3]) ^ (1 / K))
      }
      
      if (WL) {
        ab1 <- cbind(ab[, 1], ab[, 2], 1 / sqrt(sigma2), ab[, 4])
      } else {
        ab1 <- ab
      }
      
      # Draw time intensity parameter
      if (difft == 0) {
        dum <-
          Conditional(
            kk = 4,
            Mu = muI,
            Sigma = SigmaI,
            Z = ab1
          )  #time intensity
        if (MBDTI) {
          ab[, 4] <-
            DrawLambda_LNIRT(RT, ab[, 3], theta[, 2], sigma2, dum$CMU, dum$CVAR)$lambda
        } else{
          ab[, 4] <-
            DrawLambdaMBD_LNIRT(RT, ab[, 3], theta[, 2], sigma2, dum$CMU, dum$CVAR, MBDT =
                                  MBDT)
        }
      } else{
        ab[, 4] <- lambda
      }
      
      if (ident == 1) {
        ab[, 4] <- ab[, 4] - mean(ab[, 4])
      }
      
      MAB[ii, 1:K, 1:4] <- ab
      
      
      # WL: Measurement error variance for time discrimination
      # Else: Measurement error variance for time term
      if (MBDTI) {
        sigma2 <-
          SampleS2_LNIRT(
            RT = RT,
            zeta = theta[, 2],
            lambda = ab[, 4],
            phi = ab[, 3]
          )
      } else{
        sigma2 <-
          SampleS2MBD_LNIRT(
            RT = RT,
            zeta = theta[, 2],
            lambda = ab[, 4],
            phi = ab[, 3],
            MBDT = MBDT
          )
      }
      if (WL) {
        Msigma2[ii, 1:K] <- 1 / sqrt(sigma2)
      } else {
        Msigma2[ii, 1:K] <- sigma2
      }
      
      if (nopredictorp) {
        # Population mean estimate for person ability and speed
        X <- matrix(1, N, 1)
        muP <-
          SampleB(
            Y = theta,
            X = X,
            Sigma = SigmaP,
            B0 = muP0,
            V0 = SigmaP0
          )
        MmuP[ii,] <- muP$B
        muP <- muP$pred
      } else{
        muPP <- SampleBX_LNIRT(Y = theta, XPA = XPA, XPT = XPT)
        MmuP[ii,] <- muPP$B
        muP <- muPP$pred
      }
      
      
      # Covariance matrix person parameters
      SS <- crossprod(theta - muP) + SigmaP0
      SigmaP <- rwishart(2 + N, chol2inv(chol(SS)))$IW
      MSP[ii, ,] <- SigmaP
      
      X <- matrix(1, K, 1)
      if (WL) {
        ab1 <- cbind(ab[, 1], ab[, 2], 1 / sqrt(sigma2), ab[, 4])
      } else {
        ab1 <- ab
      }
      
      if (nopredictori) {
        # Population mean estimates for item parameters
        # 1: item discrimination 2: item difficulty 3: time discrimination 4: time intensity
        muI2 <-
          SampleB(
            Y = ab1,
            X = X,
            Sigma = SigmaI,
            B0 = muI0,
            V0 = SigmaI0
          )
        if (ident == 2) {
          MmuI[ii, c(1, 2, 3, 4)] <- muI2$B[c(1, 2, 3, 4)]
          muI[, c(1, 2, 3, 4)] <- muI2$pred[, c(1, 2, 3, 4)]
        } else {
          MmuI[ii, c(1, 3)] <- muI2$B[c(1, 3)]
          muI[, c(1, 3)] <- muI2$pred[, c(1, 3)]
        }
      } else {
        ## Item Predictors
        
        #### XIA and XIT should include intercept when ident=2
        set2 <- c(2, 4)
        ab2 <- ab1[, set2]
        dum <- SampleBX_LNIRT(Y = ab2, XPA = XIA, XPT = XIT)
        MmuI[ii, 2:(ncol(XIA) + 1)] <- dum$B[1:ncol(XIA)]
        MmuI[ii, (3 + ncol(XIA)):(ncol(XIA) + 2 + ncol(XIT))] <-
          dum$B[(ncol(XIA) + 1):(ncol(XIA) + ncol(XIT))]
        if (ident == 1) {
          MmuI[ii, c(2, 3 + ncol(XIA))] <- c(0, 0)
        }
        muI[, set2] <- dum$pred
        
        #mean discrimination and time discrimination
        set2 <- c(1, 3)
        SigmaI1 <- SigmaI[set2, set2]
        muI01 <- muI0[set2]
        SigmaI01 <- SigmaI0[set2, set2]
        ab2 <- ab1[, set2]
        muI2 <-
          SampleB(
            Y = ab2,
            X = X,
            Sigma = SigmaI1,
            B0 = muI01,
            V0 = SigmaI01
          )
        MmuI[ii, c(1, 2 + ncol(XIA))] <- muI2$B
        muI[, c(1, 3)] <- muI2$pred
      }
      
      
      # Covariance matrix item parameters
      muI1 <- matrix(muI,
                     ncol = 4,
                     nrow = K,
                     byrow = FALSE)
      SS <- crossprod(ab1 - muI1) + SigmaI0
      SigmaI <- rwishart(4 + K, chol2inv(chol(SS)))$IW
      MSI[ii, ,] <- SigmaI
      
      if (ii > XGresid) {
        EAPmuI <- (muI[1, 4] + (iis - 1) * EAPmuI) / iis
        EAPsigma2 <- (sigma2 + (iis - 1) * EAPsigma2) / iis
        EAPlambda <- (ab[, 4] + (iis - 1) * EAPlambda) / iis
        EAPphi <- (ab[, 3] + (iis - 1) * EAPphi) / iis
        EAPSigmaI <- (SigmaI[4, 4] + (iis - 1) * EAPSigmaI) / iis
        EAPtheta <- (theta + (iis - 1) * EAPtheta) / iis
        
        EAPmub <- (muI[1, 2] + (iis - 1) * EAPmub) / iis
        EAPsigmab <- (SigmaI[2, 2] + (iis - 1) * EAPsigmab) / iis
        if (par1) {
          EAPbeta <- (ab[, 1] * ab[, 2] + (iis - 1) * EAPbeta) / iis
        } else {
          EAPbeta <- (ab[, 2] + (iis - 1) * EAPbeta) / iis
        }
        EAPalpha <- (ab[, 1] + (iis - 1) * EAPalpha) / iis
        EAPSigmaP <- (SigmaP[1, 1] + (iis - 1) * EAPSigmaP) / iis
        
        
        ##############################
        
        if (residual) {
          if (MBDYI) {
            ## IRT Fit Evaluation
            if (par1) {
              beta1 <- ab[, 1] * ab[, 2]
            } else {
              beta1 <- ab[, 2]
            }
            if (PNO) {
              dum <-
                residualA(
                  Z = ZR,
                  Y = SR,
                  theta = theta[, 1],
                  alpha = ab[, 1],
                  beta = beta1,
                  EAPtheta = EAPtheta[, 1],
                  EAPalpha = EAPalpha,
                  EAPbeta = EAPbeta
                )
            } else {
              dum <-
                residualA(
                  Z = ZR,
                  Y = Y,
                  theta = theta[, 1],
                  alpha = ab[, 1],
                  beta = beta1,
                  EAPtheta = EAPtheta[, 1],
                  EAPalpha = EAPalpha,
                  EAPbeta = EAPbeta
                )
            }
            EAPKSA <- (dum$KS[1:K, 1] + (iis - 1) * EAPKSA) / iis
            EAPresidA <-
              (dum$presidA + (iis - 1) * EAPresidA) / iis
            
            # lZPAT <- lZPAT + dum$lZPAT
            lZPA <- lZPA + dum$lZPA
            # lZIA <- lZIA + dum$lZIA
            CF2 <-
              ifelse(dum$PFlp < 0.05, 1, 0)  #significance level = .05
            EAPCP2 <- (CF2 + (iis - 1) * EAPCP2) / iis
            
            
            EAPl0 <- ((iis - 1) * EAPl0 + dum$l0) / iis
            PFl <- ((iis - 1) * PFl + dum$PFl) / iis
            IFl <- ((iis - 1) * IFl + dum$IFl) / iis
            PFlp <- ((iis - 1) * PFlp + dum$PFlp) / iis
            IFlp <- ((iis - 1) * IFlp + dum$IFlp) / iis
          } else{
            ## IRT Fit Evaluation
            if (par1) {
              beta1 <- ab[, 1] * ab[, 2]
            } else {
              beta1 <- ab[, 2]
            }
            if (PNO) {
              dum <-
                residualAMBD(
                  Z = ZR,
                  Y = SR,
                  theta = theta[, 1],
                  alpha = ab[, 1],
                  beta = beta1,
                  EAPtheta = EAPtheta[, 1],
                  EAPalpha = EAPalpha,
                  EAPbeta = EAPbeta,
                  MBDY = MBDY
                )
            } else {
              dum <-
                residualAMBD(
                  Z = ZR,
                  Y = Y,
                  theta = theta[, 1],
                  alpha = ab[, 1],
                  beta = beta1,
                  EAPtheta = EAPtheta[, 1],
                  EAPalpha = EAPalpha,
                  EAPbeta = EAPbeta,
                  MBDY = MBDY
                )
            }
            EAPKSA <- (dum$KS[1:K, 1] + (iis - 1) * EAPKSA) / iis
            EAPresidA <-
              (dum$presidA + (iis - 1) * EAPresidA) / iis
            
            # lZPAT <- lZPAT + dum$lZPAT
            lZPA <- lZPA + dum$lZPA
            # lZIA <- lZIA + dum$lZIA
            CF2 <-
              ifelse(dum$PFlp < 0.05, 1, 0)  #significance level = .05
            EAPCP2 <- (CF2 + (iis - 1) * EAPCP2) / iis
            
            
            EAPl0 <- ((iis - 1) * EAPl0 + dum$l0) / iis
            PFl <- ((iis - 1) * PFl + dum$PFl) / iis
            IFl <- ((iis - 1) * IFl + dum$IFl) / iis
            PFlp <- ((iis - 1) * PFlp + dum$PFlp) / iis
            IFlp <- ((iis - 1) * IFlp + dum$IFlp) / iis
          }
          ##############################
          
          if (MBDTI) {
            ## Log-Normal Fit Evaluation
            if (par1) {
              lambda1 <- ab[, 3] * ab[, 4]
            } else {
              lambda1 <- ab[, 4]
            }
            dum <-
              personfitLN(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2
              )  # lZ statistic
            lZP <- lZP + dum$lZP
            lZPT <- lZPT + dum$lZPT
            CF1 <-
              ifelse(dum$lZP < 0.05, 1, 0)  #significance level = .05
            EAPCP1 <- (CF1 + (iis - 1) * EAPCP1) / iis  #speed
            EAPCP3 <- (CF1 * CF2 + (iis - 1) * EAPCP3) / iis
            
            dum <-
              itemfitLN(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2
              )
            lZI <- lZI + dum$lZI
            
            dum <-
              residualLN(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2,
                EAPtheta = EAPtheta[, 2],
                EAPlambda = EAPlambda,
                EAPphi = EAPphi,
                EAPsigma2 = EAPsigma2
              )
            EAPresid <- EAPresid + dum$presid
            EAPKS <- (dum$KS[1:K, 1] + (iis - 1) * EAPKS) / iis
          } else{
            ## Log-Normal Fit Evaluation
            if (par1) {
              lambda1 <- ab[, 3] * ab[, 4]
            } else {
              lambda1 <- ab[, 4]
            }
            dum <-
              personfitLNMBD(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2,
                MBDT = MBDT
              )  # lZ statistic
            lZP <- lZP + dum$lZP
            lZPT <- lZPT + dum$lZPT
            CF1 <-
              ifelse(dum$lZP < 0.05, 1, 0)  #significance level = .05
            EAPCP1 <- (CF1 + (iis - 1) * EAPCP1) / iis  #speed
            EAPCP3 <- (CF1 * CF2 + (iis - 1) * EAPCP3) / iis
            
            dum <-
              itemfitLNMBD(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2,
                MBDT = MBDT
              )
            lZI <- lZI + dum$lZI
            
            dum <-
              residualLNMBD(
                RT = RT,
                theta = theta[, 2],
                phi = ab[, 3],
                lambda = lambda1,
                sigma2 = sigma2,
                EAPtheta = EAPtheta[, 2],
                EAPlambda = EAPlambda,
                EAPphi = EAPphi,
                EAPsigma2 = EAPsigma2,
                MBDT = MBDT
              )
            EAPresid <-
              EAPresid + dum$presid ##for missing values set to zero
            EAPKS <- (dum$KS[1:K, 1] + (iis - 1) * EAPKS) / iis
          }
          iis <- iis + 1
        }
      }
      
      # Update progress bar
      setTxtProgressBar(pb, ii)
      
      # if (ii%%100 == 0)
      #     cat("Iteration ", ii, " ", "\n")
      # flush.console()
    }
    
    
    MT <- MT / XG # Posterior mean theta
    MT2 <- sqrt(MT2 / XG - MT ^ 2) # Posterior mean sd(theta)
    
    if (ii > XGresid & residual) {
      lZP <- lZP / (XG - XGresid)
      lZPT <- lZPT / (XG - XGresid)
      lZI <- lZI / (XG - XGresid)
      EAPresid <- EAPresid / (XG - XGresid)
      lZPAT <- lZPAT / (XG - XGresid)
      lZPA <- lZPA / (XG - XGresid)
      lZIA <- lZIA / (XG - XGresid)
    }
    
    
    kia <- kit <- 0
    if (ncol(XIA) > 0)
      kia <- ncol(XIA) - 1
    if (ncol(XIT) > 0)
      kit <- ncol(XIT) - 1
    
    # Output
    #MCMC.Samples <- list()
    #MCMC.Samples$Person.Ability <- matrix(NA, nrow = XG, ncol = N)
    #MCMC.Samples$Person.Speed <- matrix(NA, nrow = XG, ncol = N)
    MCMC.Samples$Mu.Person.Ability <- MmuP[, 1:ncol(XPA)]
    MCMC.Samples$Mu.Person.Speed <-
      MmuP[, (ncol(XPA) + 1):(ncol(XPA) + ncol(XPT))]
    MCMC.Samples$Var.Person.Ability <- MSP[, 1, 1]
    MCMC.Samples$Var.Person.Speed <- MSP[, 2, 2]
    MCMC.Samples$Cov.Person.Ability.Speed <- MSP[, 1, 2]
    MCMC.Samples$Item.Discrimination <- MAB[, , 1]
    MCMC.Samples$Item.Difficulty <- MAB[, , 2]
    MCMC.Samples$Time.Discrimination <- MAB[, , 3]
    MCMC.Samples$Time.Intensity <- MAB[, , 4]
    MCMC.Samples$Mu.Item.Discrimination <- MmuI[, 1]
    MCMC.Samples$Mu.Item.Difficulty <- MmuI[, 2:(2 + kia)]
    MCMC.Samples$Mu.Time.Discrimination <- MmuI[, 2 + kia + 1]
    MCMC.Samples$Mu.Time.Intensity <-
      MmuI[, (2 + kia + 1 + 1):(ncol(MmuI))]
    MCMC.Samples$Item.Guessing <- Mguess
    MCMC.Samples$Sigma2 <- Msigma2
    MCMC.Samples$CovMat.Item <- MSI
    
    XGburnin <- round(XG * burnin / 100, 0)
    Post.Means <- list()
    Post.Means$Person.Ability <- MT[, 1]
    Post.Means$Person.Speed <- MT[, 2]
    if (ncol(XPA) == 1)
      Post.Means$Mu.Person.Ability <- mean(MmuP[XGburnin:XG, 1])
    else
      Post.Means$Mu.Person.Ability <-
      colMeans(MmuP[XGburnin:XG, 1:ncol(XPA)])
    if (ncol(XPT) == 1)
      Post.Means$Mu.Person.Speed <-
      mean(MmuP[XGburnin:XG, (ncol(XPA) + 1):(ncol(XPA) + ncol(XPT))])
    else
      Post.Means$Mu.Person.Speed <-
      colMeans(MmuP[XGburnin:XG, (ncol(XPA) + 1):(ncol(XPA) + ncol(XPT))])
    Post.Means$Var.Person.Ability <- mean(MSP[XGburnin:XG, 1, 1])
    Post.Means$Var.Person.Speed <- mean(MSP[XGburnin:XG, 2, 2])
    Post.Means$Cov.Person.Ability.Speed <-
      mean(MSP[XGburnin:XG, 1, 2])
    Post.Means$Item.Discrimination <- colMeans(MAB[XGburnin:XG, , 1])
    Post.Means$Item.Difficulty <- colMeans(MAB[XGburnin:XG, , 2])
    Post.Means$Time.Discrimination <- colMeans(MAB[XGburnin:XG, , 3])
    Post.Means$Time.Intensity <- colMeans(MAB[XGburnin:XG, , 4])
    Post.Means$Mu.Item.Discrimination <- mean(MmuI[XGburnin:XG, 1])
    if (kia == 0)
      Post.Means$Mu.Item.Difficulty <-
      mean(MmuI[XGburnin:XG, 2:(2 + kia)])
    else
      Post.Means$Mu.Item.Difficulty <-
      colMeans(MmuI[XGburnin:XG, 2:(2 + kia)])
    Post.Means$Mu.Time.Discrimination <-
      mean(MmuI[XGburnin:XG, 2 + kia + 1])
    if (kit == 0)
      Post.Means$Mu.Time.Intensity <-
      mean(MmuI[XGburnin:XG, (2 + kia + 1 + 1):(ncol(MmuI))])
    else
      Post.Means$Mu.Time.Intensity <-
      colMeans(MmuI[XGburnin:XG, (2 + kia + 1 + 1):(ncol(MmuI))])
    Post.Means$Sigma2 <- colMeans(Msigma2[XGburnin:XG,])
    Post.Means$CovMat.Item  <-
      matrix(
        c(
          round(apply(MSI[XGburnin:XG, 1,], 2, mean), 3),
          round(apply(MSI[XGburnin:XG, 2,], 2, mean), 3),
          round(apply(MSI[XGburnin:XG, 3,], 2, mean), 3),
          round(apply(MSI[XGburnin:XG, 4,], 2, mean), 3)
        ),
        ncol = 4,
        nrow = 4,
        byrow = TRUE
      )
    if (guess)
      Post.Means$Item.Guessing <- colMeans(Mguess[XGburnin:XG, ])
    else
      Post.Means$Item.Guessing <- NULL
    
    if (!(is(data, "simLNIRT"))) {
      data <- NULL # only attach sim data for summary function
    }
    
    if (XG > XGresid & residual) {
      out <-
        list(
          Post.Means = Post.Means,
          MCMC.Samples = MCMC.Samples,
          Mtheta = MT,
          MTSD = MT2,
          MAB = MAB,
          MmuP = MmuP,
          MSP = MSP,
          MmuI = MmuI,
          MSI = MSI,
          Mguess = Mguess,
          Msigma2 = Msigma2,
          lZP = lZP,
          lZPT = lZPT,
          lZPA = lZPA,
          lZI = lZI,
          EAPresid = EAPresid,
          EAPresidA = EAPresidA,
          EAPKS = EAPKS,
          EAPKSA = EAPKSA,
          PFl = PFl,
          PFlp = PFlp,
          IFl = IFl,
          IFlp = IFlp,
          EAPl0 = EAPl0,
          RT = RT,
          Y = Y,
          EAPCP1 = EAPCP1,
          EAPCP2 = EAPCP2,
          EAPCP3 = EAPCP3,
          WL = WL,
          td = td,
          guess = guess,
          par1 = par1,
          data = data,
          XPA = XPA,
          XPT = XPT,
          XIA = XIA,
          XIT = XIT,
          XG = XG,
          burnin = burnin,
          ident = ident,
          residual = residual,
          XGresid = XGresid
        )
    }
    else {
      out <-
        list(
          Post.Means = Post.Means,
          MCMC.Samples = MCMC.Samples,
          Mtheta = MT,
          MTSD = MT2,
          MAB = MAB,
          MmuP = MmuP,
          MSP = MSP,
          MmuI = MmuI,
          MSI = MSI,
          Mguess = Mguess,
          Msigma2 = Msigma2,
          RT = RT,
          Y = Y,
          WL = WL,
          td = td,
          guess = guess,
          par1 = par1,
          data = data,
          XPA = XPA,
          XPT = XPT,
          XIA = XIA,
          XIT = XIT,
          XG = XG,
          burnin = burnin,
          ident = ident,
          residual = residual,
          XGresid = XGresid
        )
    }
    
    cat("\n\n")
    
    class(out) <- c("LNIRT", "list")
    return(out)
    
  }
