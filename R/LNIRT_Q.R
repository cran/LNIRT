## March 2015
## Variable Speed Model
## Code from LNRT (Log-Normal Response Time) plus random person effects (intercept, trend, quadratic)
## Use functions from LNIRT

#' Log-normal response time IRT modelling with variable person speed (intercept, trend, quadratic)
#'
#' @param RT
#' a Person-x-Item matrix of log-response times (time spent on solving an item).
#' @param Y
#' a Person-x-Item matrix of responses.
#' @param data
#' either a list or a simLNIRTQ object containing the response time and response matrices
#' and optionally the predictors for the item and person parameters.
#' If a simLNIRTQ object is provided, in the summary the simulated item and time parameters are shown alongside of the estimates.
#' If the required variables cannot be found in the list, or if no data object is given, then the variables are taken
#' from the environment from which LNIRTQ is called.
#' @param X
#' explanatory (time) variables for random person speed (default: (1:N.items - 1)/N.items).
#' @param XG
#' the number of MCMC iterations to perform (default: 1000).
#' @param burnin
#' the percentage of MCMC iterations to discard as burn-in period (default: 10).
#' @param XGresid
#' the number of MCMC iterations to perform before residuals are computed (default: 1000).
#' @param residual
#' compute residuals, >1000 iterations are recommended (default: false).
#'
#' @return
#' an object of class LNIRTQ.
#' @export
LNIRTQ <-
  function(Y,
           RT,
           X,
           data,
           XG = 1000,
           burnin = 10,
           XGresid = 1000,
           residual = FALSE) {
    ## Y = response matrix of dim(N=persons,K=items)
    ## RT = log-response time matrix (time spent on solving an item) of dim(N=persons,K=items)
    ## XG = number of XG iterations for the MCMC algorithm
    ## X = explanatory (time) variables for random person speed
    
    if (XG <= 0) {
      print("Error: XG must be > 0.")
      return (NULL)
    }
    if ((burnin <= 0) || (burnin >= XG)) {
      print("Error: burn-in period must be between 0% and 100%.")
      return (NULL)
    }
    if (residual && (XGresid >= XG || XGresid <= 0)) {
      print("Warning: XGresid must be < XG and > 0. Residuals will not be computed.")
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
      
      # Try to find explanatory (time) variables in the data set first
      tryCatch(
        X <- eval(substitute(X), data),
        error = function(e)
          NULL
      )
    } else {
      data <- NULL
    }
    
    Y <- as.matrix(Y)
    RT <- as.matrix(RT)
    
    if ((nrow(Y) != nrow(RT)) || (ncol(Y) != ncol(RT))) {
      print("Error: Y and RT must be of equal dimension.")
      return (NULL)
    }
    
    ## Initialise all parameters
    N <- nrow(Y) #persons
    K <- ncol(Y) #items (complete design)
    Q <- 4 #number of speed components plus ability
    
    if (missing(X)) {
      X <- 1:K
      X <- (X - 1) / K
    }
    
    cat (" \n")
    cat ("   LNIRT v", packageDescription("LNIRT")$Version, "\n", sep = "")
    cat ("   ", rep('-', 20), "\n\n", sep = "")
    # cat ("   Jean-Paul Fox \n")
    # cat ("   Konrad Klotzke \n")
    # cat ("   ", rep('-', 20), "\n\n", sep = "")
    
    #cat ("   ", rep('-', 40), "\n", sep = "")
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
    #cat ("   ", rep('-', 40), "\n\n", sep = "")
    
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
    
    ## population theta (ability - speed)
    theta <- matrix(rnorm(N), ncol = 1)
    zeta <- matrix(rnorm((Q - 1) * N), nrow = N, ncol = Q - 1)
    
    muP <- matrix(0, ncol = Q, nrow = N)
    SigmaP <-  diag(Q)
    muP0 <- matrix(0, 1, Q)
    SigmaP0 <-  diag(Q) / 10
    SigmaR <- diag(Q - 1) #covariance matrix speed components
    
    ##population item (ability - speed)
    ab <- matrix(rnorm(K * 4), ncol = 4)
    ab[, c(1, 3)] <- 1
    muI <- t(matrix(rep(c(1, 0), 2 * K), ncol = K))
    muI0 <- muI[1, ]
    SigmaI <- diag(4)
    SigmaI0 <- diag(4) / 10
    for (ii in 1:4) {
      SigmaI0[ii, ] <- SigmaI0[ii, ] * rep(c(.5, 3), 2)
    }
    
    ##storage
    MT <- array(0, dim = c(N, Q))
    MAB <- array(0, dim = c(XG, K, 4))
    MmuP <- array(0, dim = c(XG, Q))
    MmuI <- array(0, dim = c(XG, 4))
    MSP <- array(0, dim = c(XG, Q, Q))
    MSI <- array(0, dim = c(XG, 4, 4))
    sigma2 <- rep(1, K)
    Msigma2 <- matrix(0, nrow = XG, ncol = K)
    
    lZP <- 0
    lZPT <- 0
    
    EAPphi <- matrix(0, ncol = 1, nrow = K)
    EAPlambda <- matrix(0, ncol = 1, nrow = K)
    EAPzetapred <- matrix(0, ncol = K, nrow = N)
    EAPsigma2	<- matrix(0, ncol = 1, nrow = K)
    EAPresid <- matrix(0, ncol = K, nrow = N)
    EAPKS <- matrix(0, ncol = 1, nrow = K)
    
    iis <- 1
    
    D 	<- matrix(1, ncol = 1, nrow = N * K) # no missings
    D[which(is.na(Y))] <- 0 #identify missings
    D <- matrix(D, nrow = N, ncol = K)
    
    DT 	<- matrix(1, ncol = 1, nrow = N * K) # no missings
    DT[which(is.na(RT))] <- 0 #identify missings
    DT <- matrix(DT, nrow = N, ncol = K)
    
    ## Start MCMC algorithm
    
    for (ii in 1:XG) {
      #ability test
      ZR <- DrawZQ(
        alpha0 = ab[, 1],
        beta0 = ab[, 2],
        theta0 = theta,
        S = Y,
        D = D
      )
      
      ## replace mising response times
      if (sum(DT == 0) > 0) {
        RT <-
          SimulateRTQ(
            RT = RT,
            zeta = zeta,
            lambda = ab[, 4],
            phi = ab[, 3],
            sigma2 = sigma2,
            X = X,
            DT = DT
          )
      }
      
      ## Sample all the person parameters
      SigmaR <- SigmaP[2:4, 2:4]
      mu <-
        muP[, 1] + t((SigmaP[1, 2:4] %*% solve(SigmaR)) %*% t(zeta - muP[, 2:4]))
      sigma <-
        SigmaP[1, 1] - SigmaP[1, 2:4] %*% solve(SigmaR) %*% SigmaP[2:4, 1]
      theta <-
        DrawThetaQ(
          alpha0 = ab[, 1],
          beta0 = ab[, 2],
          Z = ZR,
          mu = mu,
          sigma = sigma
        )
      
      mu <-
        muP[, 2:4] + t(SigmaP[2:4, 1] %*% t(theta - muP[, 1])) / (SigmaP[1, 1])
      sigma <- SigmaR - (SigmaP[1, 2:4] %o% SigmaP[1, 2:4]) / SigmaP[1, 1]
      dum <-
        DrawZetaQ(
          RT = RT,
          phi = ab[, 3],
          lambda = ab[, 4],
          sigma2 = sigma2,
          mu = mu,
          SigmaR = sigma,
          X = X
        ) ## variable speed
      zeta <- dum$zeta
      zetapred <- dum$zetapred
      
      ### restriction average speed is equal to zero
      zeta <- zeta - matrix(
        apply(zeta, 2, mean),
        ncol = Q - 1,
        nrow = N,
        byrow = T
      )
      
      dum <- ConditionalQ(1, muI, SigmaI, ab)
      ab[, 1] <-
        abs(
          DrawPhiQ(
            T = ZR,
            lambda = -ab[, 2],
            zeta = -theta,
            sigma2 = rep(1, K),
            mu = dum$CMU,
            sigmal = dum$CVAR
          )
        )
      ab[, 1] <-
        ab[, 1] / (prod(ab[, 1]) ^ (1 / K)) #restrict item discrimination
      
      dum <- ConditionalQ(2, muI, SigmaI, ab)
      ab[, 2] <-
        -DrawLambdaQ(
          T = ZR,
          phi = -ab[, 1],
          zeta = theta,
          sigma2 = rep(1, K),
          mu = dum$CMU,
          sigmal = dum$CVAR
        )$lambda
      ab[, 2] <- ab[, 2] - mean(ab[, 2])		#restrict mean item difficulty
      dum <-
        DrawLambdaPhiQ(
          RT = RT,
          zeta = zeta,
          X = X,
          sigma2 = sigma2,
          muI = muI[, 3:4],
          SigmaI = SigmaI[3:4, 3:4]
        )
      ab[, 3] <- dum$phi
      ab[, 3] <- ab[, 3] / (prod(abs(ab[, 3])) ^ (1 / K))
      ab[, 4] <- dum$lambda
      ##ab[,4] <- ab[,4] - mean(ab[,4]) ### restriction average time-intensity equal to zero
      MAB[ii, 1:K, 1:4] <- ab
      sigma2 <-
        SampleS2Q(
          RT = RT,
          zeta = zeta,
          X = X,
          lambda = ab[, 4],
          phi = ab[, 3]
        )
      Msigma2[ii, ] <- sigma2
      
      ## 2nd level models for person and item parameters.
      X1 <- matrix(1, N, 1)
      thetazeta <- cbind(theta, zeta)
      dum <-
        SampleBQ(
          Y = thetazeta,
          X = X1,
          Sigma = SigmaP,
          B0 = muP0,
          V0 = SigmaP0
        )
      MmuP[ii, ] <- dum$B
      muP <- dum$pred
      
      ## Sample covmatrix persons (speed)
      SS <- diag(crossprod(zeta - dum$pred[, 2:Q]))
      SigmaR <- diag((SS + 1) / rchisq(Q - 1, N))
      SigmaP[2:Q, 2:Q] <- SigmaR
      
      #if(ii %% 100 == 0) cat("Iteration ",ii,"\n")
      #flush.console()
      
      ## Sample variance ability
      SigmaP <- DrawRhoQ(
        zeta = zeta,
        theta = theta,
        muP = muP,
        SigmaP = SigmaP
      )
      MSP[ii, , ] <- SigmaP
      
      X2 <- matrix(1, K, 1)
      muI2 <- SampleBQ(ab, X2, SigmaI, muI0, SigmaI0)
      MmuI[ii, c(1, 3, 4)] <- muI2$B[c(1, 3, 4)]
      muI[, c(1, 3, 4)] <- muI2$pred[, c(1, 3, 4)]
      muI1 <- matrix(muI,
                     ncol = 4,
                     nrow = K,
                     byrow = FALSE)
      SS <- crossprod(ab - muI1) + SigmaI0
      if (1 / det(SS) < 10e-10) {
        SSn <- diag(SS)
        SigmaI <- rwishart(4 + K, chol2inv(chol(SSn)))$IW
      } else{
        SigmaI <- rwishart(4 + K, chol2inv(chol(SS)))$IW
      }
      
      MSI[ii, , ] <- SigmaI
      MT[1:N, 1:4] <- MT[1:N, 1:4] + cbind(theta, zeta)
      
      if (ii > XGresid && residual) {
        EAPlambda <- (ab[, 4] + (iis - 1) * EAPlambda) / iis
        EAPzetapred <- (zetapred + (iis - 1) * zetapred) / iis
        EAPphi <- (ab[, 3] + (iis - 1) * EAPphi) / iis
        EAPsigma2 <- (sigma2 + (iis - 1) * EAPsigma2) / iis
        
        ## Log-Normal Fit Evaluation
        dum <-
          personfitLNQ(
            RT = RT,
            theta = zetapred,
            phi = ab[, 3],
            lambda = ab[, 4],
            sigma2 = sigma2
          )	# lZ statistic
        lZP <- lZP + dum$lZP
        lZPT <- lZPT + dum$lZPT
        
        dum <-
          residualLNQ(
            RT = RT,
            theta = zetapred,
            phi = ab[, 3],
            lambda = ab[, 4],
            sigma2 = sigma2,
            EAPtheta = EAPzetapred,
            EAPlambda = EAPlambda,
            EAPphi = EAPphi,
            EAPsigma2 = EAPsigma2
          )
        EAPresid <- EAPresid + dum$presid
        EAPKS <- (dum$KS[1:K, 1] + (iis - 1) * EAPKS) / iis
        
        iis <- iis + 1
      }
      
      # Update progress bar
      setTxtProgressBar(pb, ii)
      
      # if (ii%%100 == 0)
      #     cat("Iteration ", ii, " ", "\n")
      # flush.console()
      
      
    }## end iterations
    
    MT <- MT / XG
    ##MT2 <- sqrt(MT2/XG - MT^2 ) ## calculate posterior SD of person parameters.
    
    if (ii > XGresid && residual) {
      lZP <- lZP / (XG - XGresid)
      lZPT <- lZPT / (XG - XGresid)
      EAPresid <- EAPresid / (XG - XGresid)
    }
    
    if (!(any(class(data) == "simLNIRTQ"))) {
      data <- NULL # only attach sim data for summary function
    }
    
    if (XG > XGresid && residual) {
      out <-
        list(
          Mtheta = MT,
          MAB = MAB,
          MmuP = MmuP,
          MSP = MSP,
          MmuI = MmuI,
          MSI = MSI,
          Msigma2 = Msigma2,
          lZP = lZP,
          lZPT = lZPT,
          EAPresid = EAPresid,
          EAPKS = EAPKS,
          XG = XG,
          burnin = burnin,
          XGresid = XGresid,
          data = data
        )
    } else{
      out <-
        list(
          Mtheta = MT,
          MAB = MAB,
          MmuP = MmuP,
          MSP = MSP,
          MmuI = MmuI,
          MSI = MSI,
          Msigma2 = Msigma2,
          XG = XG,
          burnin = burnin,
          XGresid = XGresid,
          data = data
        )
    }
    
    cat("\n\n")
    
    class(out) <- c("LNIRTQ", "list")
    return(out)
  }
