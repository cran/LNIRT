#' Log-normal response time modelling
#' 
#' @param RT
#' a Person-x-Item matrix of log-response times (time spent on solving an item).
#' @param data
#' either a list or a simLNIRT object containing the response time matrix. 
#' If a simLNIRT object is provided, in the summary the simulated time parameters are shown alongside of the estimates.
#' If the RT variable cannot be found in the list, or if no data object is given, then the RT variable is taken
#' from the environment from which LNRT is called.
#' @param XG
#' the number of MCMC iterations to perform (default: 1000).
#' @param Discrimination
#' estimate the time-discrimination parameter (default: true).
#' @param WL
#' define the time-discrimination parameter as measurement error variance parameter (default: false).
#' 
#' @return 
#' an object of class LNRT.
#' 
#' @examples 
#' # Log-normal response time modelling
#' data <- simLNIRT(N = 500, K = 20, rho = 0.8, WL = FALSE)
#' summary(LNRT(RT = RT, data = data, XG = 500, Discrimination = TRUE, WL = FALSE))
#'  
#' @export
LNRT <- function(RT, data, XG = 1000, Discrimination = TRUE, WL = FALSE) {
    
  if (!missing(data)) {
    # Try to find RT in the data set first
    tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
  } else {
    data = NULL
  }
  
  N <- nrow(RT)
  K <- ncol(RT)
  
  if(WL) { 
    Discrimination <- TRUE   #WL <- 1  #time discrimination = 1/sqrt(error variance)
  }
  
  ## population theta (ability - speed)
  theta <- matrix(rnorm(N))
  muP <- matrix(0, N, 1)
  SigmaP <- diag(1)
  muP0 <- matrix(0, 1, 1)
  SigmaP0 <- diag(1)/100
  ingroup <- rep(1, N)
  Mingroup <- matrix(0, ncol = 2, nrow = N)
  
  if (XG > 1000) {
      flagged <- matrix(0, ncol = 1, nrow = N)
      ss <- 1
  }
  
  ## population item (ability - speed)
  ab <- matrix(rnorm(K * 2), ncol = 2)
  ab[, 1] <- 1
  muI <- t(matrix(rep(c(1, 0), K), ncol = K))
  muI0 <- muI[1, ]
  SigmaI <- diag(2)
  SigmaI0 <- diag(2) * 10
  for (ii in 1:2) {
      SigmaI0[ii, ] <- SigmaI0[ii, ] * c(0.5, 3)
  }
  
  ## storage
  MT <- MT2 <- array(0, dim = c(N))
  MAB <- array(0, dim = c(XG, K, 2))
  MmuP <- array(0, dim = c(XG, 1))
  MmuI <- array(0, dim = c(XG, 2))
  MSP <- array(0, dim = c(XG, 1, 1))
  MSI <- array(0, dim = c(XG, 2, 2))
  sigma2 <- rep(1, K)
  Msigma2 <- matrix(0, ncol = K, nrow = XG)
  lZP <- 0
  lZPT <- 0
  lZI <- 0
  EAPresid <- matrix(0, ncol = K, nrow = N)
  EAPKS <- matrix(0, ncol = 1, nrow = K)
  iis <- 1
  EAPphi <- matrix(0, ncol = 1, nrow = K)
  EAPlambda <- matrix(0, ncol = 1, nrow = K)
  EAPtheta <- matrix(0, ncol = 1, nrow = N)
  EAPsigma2 <- matrix(0, ncol = 1, nrow = K)
  
  DT <- matrix(1, ncol = 1, nrow = N * K)
  DT[which(is.na(RT))] <- 0
  DT <- matrix(DT, nrow = N, ncol = K)
  
  EAPCP <- matrix(0, ncol = 1, nrow = N)
  
  ## Start MCMC algorithm
  
  for (ii in 1:XG) {
      
      if (sum(DT == 0) > 0) {
          if (WL) {
              RT <- SimulateRT(RT = RT, zeta = theta, lambda = ab[, 2], phi = rep(1, K), sigma2 = sigma2, DT = DT)
          } else {
              RT <- SimulateRT(RT = RT, zeta = theta, lambda = ab[, 2], phi = ab[, 1], sigma2 = sigma2, DT = DT)
          }
      }
      
      theta <- DrawZeta(RT, ab[, 1], ab[, 2], sigma2, muP, SigmaP[1, 1])
      theta[1:N] <- theta[1:N] - mean(theta)
      MT[1:N] <- MT[1:N] + theta[1:N]
      MT2[1:N] <- MT2[1:N] + theta[1:N]^2
      
      if ((WL)) {
          # no time discrimination, 1/(sqrt(error variance)) = discrimination on MVN prior
          ab[, 1] <- 1
          ab1 <- cbind(1/sqrt(sigma2), ab[, 2])
          dum <- Conditional(kk = 2, Mu = muI, Sigma = SigmaI, Z = ab1)
          ab[, 2] <- DrawLambda_LNRT(RT = RT, zeta = theta, sigma2 = sigma2, mu = dum$CMU, sigma = dum$CVAR[1,1])
      } else {
          if (Discrimination) {
              dum <- DrawLambdaPhi_LNRT(RT, theta, sigma2, muI, SigmaI, ingroup)
              ab[, 1] <- dum$phi
              ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
              ab[, 2] <- dum$lambda
          } else {
              ab[, 1] <- rep(1, K)
              ab[, 2] <- DrawLambda_LNRT(RT = RT, zeta = theta, sigma2 = sigma2, mu = muI[1,2], sigma = SigmaI[2,2])
          }
          
      }
      
      MAB[ii, 1:K, 1:2] <- ab
      sigma2 <- SampleS_LNRT(RT, theta, ab[, 2], ab[, 1], ingroup)
      if (WL) {
          Msigma2[ii, 1:K] <- 1/sqrt(sigma2)
      } else {
          Msigma2[ii, 1:K] <- sigma2
      }
      
      X <- matrix(1, N, 1)
      muP <- SampleB(theta, X, SigmaP, muP0, SigmaP0)
      MmuP[ii, ] <- muP$B
      muP <- muP$pred
      
      SS <- crossprod(theta - muP) + SigmaP0
      SigmaP <- rwishart(1 + N, chol2inv(chol(SS)))$IW
      MSP[ii, , ] <- SigmaP
      
      X <- matrix(1, K, 1)
      if (WL) {
          ab1 <- cbind(1/sqrt(sigma2), ab[, 2])
      } else {
          ab1 <- ab
      }
      # muI2 <- SampleB(ab1, X, SigmaI, muI0, SigmaI0)
      # MmuI[ii, 1] <- muI2$B[1]
      # MmuI[ii, 2] <- muI2$B[2]
      # muI[, 1] <- muI2$pred[, 1]
      # muI[, 2] <- muI2$pred[, 2]
      # 
      # SS <- crossprod(ab1 - muI) + SigmaI0
      # SigmaI <- rwishart(2 + K, chol2inv(chol(SS)))$IW
      # MSI[ii, , ] <- SigmaI
      
      
      ## Adjust sampling mean item parameters 
      if((!WL) && (!Discrimination)){
        meanmuI2  <- (sum(ab[,2])/SigmaI[2,2]+muI0[2]/SigmaI0[2,2])/(1/SigmaI0[2,2]+1/(SigmaI[2,2]))		
        sdmuI2  <- sqrt(1/(1/SigmaI0[2,2]+K/(SigmaI[2,2])))	
        muI2	<- rnorm(1,mean=meanmuI2,sd=sdmuI2)
        MmuI[ii,1] <- 1
        MmuI[ii,2] <- muI2
        muI[,1] <- 1
        muI[,2] <- muI2
      }else{
        muI2 <- SampleB(ab1,X,SigmaI,muI0,SigmaI0)
        MmuI[ii,1] <- muI2$B[1]
        MmuI[ii,2] <- muI2$B[2]
        muI[,1] <- muI2$pred[,1]
        muI[,2] <- muI2$pred[,2]
      }
      
      ## Adjust sampling covariance matrix item parameters 
      if(!Discrimination){
        SS <- sum((ab1[,2]-muI[,2])**2) + SigmaI0[2,2]
        SigmaI[2,2] <- SS/rgamma(1,(K+1)/2,1/2) #change sampling
        SigmaI[1,1] <- 1
        MSI[ii,,] <- SigmaI
      }else{		
        SS <- crossprod(ab1 - muI) + SigmaI0
        SigmaI <- rwishart(2 + K, chol2inv(chol(SS)))$IW
        MSI[ii,,] <- SigmaI
      }
      
      
      if (ii > 1000) {
          EAPphi <- (ab[, 1] + (iis - 1) * EAPphi)/iis
          EAPlambda <- (ab[, 2] + (iis - 1) * EAPlambda)/iis
          EAPtheta <- (theta + (iis - 1) * EAPtheta)/iis
          EAPsigma2 <- (sigma2 + (iis - 1) * EAPsigma2)/iis
          
          dum <- personfitLN(RT = RT, theta = theta, phi = ab[, 1], lambda = ab[, 2], sigma2 = sigma2)
          lZP <- lZP + dum$lZP
          lZPT <- lZPT + dum$lZPT
          CF <- ifelse(dum$lZP < 0.05, 1, 0)  #significance level = .05
          EAPCP <- (CF + (iis - 1) * EAPCP)/iis
          
          dum <- itemfitLN(RT = RT, theta = theta, phi = ab[, 1], lambda = ab[, 2], sigma2 = sigma2)
          lZI <- lZI + dum$lZI
          
          dum <- residualLN(RT = RT, theta = theta, phi = ab[, 1], lambda = ab[, 2], sigma2 = sigma2, EAPtheta = EAPtheta, EAPlambda = EAPlambda, 
              EAPphi = EAPphi, EAPsigma2 = EAPsigma2)
          EAPresid <- EAPresid + dum$presid
          EAPKS <- (dum$KS[1:K, 1] + (iis - 1) * EAPKS)/iis
          
          iis <- iis + 1
      }
      
      if (ii%%100 == 0) 
          cat("Iteration ", ii, " ", "\n")
      flush.console()
  }
  
  MT <- MT/XG
  MT2 <- sqrt(MT2/XG - MT^2)
  if (ii > 1000) {
      lZP <- lZP/(XG - 1000)
      lZPT <- lZPT/(XG - 1000)
      lZI <- lZI/(XG - 1000)
      EAPresid <- EAPresid/(XG - 1000)
  }
  
  if (XG > 1000) {
      out <- list(Mtheta = MT, MTSD = MT2, MAB = MAB, MmuP = MmuP, MSP = MSP, MmuI = MmuI, MSI = MSI, lZP = lZP, lZPT = lZPT, Msigma2 = Msigma2, 
          theta = theta, sigma2 = sigma2, lZI = lZI, EAPresid = EAPresid, EAPKS = EAPKS, RT = RT, EAPCP = EAPCP, Discrimination = Discrimination, WL = WL, data = data)
  } else {
      out <- list(Mtheta = MT, MTSD = MT2, MAB = MAB, MmuP = MmuP, MSP = MSP, MmuI = MmuI, MSI = MSI, Msigma2 = Msigma2, theta = theta, sigma2 = sigma2, 
          RT = RT, Discrimination = Discrimination, WL = WL, data = data)
  }
  
  class(out) <- "LNRT"
  return(out)
}
