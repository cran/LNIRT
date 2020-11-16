personfitLNQ <- function(RT, theta, phi, lambda, sigma2) {
  K <- ncol(RT)
  N <- nrow(RT)
  
  diff <- matrix(lambda,
                 nrow = N,
                 ncol = K,
                 byrow = T) -
    matrix(phi,
           ncol = K,
           nrow = N,
           byrow = T) * matrix(theta, ncol = K, nrow = N)
  diff <- (RT - diff) ** 2
  lZd <- diff * matrix(1 / sigma2,
                       ncol = K,
                       nrow = N,
                       byrow = TRUE)
  lZ <- (apply(lZd, 1, sum) - K) / sqrt(2 * K)
  
  dum <- (apply(lZd, 1, mean) ^ (1 / 3) - (1 - 2 / (9 * K))) / sqrt(2 /
                                                                      (9 * K))
  lZP1 <- 1 - pnorm(dum)
  lZP2 <- 1 - pnorm(lZ)
  lZP3 <- 1 - pchisq(apply(lZd, 1, sum), df = K)
  lZPT <- apply(lZd, 1, sum)
  return(list(lZPT = lZPT, lZP = lZP3))
}

residualLNQ <-
  function(RT,
           theta,
           phi,
           lambda,
           sigma2,
           EAPtheta,
           EAPlambda,
           EAPphi,
           EAPsigma2) {
    K <- ncol(RT)
    N <- nrow(RT)
    KS <- matrix(0, ncol = 1, nrow = K)
    
    #compute fitted probabilities (approximately uniformly distributed)
    
    muik <-
      matrix(lambda,
             nrow = N,
             ncol = K,
             byrow = TRUE) - matrix(phi,
                                    ncol = K,
                                    nrow = N,
                                    byrow = TRUE) * matrix(theta, ncol = K, nrow = N)
    
    #Compute Extremeness Residuals (posterior probability greater than 2)
    
    diff <- (RT - muik) * matrix(sqrt(1 / sigma2),
                                 ncol = K,
                                 nrow = N,
                                 byrow = TRUE)
    presid <- (1 - pnorm(2, mean = diff, sd = 1)) + pnorm(-2, mean = diff, sd =
                                                            1)
    
    muik <-
      matrix(EAPlambda,
             nrow = N,
             ncol = K,
             byrow = TRUE) - matrix(EAPphi,
                                    ncol = K,
                                    nrow = N,
                                    byrow = TRUE) * matrix(EAPtheta, ncol = K, nrow = N)
    muiklong <- matrix(muik, ncol = 1, nrow = N * K)
    RTlong <- matrix(RT, ncol = 1, nrow = N * K)
    sigma2long <-
      matrix(matrix(
        sqrt(EAPsigma2),
        ncol = K,
        nrow = N,
        byrow = TRUE
      ),
      ncol = 1,
      nrow = N * K)
    errorfit <- (RTlong - muiklong) / sigma2long
    errorfit <-
      matrix(rnorm(N * K) * 1e-06 + errorfit,
             ncol = K,
             nrow = N) #to remove ties
    
    #Perform one-sample Kolmogorov Smirnov Test across Items
    
    for (kk in 1:K) {
      KS[kk, 1] <- ks.test(errorfit[, kk], y = pnorm)$p.value
    }
    
    return(list(KS = KS, presid = presid))
    
  }
