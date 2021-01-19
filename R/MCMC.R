### Common functions used in MCMC ###

Conditional <- function(kk, Mu, Sigma, Z) {
  K <- ncol(Z)
  N <- nrow(Z)
  
  if (kk == 1) {
    C <- matrix(Z[, 2:K] - Mu[, 2:K], ncol = (K - 1))
    CMEAN <-
      Mu[, 1] + Sigma[1, 2:K] %*% solve(Sigma[2:K, 2:K]) %*% t(C)
    CSD <-
      Sigma[1, 1] - Sigma[1, 2:K] %*% solve(Sigma[2:K, 2:K]) %*% Sigma[2:K, 1]
  }
  if (kk > 1) {
    if (kk < K) {
      C <- matrix(Z[1:N, 1:(kk - 1)] - Mu[, 1:(kk - 1)], ncol = (kk - 1))
      CMu1 <-
        Mu[, kk:K] + t(matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk - 1), 1:(kk - 1)]) %*% t(C))
      CSigma <-
        Sigma[kk:K, kk:K] - matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk - 1), 1:(kk - 1)]) %*% matrix((Sigma[1:(kk -
                                                                                                                                                1), kk:K]), nrow = (kk - 1))
      J <- ncol(CSigma)
      C <-
        matrix(Z[1:N, (kk + 1):K] - CMu1[1:N, 2:J], ncol = (J - 1))
      CMEAN <-
        CMu1[, 1] + CSigma[1, 2:J] %*% solve(CSigma[2:J, 2:J]) %*% t(C)
      CSD <-
        CSigma[1, 1] - CSigma[1, 2:J] %*% solve(CSigma[2:J, 2:J]) %*% CSigma[2:J, 1]
    }
    if (kk == K) {
      C <- matrix(Z[1:N, 1:(K - 1)] - Mu[, 1:(K - 1)], ncol = (K - 1))
      CMEAN <-
        Mu[, K] + t(Sigma[1:(K - 1), K]) %*% solve(Sigma[1:(K - 1), 1:(K - 1)]) %*% t(C)
      CSD <-
        Sigma[K, K] - t(Sigma[1:(K - 1), K]) %*% solve(Sigma[1:(K - 1), 1:(K - 1)]) %*% Sigma[1:(K - 1), K]
    }
  }
  return(list(CMU = CMEAN, CVAR = CSD))
}

SimulateRT <- function(RT, zeta, lambda, phi, sigma2, DT) {
  # assume MAR
  N <- nrow(RT)
  K <- ncol(RT)
  meanT <-
    t(matrix(lambda, K, N)) - t(matrix(phi, K, N)) * (zeta %*% t(rep(1, K)))
  meanT <- matrix(meanT, ncol = 1, nrow = N * K)
  sigmaL <-
    matrix(t(matrix(
      rep(sqrt(sigma2), N), nrow = K, ncol = N
    )), ncol = 1, nrow = N * K)
  RT <- matrix(RT, ncol = 1, nrow = N * K)
  RT[which(DT == 0)] <-
    rnorm(sum(DT == 0), mean = meanT[which(DT == 0)], sd = sigmaL[which(DT == 0)])
  RT <- matrix(RT, ncol = K, nrow = N)
  
  return(RT)
}


SimulateRTMBD <- function(RT, zeta, lambda, phi, sigma2, DT, MBDT) {
  #assume MAR for NAs
  #not when MBDT==0 (missing by design)
  N <- nrow(RT)
  K <- ncol(RT)
  meanT <-
    t(matrix(lambda, K, N)) - t(matrix(phi, K, N)) * (zeta %*% t(rep(1, K)))
  meanT <- matrix(meanT, ncol = 1, nrow = N * K)
  sigmaL <-
    matrix(t(matrix(
      rep(sqrt(sigma2), N), nrow = K, ncol = N
    )), ncol = 1, nrow = N * K)
  RT <- matrix(RT, ncol = 1, nrow = N * K)
  RT[which(DT == 0 &
             MBDT == 1)] <-
    rnorm(sum(DT == 0 &
                MBDT == 1), mean = meanT[which(DT == 0 &
                                                 MBDT == 1)], sd = sigmaL[which(DT == 0 & MBDT == 1)])
  RT <- matrix(RT, ncol = K, nrow = N)
  
  return(RT)
}


SimulateY <- function(Y, theta, alpha0, beta0, guess0, D) {
  # with guessing
  N <- nrow(Y)
  K <- ncol(Y)
  
  G <- matrix(0, ncol = K, nrow = N)
  
  for (kk in 1:K) {
    G[, kk] <- rbinom(N, size = 1, prob = guess0[kk])
  }
  Y[(D == 0) & (G == 1)] <- 1  #missing: guessed correctly
  
  par <-
    theta %*% matrix(alpha0, nrow = 1, ncol = K) - t(matrix(beta0, nrow = K, ncol = N))
  probs <- matrix(pnorm(par), ncol = K, nrow = N)
  Yn <- matrix(runif(N * K), nrow = N, ncol = K)
  Yn <- ifelse(Yn < probs, 1, 0)
  Y[(D == 0) &
      (G == 0)] <- Yn[(D == 0) & (G == 0)]  #missing: response generated
  
  return(Y)
}

SimulateYMBD <- function(Y, theta, alpha0, beta0, guess0, D, MBDY) {
  #with guessing
  #with missing by design
  N <- nrow(Y)
  K <- ncol(Y)
  
  G <- matrix(0, ncol = K, nrow = N)
  
  for (kk in 1:K) {
    G[, kk] <- rbinom(N, size = 1, prob = guess0[kk])
  }
  Y[(D == 0 &
       MBDY == 1) & (G == 1)] <-
    1 #missing not by design: guessed correctly
  
  par <-
    theta %*% matrix(alpha0, nrow = 1, ncol = K) - t(matrix(beta0, nrow = K, ncol =
                                                              N))
  probs <- matrix(pnorm(par), ncol = K, nrow = N)
  Yn  <- matrix(runif(N * K), nrow = N, ncol = K)
  Yn <- ifelse(Yn < probs, 1, 0)
  Y[(D == 0 &
       MBDY == 1) &
      (G == 0)] <-
    Yn[(D == 0 & MBDY == 1) & (G == 0)] #missing: response generated
  
  return(Y)
}

DrawZeta <- function(RT, phi, lambda, sigma2, mu, sigmaz) {
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- matrix(lambda, nrow = K, ncol = N) - t(RT)
  X <- matrix(phi, K, 1)
  sigma2inv <- diag(1 / sigma2[1:K])
  vartheta <- (1 / ((t(phi) %*% sigma2inv) %*% phi + 1 / sigmaz))[1, 1]
  meantheta <-
    matrix(((t(phi) %*% sigma2inv) %*% Z + t(mu / sigmaz)) * vartheta, ncol = 1, nrow = N)
  zeta <-
    matrix(rnorm(N, mean = meantheta, sd = sqrt(vartheta)),
           ncol = 1,
           nrow = N)
  
  return(zeta)
}

DrawZetaMBD <- function(RT, phi, lambda, sigma2, mu, sigmaz, MBDT) {
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- (matrix(lambda, nrow = K, ncol = N) - t(RT)) * t(MBDT)
  X <- matrix(phi, K, 1)
  sigma2inv <- diag(1 / sigma2[1:K])
  
  mphi <- matrix(phi,
                 ncol = K,
                 nrow = N,
                 byrow = T) * MBDT
  vartheta <- 1 / (mphi %*% sigma2inv %*% phi + 1 / sigmaz)
  meantheta <-
    matrix(((t(phi) %*% sigma2inv) %*% Z + t(mu / sigmaz)) * t(vartheta), ncol =
             1, nrow = N)
  zeta <-
    matrix(rnorm(N, mean = meantheta, sd = sqrt(vartheta)),
           ncol = 1,
           nrow = N)
  
  return(zeta)
}


SampleB <- function(Y, X, Sigma, B0, V0) {
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x% solve(crossprod(X))) + solve(V0))
  Btilde <-
    Bvar %*% (
      solve(Sigma %x% diag(m)) %*% matrix(crossprod(X, Y), ncol = 1) + solve(V0) %*% matrix(B0, ncol = 1)
    )
  B <-
    Btilde + chol(Bvar) %*% matrix(rnorm(length(Btilde)), ncol = 1)
  pred <- X %*% matrix(B, ncol = k)
  return(list(B = B, pred = pred))
}

rwishart <- function(nu, V) {
  m = nrow(V)
  df = (nu + nu - m + 1) - (nu - m + 1):nu
  if (m > 1) {
    RT = diag(sqrt(rchisq(c(rep(
      1, m
    )), df)))
    RT[lower.tri(RT)] = rnorm((m * (m + 1) / 2 - m))
  } else {
    RT = sqrt(rchisq(1, df))
  }
  U = chol(V)
  C = t(RT) %*% U
  CI = backsolve(C, diag(m))
  return(list(IW = crossprod(t(CI))))
}


### MCMC functions used in LNRT ###

DrawLambdaPhi_LNRT <-
  function(RT, theta, sigma2, muI, sigmaI, ingroup) {
    # library(MASS)
    K <- ncol(RT)
    N <- nrow(RT)
    invSigmaI <- solve(sigmaI)
    H <- matrix(c(-theta, rep(1, N)), ncol = 2, nrow = N) * ingroup
    varest <-
      solve(kronecker(diag(1 / sigma2[1:K]), (t(H) %*% H)) + kronecker(diag(1, K), invSigmaI))
    
    meanest <-
      t((t(H) %*% RT) / (t(matrix(
        sigma2, nrow = K, ncol = 2
      ))) + matrix(t(muI[1,] %*% invSigmaI), ncol = K, nrow = 2))
    meanest <-
      apply((matrix(rep(meanest, K), ncol = 2 * K) %*% varest) * t(kronecker(diag(K), c(1, 1))), 2, sum)
    
    lambdaphi <- MASS::mvrnorm(1, mu = meanest, Sigma = varest)
    lambdaphi <- matrix(lambdaphi,
                        ncol = 2,
                        nrow = K,
                        byrow = RT)
    
    set <- which(lambdaphi[, 1] < 0.3)
    if (length(set) >= 1) {
      lambdaphi[set, 1] <- 0.3
    }
    
    return(list(phi = lambdaphi[, 1], lambda = lambdaphi[, 2]))
  }

SampleS_LNRT <- function(RT, zeta, lambda, phi, ingroup) {
  K <- ncol(RT)
  N <- nrow(RT)
  Nr <- sum(ingroup)
  ss0 <- 10
  ingroup <- matrix(ingroup, ncol = K, nrow = N)
  Z <-
    (RT + t(matrix(phi, K, N)) * matrix(zeta, N, K) - t(matrix(lambda, K, N))) * ingroup
  sigma2 <- (apply(Z * Z, 2, sum) + ss0) / rchisq(K, Nr)
  return(sigma2)
}

DrawLambda_LNRT <- function(RT, zeta, sigma2, mu, sigma) {
  # prior mu,sigma
  N <- nrow(RT)
  K <- ncol(RT)
  zetaT <- matrix(zeta,
                  ncol = K,
                  nrow = N,
                  byrow = F)
  RTs <- matrix(RT + zetaT, ncol = K, nrow = N)
  XX <- matrix(1, ncol = K, nrow = N)
  pvar <- diag(t(XX) %*% XX) / sigma2 + 1 / sigma
  betahat <- diag(t(XX) %*% RT) / sigma2
  mu <- (betahat + mu / sigma) / pvar
  beta <- rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(beta)
}

SampleBX_LNRT <- function(Y, XPT) {
  N <- nrow(Y)
  kt <- ncol(XPT)
  V0 <- diag(kt) #inverse prior covariance matrix
  B0 <- matrix(0, ncol = 1, nrow = kt)
  V0B0 <- matrix(V0 %*% matrix(B0, ncol = 1), ncol = 1)
  
  Bvart <- solve(crossprod(XPT) + V0)
  Btildet <- Bvart %*% (crossprod((XPT), Y) + V0B0)
  Bt <-
    matrix(MASS::mvrnorm(1, mu = Btildet, Sigma = Bvart), nrow = 1)
  B <- Bt
  pred <- XPT %*% matrix(Bt, nrow = kt, ncol = 1)
  
  return(list(B = B, pred = pred))
}

### MCMC functions used in LNIRT ###

DrawS_LNIRT <- function(alpha0, beta0, guess0, theta0, Y) {
  N <- nrow(Y)
  K <- ncol(Y)
  eta <-
    t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow = N) - t(matrix(beta0, nrow = K, ncol = N))
  eta <- matrix(pnorm(eta), ncol = K, nrow = N)
  probS <-
    eta / (eta + t(matrix(
      guess0, nrow = K, ncol = N
    )) * (matrix(1, ncol = K, nrow = N) - eta))
  S <- matrix(runif(N * K), ncol = K, nrow = N)
  S <- matrix(ifelse(S > probS, 0, 1), ncol = K, nrow = N)
  S <- S * Y
  return(S)
}

DrawSMBD_LNIRT <- function(alpha0, beta0, guess0, theta0, Y, MBDY) {
  N  	<- nrow(Y)
  K		<- ncol(Y)
  eta 	<-
    t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow =
                                                     N) - t(matrix(beta0, nrow = K, ncol = N))
  eta 	<- matrix(pnorm(eta), ncol = K, nrow = N)
  probS 	<-
    eta / (eta + t(matrix(
      guess0, nrow = K, ncol = N
    )) * (matrix(1, ncol = K, nrow = N) - eta))
  S		<- matrix(runif(N * K), ncol = K, nrow = N)
  S		<- matrix(ifelse(S > probS, 0, 1), ncol = K, nrow = N)
  S[which(MBDY == 1)] <-
    S[which(MBDY == 1)] * Y[which(MBDY == 1)] #not missing by design
  return(S)
}


DrawZ_LNIRT <- function(alpha0, beta0, theta0, S, D) {
  N <- nrow(S)
  K <- ncol(S)
  eta <-
    t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow = N) - t(matrix(beta0, ncol = N, nrow = K))
  BB <- matrix(pnorm(-eta), ncol = K, nrow = N)
  BB[which(BB < 1e-05)] <- 1e-05
  BB[which(BB > (1 - 1e-05))] <- (1 - 1e-05)
  u <- matrix(runif(N * K), ncol = K, nrow = N)
  tt <-
    matrix((BB * (1 - S) + (1 - BB) * S) * u + BB * S, ncol = K, nrow = N)
  Z <- matrix(0, ncol = 1, nrow = N * K)
  tt <- matrix(tt, ncol = 1, nrow = N * K)
  eta <- matrix(eta, ncol = 1, nrow = N * K)
  Z[which(D == 1)] <-
    qnorm(tt)[which(D == 1)] + eta[which(D == 1)]
  Z[which(D == 0)] <- rnorm(sum(D == 0)) + eta[which(D == 0)]
  Z <- matrix(Z, ncol = K, nrow = N)
  
  return(Z)
}

DrawZMBD_LNIRT   <- function(alpha0, beta0, theta0, S, D, MBDY) {
  ##Y[MBDY==0] <- 0 #avoid problems with NAs
  N 		<- nrow(S)
  K 		<- ncol(S)
  eta	<-
    t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow =
                                                     N) -
    t(matrix(beta0, ncol = N, nrow = K))
  BB		<- matrix(pnorm(-eta), ncol = K, nrow = N)
  BB[which(BB < .00001)] <- .00001
  BB[which(BB > (1 - .00001))] <- (1 - .00001)
  u		<- matrix(runif(N * K), ncol = K, nrow = N)
  tt		<-
    matrix((BB * (1 - S) + (1 - BB) * S) * u + BB * S, ncol = K, nrow = N)
  Z  <- matrix(0, ncol = 1, nrow = N * K)
  tt <- matrix(tt, ncol = 1, nrow = N * K)
  eta <- matrix(eta, ncol = 1, nrow = N * K)
  Z[which(D == 1)] <-  qnorm(tt)[which(D == 1)] + eta[which(D == 1)]
  Z[which(D == 0)] <-  rnorm(sum(D == 0)) + eta[which(D == 0)]
  Z <- matrix(Z, ncol = K, nrow = N)
  
  Z[which(MBDY == 0)] <-  0 ## missing by design (no imputation)
  
  return(Z)
}

DrawTheta_LNIRT <- function(alpha0, beta0, Z, mu, sigma) {
  N <- nrow(Z)
  K <- ncol(Z)
  pvar <- (sum(alpha0 ^ 2) + 1 / as.vector(sigma))
  thetahat <-
    (((Z + t(
      matrix(beta0, ncol = N, nrow = K)
    )) %*% matrix(
      alpha0, ncol = 1, nrow = K
    )))
  mu <- (thetahat + as.vector(mu) / as.vector(sigma)) / pvar
  theta <- rnorm(N, mean = mu, sd = sqrt(1 / pvar))
  
  return(theta)
}


DrawThetaMBD_LNIRT <- function(alpha0, beta0, Z, mu, sigma, MBDY) {
  ## account for missing by design
  N			<- nrow(Z)
  K 			<- ncol(Z)
  pvar 		<-
    apply((matrix(
      alpha0 ^ 2,
      ncol = K,
      nrow = N,
      byrow = T
    ) * MBDY), 1, sum) + 1 / as.vector(sigma) #correct for MBDY
  thetahat	<-
    (((Z + t(
      matrix(beta0, ncol = N, nrow = K)
    )) %*% matrix(
      alpha0, ncol = 1, nrow = K
    )))
  mu			<- (thetahat + as.vector(mu) / as.vector(sigma)) / pvar
  theta		<- rnorm(N, mean = mu, sd = sqrt(1 / pvar))
  
  return(theta)
}


DrawC_LNIRT <- function(S, Y) {
  N <- nrow(Y)
  K <- ncol(Y)
  Q1 <-
    20 + apply((S == 0) * (Y == 1), 2, sum)  #answer unknown and guess correctly
  Q2 <- 80 + apply(S == 0, 2, sum)  #answer unknown
  guess <- rbeta(K, Q1, Q2)
  
  return(guess)
}


DrawCMBD_LNIRT <- function(S, Y, MBDY) {
  N		<- nrow(Y)
  K		<- ncol(Y)
  S[MBDY == 0]	<- 9 #missing by design
  Y[MBDY == 0]	<- 9 #missing by design
  Q1		<-
    20 + apply((S == 0) * (Y == 1), 2, sum) #answer unknown and guess correctly
  Q2		<- 80 + apply(S == 0, 2, sum) #answer unknown
  guess	<- rbeta(K, Q1, Q2)
  return(guess)
}


DrawBeta_LNIRT <- function(theta, alpha, Z, mu, sigma) {
  # prior mu,sigma
  N <- nrow(Z)
  K <- ncol(Z)
  alphaT <- matrix(-alpha,
                   ncol = K,
                   nrow = N,
                   byrow = T)
  thetaT <- matrix(theta,
                   ncol = K,
                   nrow = N,
                   byrow = F)
  Z <- matrix(Z - alphaT * thetaT, ncol = K, nrow = N)
  
  XX <- alphaT
  pvar <- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  betahat <- diag(t(XX) %*% Z)
  mu <- (betahat + mu / sigma[1, 1]) / pvar
  beta <- rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(beta)
}

DrawBetaMBD_LNIRT <- function(theta, alpha, Z, mu, sigma, MBDY) {
  #missing by design
  #par1=1
  #prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z)
  thetaT	<-
    matrix(theta,
           ncol = K,
           nrow = N,
           byrow = F) * matrix(alpha,
                               ncol = K,
                               nrow = N,
                               byrow = T)
  Z 		<- matrix(Z - thetaT, ncol = K, nrow = N) * MBDY
  XX		<- matrix(-alpha,
                ncol = K,
                nrow = N,
                byrow = T) * MBDY
  pvar	<- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  betahat <- diag(t(XX) %*% Z)
  mu		<- (betahat + mu / sigma[1, 1]) / pvar
  beta	<-  rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  #beta[which(beta > 4)] <- 4 ##restrict upperbound
  #beta[which(beta < -4)] <- -4 ##restrict lowerbound
  
  return(beta)
}

DrawBetaMBD2_LNIRT <- function(theta, alpha, Z, mu, sigma, MBDY) {
  #missing by design
  #prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z)
  thetaT	<-
    matrix(theta,
           ncol = K,
           nrow = N,
           byrow = F) * matrix(alpha,
                               ncol = K,
                               nrow = N,
                               byrow = T)
  Z 		<- matrix(Z - thetaT, ncol = K, nrow = N) * MBDY
  XX		<- matrix(-1,
                ncol = K,
                nrow = N,
                byrow = T) * MBDY
  pvar	<- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  betahat <- diag(t(XX) %*% Z)
  mu		<- (betahat + mu / sigma[1, 1]) / pvar
  beta	<-  rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(beta)
}


DrawLambda_LNIRT <- function(RT, phi, zeta, sigma2, mu, sigmal) {
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- RT + t(matrix(phi, K, N)) * matrix(zeta, N, K)
  X <- matrix(1, N, 1)
  Sigma <- diag(sigma2)
  Sigma0 <- diag(K) * as.vector(sigmal)
  lambda <- SampleB(Z, X, Sigma, as.vector(mu), Sigma0)$B
  return(list(lambda = lambda))
}


DrawLambdaMBD_LNIRT <- function(RT, phi, zeta, sigma2, mu, sigmal, MBDT) {
  K <- ncol(RT)
  N <- nrow(RT)
  mu <- as.vector(mu)
  sigmal <- sigmal[1, 1]
  RTd <- (RT + t(matrix(phi, K, N)) * matrix(zeta, N, K)) * MBDT
  XX 	<- matrix(1, ncol = K, nrow = N) * MBDT
  pvar	<- diag(t(XX) %*% XX) / sigma2 + 1 / sigmal
  lambdahat <- diag(t(XX) %*% (RTd)) / sigma2
  mu		<- (lambdahat + mu / sigmal) / pvar
  lambda	<- rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(lambda)
}

DrawAlpha_LNIRT <- function(theta, beta, Z, mu, sigma) {
  # prior mu,sigma
  N <- nrow(Z)
  K <- ncol(Z)
  betaT <- matrix(beta,
                  ncol = K,
                  nrow = N,
                  byrow = T)
  XX <- matrix(c(theta - betaT), nrow = N, ncol = K)
  pvar <- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  alphahat <- diag(t(XX) %*% Z)
  mu <- (alphahat + mu / sigma[1, 1]) / pvar
  alpha <- rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(alpha)
}

DrawAlphaMBD_LNIRT <- function(theta, beta, Z, mu, sigma, MBDY) {
  #par1=1
  #prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z)
  betaT	<- matrix(beta,
                  ncol = K,
                  nrow = N,
                  byrow = T)
  XX		<-
    matrix(c(theta - betaT), nrow = N, ncol = K) * MBDY #correction missing by design
  pvar	<- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  alphahat <- diag(t(XX) %*% Z)
  mu		<- (alphahat + mu / sigma[1, 1]) / pvar
  alpha	<-  rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  ##alpha[which(alpha > 4)] <- 4 ##restrict upperbound
  
  return(alpha)
}

DrawAlphaMBD2_LNIRT <- function(theta, beta, Z, mu, sigma, MBDY) {
  #include missing by design
  #prior mu,sigma
  N		<- nrow(Z)
  K 		<- ncol(Z)
  betaT	<- matrix(beta,
                  ncol = K,
                  nrow = N,
                  byrow = T) * MBDY
  XX		<- matrix(theta, nrow = N, ncol = K) * MBDY
  pvar	<- diag(t(XX) %*% XX) + 1 / sigma[1, 1]
  alphahat <- diag(t(XX) %*% (Z + betaT))
  mu		<- (alphahat + mu / sigma[1, 1]) / pvar
  alpha	<-  rnorm(K, mean = mu, sd = sqrt(1 / pvar))
  
  return(alpha)
}


DrawPhi_LNIRT <- function(RT, lambda, zeta, sigma2, mu, sigmal) {
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- -RT + t(matrix(lambda, K, N))
  X <- matrix(zeta, N, 1)
  Sigma <- diag(sigma2)
  Sigma0 <- diag(K) * as.vector(sigmal)
  phi <- SampleB(Z, X, Sigma, as.vector(mu), Sigma0)$B
  return(phi)
}

DrawPhiMBD_LNIRT <-
  function(RT, lambda, zeta, sigma2, mu, sigmal, MBDT = MBDT) {
    K <- ncol(RT)
    N <- nrow(RT)
    mu <- as.vector(mu)
    sigmal <- sigmal[1, 1]
    RTd	<- (matrix(
      lambda,
      ncol = K,
      nrow = N,
      byrow = T
    ) - RT) * MBDT
    XX		<- matrix(zeta, nrow = N, ncol = K) * MBDT
    pvar	<- diag(t(XX) %*% XX) / sigma2 + 1 / sigmal
    phihat <- diag(t(XX) %*% (RTd)) / sigma2
    mu		<- (phihat + mu / sigmal) / pvar
    phi	<-  rnorm(K, mean = mu, sd = sqrt(1 / pvar))
    
    return(phi)
  }


SampleS2_LNIRT <- function(RT, zeta, lambda, phi) {
  K <- ncol(RT)
  N <- nrow(RT)
  ss0 <- 10
  Z <-
    RT + t(matrix(phi, K, N)) * matrix(zeta, N, K) - t(matrix(lambda, K, N))
  sigma2 <- (apply(Z * Z, 2, sum) + ss0) / rchisq(K, N)
  return(sigma2)
}

SampleS2MBD_LNIRT <- function(RT, zeta, lambda, phi, MBDT) {
  K <- ncol(RT)
  N <- nrow(RT)
  Nt <- apply(MBDT, 2, sum)
  ss0 <- 10
  Z <-
    (RT + t(matrix(phi, K, N)) * matrix(zeta, N, K) - t(matrix(lambda, K, N))) *
    MBDT
  sigma2 <- (apply(Z * Z, 2, sum) + ss0) / rchisq(K, Nt)
  return(sigma2)
}


SampleBX_LNIRT <- function(Y, XPA, XPT) {
  N <- nrow(Y)
  ka <- ncol(XPA)
  kt <- ncol(XPT)
  V0 <- diag(ka + kt) #inverse prior covariance matrix
  B0 <- matrix(0, ncol = 1, nrow = ka + kt)
  V0B0 <- matrix(V0 %*% matrix(B0, ncol = 1), ncol = 1)
  
  Bvara <- solve(crossprod(XPA) + V0[1:ka, 1:ka])
  Bvart <- solve(crossprod(XPT) + V0[(ka + 1):(ka + kt), (ka + 1):(ka +
                                                                     kt)])
  Btildea <- Bvara %*% (crossprod((XPA), Y[, 1]) + V0B0[1:ka])
  Btildet <-
    Bvart %*% (crossprod((XPT), Y[, 2]) + V0B0[(ka + 1):(ka + kt)])
  Ba <- matrix(MASS::mvrnorm(1, mu = Btildea, Sigma = Bvara), nrow = 1)
  Bt <- matrix(MASS::mvrnorm(1, mu = Btildet, Sigma = Bvart), nrow = 1)
  B <- matrix(cbind(Ba, Bt), ncol = 1)
  pred <-
    matrix(c(
      XPA %*% matrix(Ba, nrow = ka, ncol = 1),
      XPT %*% matrix(Bt, nrow = kt, ncol = 1)
    ),
    ncol = 2,
    nrow = N)
  
  return(list(B = B, pred = pred))
}
