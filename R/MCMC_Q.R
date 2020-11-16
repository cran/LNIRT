DrawZetaQ <- function(RT,
                      phi,
                      lambda,
                      sigma2,
                      mu,
                      sigmaz,
                      SigmaR,
                      X) {
  ## returns random draw from the posterior of the speed parameters
  ## where zeta[i,] ~ N(mu0[i]+mu1[i]*X+mu2[i]*X2,sigmaz)
  ## mu mean random speed components (3*1)
  ## SigmaR covar random speed components
  ## X matrix of measurement occasions
  
  K <- ncol(RT)
  N <- nrow(RT)
  Q <- 3 #intercept, slope, quadratic
  zeta <- matrix(0, ncol = Q, nrow = N)
  zetapred <- matrix(0, ncol = K, nrow = N)
  
  Z <- matrix(lambda,
              nrow = N,
              ncol = K,
              byrow = T) - RT
  muzeta <- matrix(mu %*% solve(SigmaR), nrow = N, ncol = 3)
  
  Xn <- cbind(phi, phi * X, phi * X ** 2) #equal X matrix over persons
  invsigma2K <- diag(1 / sigma2)
  invSigmax <- t(Xn) %*% invsigma2K %*% Xn
  Sigma <- solve(invSigmax + solve(SigmaR))
  zetahat <- Sigma %*% (t(Xn) %*% invsigma2K %*% t(Z) + t(muzeta))
  for (ii in 1:N) {
    zeta[ii, ] <- MASS::mvrnorm(n = 1, mu = zetahat[, ii], Sigma = Sigma)
    zetapred[ii, ] <- Xn %*% matrix(zeta[ii, ], ncol = 1)
  }
  return(list(zeta = zeta, zetapred = zetapred))
}

DrawZQ   <- function(alpha0, beta0, theta0, S, D) {
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
  ##Z		<- matrix(qnorm(tt),ncol = K, nrow = N) + matrix(eta,ncol = K, nrow = N)
  
  #deal with missings MAR#
  Z  <- matrix(0, ncol = 1, nrow = N * K)
  tt <- matrix(tt, ncol = 1, nrow = N * K)
  eta <- matrix(eta, ncol = 1, nrow = N * K)
  Z[which(D == 1)] <-  qnorm(tt)[which(D == 1)] + eta[which(D == 1)]
  Z[which(D == 0)] <-  rnorm(sum(D == 0)) + eta[which(D == 0)]
  Z <- matrix(Z, ncol = K, nrow = N)
  
  return(Z)
}

DrawThetaQ <- function(alpha0, beta0, Z, mu, sigma) {
  #prior theta MVN(muP,SigmaP)
  N			<- nrow(Z)
  K 			<- ncol(Z)
  pvar		<- (sum(alpha0 ^ 2) + 1 / as.vector(sigma))
  thetahat	<-
    (((Z + t(
      matrix(beta0, ncol = N, nrow = K)
    )) %*% matrix(
      alpha0, ncol = 1, nrow = K
    )))
  mu			<- (thetahat + as.vector(mu) / as.vector(sigma)) / pvar
  theta		<- rnorm(N, mean = mu, sd = sqrt(1 / pvar))
  #theta		<- (theta - mean(theta))#/sqrt(var(theta))
  return(theta)
}

DrawCQ <- function(S, Y) {
  # Informative beta prior: c ~ B(2,6)
  N		<- nrow(Y)
  K		<- ncol(Y)
  Q1 		<- 2 + apply((matrix(1, ncol = K, nrow = N) - S) * Y, 2, sum)
  Q2		<- 6 + apply((matrix(1, ncol = K, nrow = N) - S), 2, sum) - Q1
  guess	<- rbeta(K, Q1, Q2)
  return(guess)
}

DrawLambdaQ <- function(T, phi, zeta, sigma2, mu, sigmal) {
  ## write it as a multivariate regression problem and use SampleBQ function
  K <- ncol(T)
  N <- nrow(T)
  Z <- T + t(matrix(phi, K, N)) * matrix(zeta, N, K)
  X <-
    matrix(1, N, 1) ## will only fit an intercept term, e.g., lambda
  Sigma <- diag(sigma2) # diagonal matrix with sigma2
  Sigma0 <- diag(K) * as.vector(sigmal)
  lambda <- SampleBQ(Z, X, Sigma, as.vector(mu), Sigma0)$B
  return(list(lambda = lambda))
}

DrawPhiQ <- function(T, lambda, zeta, sigma2, mu, sigmal) {
  K <- ncol(T)
  N <- nrow(T)
  Z <- -T + t(matrix(lambda, K, N))
  X <- matrix(zeta, N, 1)
  Sigma <- diag(sigma2) # diagonal matrix with sigma2
  Sigma0 <- diag(K) * as.vector(sigmal)
  phi <- SampleBQ(Z, X, Sigma, as.vector(mu), Sigma0)$B
  return(phi)
}

SampleS2Q <- function(RT,
                      zeta = zeta,
                      X = X,
                      lambda,
                      phi) {
  ## calculates sum of squares for each item K and returns a draw from the
  ## posterior of the item-specific residual variance.
  ## prior: ss0 with degrees of freedom = 1.
  K <- ncol(RT)
  N <- nrow(RT)
  Xn <- cbind(phi, phi * X, phi * X ** 2)
  speed <- zeta %*% t(Xn)
  ss0 <- 10
  Z <- (RT + speed - t(matrix(lambda, K, N)))
  sigma2 <- (apply(Z * Z, 2, sum) + ss0) / rchisq(K, N)
  return(sigma2)
}

SampleBQ <- function(Y, X, Sigma, B0, V0) {
  ## multivariate normal regression of Y(N,k) on X(N,p)
  ## returns: regression coefficients (as vector), multivariate predictor
  ## model: Y = XB + E, E(0,Sigma)
  ## prior: vec(B) ~ N(vec(B0),V0)
  ## author: Rinke Klein Entink, adapted from Peter Rossi
  
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x% solve(crossprod(X))) + solve(V0))
  Btilde <-
    Bvar %*% (
      solve(Sigma %x% diag(m)) %*% matrix(crossprod(X, Y), ncol = 1) + solve(V0) %*%
        matrix(B0, ncol = 1)
    )
  B <- Btilde + chol(Bvar) %*% matrix(rnorm(length(Btilde)), ncol = 1)
  pred <- X %*% matrix(B, ncol = k)
  return(list(B = B, pred = pred))
}

ConditionalQ <- function(kk, Mu, Sigma, Z) {
  ## Returns the univariate conditional distribution of outcome variable Z[kk] from a K-dimensional MVN model
  ## Z[1,...,K] = mu[1,...,K] + E, E~MVN(0,Sigma)
  
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
        Mu[, kk:K] + t(matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk -
                                                                                                1), 1:(kk - 1)]) %*% t(C))
      CSigma <-
        Sigma[kk:K, kk:K] - matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk -
                                                                                                     1), 1:(kk - 1)]) %*% matrix((Sigma[1:(kk - 1), kk:K]), nrow = (kk - 1))
      J <- ncol(CSigma)
      C <- matrix(Z[1:N, (kk + 1):K] - CMu1[1:N, 2:J], ncol = (J - 1))
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
        Sigma[K, K] - t(Sigma[1:(K - 1), K]) %*% solve(Sigma[1:(K - 1), 1:(K -
                                                                             1)]) %*% Sigma[1:(K - 1), K]
    }
  }
  return(list(CMU = CMEAN, CVAR = CSD))
}

DrawLambdaPhiQ <- function(RT, zeta, X, sigma2, muI, SigmaI) {
  ## sample item parameters through multivariate matrix design ##
  ## Speed of dimension N by K
  K <- ncol(RT)
  N <- nrow(RT)
  phi <- matrix(NA, ncol = 1, nrow = K)
  lambda <- matrix(NA, ncol = 1, nrow = K)
  invSigmaI <- solve(SigmaI)
  Xn <- cbind(1, X, X ** 2)
  speed <- zeta %*% t(Xn)
  for (kk in 1:K) {
    H <- matrix(c(-speed[, kk], rep(1, N)), ncol = 2, nrow = N)
    varest <- solve(((1 / sigma2[kk]) * (t(H) %*% H)) + invSigmaI)
    meanest <-
      varest %*% ((t(H) %*% RT[, kk]) / sigma2[kk] + t(muI[1, ] %*% invSigmaI))
    lambdaphi <- MASS::mvrnorm(1, mu = t(meanest), Sigma = varest)
    phi[kk] <- lambdaphi[1]
    lambda[kk] <- lambdaphi[2]
  }
  
  return(list(phi = phi, lambda = lambda))
}

DrawRhoQ <- function(zeta, theta, muP, SigmaP) {
  ## Covariance between ability and Q-1 speed components
  ## priors : murho <- rep(0,Q-1),varrho <- diag(Q-1)*10
  
  Q <- ncol(SigmaP)
  SigmaR <- SigmaP[2:Q, 2:Q]
  N <- nrow(zeta)
  
  ##When SigmaR is a diagonal matrix, compute inverse as
  InvSigmaR <- diag(1 / diag(SigmaR))
  sigma <- SigmaP[1, 1] - SigmaP[1, 2:4] %*% (InvSigmaR) %*% SigmaP[2:4, 1]
  XX <- ((zeta - muP[, 2:Q]) %*% (InvSigmaR))
  
  ##else --> InvSigmaR <- solve(SigmaR)
  ##sigma <- SigmaP[1,1]-SigmaP[1,2:4]%*%solve(SigmaR)%*%SigmaP[2:4,1]
  ##XX <- ((zeta - muP[,2:Q])%*%solve(SigmaR))
  Sigmarho <- (t(XX) %*% XX / sigma[1] + diag(Q - 1) / 100)
  
  if (1 / det(Sigmarho) < 10e-14) {
    #only covariance between ability and intercept speed
    InvSigmaR <- (1 / SigmaP[2, 2])
    sigma <- SigmaP[1, 1] - SigmaP[1, 2] %*% (InvSigmaR) %*% SigmaP[2, 1]
    XX <- ((zeta[, 1] - muP[, 2]) * (InvSigmaR))
    Sigmarho <- 1 / (t(XX) %*% XX / sigma[1] + 1 / 100)
    if (Sigmarho < 10e-10) {
      Sigmarho <- .00001
    }
    rhohat <- Sigmarho %*% (t(XX) %*% (theta - muP[, 1]) / sigma[1])
    rho1 <- rnorm(1, mean = rhohat, sd = sqrt(Sigmarho))
    rho <- c(rho1, rep(0, Q - 2))
    SigmaP[1, 2:Q] <- SigmaP[2:Q, 1] <- rho
    
    SigmaR <- SigmaP[2:Q, 2:Q]
    InvSigmaR <- diag(1 / diag(SigmaR))
    XX <- ((zeta - muP[, 2:Q]) %*% (InvSigmaR))
    SS <- sum(((theta - muP[, 1]) - XX %*% rho) ** 2)
    sigman <- (SS + 10) / rgamma(1, shape = N / 2, rate = 1 / 2)
    SigmaP[1, 1] <-
      sigman + SigmaP[1, 2:4] %*% InvSigmaR %*% SigmaP[2:4, 1]
  } else{
    Sigmarho <- solve(Sigmarho)
    rhohat <- Sigmarho %*% (t(XX) %*% (theta - muP[, 1]) / sigma[1])
    rho <- MASS::mvrnorm(1, mu = t(rhohat), Sigma = Sigmarho)
    SigmaP[1, 2:Q] <- SigmaP[2:Q, 1] <- rho
    
    SS <- sum(((theta - muP[, 1]) - XX %*% rho) ** 2)
    sigman <- (SS + 10) / rgamma(1, shape = N / 2, rate = 1 / 2)
    SigmaP[1, 1] <-
      sigman + SigmaP[1, 2:4] %*% InvSigmaR %*% SigmaP[2:4, 1]
    ##SigmaP[1,1] <- sigman + SigmaP[1,2:4]%*%solve(SigmaR)%*%SigmaP[2:4,1]
  }
  
  return(SigmaP)
}

SimulateRTQ <- function(RT, zeta, lambda, phi, sigma2, X, DT) {
  #assume MAR
  N <- nrow(RT)
  K <- ncol(RT)
  meanT <-
    t(matrix(lambda, K, N)) - t(matrix(phi, K, N)) * (zeta %*% t(cbind(1, X, X **
                                                                         2)))
  meanT <- matrix(meanT, ncol = 1, nrow = N * K)
  sigmaL <-
    matrix(t(matrix(
      rep(sqrt(sigma2), N), nrow = K, ncol = N
    )), ncol = 1, nrow = N * K)
  RT <- matrix(RT, ncol = 1, nrow = N * K)
  RT[which(DT == 0)] <-
    rnorm(sum(DT == 0), mean = meanT[which(DT == 0)], sd = sigmaL[which(DT ==
                                                                          0)])
  RT <- matrix(RT, ncol = K, nrow = N)
  
  return(RT)
}

DrawLambdaQ1 <- function(RT, zeta, X, sigma2, muI, SigmaI) {
  ## sample time-intensity parameters through multivariate matrix design ##
  ## Speed of dimension N by K
  K <- ncol(RT)
  N <- nrow(RT)
  Xn <- cbind(1, X, X ** 2)
  speed <- zeta %*% t(Xn)
  Y <- RT + speed
  var1 <- 1 / (N / sigma2 + rep(1 / SigmaI, K))
  lambda <-
    rnorm(K,
          mean = var1 * (apply(Y, 2, sum) / sigma2 + rep(muI / SigmaI, K)),
          sd = sqrt(var1))
  
  return(list(lambda = lambda))
}

DrawLambdaPhiQ1 <- function(RT, zeta, X, sigma2, muI, SigmaI) {
  ## sample item parameters through multivariate matrix design ##
  ## Speed of dimension N by K
  K <- ncol(RT)
  N <- nrow(RT)
  phi <- matrix(NA, ncol = 1, nrow = K)
  lambda <- matrix(NA, ncol = 1, nrow = K)
  invSigmaI <- solve(SigmaI)
  Xn <- cbind(1, X, X ** 2)
  speed <- zeta %*% t(Xn)
  for (kk in 1:K) {
    H <- matrix(c(-speed[, kk], rep(1, N)), ncol = 2, nrow = N)
    varest <- solve(((1 / sigma2[kk]) * (t(H) %*% H)) + invSigmaI)
    meanest <-
      varest %*% ((t(H) %*% RT[, kk]) / sigma2[kk] + t(muI[1, ] %*% invSigmaI))
    lambdaphi <- MASS::mvrnorm(1, mu = t(meanest), Sigma = varest)
    phi[kk] <- lambdaphi[1]
    lambda[kk] <- lambdaphi[2]
  }
  
  return(list(phi = phi, lambda = lambda))
}

SampleS2Q1 <- function(RT,
                       zeta = zeta,
                       X = X,
                       lambda,
                       phi) {
  ## calculates sum of squares for each item K and returns a draw from the
  ## posterior of the item-specific residual variance.
  ## prior: ss0 with degrees of freedom = 1.
  K <- ncol(RT)
  N <- nrow(RT)
  Xn <- cbind(phi, phi * X, phi * X ** 2)
  speed <- zeta %*% t(Xn)
  ss0 <- 10
  Z <- (RT + speed - t(matrix(lambda, K, N)))
  sigma2 <- (apply(Z * Z, 2, sum) + ss0) / rchisq(K, N)
  return(sigma2)
}

DrawZetaQ1 <- function(RT,
                       phi,
                       lambda,
                       sigma2,
                       mu,
                       sigmaz,
                       SigmaR,
                       X) {
  ## returns random draw from the posterior of the speed parameters
  ## where zeta[i,] ~ N(mu0[i]+mu1[i]*X+mu2[i]*X2,sigmaz)
  ## mu mean random speed components (3*1)
  ## SigmaR covar random speed components
  ## X matrix of measurement occasions
  
  K <- ncol(RT)
  N <- nrow(RT)
  Q <- 3 #intercept, slope, quadratic
  zeta <- matrix(0, ncol = Q, nrow = N)
  zetapred <- matrix(0, ncol = K, nrow = N)
  
  Z <- matrix(lambda,
              nrow = N,
              ncol = K,
              byrow = T) - RT
  muzeta <- matrix(solve(SigmaR) %*% mu, ncol = N, nrow = 3)
  
  Xn <- cbind(phi, phi * X, phi * X ** 2) #equal X matrix over persons
  invsigma2K <- diag(1 / sigma2)
  invSigmax <- t(Xn) %*% invsigma2K %*% Xn
  Sigma <- solve(invSigmax + solve(SigmaR))
  zetahat <- Sigma %*% (t(Xn) %*% invsigma2K %*% t(Z) + muzeta)
  for (ii in 1:N) {
    zeta[ii, ] <- MASS::mvrnorm(n = 1, mu = zetahat[, ii], Sigma = Sigma)
    zetapred[ii, ] <- Xn %*% matrix(zeta[ii, ], ncol = 1)
  }
  return(list(zeta = zeta, zetapred = zetapred))
}

SampleBQ1 <- function(Y, X, Sigma, B0, V0) {
  ## multivariate normal regression of Y(N,k) on X(N,p)
  ## returns: regression coefficients (as vector), multivariate predictor
  ## model: Y = XB + E, E(0,Sigma)
  ## prior: vec(B) ~ N(vec(B0),V0)
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x% solve(crossprod(X))) + solve(V0))
  Btilde <-
    Bvar %*% (
      solve(Sigma %x% diag(m)) %*% matrix(crossprod(X, Y), ncol = 1) + solve(V0) %*%
        matrix(B0, ncol = 1)
    )
  B <- Btilde + chol(Bvar) %*% matrix(rnorm(length(Btilde)), ncol = 1)
  pred <- X %*% matrix(B, ncol = k)
  return(list(B = B, pred = pred))
}
