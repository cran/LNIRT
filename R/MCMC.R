### Common functions used in MCMC ###

Conditional <- function(kk, Mu, Sigma, Z) {

    K <- ncol(Z)
    N <- nrow(Z)
    
    if (kk == 1) {
        C <- matrix(Z[, 2:K] - Mu[, 2:K], ncol = (K - 1))
        CMEAN <- Mu[, 1] + Sigma[1, 2:K] %*% solve(Sigma[2:K, 2:K]) %*% t(C)
        CSD <- Sigma[1, 1] - Sigma[1, 2:K] %*% solve(Sigma[2:K, 2:K]) %*% Sigma[2:K, 1]
    }
    if (kk > 1) {
        if (kk < K) {
            C <- matrix(Z[1:N, 1:(kk - 1)] - Mu[, 1:(kk - 1)], ncol = (kk - 1))
            CMu1 <- Mu[, kk:K] + t(matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk - 1), 1:(kk - 1)]) %*% t(C))
            CSigma <- Sigma[kk:K, kk:K] - matrix(t(Sigma[1:(kk - 1), kk:K]), ncol = (kk - 1)) %*% solve(Sigma[1:(kk - 1), 1:(kk - 1)]) %*% matrix((Sigma[1:(kk - 
                1), kk:K]), nrow = (kk - 1))
            J <- ncol(CSigma)
            C <- matrix(Z[1:N, (kk + 1):K] - CMu1[1:N, 2:J], ncol = (J - 1))
            CMEAN <- CMu1[, 1] + CSigma[1, 2:J] %*% solve(CSigma[2:J, 2:J]) %*% t(C)
            CSD <- CSigma[1, 1] - CSigma[1, 2:J] %*% solve(CSigma[2:J, 2:J]) %*% CSigma[2:J, 1]
        }
        if (kk == K) {
            C <- matrix(Z[1:N, 1:(K - 1)] - Mu[, 1:(K - 1)], ncol = (K - 1))
            CMEAN <- Mu[, K] + t(Sigma[1:(K - 1), K]) %*% solve(Sigma[1:(K - 1), 1:(K - 1)]) %*% t(C)
            CSD <- Sigma[K, K] - t(Sigma[1:(K - 1), K]) %*% solve(Sigma[1:(K - 1), 1:(K - 1)]) %*% Sigma[1:(K - 1), K]
        }
    }
    return(list(CMU = CMEAN, CVAR = CSD))
}

SimulateRT <- function(RT, zeta, lambda, phi, sigma2, DT) {
    # assume MAR
    N <- nrow(RT)
    K <- ncol(RT)
    meanT <- t(matrix(lambda, K, N)) - t(matrix(phi, K, N)) * (zeta %*% t(rep(1, K)))
    meanT <- matrix(meanT, ncol = 1, nrow = N * K)
    sigmaL <- matrix(t(matrix(rep(sqrt(sigma2), N), nrow = K, ncol = N)), ncol = 1, nrow = N * K)
    RT <- matrix(RT, ncol = 1, nrow = N * K)
    RT[which(DT == 0)] <- rnorm(sum(DT == 0), mean = meanT[which(DT == 0)], sd = sigmaL[which(DT == 0)])
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
    
    par <- theta %*% matrix(alpha0, nrow = 1, ncol = K) - t(matrix(beta0, nrow = K, ncol = N))
    probs <- matrix(pnorm(par), ncol = K, nrow = N)
    Yn <- matrix(runif(N * K), nrow = N, ncol = K)
    Yn <- ifelse(Yn < probs, 1, 0)
    Y[(D == 0) & (G == 0)] <- Yn[(D == 0) & (G == 0)]  #missing: response generated
    
    return(Y)
}

DrawZeta <- function(RT, phi, lambda, sigma2, mu, sigmaz) {
  K <- ncol(RT)
  N <- nrow(RT)
  Z <- matrix(lambda, nrow = K, ncol = N) - t(RT)
  X <- matrix(phi, K, 1)
  sigma2inv <- diag(1/sigma2[1:K])
  vartheta <- (1/((t(phi) %*% sigma2inv) %*% phi + 1/sigmaz))[1, 1]
  meantheta <- matrix(((t(phi) %*% sigma2inv) %*% Z + t(mu/sigmaz)) * vartheta, ncol = 1, nrow = N)
  zeta <- matrix(rnorm(N, mean = meantheta, sd = sqrt(vartheta)), ncol = 1, nrow = N)
  
  return(zeta)
}

SampleB <- function(Y, X, Sigma, B0, V0) {
  m <- ncol(X)
  k <- ncol(Y)
  Bvar <- solve(solve(Sigma %x% solve(crossprod(X))) + solve(V0))
  Btilde <- Bvar %*% (solve(Sigma %x% diag(m)) %*% matrix(crossprod(X, Y), ncol = 1) + solve(V0) %*% matrix(B0, ncol = 1))
  B <- Btilde + chol(Bvar) %*% matrix(rnorm(length(Btilde)), ncol = 1)
  pred <- X %*% matrix(B, ncol = k)
  return(list(B = B, pred = pred))
}

rwishart <- function(nu, V) {
  m = nrow(V)
  df = (nu + nu - m + 1) - (nu - m + 1):nu
  if (m > 1) {
    RT = diag(sqrt(rchisq(c(rep(1, m)), df)))
    RT[lower.tri(RT)] = rnorm((m * (m + 1)/2 - m))
  } else {
    RT = sqrt(rchisq(1, df))
  }
  U = chol(V)
  C = t(RT) %*% U
  CI = backsolve(C, diag(m))
  return(list(IW = crossprod(t(CI))))
}


### MCMC functions used in LNRT ###

DrawLambdaPhi_LNRT <- function(RT, theta, sigma2, muI, sigmaI, ingroup) {
    # library(MASS)
    K <- ncol(RT)
    N <- nrow(RT)
    invSigmaI <- solve(sigmaI)
    H <- matrix(c(-theta, rep(1, N)), ncol = 2, nrow = N) * ingroup
    varest <- solve(kronecker(diag(1/sigma2[1:K]), (t(H) %*% H)) + kronecker(diag(1, K), invSigmaI))
    
    meanest <- t((t(H) %*% RT)/(t(matrix(sigma2, nrow = K, ncol = 2))) + matrix(t(muI[1, ] %*% invSigmaI), ncol = K, nrow = 2))
    meanest <- apply((matrix(rep(meanest, K), ncol = 2 * K) %*% varest) * t(kronecker(diag(K), c(1, 1))), 2, sum)

    lambdaphi <- MASS::mvrnorm(1, mu = meanest, Sigma = varest)
    lambdaphi <- matrix(lambdaphi, ncol = 2, nrow = K, byrow = RT)
    
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
    Z <- (RT + t(matrix(phi, K, N)) * matrix(zeta, N, K) - t(matrix(lambda, K, N))) * ingroup
    sigma2 <- (apply(Z * Z, 2, sum) + ss0)/rchisq(K, Nr)
    return(sigma2)
}

DrawLambda_LNRT <- function(RT, zeta, sigma2, mu, sigma) {
    
    # prior mu,sigma
    N <- nrow(RT)
    K <- ncol(RT)
    zetaT <- matrix(zeta, ncol = K, nrow = N, byrow = F)
    RTs <- matrix(RT + zetaT, ncol = K, nrow = N)
    XX <- matrix(1, ncol = K, nrow = N)
    pvar <- diag(t(XX) %*% XX)/sigma2 + 1/sigma[1, 1]
    betahat <- diag(t(XX) %*% RT)/sigma2
    mu <- (betahat + mu/sigma[1, 1])/pvar
    beta <- rnorm(K, mean = mu, sd = sqrt(1/pvar))
    
    return(beta)
}



### MCMC functions used in LNIRT ###

DrawS_LNIRT <- function(alpha0, beta0, guess0, theta0, Y) {
    N <- nrow(Y)
    K <- ncol(Y)
    eta <- t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow = N) - t(matrix(beta0, nrow = K, ncol = N))
    eta <- matrix(pnorm(eta), ncol = K, nrow = N)
    probS <- eta/(eta + t(matrix(guess0, nrow = K, ncol = N)) * (matrix(1, ncol = K, nrow = N) - eta))
    S <- matrix(runif(N * K), ncol = K, nrow = N)
    S <- matrix(ifelse(S > probS, 0, 1), ncol = K, nrow = N)
    S <- S * Y
    return(S)
}


DrawZ_LNIRT <- function(alpha0, beta0, theta0, S, D) {
    N <- nrow(S)
    K <- ncol(S)
    eta <- t(matrix(alpha0, ncol = N, nrow = K)) * matrix(theta0, ncol = K, nrow = N) - t(matrix(beta0, ncol = N, nrow = K))
    BB <- matrix(pnorm(-eta), ncol = K, nrow = N)
    BB[which(BB < 1e-05)] <- 1e-05
    BB[which(BB > (1 - 1e-05))] <- (1 - 1e-05)
    u <- matrix(runif(N * K), ncol = K, nrow = N)
    tt <- matrix((BB * (1 - S) + (1 - BB) * S) * u + BB * S, ncol = K, nrow = N)
    
    
    
    Z <- matrix(0, ncol = 1, nrow = N * K)
    tt <- matrix(tt, ncol = 1, nrow = N * K)
    eta <- matrix(eta, ncol = 1, nrow = N * K)
    Z[which(D == 1)] <- qnorm(tt)[which(D == 1)] + eta[which(D == 1)]
    Z[which(D == 0)] <- rnorm(sum(D == 0)) + eta[which(D == 0)]
    Z <- matrix(Z, ncol = K, nrow = N)
    
    return(Z)
}


DrawTheta_LNIRT <- function(alpha0, beta0, Z, mu, sigma) {

    N <- nrow(Z)
    K <- ncol(Z)
    pvar <- (sum(alpha0^2) + 1/as.vector(sigma))
    thetahat <- (((Z + t(matrix(beta0, ncol = N, nrow = K))) %*% matrix(alpha0, ncol = 1, nrow = K)))
    mu <- (thetahat + as.vector(mu)/as.vector(sigma))/pvar
    theta <- rnorm(N, mean = mu, sd = sqrt(1/pvar))
    
    return(theta)
}


DrawC_LNIRT <- function(S, Y) {
    
    N <- nrow(Y)
    K <- ncol(Y)
    Q1 <- 20 + apply((S == 0) * (Y == 1), 2, sum)  #answer unknown and guess correctly
    Q2 <- 80 + apply(S == 0, 2, sum)  #answer unknown
    guess <- rbeta(K, Q1, Q2)
    
    return(guess)
}

DrawBeta_LNIRT <- function(theta, alpha, Z, mu, sigma) {

    # prior mu,sigma
    N <- nrow(Z)
    K <- ncol(Z)
    alphaT <- matrix(-alpha, ncol = K, nrow = N, byrow = T)
    thetaT <- matrix(theta, ncol = K, nrow = N, byrow = F)
    Z <- matrix(Z - alphaT * thetaT, ncol = K, nrow = N)
    
    XX <- alphaT
    pvar <- diag(t(XX) %*% XX) + 1/sigma[1, 1]
    betahat <- diag(t(XX) %*% Z)
    mu <- (betahat + mu/sigma[1, 1])/pvar
    beta <- rnorm(K, mean = mu, sd = sqrt(1/pvar))

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

DrawAlpha_LNIRT <- function(theta, beta, Z, mu, sigma) {
    
    # prior mu,sigma
    N <- nrow(Z)
    K <- ncol(Z)
    betaT <- matrix(beta, ncol = K, nrow = N, byrow = T)
    XX <- matrix(c(theta - betaT), nrow = N, ncol = K)
    pvar <- diag(t(XX) %*% XX) + 1/sigma[1, 1]
    alphahat <- diag(t(XX) %*% Z)
    mu <- (alphahat + mu/sigma[1, 1])/pvar
    alpha <- rnorm(K, mean = mu, sd = sqrt(1/pvar))
    
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


SampleS2_LNIRT <- function(RT, zeta, lambda, phi) {
    K <- ncol(RT)
    N <- nrow(RT)
    ss0 <- 10
    Z <- RT + t(matrix(phi, K, N)) * matrix(zeta, N, K) - t(matrix(lambda, K, N))
    sigma2 <- (apply(Z * Z, 2, sum) + ss0)/rchisq(K, N)
    return(sigma2)
}




