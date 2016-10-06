#' Simulate data for log-normal response time IRT modelling
#' 
#' @param N
#' the number of persons.
#' @param K
#' the number of items.
#' @param rho
#' the correlation between the person ability and person speed parameter.
#' @param td
#' set time-discrimination to one (default: false).
#' @param WL
#' define the time-discrimination parameter as measurement error variance parameter (default: false).
#' 
#' @return 
#' an object of class simLNIRT.
#' 
#' @export
simLNIRT <- function(N, K, rho, td, WL) {
    if (missing(td)) {
        td <- 0
    } else {
        td <- 1  #time discrimination <- 1
    }
    
    if (missing(WL)) {
        WL <- 0
    } else {
        WL <- 1  #time discrimination = 1/sqrt(sigma2)
    }
    
    mutheta <- rep(0, 2)
    covtheta <- diag(2)
    covtheta[1, 2] <- covtheta[2, 1] <- rho
    theta <- MASS::mvrnorm(N, mutheta, covtheta, empirical = TRUE)
    
    covitem <- diag(4)
    for (ii in 1:4) {
        covitem[ii, ] <- covitem[ii, ] * rep(c(0.05, 1), 2)
    }
    muitem <- rep(c(1, 0), 2)
    
    sigma2 <- rlnorm(K, meanlog = 0, sdlog = 0.3)
    
    if (td == 1 & WL == 0) {
        ab <- MASS::mvrnorm(K, muitem, covitem)
        ab[, 3] <- rep(1, K)
        abn <- MASS::mvrnorm(K, muitem[-3], covitem[-3, -3])
        ab[, c(2, 4)] <- abn[, c(2, 3)] - t(matrix(colMeans(abn[, c(2, 3)]), 2, K))
        ab[, 1] <- abs(ab[, 1])
        ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
    }
    if (WL == 1) {
        muitem <- c(1, 0, 1, 0)
        ab <- MASS::mvrnorm(K, muitem, covitem)
        sigma2 <- 1/ab[, 3]^2
        ab[, 1] <- abs(ab[, 1])
        ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
    }
    
    if (td == 0 & WL == 0) {
        ab <- MASS::mvrnorm(K, muitem, covitem)
        ab[, c(2, 4)] <- ab[, c(2, 4)] - t(matrix(colMeans(ab[, c(2, 4)]), 2, K))
        ab[, 1] <- abs(ab[, 1])
        ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
        ab[, 3] <- abs(ab[, 3])
        ab[, 3] <- ab[, 3]/(prod(ab[, 3])^(1/K))
    }
    
    
    par <- theta[, 1] %*% matrix(ab[, 1], nrow = 1, ncol = K) - t(matrix(ab[, 2], nrow = K, ncol = N))
    probs <- matrix(pnorm(par), ncol = K, nrow = N)
    Y <- matrix(runif(N * K), nrow = N, ncol = K)
    Y <- ifelse(Y < probs, 1, 0)
    
    # quess <- rbeta(K,20,90)
    quess <- rep(0.1, K)
    S <- matrix(0, ncol = K, nrow = N)  #S=1:guess response
    Y1g <- Yg <- matrix(0, ncol = K, nrow = N)
    for (kk in 1:K) {
        S[, kk] <- rbinom(N, 1, quess[kk])
    }
    Yg[S == 1] <- 1
    Yg[S == 0 & Y == 1] <- 1
    
    # parameterisatie a*(theta-b) # ab[,2] = ab[,2]/ab[,1]
    par <- matrix(ab[, 1], ncol = K, nrow = N, byrow = T) * (matrix(theta[, 1], ncol = K, nrow = N) - matrix(ab[, 2], ncol = K, nrow = N, byrow = T))
    probs <- matrix(pnorm(par), ncol = K, nrow = N)
    Y1 <- matrix(runif(N * K), nrow = N, ncol = K)
    Y1 <- ifelse(Y1 < probs, 1, 0)
    Y1g[S == 1] <- 1
    Y1g[S == 0 & Y1 == 1] <- 1
    
    # response time
    RT <- RT1 <- matrix(0, ncol = K, nrow = N)
    for (kk in 1:K) {
        time <- matrix(rnorm(N, sd = sqrt(sigma2[kk])), nrow = N, ncol = 1)
        RT1[1:N, kk] <- (ab[kk, 4] - theta[, 2]) + time[1:N]
        RT[1:N, kk] <- ab[kk, 4] - ab[kk, 3] * theta[, 2] + time[1:N]
    }
    
    out <- list(Y = Y, Yg = Yg, Y1 = Y1, Y1g = Y1g, RT = RT, RT1 = RT1, theta = theta, ab = ab, sigma2 = sigma2, quess = quess)
    class(out) <- "simLNIRT"
    return(out)
}

