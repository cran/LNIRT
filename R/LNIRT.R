#' Log-normal response time IRT modelling
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats ks.test pchisq pgamma pnorm qnorm rbeta
#' rbinom rchisq rgamma rlnorm rnorm runif var
#' @importFrom utils flush.console
#' 
#' @param RT
#' a Person-x-Item matrix of log-response times (time spent on solving an item).
#' @param Y
#' a Person-x-Item matrix of responses.
#' @param data
#' either a list or a simLNIRT object containing the response time and response matrices. 
#' If a simLNIRT object is provided, in the summary the simulated item and time parameters are shown alongside of the estimates.
#' If the required variables cannot be found in the list, or if no data object is given, then the variables are taken
#' from the environment from which LNIRT is called.
#' @param XG
#' the number of MCMC iterations to perform (default: 1000).
#' @param guess 
#' include guessing parameters in the IRT model (default: false).
#' @param par1
#' use alternative parameterization (default: false).
#' @param residual
#' compute residuals (default: false).
#' @param WL
#' define the time-discrimination parameter as measurement error variance parameter (default: false).
#' @param td
#' set time-discrimination to one (default: false).
#' @param alpha
#' an optional vector of pre-defined item-discrimination parameters.
#' @param beta
#' an optional vector of pre-defined item-difficulty parameters.
#' 
#' @return 
#' an object of class LNIRT.
#' 
#' @examples 
#' # Log-normal response time IRT modelling
#' data <- simLNIRT(N = 500, K = 20, rho = 0.8, WL = FALSE)
#' summary(LNIRT(RT = RT1, Y = Y, data = data, XG = 500, WL = FALSE))
#'  
#' @export
LNIRT <- function(RT, Y, data, XG = 1000, guess = FALSE, par1 = FALSE, residual = FALSE, WL = FALSE, td = FALSE, alpha, beta) {

    ## ident = 1: Identification : fix mean item difficulty(intensity) and product item (time) discrimination responses and response times ident =
    ## 2: Identification : fix mean ability and speed and product item discrimination responses and response times ident <- 1 (default)
    ident <- 2  #(to investigate person fit using latent scores)
    
    if (!missing(data)) {
      # Try to find RT and Y in the data set first
      tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
      tryCatch(Y <- eval(substitute(Y), data), error=function(e) NULL)
    } else {
      data = NULL
    }
    
    
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
    
    N <- nrow(Y)
    K <- ncol(Y)
    
    ## population theta (ability - speed)
    theta <- matrix(rnorm(N * 2), ncol = 2) # 1: ability, 2: speed
    muP <- matrix(0, N, 2) # Mean estimates for person parameters 
    SigmaP <- diag(2) # Person covariance matrix
    muP0 <- matrix(0, 1, 2)
    SigmaP0 <- diag(2)
    
    ## population item (ability - speed)
    ab <- matrix(rnorm(K * 4), ncol = 4) # 1: item discrimination 2: item difficulty 3: time discriminiation 4: time intensity
    ab[, c(1, 3)] <- 1
    muI <- t(matrix(rep(c(1, 0), 2 * K), ncol = K)) # Mean estimates for item parameters 
    muI0 <- muI[1, ]
    SigmaI <- diag(4) # Item covariance matrix
    SigmaI0 <- diag(4)/10
    for (ii in 1:4) {
        SigmaI0[ii, ] <- SigmaI0[ii, ] * rep(c(0.5, 3), 2)
    }
    ifelse(PNO, guess0 <- rep(0.2, K), guess0 <- rep(0, K))
    
    ## storage
    MT <- MT2 <- array(0, dim = c(N, 2))
    MAB <- array(0, dim = c(XG, K, 4))
    Mguess <- matrix(0, XG, K)
    MmuP <- array(0, dim = c(XG, 2))
    MmuI <- array(0, dim = c(XG, 4))
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
    
    if (XG > 1000) {
        MPF <- matrix(0, ncol = N, nrow = XG - 1000)
        MPFb <- matrix(0, ncol = N, nrow = XG - 1000)
        MPFp <- matrix(0, ncol = N, nrow = XG - 1000)
        MPFbp <- matrix(0, ncol = N, nrow = XG - 1000)
        
        MPFTb <- matrix(0, ncol = N, nrow = XG - 1000)
        MPFTbp <- matrix(0, ncol = N, nrow = XG - 1000)
    }
    EAPbeta <- rep(0, K)
    EAPalpha <- rep(0, K)
    EAPmub <- 0
    EAPsigmab <- 1
    EAPSigmaP <- 1
    EAPmuI <- 0
    EAPSigmaI <- 1
    iis <- 1
    
    # Indicate which responses are missing
    D <- matrix(1, ncol = 1, nrow = N * K)
    D[which(is.na(Y))] <- 0
    D <- matrix(D, nrow = N, ncol = K)
    
    # Indicate which RT's are missing
    DT <- matrix(1, ncol = 1, nrow = N * K)
    DT[which(is.na(RT))] <- 0
    DT <- matrix(DT, nrow = N, ncol = K)
    
    EAPCP1 <- matrix(0, ncol = 1, nrow = N)
    EAPCP2 <- matrix(0, ncol = 1, nrow = N)
    EAPCP3 <- matrix(0, ncol = 1, nrow = N)
    
    ## Start MCMC algorithm
    
    for (ii in 1:XG) { # For each iteration
        if (sum(DT == 0) > 0) { # Simulate missing RT's 
            if (WL) {
                RT <- SimulateRT(RT = RT, zeta = theta[, 2], lambda = ab[, 4], phi = rep(1, K), sigma2 = sigma2, DT = DT)
            } else {
                RT <- SimulateRT(RT = RT, zeta = theta[, 2], lambda = ab[, 4], phi = ab[, 3], sigma2 = sigma2, DT = DT)
            }
        }
        # ability test
        if (PNO) { # If guessing
            if (sum(D == 0) > 0) { # Simulate missing responses 
                Y <- SimulateY(Y = Y, theta = theta[, 1], alpha0 = ab[, 1], beta0 = ab[, 2], guess0 = guess0, D = D)
            }
            if (par1) {
                SR <- DrawS_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2] * ab[, 1], guess0 = guess0, theta0 = theta[, 1], Y = Y)
                ZR <- DrawZ_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2] * ab[, 1], theta0 = theta[, 1], S = SR, D = D)
            } else {
                SR <- DrawS_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2], guess0 = guess0, theta0 = theta[, 1], Y = Y)
                ZR <- DrawZ_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2], theta0 = theta[, 1], S = SR, D = D)
            }
        } else {
            if (par1) {
                ZR <- DrawZ_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2] * ab[, 1], theta0 = theta[, 1], S = Y, D = D)
            } else {
                ZR <- DrawZ_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2], theta0 = theta[, 1], S = Y, D = D)
            }
        }
        
        # Draw ability parameter 
        dum <- Conditional(1, muP, SigmaP, theta)
        if (par1) {
            theta[, 1] <- DrawTheta_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2] * ab[, 1], Z = ZR, mu = dum$CMU, sigma = dum$CVAR)
        } else {
            theta[, 1] <- DrawTheta_LNIRT(alpha0 = ab[, 1], beta0 = ab[, 2], Z = ZR, mu = dum$CMU, sigma = dum$CVAR)
        }
        if (ident == 2) {
            # rescale for identification
            theta[, 1] <- theta[, 1] - mean(theta[, 1])
        }
        
        # Draw speed parameter 
        dum <- Conditional(2, muP, SigmaP, theta)
        if (par1) {
            theta[, 2] <- DrawZeta(RT = RT, phi = ab[, 3], lambda = ab[, 3] * ab[, 4], sigma2 = sigma2, mu = as.vector(dum$CMU), sigmaz = as.vector(dum$CVAR))  ## speed 
        } else {
            theta[, 2] <- DrawZeta(RT = RT, phi = ab[, 3], lambda = ab[, 4], sigma2 = sigma2, mu = as.vector(dum$CMU), sigmaz = as.vector(dum$CVAR))  ## speed 
        }
        
        # Rescale for identification
        if (ident == 2) {
            theta[, 2] <- theta[, 2] - mean(theta[, 2])
        }
        
        MT[1:N, 1:2] <- MT[1:N, 1:2] + theta # for theta
        MT2[1:N, 1:2] <- MT2[1:N, 1:2] + theta^2 # for sd(theta) 
        
        # Draw guessing parameter from beta distribution
        if (PNO) {
            guess0 <- DrawC_LNIRT(S = SR, Y = Y)
        }
        Mguess[ii, ] <- guess0
        
        # ab - 1: item discrimination 2: item difficulty 3: time discriminiation 4: time intensity
        if (WL) {
            ab1 <- cbind(ab[, 1], ab[, 2], 1/sqrt(sigma2), ab[, 4])
        } else {
            ab1 <- ab
        }
        
        # Draw item discrimination parameter
        if (discr) {
            dum <- Conditional(kk = 1, Mu = muI, Sigma = SigmaI, Z = ab1)  #discrimination
            if (par1 == 1) {
                ab[, 1] <- abs(DrawAlpha_LNIRT(theta = theta[, 1], beta = ab[, 2], Z = ZR, mu = dum$CMU, sigma = dum$CVAR))
            } else {
                ab[, 1] <- abs(DrawPhi_LNIRT(RT = ZR, lambda = -ab[, 2], zeta = -theta[, 1], sigma2 = rep(1, K), mu = dum$CMU, sigmal = dum$CVAR))
            }
        } else {
            ab[, 1] <- alpha
        }
        ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
        
        if (WL) {
            ab1 <- cbind(ab[, 1], ab[, 2], 1/sqrt(sigma2), ab[, 4])
        } else {
            ab1 <- ab
        }
        
        # Draw item difficulty parameter
        if (diffc) {
            dum <- Conditional(kk = 2, Mu = muI, Sigma = SigmaI, Z = ab1)  #difficulty
            if (par1 == 1) {
                ab[, 2] <- DrawBeta_LNIRT(theta = theta[, 1], alpha = ab[, 1], Z = ZR, mu = dum$CMU, sigma = dum$CVAR)
            } else {
                ab[, 2] <- -DrawLambda_LNIRT(ZR, -ab[, 1], theta[, 1], rep(1, K), dum$CMU, dum$CVAR)$lambda
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
        if ((WL) | (td)) {
            # no time discrimination, 1/(sqrt(error variance)) = discrimination on MVN prior
            ab[, 3] <- 1
        } else {
            dum <- Conditional(3, muI, SigmaI, ab)  #time discrimination
            # if(par1==1){ ab[,3] <- abs(DrawPhi_LNIRT(RT,ab[,3]*ab[,4],theta[,2],sigma2,dum$CMU,dum$CVAR)) }else{
            ab[, 3] <- abs(DrawPhi_LNIRT(RT, ab[, 4], theta[, 2], sigma2, dum$CMU, dum$CVAR))
            # }
            ab[, 3] <- ab[, 3]/(prod(ab[, 3])^(1/K))
        }
        
        if (WL) {
            ab1 <- cbind(ab[, 1], ab[, 2], 1/sqrt(sigma2), ab[, 4])
        } else {
            ab1 <- ab
        }
        
        # Draw time intensity parameter
        dum <- Conditional(kk = 4, Mu = muI, Sigma = SigmaI, Z = ab1)  #time intensity
        ab[, 4] <- DrawLambda_LNIRT(RT, ab[, 3], theta[, 2], sigma2, dum$CMU, dum$CVAR)$lambda
        # if(par1==1){ ab[,4] <- ab[,4]/ab[,3] } rescale for identification
        if (ident == 1) {
            ab[, 4] <- ab[, 4] - mean(ab[, 4])
        }
        
        MAB[ii, 1:K, 1:4] <- ab
        
        # WL: Measurement error variance for time discrimination
        # Else: Measurement error variance for time term 
        sigma2 <- SampleS2_LNIRT(RT = RT, zeta = theta[, 2], lambda = ab[, 4], phi = ab[, 3])
        if (WL) {
            Msigma2[ii, 1:K] <- 1/sqrt(sigma2)
        } else {
            Msigma2[ii, 1:K] <- sigma2
        }
        
        # Population mean estimate for person ability and speed
        X <- matrix(1, N, 1)
        muP <- SampleB(theta, X, SigmaP, muP0, SigmaP0)
        MmuP[ii, ] <- muP$B
        muP <- muP$pred
        
        # Covariance matrix person parameters
        SS <- crossprod(theta - muP) + SigmaP0
        SigmaP <- rwishart(2 + N, chol2inv(chol(SS)))$IW
        MSP[ii, , ] <- SigmaP
        
        X <- matrix(1, K, 1)
        if (WL) {
            ab1 <- cbind(ab[, 1], ab[, 2], 1/sqrt(sigma2), ab[, 4])
        } else {
            ab1 <- ab
        }
        
        # Population mean estimates for item parameters 
        # 1: item discrimination 2: item difficulty 3: time discrimination 4: time intensity
        muI2 <- SampleB(Y = ab1, X = X, Sigma = SigmaI, B0 = muI0, V0 = SigmaI0)
        if (ident == 2) {
            MmuI[ii, c(1, 2, 3, 4)] <- muI2$B[c(1, 2, 3, 4)]
            muI[, c(1, 2, 3, 4)] <- muI2$pred[, c(1, 2, 3, 4)]
        } else {
            MmuI[ii, c(1, 3)] <- muI2$B[c(1, 3)]
            muI[, c(1, 3)] <- muI2$pred[, c(1, 3)]
        }
        
        # Covariance matrix item parameters
        muI1 <- matrix(muI, ncol = 4, nrow = K, byrow = FALSE)
        SS <- crossprod(ab1 - muI1) + SigmaI0
        SigmaI <- rwishart(4 + K, chol2inv(chol(SS)))$IW
        MSI[ii, , ] <- SigmaI
        
        if (ii > 1000) {
            EAPmuI <- (muI[1, 4] + (iis - 1) * EAPmuI)/iis
            EAPsigma2 <- (sigma2 + (iis - 1) * EAPsigma2)/iis
            EAPlambda <- (ab[, 4] + (iis - 1) * EAPlambda)/iis
            EAPphi <- (ab[, 3] + (iis - 1) * EAPphi)/iis
            EAPSigmaI <- (SigmaI[4, 4] + (iis - 1) * EAPSigmaI)/iis
            EAPtheta <- (theta + (iis - 1) * EAPtheta)/iis
            
            EAPmub <- (muI[1, 2] + (iis - 1) * EAPmub)/iis
            EAPsigmab <- (SigmaI[2, 2] + (iis - 1) * EAPsigmab)/iis
            if (par1) {
                EAPbeta <- (ab[, 1] * ab[, 2] + (iis - 1) * EAPbeta)/iis
            } else {
                EAPbeta <- (ab[, 2] + (iis - 1) * EAPbeta)/iis
            }
            EAPalpha <- (ab[, 1] + (iis - 1) * EAPalpha)/iis
            EAPSigmaP <- (SigmaP[1, 1] + (iis - 1) * EAPSigmaP)/iis
            
            
            ############################## 
            
            if (residual) {
                
                ## IRT Fit Evaluation
                if (par1) {
                  beta1 <- ab[, 1] * ab[, 2]
                } else {
                  beta1 <- ab[, 2]
                }
                if (PNO) {
                  dum <- residualA(Z = ZR, Y = SR, theta = theta[, 1], alpha = ab[, 1], beta = beta1, EAPtheta = EAPtheta[, 1], EAPalpha = EAPalpha, 
                    EAPbeta = EAPbeta)
                } else {
                  dum <- residualA(Z = ZR, Y = Y, theta = theta[, 1], alpha = ab[, 1], beta = beta1, EAPtheta = EAPtheta[, 1], EAPalpha = EAPalpha, 
                    EAPbeta = EAPbeta)
                }
                EAPKSA <- (dum$KS[1:K, 1] + (iis - 1) * EAPKSA)/iis
                EAPresidA <- (dum$presidA + (iis - 1) * EAPresidA)/iis
                
                # lZPAT <- lZPAT + dum$lZPAT
                lZPA <- lZPA + dum$lZPA
                # lZIA <- lZIA + dum$lZIA
                CF2 <- ifelse(dum$PFlp < 0.05, 1, 0)  #significance level = .05
                EAPCP2 <- (CF2 + (iis - 1) * EAPCP2)/iis
                
                
                EAPl0 <- ((iis - 1) * EAPl0 + dum$l0)/iis
                PFl <- ((iis - 1) * PFl + dum$PFl)/iis
                IFl <- ((iis - 1) * IFl + dum$IFl)/iis
                PFlp <- ((iis - 1) * PFlp + dum$PFlp)/iis
                IFlp <- ((iis - 1) * IFlp + dum$IFlp)/iis
                
                ############################## 
                
                ## Log-Normal Fit Evaluation
                if (par1) {
                  lambda1 <- ab[, 3] * ab[, 4]
                } else {
                  lambda1 <- ab[, 4]
                }
                dum <- personfitLN(RT = RT, theta = theta[, 2], phi = ab[, 3], lambda = lambda1, sigma2 = sigma2)  # lZ statistic
                lZP <- lZP + dum$lZP
                lZPT <- lZPT + dum$lZPT
                CF1 <- ifelse(dum$lZP < 0.05, 1, 0)  #significance level = .05
                EAPCP1 <- (CF1 + (iis - 1) * EAPCP1)/iis  #speed
                EAPCP3 <- (CF1 * CF2 + (iis - 1) * EAPCP3)/iis
                
                dum <- itemfitLN(RT = RT, theta = theta[, 2], phi = ab[, 3], lambda = lambda1, sigma2 = sigma2)
                lZI <- lZI + dum$lZI
                
                dum <- residualLN(RT = RT, theta = theta[, 2], phi = ab[, 3], lambda = lambda1, sigma2 = sigma2, EAPtheta = EAPtheta[, 2], EAPlambda = EAPlambda, 
                  EAPphi = EAPphi, EAPsigma2 = EAPsigma2)
                EAPresid <- EAPresid + dum$presid
                EAPKS <- (dum$KS[1:K, 1] + (iis - 1) * EAPKS)/iis
                
                iis <- iis + 1
            }
            
        }
        
        if (ii%%100 == 0) 
            cat("Iteration ", ii, " ", "\n")
        flush.console()
        
    }
    
    
    MT <- MT/XG # Posterior mean theta
    MT2 <- sqrt(MT2/XG - MT^2) # Posterior mean sd(theta)
    
    if (ii > 1000) {
        if (residual) {
            lZP <- lZP/(XG - 1000)
            lZPT <- lZPT/(XG - 1000)
            lZI <- lZI/(XG - 1000)
            EAPresid <- EAPresid/(XG - 1000)
            lZPAT <- lZPAT/(XG - 1000)
            lZPA <- lZPA/(XG - 1000)
            lZIA <- lZIA/(XG - 1000)
        }
    }

    if (XG > 1000) {
        if (residual) {
            out <- list(Mtheta = MT, MTSD = MT2, MAB = MAB, MmuP = MmuP, MSP = MSP, MmuI = MmuI, MSI = MSI, Mguess = Mguess, Msigma2 = Msigma2, 
                lZP = lZP, lZPT = lZPT, lZPA = lZPA, lZI = lZI, EAPresid = EAPresid, EAPresidA = EAPresidA, EAPKS = EAPKS, EAPKSA = EAPKSA, PFl = PFl, 
                PFlp = PFlp, IFl = IFl, IFlp = IFlp, EAPl0 = EAPl0, RT = RT, Y = Y, EAPCP1 = EAPCP1, EAPCP2 = EAPCP2, EAPCP3 = EAPCP3, WL = WL, td = td, guess = guess, par1 = par1, data = data)
        } else {
            out <- list(Mtheta = MT, MTSD = MT2, MAB = MAB, MmuP = MmuP, MSP = MSP, MmuI = MmuI, MSI = MSI, Mguess = Mguess, Msigma2 = Msigma2, 
                RT = RT, Y = Y, WL = WL, td = td, guess = guess, par1 = par1, data = data)
        }
    } else {
        out <- list(Mtheta = MT, MTSD = MT2, MAB = MAB, MmuP = MmuP, MSP = MSP, MmuI = MmuI, MSI = MSI, Mguess = Mguess, Msigma2 = Msigma2, RT = RT, 
            Y = Y, WL = WL, td = td, guess = guess, par1 = par1, data = data)
    }
    
    class(out) <- "LNIRT"
    return(out)
    
}
