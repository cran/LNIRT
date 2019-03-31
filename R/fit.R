### Person fit ###

personfit <- function(RT, theta, phi, lambda, sigma2, nug) {
    
    if (missing(nug)) {
        nugpf <- FALSE
    } else {
        nugpf <- TRUE
    }
    
    if (missing(sigma2)) {
        sigm <- FALSE
    } else {
        sigm <- TRUE
    }
    
    if (nugpf) {
        
        return(personfitLG(RT, theta, phi, lambda, nug))
        
    }
    
    if (sigm) {
        return(personfitLN(RT, theta, phi, lambda, sigma2))
    }
    
}

personfitLN <- function(RT, theta, phi, lambda, sigma2) {
    
    K <- ncol(RT)
    N <- nrow(RT)
    
    diff <- matrix(lambda, nrow = N, ncol = K, byrow = T) - matrix(phi, ncol = K, nrow = N, byrow = T) * matrix(theta, ncol = K, nrow = N)
    diff <- (RT - diff)^2
    lZd <- diff * matrix(1/sigma2, ncol = K, nrow = N, byrow = TRUE)
    lZ <- (apply(lZd, 1, sum) - K)/sqrt(2 * K)
    
    dum <- (apply(lZd, 1, mean)^(1/3) - (1 - 2/(9 * K)))/sqrt(2/(9 * K))
    lZP1 <- 1 - pnorm(dum)
    lZP2 <- 1 - pnorm(lZ)
    lZP3 <- 1 - pchisq(apply(lZd, 1, sum), df = K)
    lZPT <- apply(lZd, 1, sum)
    return(list(lZPT = lZPT, lZP = lZP3))
    
    
}

personfitLNMBD <- function(RT,theta,phi,lambda,sigma2,MBDT){
  #correct missing by design
  K <- ncol(RT)
  N <- nrow(RT)
  
  RTd <- RT*MBDT	
  diff <- (matrix(lambda,nrow=N,ncol=K,byrow=T) - 
    matrix(phi,ncol=K,nrow=N,byrow=T)*matrix(theta,ncol=K,nrow=N))*MBDT
  diff <- (RTd - diff)**2

  Kt <- apply(MBDT,1,sum)	
  lZd <- diff*matrix(1/sigma2,ncol=K,nrow=N,byrow=TRUE)*MBDT
  lZ <- (apply(lZd,1,sum) - Kt)/sqrt(2*Kt)  
  dum <- (apply(lZd,1,mean)^(1/3) - (1-2/(9*Kt)))/sqrt(2/(9*Kt))

  lZP1 <- 1 - pnorm(dum)    
  lZP2 <- 1 - pnorm(lZ)     

  lZP3 <- 1 - pchisq(apply(lZd,1,sum),df=Kt) 
  lZPT <- apply(lZd,1,sum)

  return(list(lZPT=lZPT,lZP=lZP3))
}   


personfitLG <- function(RT, theta, phi, lambda, nug) {
    
    K <- ncol(RT)
    N <- nrow(RT)
    logspeed <- t(matrix(lambda, nrow = K, ncol = N)) - matrix(theta, ncol = K, nrow = N)
    logspeed1 <- matrix(t(logspeed), ncol = 1, nrow = N * K)
    lambdat1 <- matrix(t(exp(logspeed)), ncol = 1, nrow = N * K)
    nug1 <- matrix((matrix(nug, nrow = K, ncol = N)), ncol = 1, nrow = N * K)
    
    ord <- matrix(t(matrix(rep(1:N, K), ncol = K)), ncol = 1)
    RT1 <- matrix(t(matrix(RT, ncol = K, nrow = N)), ncol = 1, nrow = N * K)
    
    
    lGP1 <- (nug1 - 1) * (log(nug1) - (logspeed1) + (log(RT1) - digamma(nug1))) - nug1 * (RT1/lambdat1 - 1)
    lGP1 <- tapply(lGP1^2, ord, sum)
    
    
    RTn <- rgamma(N * K, shape = nug1, rate = nug1/lambdat1)
    lGPn <- (nug1 - 1) * (log(nug1) - (logspeed1) + (log(RTn) - digamma(nug1))) - nug1 * (RTn/lambdat1 - 1)
    lGPn <- tapply(lGPn^2, ord, sum)
    
    lGP <- rep(0, N)
    lGP2 <- rep(0, N)
    lGP[which(lGPn > lGP1)] <- 1
    
    
    # gamma parameters for each person
    lambdat1t <- tapply(lambdat1, ord, sum)
    lambdat1n <- tapply(lambdat1^2/nug1, ord, sum)
    nugS <- lambdat1t^2/lambdat1n
    lambdat1S <- lambdat1n/lambdat1t
    RT1S <- tapply(RT1, ord, sum)
    lGP2 <- pgamma(RT1S, shape = nugS, scale = lambdat1S)
    lGP2 <- pgamma(RT1S/lambdat1S, shape = nugS, scale = 1)
    
    return(list(lGP = lGP, lGP2 = lGP2))
}



### Item fit ###

itemfitLN <- function(RT, theta, phi, lambda, sigma2) {
    
    K <- ncol(RT)
    N <- nrow(RT)
    
    diff <- matrix(lambda, nrow = N, ncol = K, byrow = T) - matrix(phi, ncol = K, nrow = N, byrow = T) * matrix(theta, ncol = K, nrow = N)
    diff <- (RT - diff)^2
    lZd <- diff * matrix(1/sigma2, ncol = K, nrow = N, byrow = TRUE)
    lI <- (apply(lZd, 2, sum) - N)/sqrt(2 * N)
    lZI <- 1 - pnorm(lI)
    
    lZI <- 1 - pchisq(apply(lZd, 2, sum), df = N)
    
    return(list(lZI = lZI))
}

itemfitLNMBD <- function(RT,theta,phi,lambda,sigma2,MBDT=MBDT){
  
  K <- ncol(RT)
  N <- nrow(RT)
  Nt <- apply(MBDT,2,sum)
  RTd <- RT*MBDT	
  diff <- (matrix(lambda,nrow=N,ncol=K,byrow=T) - 
    matrix(phi,ncol=K,nrow=N,byrow=T)*matrix(theta,ncol=K,nrow=N))*MBDT
  diff <- (RTd - diff)**2
  lZd <- diff*matrix(1/sigma2,ncol=K,nrow=N,byrow=TRUE)*MBDT
  lI <- (apply(lZd,2,sum) - Nt)/sqrt(2*Nt)  
  lZI <- 1 - pnorm(lI)     

  lZI <- 1 - pchisq(apply(lZd,2,sum),df=Nt) 

  return(list(lZI=lZI))
}   


## Residuals ###

residualLN <- function(RT, theta, phi, lambda, sigma2, EAPtheta, EAPlambda, EAPphi, EAPsigma2) {
    
    K <- ncol(RT)
    N <- nrow(RT)
    KS <- matrix(0, ncol = 1, nrow = K)
    
    # compute fitted probabilities (approximately uniformly distributed)
    
    muik <- matrix(lambda, nrow = N, ncol = K, byrow = TRUE) - matrix(phi, ncol = K, nrow = N, byrow = TRUE) * matrix(theta, ncol = K, nrow = N)
    
    # Compute Extremeness Residuals (posterior probability greater than 2)
    
    diff <- (RT - muik) * matrix(sqrt(1/sigma2), ncol = K, nrow = N, byrow = TRUE)
    presid <- (1 - pnorm(2, mean = diff, sd = 1)) + pnorm(-2, mean = diff, sd = 1)
    
    muik <- matrix(EAPlambda, nrow = N, ncol = K, byrow = TRUE) - matrix(EAPphi, ncol = K, nrow = N, byrow = TRUE) * matrix(EAPtheta, ncol = K, 
        nrow = N)
    muiklong <- matrix(muik, ncol = 1, nrow = N * K)
    RTlong <- matrix(RT, ncol = 1, nrow = N * K)
    sigma2long <- matrix(matrix(sqrt(EAPsigma2), ncol = K, nrow = N, byrow = TRUE), ncol = 1, nrow = N * K)
    errorfit <- (RTlong - muiklong)/sigma2long
    errorfit <- matrix(rnorm(N * K) * 1e-06 + errorfit, ncol = K, nrow = N)  #to remove ties
    
    # Perform one-sample Kolmogorov Smirnov Test across Items
    
    for (kk in 1:K) {
        KS[kk, 1] <- ks.test(errorfit[, kk], y = pnorm)$p.value
    }
    
    return(list(KS = KS, presid = presid))
    
}

residualLNMBD <- function(RT,theta,phi,lambda,sigma2,EAPtheta,EAPlambda,EAPphi,EAPsigma2,MBDT=MBDT){

K <- ncol(RT)
N <- nrow(RT)
KS <- matrix(0,ncol=1,nrow=K)

#compute fitted probabilities (approximately uniformly distributed)	

muik <- matrix(lambda,nrow=N,ncol=K,byrow=TRUE) - matrix(phi,ncol=K,nrow=N,byrow=TRUE)*matrix(theta,ncol=K,nrow=N) 
muik <- muik * MBDT

#Compute Extremeness Residuals (posterior probability greater than 2)
RTd <- RT*MBDT
diff <- (RTd-muik)*matrix(sqrt(1/sigma2),ncol=K,nrow=N,byrow=TRUE)
presid <- (1-pnorm(2,mean=diff,sd=1)) + pnorm(-2,mean=diff,sd=1)
presid[MBDT==0] <- 0 ##fix to zero

muik <-  matrix(EAPlambda,nrow=N,ncol=K,byrow=TRUE) - matrix(EAPphi,ncol=K,nrow=N,byrow=TRUE)*matrix(EAPtheta,ncol=K,nrow=N) 
muik <- muik*MBDT
muiklong <- matrix(muik,ncol=1,nrow=N*K)
RTd <- RT*MBDT
RTlong <- matrix(RTd,ncol=1,nrow=N*K)
sigma2long <- matrix(matrix(sqrt(EAPsigma2),ncol=K,nrow=N,byrow=TRUE),ncol=1,nrow=N*K)
errorfit <- (RTlong - muiklong)/sigma2long
errorfit <- matrix(rnorm(N*K)*1e-06 + errorfit,ncol=K,nrow=N) #to remove ties

#Perform one-sample Kolmogorov Smirnov Test across Items

for(kk in 1:K){
	KS[kk,1] <- ks.test(errorfit[MBDT[,kk]==1,kk],y=pnorm)$p.value ##remove missings by design
}

return(list(KS=KS,presid=presid))

}


residualA <- function(Z, Y, theta, alpha, beta, EAPtheta, EAPalpha, EAPbeta) {
    
    K <- ncol(Y)
    N <- nrow(Y)
    KS <- matrix(0, ncol = 1, nrow = K)
    
    # compute fitted probabilities (approximately uniformly distributed)
    ## PF UNDER CONSTRUCTION
	#compute fitted probabilities (approximately uniformly distributed)	
	muik <- t(matrix(EAPalpha,ncol=N,nrow=K))*matrix(EAPtheta,ncol=K,nrow=N)-t(matrix(EAPbeta,ncol=N,nrow=K))
	diff <- matrix((Z-muik)**2,ncol=K,nrow=N)
	lZPAT <- apply(diff,1,sum) 
	lZPA  <- 1-pchisq(lZPAT,K) #person fit
    
    # Compute Extremeness Residuals (posterior probability greater than 2)
    presidA <- matrix((pnorm(-2)/pnorm(muik)), ncol = K, nrow = N) * Y + matrix((pnorm(-2)/(1 - pnorm(muik))), ncol = K, nrow = N) * (1 - Y)
    
    pmuik <- pnorm(muik)
    # l0 is the natural logarithm of the likelihood :
    l0 <- matrix((Y * log(pmuik)) + ((1 - Y) * log(1 - pmuik)), nrow = N, ncol = K)
    lP0 <- rowSums(l0, na.rm = T)
    lI0 <- colSums(l0, na.rm = T)
    
    # Expected l0 :
    El0 <- matrix(pmuik * log(pmuik) + (1 - pmuik) * log(1 - pmuik), nrow = N, ncol = K)
    ElP0 <- rowSums(El0, na.rm = T)
    ElI0 <- colSums(El0, na.rm = T)
    
    # conditional variance :
    Vl <- matrix(pmuik * (1 - pmuik) * ((log(pmuik/(1 - pmuik)))^2), ncol = K, nrow = N)
    VlP0 <- rowSums(Vl, na.rm = T)
    VlI0 <- colSums(Vl, na.rm = T)
    
    # The person Fit , Item Fit:
    PFl <- -(lP0 - ElP0)/sqrt(VlP0)
    IFl <- -(lI0 - ElI0)/sqrt(VlI0)
    
    PFlp <- 1 - pnorm(PFl)
    IFlp <- 1 - pnorm(IFl)
    
    muik <- t(matrix(alpha, ncol = N, nrow = K)) * matrix(theta, ncol = K, nrow = N) - t(matrix(beta, ncol = N, nrow = K))
    muiklong <- matrix(muik, ncol = 1, nrow = N * K)
    Zlong <- matrix(Z, ncol = 1, nrow = N * K)
    errorfit <- (Zlong - muiklong)
    errorfit <- matrix(rnorm(N * K) * 1e-06 + errorfit, ncol = K, nrow = N)  #to remove ties
    
    # Perform one-sample Kolmogorov Smirnov Test across Items
    
    for (kk in 1:K) {
        KS[kk, 1] <- ks.test(errorfit[, kk], y = pnorm)$p.value
    }
    
    # return(list(KS=KS,presidA=presidA,lZPAT=lZPAT,lZPA=lZPA,lZIA=lZIA,PFl=PFl,IFl=IFl,PFlp=PFlp,IFlp=IFlp,l0=l0))
    
    return(list(KS = KS, presidA = presidA, PFl = PFl, IFl = IFl, PFlp = PFlp, IFlp = IFlp, l0 = l0, lZPA = lZPA))
    
}


residualAMBD <- function(Z,Y,theta,alpha,beta,EAPtheta,EAPalpha,EAPbeta,MBDY){

## UNDER CONSTRUCTION ##
## lZPA not working 

 K <- ncol(Y)
 N <- nrow(Y)
 KS <- matrix(0,ncol=1,nrow=K)

#compute fitted probabilities (approximately uniformly distributed)	

## PF UNDER CONSTRUCTION
#compute fitted probabilities (approximately uniformly distributed)	
 muik <- t(matrix(EAPalpha,ncol=N,nrow=K))*matrix(EAPtheta,ncol=K,nrow=N)-t(matrix(EAPbeta,ncol=N,nrow=K))
 diff <- matrix((Z-muik)**2,ncol=K,nrow=N)*MBDY
 Kt <- apply(MBDY,1,sum)
 lZPAT <- apply(diff,1,sum) 
 lZPA  <- 1-pchisq(lZPAT,Kt) #person fit

#Compute Extremeness Residuals (posterior probability greater than 2)
presidA <- matrix((pnorm(-2)/pnorm(muik)),ncol=K,nrow=N)*Y + matrix((pnorm(-2)/(1-pnorm(muik))),ncol=K,nrow=N)*(1-Y)
presidA[MBDY==0] <- 0 ##fix to zero

  pmuik <- pnorm(muik)
  # l0 is the natural logarithm of the likelihood :
  Y[MBDY==0] <- 0	
  l0 <- matrix((Y*log(pmuik)) + ((1-Y)*log(1-pmuik)),nrow=N,ncol=K)*MBDY
  lP0 <- rowSums(l0,na.rm=T)
  lI0 <- colSums(l0,na.rm=T)

  # Expected l0 :
  El0 <- matrix(pmuik*log(pmuik) + (1-pmuik)*log(1-pmuik),nrow=N,ncol=K)*MBDY
  ElP0 <- rowSums(El0,na.rm=T)
  ElI0 <- colSums(El0,na.rm=T)

  # conditional variance :
  Vl <- matrix(pmuik*(1-pmuik)*((log(pmuik/(1-pmuik)))^2),ncol=K,nrow=N)*MBDY
  VlP0 <- rowSums(Vl,na.rm=T)
  VlI0 <- colSums(Vl,na.rm=T)
  
  # The person Fit , Item Fit:
  PFl <- -(lP0-ElP0)/sqrt(VlP0)
  IFl <- -(lI0-ElI0)/sqrt(VlI0)
  
  PFlp <- 1-pnorm(PFl)
  IFlp <- 1-pnorm(IFl)

 muik <- (t(matrix(alpha,ncol=N,nrow=K))*matrix(theta,ncol=K,nrow=N)-t(matrix(beta,ncol=N,nrow=K)))
 muiklong <- matrix(muik,ncol=1,nrow=N*K)
 Zlong <- matrix(Z,ncol=1,nrow=N*K)
 errorfit <- (Zlong - muiklong)
 errorfit <- matrix(rnorm(N*K)*1e-06 + errorfit,ncol=K,nrow=N) #to remove ties

#Perform one-sample Kolmogorov Smirnov Test across Items

for(kk in 1:K){
	KS[kk,1] <- ks.test(errorfit[MBDY[,kk]==1,kk],y=pnorm)$p.value
}

#return(list(KS=KS,presidA=presidA,lZPAT=lZPAT,lZPA=lZPA,lZIA=lZIA,PFl=PFl,IFl=IFl,PFlp=PFlp,IFlp=IFlp,l0=l0))

return(list(KS=KS,presidA=presidA,PFl=PFl,IFl=IFl,PFlp=PFlp,IFlp=IFlp,l0=l0,lZPA=lZPA))

}
 
