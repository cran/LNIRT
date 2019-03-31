## Fox March 2015
##  
## Code from LNRT (Log-Normal Response Time) plus random person effects (intercept, trend, quadratic)
## Use functions from LNIRT

#' Log-normal response time modelling with variable person speed (intercept, trend, quadratic)
#' 
#' @param RT
#' a Person-x-Item matrix of log-response times (time spent on solving an item).
#' @param X
#' explanatory (time) variables for random person speed (default: (1:N.items - 1)/N.items). 
#' @param XG
#' the number of MCMC iterations to perform (default: 1000).
#' 
#' @return 
#' an object of class LNRTQ. 
#' @export
LNRTQ <- function(RT, X, XG = 1000){
  ## RT = log-response time matrix (time spent on solving an item) of dim(N=persons,K=items) 
  ## X = explanatory (time) variables for random person speed 		
  ## XG = number of iterations for the MCMC algorithm

  N <- nrow(RT) #persons
  K <- ncol(RT) #items (complete design)	
  Q <- 3 #(intercept, slope, quadratic) speed components	
	
  if(missing(X)) {
   X <- 1:K
   X <- (X-1)/K
  }
  
## population zeta (speed)
  muP <- matrix(0,ncol=1,nrow=Q)
  SigmaP <-  diag(Q)
  muP0 <- matrix(0,1,Q)
  SigmaP0 <-  diag(Q)/100
  SigmaR <-  diag(Q)

  ab <- matrix(rnorm(K*2),ncol=2)
  ab[,1] <- 1
  muI <- t(matrix(rep(c(1,0),K),ncol=K))
  muI0 <- muI[1,]
  SigmaI <- diag(2)
  SigmaI0 <- diag(2)*10
  for(ii in 1:2) {
	SigmaI0[ii,] <- SigmaI0[ii,]*c(.5,3)
  }

##storage
  MZ <- MZ2 <-list() 	
  for(jj in 1:Q){
	  MZ[[jj]] <- MZ2[[jj]] <- matrix(0,ncol=1,nrow=N)
  }	
  MAB <- array(0,dim=c(XG,K,2))
  MmuP <- matrix(0,ncol=Q,nrow=XG)
  MmuI <- array(0,dim=c(XG,2))
  MSP <- array(0,dim=c(XG,Q,Q))
  MSI <- array(0,dim=c(XG,2,2))
  sigma2 <- rep(1,K)
  Msigma2 <- matrix(0,ncol=K,nrow=XG)

  ## Start MCMC algorithm
  
  for (ii in 1:XG){
     dum <- DrawZetaQ1(RT=RT,phi=ab[,1],lambda=ab[,2],sigma2=sigma2,mu=muP,SigmaR=SigmaR,X=X) ## variable speed 
	speed <- dum$zetapred #product of X times zeta
	zeta <- dum$zeta

	### restriction average speed is equal to zero 
	zeta <- zeta - matrix(apply(zeta,2,mean),ncol=Q,nrow=N,byrow=T)
	for(jj in 1:Q){
	    ##zeta[1:N,jj] <- zeta[1:N,jj] - mean(zeta[,jj]) ### restriction average speed is equal to zero 
	    MZ[[jj]][1:N] <- MZ[[jj]][1:N] + zeta[1:N,jj]
	    MZ2[[jj]][1:N] <- MZ2[[jj]][1:N] + zeta[1:N,jj]^2 
	}
    	
    	## Sample the item parameters of all models and rescale for identification
	    dum <- DrawLambdaPhiQ1(RT=RT,zeta=zeta,X=X,sigma2=sigma2,muI=muI,SigmaI=SigmaI)
	    ab[,1] <- dum$phi	
	    ab[,1] <- ab[,1]/(prod(ab[,1])^(1/K))
	    ab[,2] <- dum$lambda
	    ##ab[,2] <- ab[,2] - mean(ab[,2]) ### restriction average time-intensity equal to zero

    MAB[ii,1:K,1:2] <- ab
    sigma2 <- SampleS2Q1(RT=RT,zeta=zeta,X=X,lambda=ab[,2],phi=ab[,1])
    Msigma2[ii,] <- sigma2 
    
    ## 2nd level model for person
    X1 <- matrix(1,N,1)
    dum <- SampleBQ1(Y=zeta,X=X1,Sigma=SigmaR,B0=muP0,V0=SigmaP0)
    MmuP[ii,] <- dum$B
    muP <- dum$B

	## Fix Covariance Terms	
	SS <- diag(crossprod(zeta - dum$pred))
	SigmaR <- diag((SS + 1)/rchisq(Q,N))
	MSP[ii,,] <- SigmaR
		## Sample covmatrix persons
		##SS <- crossprod(zeta - dum$pred) + SigmaP0
		##SigmaR <- rwishart(Q + N, chol2inv(chol(SS)))$IW
		##MSP[ii,,] <- SigmaR

    ## 2nd level model for person
	    X1 <- matrix(1,K,1)
	    muI2 <- SampleBQ1(Y=ab,X=X1,Sigma=SigmaI,B0=muI0,V0=SigmaI0)
	    MmuI[ii,1] <- muI2$B[1]
	    MmuI[ii,2] <- muI2$B[2]
	    muI[,1] <- muI2$pred[,1]
	    muI[,2] <- muI2$pred[,2]
	
	    SS <- crossprod(ab - muI) + SigmaI0
	    SigmaI <- rwishart(2 + K, chol2inv(chol(SS)))$IW
	    MSI[ii,,] <- SigmaI

    	if(ii %% 100 == 0) cat("Iteration ",ii,"\n")
	flush.console()
  }

	## end iterations

  for(jj in 1:Q){
	  MZ[[jj]] <- MZ[[jj]]/XG
	  MZ2[[jj]] <- sqrt(MZ2[[jj]]/XG - MZ[[jj]]^2 ) ## calculate posterior SD of person parameters. 
  }	

  out <- list(Mzeta=MZ,MZSD =MZ2, MAB=MAB,MmuP = MmuP,MSP= MSP, MmuI=MmuI,
				MSI = MSI,Msigma2 = Msigma2)
  class(out) <- c("LNRTQ", "list")
  return(out)
}


