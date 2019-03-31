#' Simulate data for log-normal response time IRT modelling with variable person speed (intercept, trend, quadratic)
#' 
#' @param N
#' the number of persons.
#' @param K
#' the number of items.
#' @param ...
#' optional arguments.
#' 
#' @return 
#' an object of class simLNIRTQ.
#' 
#' @export
simLNIRTQ <- function(N,K,...){

##variable speed 
##zeta = speed components
##speed = realized speed per item 

  # N respondents, 
  # K items
  Q <- 4 #ability and speed components	

	dots <- list(...)


# define time scale
  X <- 1:K
  X <- (X-1)/K
  RT <- speed <- matrix(0,ncol=K,nrow=N)

  Sigmazeta <- matrix(0,ncol=Q-1,nrow=Q-1)	
  diag(Sigmazeta) <- (c(1,.50,.50))	
  
	if(hasArg(zeta)){
		zeta <- matrix(dots$zeta,nrow=N,ncol=Q-1)
	}else{
		zeta <- mvrnorm(N,mu=rep(0,ncol(Sigmazeta)),Sigma=Sigmazeta)
	}

#sample theta given speed components
SigmaR <- matrix(c(1,.7,.2,.1,.7,1,0,0,.2,0,.25,0,.1,0,0,.5),ncol=Q)
#dd[1,1] - dd[1,2:4]%*%solve(dd[2:4,2:4])%*%dd[2:4,1]
#covas = cov(theta,zeta)
covas <- matrix(c(.7,.2,.1),ncol=Q-1)

	if(hasArg(theta)){
		theta <- matrix(dots$theta,nrow=N,ncol=1)
	}else{
		mutheta <- 0 + t((covas %*%solve(Sigmazeta))%*%t(zeta))
		sdtheta <- sqrt(1 - covas%*%solve(Sigmazeta)%*%t(covas))
		theta <- rnorm(N,mean=mutheta,sd=sdtheta) 
	}

## item parameters 
covitem <- diag(Q) 
for(ii in 1:4) {
	covitem[ii,] <- covitem[ii,]*rep(c(.05,1),2)
}
muitem <- rep(c(1,0),2)
	if(hasArg(ab)){
  	  ab <- matrix(ab,ncol=4,nrow=K)	
	}else{
	  ab <- mvrnorm(K, muitem, covitem)
	  ab[,c(2,4)] <- ab[,c(2,4)] - t(matrix(colMeans(ab[,c(2,4)]),2,K))
	  ab[,1] <- abs(ab[,1])
	  ab[,1] <- ab[,1]/(prod(ab[,1])^(1/K))
	  ab[,3] <- abs(ab[,3])
	  ab[,3] <- ab[,3]/(prod(ab[,3])^(1/K))
	}
  
  # itemcorrect
  par    <-  theta %*% matrix(ab[,1],nrow=1,ncol=K)-t( matrix(ab[,2],nrow=K,ncol=N))
  probs <-	matrix(pnorm(par),ncol = K,nrow = N)	
  Y	  <-	matrix(runif(N*K),nrow = N, ncol = K)
  Y	  <-	ifelse(Y < probs,1,0)
  
  Xn <- cbind(ab[,3],ab[,3]*X,ab[,3]*X**2) #equal X matrix over persons
  speed <- zeta%*%t(Xn)

  sigma2 <- rep(.5,K)
  # response time
  for (kk in 1:K) {
    RT[1:N,kk] <- ab[kk,4] - speed[1:N,kk] + rnorm(N,mean=0,sd=sqrt(sigma2[kk]))
  }
  
  out <- list(Y=Y,RT=RT,X=X,theta=theta,speed=speed,ab=ab,zeta=zeta,SigmaR=SigmaR,sigma2=sigma2)
  class(out) <- c("simLNIRTQ", "list")
  return(out)
}
