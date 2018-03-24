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
#' @param kpa
#' the number of predictors for the person ability parameters (optional).
#' @param kpt
#' the number of predictors for the person speed parameters (optional).
#' @param kia
#' the number of predictors for the item-difficulty parameters (optional).
#' @param kit
#' the number of predictors for the item time intensity parameters (optional).
#' 
#' @return 
#' an object of class simLNIRT.
#' 
#' @export
simLNIRT <- function(N, K, rho, td = FALSE, WL = FALSE, kpa, kpt, kia, kit) {
  #WL = TRUE: time discrimination = 1/sqrt(sigma2)
  #td = TRUE: time discrimination <- 1
  
  if (missing(kpa)) {
    kpa <- 0
    XPA <- NULL
    Ba <- NULL
    #Ba <- 0
  }
  else if (kpa == 0) {
    XPA <- NULL
    Ba <- NULL
  }
  else {
    XPA <- matrix(rnorm(kpa*N,sd=.5),ncol=kpa,nrow=N) #predictors ability
    XPA <- matrix(XPA,ncol=ncol(XPA),nrow=N) - matrix(apply(XPA,2,mean),ncol=ncol(XPA),nrow=N,byrow=T)
    Ba <- matrix(rnorm(kpa),ncol=1,nrow=kpa)
  }
  
  if(missing(kpt)) {
    kpt <- 0
    XPT <- NULL
    Bt <- NULL
    #Bt <- 0
  }
  else if (kpt == 0) {
    XPT <- NULL
    Bt <- NULL
  }
  else {
    XPT <- matrix(rnorm(kpt*N,sd=.5),ncol=kpt,nrow=N) #predictors speed
    XPT <- matrix(XPT,ncol=ncol(XPT),nrow=N) - matrix(apply(XPT,2,mean),ncol=ncol(XPT),nrow=N,byrow=T)
    Bt <- matrix(rnorm(kpt),ncol=1,nrow=kpt)
  }
  
  if(kpa!=0){
    XBa <- XPA%*%Ba
  }else{
    XBa <- rep(0,N)
  }
  
  if(kpt!=0){
    XBt <- XPT%*%Bt
  }else{
    XBt <- rep(0,N)
  }
  
  mutheta <- rep(0, 2)
  covtheta <- diag(2)
  covtheta[1, 2] <- covtheta[2, 1] <- rho
  theta <- MASS::mvrnorm(N, mutheta, covtheta, empirical = TRUE)
  
  theta[,1] <- theta[,1]+XBa	#add predictor effects ability
  theta[,2] <- theta[,2]+XBt	#add predictor effects speed
  
  
  if (missing(kia)) {
    kia <- 0
    XIA <- NULL
    Bia <- NULL
  }
  else if (kia == 0) {
    XIA <- NULL
    Bia <- NULL
  }
  else {
    XIA <- matrix(rnorm(kia*K,sd=.25),ncol=kia,nrow=N) #predictors ability
    XIA <- matrix(XIA,ncol=ncol(XIA),nrow=K) - matrix(apply(XIA,2,mean),ncol=ncol(XIA),nrow=K,byrow=T)
    Bia <- matrix(rnorm(kia,sd=.75),ncol=1,nrow=kia)
  }
  
  if(missing(kit)) {
    kit <- 0
    XIT <- NULL
    Bit <- NULL
  }
  else if (kit == 0) {
    XIT <- NULL
    Bit <- NULL
  }
  else {
    XIT <- matrix(rnorm(kit*K,sd=.25),ncol=kit,nrow=K) #predictors speed
    XIT <- matrix(XIT,ncol=ncol(XIT),nrow=K) - matrix(apply(XIT,2,mean),ncol=ncol(XIT),nrow=K,byrow=T)
    Bit <- matrix(rnorm(kit,sd=.75),ncol=1,nrow=kit)
  }
  
  if(kia!=0){
    XBIa <- XIA%*%Bia
  }else{
    XBIa <- rep(0,K)
  }
  
  if(kit!=0){
    XBIt <- XIT%*%Bit
  }else{
    XBIt <- rep(0,K)
  }
  
  if(kia!=0 & kit != 0){
    covitem <- diag(c(.05,.3,.05,.3))
  }
  if(kia!=0 & kit == 0){
    covitem <- diag(c(.05,.3,.05,1))	
  }
  if(kia==0 & kit != 0){
    covitem <- diag(c(.05,1,.05,.3))	
  }
  if(kia==0 & kit == 0){
    covitem <- diag(c(.05,1,.05,1))	
  }
  
  muitem <- rep(c(1, 0), 2)
  sigma2 <- rlnorm(K, meanlog = 0, sdlog = 0.3)
  
  if (td && !WL) {
    ab <- MASS::mvrnorm(K, muitem, covitem)
    ab[, 3] <- rep(1, K)
    abn <- MASS::mvrnorm(K, muitem[-3], covitem[-3, -3])
    ab[, c(2, 4)] <- abn[, c(2, 3)] - t(matrix(colMeans(abn[, c(2, 3)]), 2, K))
    ab[, 1] <- abs(ab[, 1])
    ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
  }
  if (WL) {
    muitem <- c(1, 0, 1, 0)
    ab <- MASS::mvrnorm(K, muitem, covitem)
    sigma2 <- 1/ab[, 3]^2
    ab[, 1] <- abs(ab[, 1])
    ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
  }
  if (!td && !WL) {
    ab <- MASS::mvrnorm(K, muitem, covitem)
    ab[, c(2, 4)] <- ab[, c(2, 4)] - t(matrix(colMeans(ab[, c(2, 4)]), 2, K))
    ab[, 1] <- abs(ab[, 1])
    ab[, 1] <- ab[, 1]/(prod(ab[, 1])^(1/K))
    ab[, 3] <- abs(ab[, 3])
    ab[, 3] <- ab[, 3]/(prod(ab[, 3])^(1/K))
  }
  
  
  #add item predictors
  ab[,2] <- ab[,2] + XBIa #difficulty
  ab[,4] <- ab[,4] + XBIt	#time intensity
  
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
  
  out <- list(Y = Y, Yg = Yg, Y1 = Y1, Y1g = Y1g, RT = RT, RT1 = RT1, theta = theta, ab = ab, sigma2 = sigma2, quess = quess,
              XPA = XPA, XPT = XPT, Ba = Ba, Bt = Bt, XIA = XIA, XIT = XIT, Bia = Bia, Bit = Bit)
  class(out) <- "simLNIRT"
  return(out)
}

