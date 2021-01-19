\dontrun{
  
  ###
  ### EXAMPLE APPLICATION AMSTERDAM CHESS DATA van der Maas and Wagenmakers (2005).
  ###
  
  library(LNIRT)
  data(AmsterdamChess)
  head(AmsterdamChess)
  
  N <- nrow(AmsterdamChess)
  Y <-as.matrix(AmsterdamChess[c(which(colnames(AmsterdamChess)=="Y1")
                                 :which(colnames(AmsterdamChess)=="Y40"))])
  K <- ncol(Y)
  ## replace missing 9 with NA
  Y[Y==9] <- NA
  ## Test takers with NAs
  ## Y[147,]
  ## Y[201,]
  ## Y[209,]
  
  RT <- as.matrix(AmsterdamChess[c(which(colnames(AmsterdamChess)=="RT1")
                                   :which(colnames(AmsterdamChess)=="RT40"))])
  ## replace missing 10000.000 with NA
  RT[RT==10000.000] <- NA
  RT<-log(RT) #logarithm of RTs
  
  # Define Time Scale
  X <- 1:K
  X <- (X - 1)/K
  
  set.seed(12345) ## used to obtain the results reported in the paper ##
  outchess <- LNIRTQ(Y=Y,RT=RT,X=X,XG=10000)
  summary(outchess)
  
  ## Check MCMC convergence
  ##
  ## check several MCMC chains 
  ## 
  
  library(mcmcse)
  ess(outchess1$MAB[1001:10000,1,1]) ## effective sample size
  mcse(outchess1$MAB[1001:10000,1,1], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(outchess$MAB[1001:10000,1,2]) ## effective sample size
  mcse(outchess$MAB[1001:10000,1,2], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(outchess$MAB[1001:10000,1,3]) ## effective sample size
  mcse(outchess$MAB[1001:10000,1,3], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(outchess$MAB[1001:10000,1,4]) ## effective sample size
  mcse(outchess$MAB[1001:10000,1,4], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(outchess$MSP[1001:10000,1,4]) ## effective sample size
  mcse(outchess$MSP[1001:10000,1,4], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(outchess$MSI[1001:10000,1,4]) ## effective sample size
  mcse(outchess$MSI[1001:10000,1,4], size = 100
       , g = NULL,method = "bm", warn = FALSE) #standard error
  
  ## Convergence Checks
  library(coda) 
  summary(as.mcmc(outchess$MAB[1001:10000,1,1]))
  summary(as.mcmc(outchess$MAB[1001:10000,1,4]))
  summary(as.mcmc(outchess$MSI[1001:10000,1,1]))
  summary(as.mcmc(outchess$MSI[1001:10000,1,4]))
  summary(as.mcmc(outchess$MSP[1001:10000,1,4]))
  summary(as.mcmc(outchess$MSP[1001:10000,1,1]))
  
  ## check some chains on convergence
  geweke.diag(as.mcmc(outchess$MAB[1001:10000,1,1]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(outchess$MAB[1001:10000,1,1]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(outchess$MAB[1001:10000,1,1]), eps=0.1, pvalue=0.05)
  
  geweke.diag(as.mcmc(outchess$MSI[1001:10000,1,1]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(outchess$MSI[1001:10000,1,1]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(outchess$MSI[1001:10000,1,1]), eps=0.1, pvalue=0.05)
  
  geweke.diag(as.mcmc(outchess$MSP[1001:10000,3,3]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(outchess$MSP[1001:10000,3,3]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(outchess$MSP[1001:10000,3,3]), eps=0.1, pvalue=0.05)
  
  ## complete missing data
  outchess$Mtheta[147,]
  
  ######################################################################
  ### THIS PART IS NOT DISCUSSED IN THE PAPER                        ###
  ######################################################################
  
  # PLOT PERSON PARAMETERS
  # ABILITY VS SPEED
  
  par(mar=c(5, 5, 2,4), xpd=F)
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(2,2), heights=c(2,2))
  
  plot(outchess$Mtheta[,1],outchess$Mtheta[,2]
       ,xlab=expression(paste("Ability"~~(theta)))
       ,ylab=expression(paste("Intercept"~~(zeta[0]))),
       xlim=c(-2,2),ylim=c(-1,1),bty="l",cex.lab=1.5,cex.axis=1.25)
  abline(lm(outchess$Mtheta[,2]~outchess$Mtheta[,1])) 
  abline(h = 0,lty = 2)
  
  plot(outchess$Mtheta[,1],outchess$Mtheta[,3]
       ,xlab=expression(paste("Ability"~~(theta)))
       ,ylab=expression(paste("Linear slope"~~(zeta[1]))),
       xlim=c(-2,2),ylim=c(-1,1),bty="l",cex.lab=1.5,cex.axis=1.25)
  abline(lm(outchess$Mtheta[,3]~outchess$Mtheta[,1]))
  abline(h = 0,lty = 2)
  
  plot(outchess$Mtheta[,1],outchess$Mtheta[,4]
       ,xlab=expression(paste("Ability"~~(theta)))
       ,ylab=expression(paste("Quadratic slope"~~(zeta[2]))),
       xlim=c(-2,2),ylim=c(-1,1),bty="l",cex.lab=1.5,cex.axis=1.25)
  abline(lm(outchess$Mtheta[,4]~outchess$Mtheta[,1]))
  abline(h = 0,lty = 2)
  
  plot(outchess$Mtheta[,3],outchess$Mtheta[,4]
       ,xlab=expression(paste("Linear slope"~~(zeta[1])))
       ,ylab=expression("Quadratic slope"~~paste((zeta[2]))),
       xlim=c(-0.5,1),ylim=c(-0.5,0.5),bty="l",cex.lab=1.5,cex.axis=1.25)
  abline(lm(outchess$Mtheta[,4]~outchess$Mtheta[,3]))
  abline(h = 0,lty = 2)
  
  ### include residual analysis ###
  
  set.seed(12345)
  outchessr <- LNIRTQ(Y=Y,RT=RT,X=X,XG=10000,burnin=10,XGresid=1000,resid=TRUE)
  summary(outchessr)
  
  ## plot of person-fit scores for RT patterns
  plot(outchessr$lZPT,outchessr$lZP)
  
}