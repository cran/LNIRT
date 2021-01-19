\dontrun{
  
  ###
  ### EXAMPLE APPLICATION CREDENTIAL FORM1 DATA CIZEK and WOLLACK (2016).
  ###
  
  library(LNIRT)
  data(CredentialForm1)
  
  ### DATA OBJECTS FOR LNIRT 
  ### RA Data
  Y <- as.matrix(CredentialForm1[c(which(colnames(CredentialForm1)=="iraw.1")
                                  :which(colnames(CredentialForm1)=="iraw.170"))])
  N <- nrow(Y) 
  
  ### RT Data 
  RT<-as.matrix(CredentialForm1[c(which(colnames(CredentialForm1)=="idur.1")
                                  :which(colnames(CredentialForm1)=="idur.170"))])
  RT[RT==0]<-NA ## zero RTs are coded as missing values
  RT<-log(RT) ## logarithmic transformation of RT
  
  ## RUN LNIRT MODEL 0
  set.seed(12345) ## used to obtain the results reported in the paper ##
  out0 <- LNIRT(RT=RT,Y=Y,XG=5000,burnin=10,ident=2)
  summary(out0)
  
  ## Check MCMC convergence
  
  library(mcmcse)
  ##
  ## check several MCMC chains 
  ## 
  ## effective sample size and effective sample size
  ess(out0$MCMC.Samples$Cov.Person.Ability.Speed[1001:5000]) ## effective sample size
  mcse(out0$MCMC.Samples$Cov.Person.Ability.Speed[1001:5000]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Var.Person.Ability[1001:5000]) ## effective sample size
  mcse(out0$MCMC.Samples$Var.Person.Ability[1001:5000]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Var.Person.Speed[1001:5000]) ## effective sample size
  mcse(out0$MCMC.Samples$Var.Person.Speed[1001:5000]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Item.Discrimination[1001:5000,155]) ## effective sample size
  mcse(out0$MCMC.Samples$Item.Discrimination[1001:5000,155]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Time.Discrimination[1001:5000,1]) ## effective sample size
  mcse(out0$MCMC.Samples$Time.Discrimination[1001:5000,1]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Person.Ability[1001:5000,1]) ## effective sample size
  mcse(out0$MCMC.Samples$Person.Ability[1001:5000,1]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ess(out0$MCMC.Samples$Person.Speed[1001:5000,1]) ## effective sample size
  mcse(out0$MCMC.Samples$Person.Speed[1001:5000,1]
       , size = 100, g = NULL,method = "bm", warn = FALSE) #standard error
  
  ## Convergence Checks
  library(coda) 
  summary(as.mcmc(out0$MCMC.Samples$Cov.Person.Ability.Speed[1001:5000]))
  summary(as.mcmc(out0$MCMC.Samples$Var.Person.Ability[1001:5000]))
  summary(as.mcmc(out0$MCMC.Samples$Var.Person.Speed[1001:5000]))
  summary(as.mcmc(out0$MCMC.Samples$Item.Discrimination[1001:5000,155]))
  summary(as.mcmc(out0$MCMC.Samples$Time.Discrimination[1001:5000,1]))
  summary(as.mcmc(out0$MCMC.Samples$Person.Ability[1001:5000,1]))
  summary(as.mcmc(out0$MCMC.Samples$Person.Speed[1001:5000,1]))
  
  ## check some chains on convergence
  geweke.diag(as.mcmc(out0$MCMC.Samples$Cov.Person.Ability.Speed[500:5000]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(out0$MCMC.Samples$Cov.Person.Ability.Speed[500:5000]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(out0$MCMC.Samples$Cov.Person.Ability.Speed[500:5000], eps=0.1, pvalue=0.05))
  
  geweke.diag(as.mcmc(out0$MCMC.Samples$Item.Discrimination[500:5000,155]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(out0$MCMC.Samples$Item.Discrimination[500:5000,155]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(out0$MCMC.Samples$Item.Discrimination[500:5000,155]), eps=0.1, pvalue=0.05)
  
  geweke.diag(as.mcmc(out0$MCMC.Samples$Person.Ability[500:5000,1]), frac1=0.1, frac2=0.5)
  geweke.plot(as.mcmc(out0$MCMC.Samples$Person.Ability[500:5000,1]), frac1=0.1, frac2=0.5)
  heidel.diag(as.mcmc(out0$MCMC.Samples$Person.Ability[500:5000,1]), eps=0.1, pvalue=0.05)
  
  ## Item parameter estimates
  
  min(apply(out0$MAB[500:5000,,1],2,mean))
  max(apply(out0$MAB[500:5000,,1],2,mean))
  
  min(apply(out0$MAB[500:5000,,2],2,mean))
  max(apply(out0$MAB[500:5000,,2],2,mean))
  
  min(apply(out0$MAB[500:5000,,3],2,mean))
  max(apply(out0$MAB[500:5000,,3],2,mean))
  
  min(apply(out0$MAB[500:5000,,4],2,mean))
  max(apply(out0$MAB[500:5000,,4],2,mean))
  
  plot(apply(out0$MAB[500:5000,,4],2,mean),(apply(RT,2,mean,na.rm=TRUE)))
  
  ### Explanatory Variables Test-takers
  XFT <- data.frame(CredentialForm1[1:10],stringsAsFactors=TRUE) #Background Variables 
  XFT$Tot_time <- (XFT$Tot_time-mean(XFT$Tot_time))/sqrt(var(XFT$Tot_time))
  
  ## DUMMY CODING FOR CATEGORICAL PREDICTORS
  ## Pretest Groups
  XFT$Pgroup <- matrix(0,ncol=2,nrow=N)
  XFT$Pgroup[XFT$Pretest==6,1] <- -1  
  XFT$Pgroup[XFT$Pretest==6,2] <- -1
  XFT$Pgroup[XFT$Pretest==7,1] <- 1
  XFT$Pgroup[XFT$Pretest==8,2] <- 1
  
  ## Countries 
  XFT$Cgroup <- matrix(0,ncol=3,nrow=N)
  XFT$Cgroup[XFT$Country=="USA",1] <- 1 
  XFT$Cgroup[XFT$Country=="Philippines",2] <- 1
  XFT$Cgroup[XFT$Country=="India",3] <- 1
  XFT$Cgroup[c(XFT$Country!="USA" & XFT$Country!="India" & XFT$Country!="Philippines"),1:3] <- -1
  
  XA <- matrix(unlist(XFT[,c("Pgroup","Tot_time")]),ncol=3,nrow=N)
  XT <- matrix(unlist(XFT[,c("Pgroup")]),ncol=2,nrow=N)
  
  ## RUN LNIRT MODEL 1 (Pretest and total test time)
  ## Include residual analysis 
  set.seed(12345) ## used to obtain the results reported in the paper ##
  out1 <- LNIRT(RT=RT,Y=Y,XG=5000,XPA=XA,XPT=XT,residual=TRUE)
  summary(out1)
  
  
  ######################################################################
  ### THIS PART IS NOT DISCUSSED IN THE PAPER                        ###
  ######################################################################
  
  ## RUN LNIRT MODEL 2 (Pretest and Country)
  XA <- matrix(unlist(XFT[,c("Pgroup","Cgroup")]),ncol=5,nrow=N)
  XT <- matrix(unlist(XFT[,c("Pgroup","Cgroup")]),ncol=5,nrow=N)
  
  set.seed(12345) ## 
  out2 <- LNIRT(RT=RT,Y=Y,XG=5000,XPA=XA,XPT=XT)
  summary(out2)
  
  XA <- matrix(unlist(XFT[,c("Pgroup","Cgroup","Tot_time")]),ncol=6,nrow=N)
  XT <- matrix(unlist(XFT[,c("Pgroup","Cgroup")]),ncol=5,nrow=N)
  
  ## RUN LNIRT MODEL 3 
  set.seed(12345) ## 
  out3 <- LNIRT(RT=RT,Y=Y,XG=5000,XPA=XA,XPT=XT)
  summary(out3)
  
  #########################################################################
  #########################################################################
  #########################################################################
  
  
  ######################################################################
  ### THIS PART IS DISCUSSED IN THE PAPER                            ###
  ######################################################################
  
  ## Subsection "Planned Missing By Design"
  ## Include pretest item data
  MBDM<-matrix(rep(0,1636*200),nrow=1636,ncol=200)
  MBDM[XFT$Pretest==6,171:180]<-1
  MBDM[XFT$Pretest==7,181:190]<-1
  MBDM[XFT$Pretest==8,191:200]<-1
  MBDM[,1:170]<-1
  
  Yt <- CredentialForm1[c(which(colnames(CredentialForm1)=="iraw.1")
                          :which(colnames(CredentialForm1)=="iraw.200"))]
  ## transform pretest data to numeric
  Yt[,171:200] <- unlist(lapply(Yt[,171:200]
                                ,function(x) as.numeric(x))) #warnings about NA can be ignored
  Yt <- as.matrix(Yt,ncol=200,nrow=1636)
  
  RTt <- (CredentialForm1[as.numeric(c(which(colnames(CredentialForm1)=="idur.1")
                                       :which(colnames(CredentialForm1)=="idur.200")))])
  RTt[,171:200] <- unlist(lapply(RTt[,171:200]
                   , function(x) as.numeric(as.character(x)))) #warnings about NA can be ignored
  RTt[RTt==0] <- NA ## zero RTs are coded as missing values
  RTt <- log(RTt) ## logarithmic transformation of RT
  RTt <- as.matrix(RTt,ncol=200,nrow=1636)
  
  # To fit the model, item discrimination parameters are restricted to one.  
  alpha1<-rep(1,200) ### Pre-defined item discrimination parameters
  ## RUN LNIRT MODEL 4
  set.seed(12345) ## used to obtain the results reported in the paper ##
  out4 <- LNIRT(RT=RTt,Y=Yt,XG=5000,alpha=alpha1,MBDY=MBDM,MBDT=MBDM)
  summary(out4)
  
  ### Subsection "Model-Fit Analysis"
  ### Return to output of out1
  
  #report fit results
  summary(out1)
  
  ## estimated average residual variance
  mean(out1$Msigma2[500:5000,])
  
  #recoding of number of zero attempts 
  XFT$Attempt[XFT$Attempt==0] <- 1
  
  ## explain heterogeneity in person-fit statistics RA and RT
  summary(lm(out1$PFl ~ as.factor(XFT$Attempt)+(XFT$Cgroup)+(XFT$Pgroup)))
  summary(lm(out1$lZPT ~ as.factor(XFT$Attempt)+(XFT$Cgroup)+(XFT$Pgroup)))
  
  ### overview plot of person fit RA versus person-fit RT per country 
  dev.new()
  plot(out1$PFl,out1$lZPT,xlab="Person-fit Statistic RA",ylab="Person-fit Statistic RT",
       col="black",cex=.5,bty="l",xlim=c(-3,3)
       , ylim=c(0,500),cex.main=.8,cex.axis=.7,cex.lab=.8,pch=15)
  ## US
  set1 <- which(XFT$Country=="USA")
  points(out1$PFl[set1],out1$lZPT[set1],col="blue",pch=10,cex=.5)
  ## India
  set2 <- which(XFT$Country=="India")
  points(out1$PFl[set2],out1$lZPT[set2],col="red",pch=13,cex=.5)
  ## Philippines
  set3 <- which(XFT$Country=="Philippines")
  points(out1$PFl[set3],out1$lZPT[set3],col="green",pch=16,cex=.5)
  abline(h = qchisq(.95, df= 170),lty = 2,col="red")
  abline(v = qnorm(.95),lty = 2,col="red")
  legend(-3,500,c("India","US","Philippines","Other"), 
         col=c("red","blue","green","black"),pch = c(13,10,16,15), bg = "gray95",cex=.7)
  
  ###################################################################################
  ###################################################################################
  ###################################################################################
  
}