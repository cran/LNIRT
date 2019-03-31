## Fox 2019

#' Summary Function for LNIRTQ
#' 
#' @param out
#' a LNIRTQ object (the fitted model)
#' @param data
#' a simLNIRTQ object (the simulated data, optional)
#' @export
summaryIRTQ <- function(out, data){  
  
  if(missing(data)){
    simv <- FALSE
  } 
  else{
    simv <- TRUE
  }  
  
  if("Mguess" %in% names(out)){
    guess <- TRUE 
  }
  else{
    guess <- FALSE
  }
  
  if("Mnug" %in% names(out)){
    gammamodel <- TRUE 
  }
  else{
    gammamodel <- FALSE
  }
  
  
  if("sigma2" %in% names(out)){
    ln <- TRUE 
  }
  else{
    ln <- FALSE
  }
  
  
  if(ncol(out$MmuP)==4){
    Qmodel <- TRUE 
  }
  else{
    Qmodel<- FALSE
  }
    
  N <- nrow(out$Mtheta)
  K <- ncol(out$MAB[,,1])
  XG <- length(out$MAB[,1,1])
  bi <- round(.10*XG,0)  #burnin
  
  #######################################
  
  ##item parameter estimates
  
  idiscr <- apply(out$MAB[bi:XG,,1],2,mean)   
  idiff <- apply(out$MAB[bi:XG,,2],2,mean) 
  tdiscr <- apply(out$MAB[bi:XG,,3],2,mean)
  tintens <- apply(out$MAB[bi:XG,,4],2,mean)
    
  seidiscr <- round(sqrt(apply(out$MAB[bi:XG,,1],2,var)),3) # 3 = jumlah decimal
  seidiff <- round(sqrt(apply(out$MAB[bi:XG,,2],2,var)),3)
  setdiscr <- round(sqrt(apply(out$MAB[bi:XG,,3],2,var)),3)
  setintens <- round(sqrt(apply(out$MAB[bi:XG,,4],2,var)),3) 
  
  if(guess){    
  guess <- round(apply(out$Mguess,2,mean),3)
  seguess <- round(sqrt(apply(out$Mguess,2,var)),3)
  }
  
  
  #######################################
  ## ITEM POPULATION PARAMETER ESTIMATES
  
  pdiscr2 <- round(apply(out$MmuI[bi:XG,],2,mean),3)    
  sepdiscr2 <- round(sqrt(apply(out$MmuI[bi:XG,],2,var)),3)
  
  
  pSdiscr2 <- c(round(apply(out$MSI[bi:XG,,1],2,mean),3),round(apply(out$MSI[bi:XG,,2],2,mean),3),round(apply(out$MSI[bi:XG,,3],2,mean),3),round(apply(out$MSI[bi:XG,,4],2,mean),3))
  sepSdiscr2 <- c(round(sqrt(apply(out$MSI[bi:XG,,1],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,,2],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,,3],2,var)),3),round(sqrt(apply(out$MSI[bi:XG,,4],2,var)),3))
  
  ppers2 <- round(apply(out$MmuP[bi:XG,],2,mean),3)
  seppers2 <- round(sqrt(apply(out$MmuP[bi:XG,],2,var)),3)
  
    
  if(ln){
  msp1 <- round(apply(out$MSP[bi:XG,,1],2,mean),3)
  msp2 <- round(apply(out$MSP[bi:XG,,2],2,mean),3)
  semsp1 <- round(apply(out$MSP[bi:XG,,1],2,sd),3)
  semsp2 <- round(apply(out$MSP[bi:XG,,2],2,sd),3)
  
  mat_mean1<-matrix(c(msp1,msp2),2,2)
  rownames(mat_mean1)<- c("Theta","Zeta")
  colnames(mat_mean1)<- c("Theta","Zeta")
  
  mat_sd1<-matrix(c(semsp1,semsp2),2,2)
  rownames(mat_sd1)<- c("Theta","Zeta")
  colnames(mat_sd1)<- c("Theta","Zeta")
  
  Cor1<-cov2cor(mat_mean1)
  }
  
  
  
  if(Qmodel){
    
    L1<-round(apply(out$MSP[bi:XG,,1],2,mean),3)
    seL1<-round(apply(out$MSP[bi:XG,,1],2,sd),3)
    L2<-round(apply(out$MSP[bi:XG,,2],2,mean),3)
    seL2<-round(apply(out$MSP[bi:XG,,2],2,sd),3)
    L3<-round(apply(out$MSP[bi:XG,,3],2,mean),3)
    seL3<-round(apply(out$MSP[bi:XG,,3],2,sd),3)
    L4<-round(apply(out$MSP[bi:XG,,4],2,mean),3)
    seL4<-round(apply(out$MSP[bi:XG,,4],2,sd),3)
    
    if(simv){
      
      cov1<-matrix(c(L1,L2,L3,L4),4,4)
      rownames(cov1)<- c("Theta","Intercept","Slope1","Slope2")
      colnames(cov1)<- c("Theta","Intercept","Slope1","Slope2")
      
      sigR<-data$SigmaR  
      rownames(sigR)<- c("Theta","Intercept","Slope1","Slope2")
      colnames(sigR)<- c("Theta","Intercept","Slope1","Slope2")
      
      cov2<-matrix(c(seL1,seL2,seL3,seL4),4,4)
      rownames(cov2)<- c("Theta","Intercept","Slope1","Slope2")
      colnames(cov2)<- c("Theta","Intercept","Slope1","Slope2")
      
      Cor1<-cov2cor(cov1)
    }
    else{
      
      cov1<-matrix(c(L1,L2,L3,L4),4,4)
      rownames(cov1)<- c("Theta","Intercept","Slope1","Slope2")
      colnames(cov1)<- c("Theta","Intercept","Slope1","Slope2")
      
      cov2<-matrix(c(seL1,seL2,seL3,seL4),4,4)
      rownames(cov2)<- c("Theta","Intercept","Slope1","Slope2")
      colnames(cov2)<- c("Theta","Intercept","Slope1","Slope2")
      
      Cor1<-cov2cor(cov1)
    }
  }
  
  
  #######################################
  
  ## MEASUREMENT ERROR PARAMETER ESTIMATES
  
  ##Shape parameter Gamma
  if(gammamodel){
    estnug <- mean(out$Mnug[bi:XG])
    seestnug <- round(sqrt(var(out$Mnug[bi:XG])),3)
  }
  if(Qmodel){
    
    estsigma2 <- round(apply(out$Msigma2[bi:XG,],2,mean),3)
    seestsigma2 <- round(sqrt(apply(out$Msigma2[bi:XG,],2,var)),3)
    
  }
  if(ln){
    estsigma2 <- round(out$sigma2,3)
    seestsigma2 <- round(sqrt(out$sigma2),3)
  }
    
  
  #######################################
  
  ## PERSON FIT OUTCOMES
  
  if(gammamodel){
    cat("\n", "Gamma RT Modeling, 2013, J.P. Fox")
  }
  if(ln){
    cat("\n", "Log-Normal RT Modeling, 2013, J.P. Fox")
  }
  if(Qmodel){
    cat("\n", "LNIRT-Q Modeling, 2013, J.P. Fox")
  }
  cat("\n", "Summary of results")
  
  
  if(simv){
    cat("\n\n\t", "Item Discrimination parameter", "\t\t", "Item Difficulty parameter","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t", "Sim","\t\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
  }
  else{
    cat("\n\n\t", "Item Discrimination parameter", "\t\t", "Item Difficulty parameter","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t\t","item", "\t", "EAP", "\t", "SD","\n")
  }
  for(ii in 1:K){
    # Item Discrimination parameter
    cat("\n\t",ii,"\t") 
    printF(idiscr[ii], w=6, d= 3) # EAP
    cat("\t")
    printF(seidiscr[ii],w=6, d=3) # SD
    if(simv){
      cat("\t")
      printF(data$ab[ii,1],w=6, d=3) # SIM
      cat("\t\t",ii,"\t")
    }
    else{
      cat("\t\t",ii,"\t")
    }
    # Item Difficulty parameter
    printF(idiff[ii], w=6,d=3) # EAP
    cat("\t") 
    printF(seidiff[ii], w=6, d=3) # SD
    if(simv){
      cat("\t")
      printF(data$ab[ii,2],w=6, d=3) #SIM
      
    }
    else{
      cat("\t")
    }
  }
  
  
  if(simv){
    cat("\n\n\t","Time Discrimination","\t","Time Intensity","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t", "Sim","\t\t","item", "\t", "EAP", "\t", "SD","\t","Sim","\n")
  }
  else{
    cat("\n\n\t","Time Discrimination","\t","Time Intensity","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\t","item", "\t", "EAP", "\t", "SD","\n")
  }
  for(ii in 1:K){  
    cat("\n\t", ii,"\t") 
    # Time Discrimination
    printF(tdiscr[ii], w=6,d=3) # EAP
    cat("\t") 
    printF(setdiscr[ii], w=6, d=3) # SD
    if(simv){
      cat("\t")
      printF(data$ab[ii,3],w=6, d=3) #SIM      
    }
    else{
      cat("\t")
    }  
    cat("\t", ii,"\t")    # column "item" 
    # Time Intensity
    printF(tintens[ii], w=6,d=3) # EAP
    cat("\t") 
    printF(setintens[ii], w=6, d=3) # SD
    if(simv){
      cat("\t")
      printF(data$ab[ii,4],w=6, d=3) #SIM      
    }
    else{
      cat("\t")
    }  
    
  }  
    
  ## Measurement Error Variance
  
  if(sum(guess)>0){  
    cat("\n\n\t","Guessing Parameter","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\n")  
    for(ii in 1:K){
    cat("\n\t", ii,"\t") 
    # Guessing Parameter
    printF(guess[ii], w=6,d=3) # EAP
    cat("\t") 
    printF(seguess[ii], w=6, d=3) # SD
    cat("\t")       
    
  }
  }    
     
    cat("\n\n\t","Measurement Error Variance","\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD","\n")
    for(ii in 1:K){
      cat("\n\t", ii,"\t")     # column "item"  
    #cat("\t")
    printF(estsigma2[ii],w=6, d=3) # EAP
    cat("\t")
    printF(seestsigma2[ii],w=6, d=3) #SD       
    }
  
  
  if(round(pdiscr2[2],3)==0 && round(pdiscr2[4],3)==0){
    cat("\n\n\t", "Mean and Covariance matrix Items", "\n")
    cat("\n\t", "--- Population Mean Item ---", "\n")
    cat("\t","mu_a","\t", "SD","\t","mu_phi","\t","SD","\n\t")
    for(jj in c(1,3)){
      printF(pdiscr2[jj], w=6,d=3)
      cat("\t")
      printF(sepdiscr2[jj], w=6,d=3)
      cat("\t")
    }
  }
  else{
    cat("\n\n\t", "Mean and Covariance matrix Items", "\n")
    cat("\n\t", "--- Population Mean Item ---", "\n")
    cat("\t","mu_a","\t", "SD","\t","mu_b","\t","SD","\t","mu_phi","\t", "SD","\t","mu_lambda","\t","SD","\n\t")
    for(jj in c(1,2,3,4)){
      printF(pdiscr2[jj], w=6,d=3)
      cat("\t")
      printF(sepdiscr2[jj], w=6,d=3)
      cat("\t")
    }    
  } 
  
  
  cat("\n\n\t", "--- Covariance matrix Items ---", "\n")
  cat("\t\t\t\t\t\t\t", "a","\t\t\t\t\t\t\t","b","\t\t\t\t\t\t\t","phi","\t\t\t\t\t\t\t","lambda","\n\t")
    
  
  cat("a")
  for(ii in c(1,2,3,4)){
    printF(pSdiscr2[ii], w=10,d=3)
    cat("\t") 
  }  
  
  cat("\t","\n","b")
  cat("\t","      ","\t")
  for(ii in c(6,7,8)){   
    printF(pSdiscr2[ii], w=10,d=3)
    cat("\t") 
    
  }  
  
  cat("\t","\n","phi")
  cat("\t","      ","\t","      ","\t")
  for(ii in c(11,12)){
    printF(pSdiscr2[ii], w=10,d=3)
    cat("\t") 
    
  }
  
  cat("\t","\n","lambda")
  cat("\t","      ","\t","      ","\t","      ","\t")
  for(ii in c(16)){  
    printF(pSdiscr2[ii], w=9,d=3)
    cat("\t") 
    
  } 
  
  cat("\n\n\t", "--- Standard Error of Estimated Covariance ---", "\n")
  cat("\t\t\t\t\t\t\t", "a","\t\t\t\t\t\t\t","b","\t\t\t\t\t\t\t","phi","\t\t\t\t\t\t\t","lambda","\n\t")
    
  cat("a")
  for(ii in c(1,2,3,4)){
    printF(sepSdiscr2[ii], w=10,d=3)  
    cat("\t")   
  }
  
  cat("\t","\n","b")
  cat("\t","      ","\t")
  for(ii in c(6,7,8)){   
    printF(sepSdiscr2[ii], w=10,d=3)
    cat("\t")     
  }
    
  cat("\t","\n","phi")
  cat("\t","      ","\t","      ","\t")
  for(ii in c(11,12)){
    printF(sepSdiscr2[ii], w=10,d=3)
    cat("\t") 
  }
    
  cat("\t","\n","lambda")
  cat("\t","      ","\t","      ","\t","      ","\t")
  for(ii in c(16)){  
    printF(sepSdiscr2[ii], w=9,d=3)
    cat("\t\n") 
  }  
  
  cat("\n","Item Matrix Correlation","\n")
  print(cov2cor(matrix(pSdiscr2,4,4,byrow=TRUE)))
  
  cat("\n\n\t", "Mean and Covariance matrix Persons", "\n")
  cat("\n\t", "--- Population Mean Person ---", "\n")
  
  if(Qmodel){
    cat("\t","Ability", "\t", "SD","\t","Intercept","\t ","SD","\t","Slope1","\t ","SD","\t","Slope2","\t ","SD","\n\t")
    for(jj in c(1,2,3,4)){
      printF(ppers2[jj], w=6,d=3)
      cat("\t")
      printF(seppers2[jj], w=6,d=3)
      cat("\t")
    } 
  }
  else{
  cat("\t","mu_Y", "\t", "SD","\t","mu_zeta","\t ","SD","\t","\n\t")
  for(jj in c(1,2)){
    printF(ppers2[jj], w=6,d=3)
    cat("\t")
    printF(seppers2[jj], w=6,d=3)
    cat("\t")
  }
  }
    
  if(ln){
    cat("\n\n\t", "--- Covariance matrix Person ---", "\n")
    print(mat_mean1)
    cat("\n","Standard Error of Estimated Covariance","\n")
    print(mat_sd1)
    cat("\n","Person Matrix Correlation","\n")
    print(Cor1)
  }
  
  if(gammamodel){
    cat("\n\n\t", "--- Shape Parameter Gamma ---", "\n")
    cat("\t","EAP", "\t", "SD","\n\t")
    printF(estnug, w=6,d=3)
    cat("\t")
    printF(seestnug[jj], w=6,d=3)
    cat("\t")
  }
  
  
  if(Qmodel){
    cat("\n\n\t", "--- Covariance matrix Person ---", "\n")
        if(simv){      
      cat("\n\n\t","True Value")
      cat("\n","Covariance Matrix","\n")
      print(sigR)
      
      cat("\n","Estimated Value")
      cat("\n","Covariance Matrix","\n")
      print(cov1)
      
      cat("\n","Standard Error of Estimated Covariance","\n")
      print(cov2)
      
      cat("\n","Person Matrix Correlation","\n")
      print(Cor1)
    }
    else{
        
      cat("\n","Estimated Value")
      cat("\n","Covariance Matrix","\n")
      print(cov1)
      
      cat("\n","Standard Error of Estimated Covariance","\n")
      print(cov2)
      
      cat("\n","Person Matrix Correlation","\n")
      print(Cor1)
    }
  }
  
}



###################################











