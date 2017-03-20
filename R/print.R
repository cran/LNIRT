#' @method print summary.LNIRT
#' @export
print.summary.LNIRT <- function(x, ...)
{
  if (x$gammamodel) {
    cat("\n", "Gamma RT-IRT Modeling, 2013, J.P. Fox")
  } else {
    cat("\n", "Log-Normal RT-IRT Modeling, 2013, J.-P. Fox")
  }
  cat("\n", "Summary of results")
  
  if (x$simv) {
    cat("\n\n\t", "Item Discrimination parameter", "\t", "Item Difficulty parameter", "\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\n")
  } else {
    cat("\n\n\t", "Item Discrimination parameter", "\t", "Item Difficulty parameter", "\n")
    cat("\t", "item", "\t", "EAP", "\t", "SD", "\t\t", "item", "\t", "EAP", "\t", "SD", "\n")
  }
  for (ii in 1:x$K) {
    cat("\t")
    if (x$simv) {
      cat("\n\t", ii, "\t")
      printF(x$idiscr[ii], w = 6, d = 3)  # EAP
      cat("\t")
      printF(x$seidiscr[ii], w = 6, d = 3)  # SD
      cat("\t")
      printF(x$data$ab[ii, 1], w = 6, d = 3)  # SIM
      cat("\t", ii, "\t")
      printF(x$idiff[ii], w = 6, d = 3)
      cat("\t")
      printF(x$seidiff[ii], w = 6, d = 3)
      cat("\t")
      printF(x$data$ab[ii, 2], w = 6, d = 3)
    } else {
      cat("\n\t", ii, "\t")
      printF(x$idiscr[ii], w = 6, d = 3)  # EAP
      cat("\t")
      printF(x$seidiscr[ii], w = 6, d = 3)  # SD
      cat("\t\t", ii, "\t")
      printF(x$idiff[ii], w = 6, d = 3)
      cat("\t")
      printF(x$seidiff[ii], w = 6, d = 3)
    }
  }
  
  if (x$simv) {
    if (x$WL) {
      cat("\n\n\t", "Time Discrimination (Measurement Error)", " ", "Time Intensity", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\n")
    } else {
      cat("\n\n\t", "Time Discrimination", "\t\t", "Time Intensity", "\t\t", "Measurement Error Variance", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "EAP", "\t", "SD", 
          "\t", "Sim", "\n")
    }
  } else {
    if (x$WL) {
      cat("\n\n\t", "Time Discrimination (Measurement Error)", "\t", "Time Intensity", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\t\t\t", "item", "\t", "EAP", "\t", "SD", "\t", "\n")
    } else {
      cat("\n\n\t", "Time Discrimination", "\t", "Time Intensity", "\t", "Measurement Error Variance", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "EAP", "\t", "SD", "\n")
    }
  }
  for (ii in 1:x$K) {
    cat("\t")
    if (x$simv) {
      cat("\n\t", ii, "\t")
      if (x$WL) {
        printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
        cat("\t")
        printF(1/sqrt(x$data$sigma2[ii]), w = 6, d = 3)  #SIM      
      } else {
        printF(x$tdiscr[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$setdiscr[ii], w = 6, d = 3)  # SD
        cat("\t")
        printF(x$data$ab[ii, 3], w = 6, d = 3)  #SIM   
      }
      cat("\t", ii, "\t")
      printF(x$tintens[ii], w = 6, d = 3)  # EAP
      cat("\t")
      printF(x$setintens[ii], w = 6, d = 3)  # SD
      cat("\t")
      printF(x$data$ab[ii, 4], w = 6, d = 3)  # SIM      
      cat("\t")
      if (x$WL) {
        cat(" ")
      } else {
        printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
        cat("\t")
        printF(x$data$sigma2[ii], w = 6, d = 3)  #SIM      
      }
    } else {
      cat("\n\t", ii, "\t")
      if (x$WL) {
        printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
        cat("\t\t")
      } else {
        printF(x$tdiscr[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$setdiscr[ii], w = 6, d = 3)  # SD
      }
      cat("\t", ii, "\t")
      printF(x$tintens[ii], w = 6, d = 3)  # EAP
      cat("\t")
      printF(x$setintens[ii], w = 6, d = 3)  # SD
      cat("\t")
      if (x$WL) {
        cat("")
      } else {
        printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seestsigma2[ii], w = 6, d = 3)  #SD     
      }
    }
  }
  
  if (x$guess) {
    
    if (x$simv) {
      cat("\n\n\t", "Guessing Parameter", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\n")
      for (ii in 1:x$K) {
        cat("\n\t", ii, "\t")
        # Guessing Parameter
        printF(x$iguess[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seiguess[ii], w = 6, d = 3)  # SD
        cat("\t")
        printF(x$data$quess[ii], w = 6, d = 3)  # true value
      }
    } else {
      cat("\n\n\t", "Guessing Parameter", "\n")
      cat("\t", "item", "\t", "EAP", "\t", "SD", "\n")
      for (ii in 1:x$K) {
        cat("\n\t", ii, "\t")
        # Guessing Parameter
        printF(x$iguess[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$seiguess[ii], w = 6, d = 3)  # SD
      }
    }
  }
  
  if(x$nopredictori){
    if(round(x$pdiscr2[2],3)==0 && round(x$pdiscr2[4],3)==0 && (!x$WL)){
      cat("\n\n\t", "Mean and Covariance matrix Items (mu_a,mu_phi)", "\n")
      cat("\n\t", "--- Population Mean Item ---", "\n")
      cat("\t","mu_a","\t", "SD","\t","mu_phi"," ","SD","\n")
      for(jj in c(1,3)){
        cat("\t")
        printF(x$pdiscr2[jj], w=6,d=3)
        cat("\t")
        printF(x$sepdiscr2[jj], w=6,d=3)
      }
    }else{
      cat("\n\n\t", "Mean and Covariance matrix Items (mu_a,mu_b,mu_phi,mu_lambda)", "\n")
      cat("\n\t", "--- Population Mean Item ---","\n")
      cat("\t","mu_a","\t", "SD","\t","mu_b","\t","SD","\t","mu_phi","", "SD","\t","mu_lam","","SD","\n\t")
      for(jj in c(1,2,3,4)){
        printF(x$pdiscr2[jj], w=6,d=3)
        cat("\t")
        printF(x$sepdiscr2[jj], w=6,d=3)
        cat("\t")
      }    
    }
    if(x$WL){
      cat("\n\n\t", "--- Covariance matrix Items (a,b,error variance,lambda)---", "\n")
      cat("\t\t","SigmaI","\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
    }else{
      cat("\n\n\t", "--- Covariance matrix Items (a,b,phi,lambda)---", "\n")
      cat("\t\t","SigmaI","\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
    }
    for(ii in 1:4){
      cat("\t") 
      printF(x$pSdiscr2[ii,], w=6,d=3)
      cat("\t")
      if(max(x$sepSdiscr2)<1e02){
        printF(x$sepSdiscr2[ii,], w=6,d=3)
      }else{
        printF(x$sepSdiscr2[ii,], w=8,d=3)
      }
      cat("\t") 
      printF(x$SigmaIcor[ii,], w=6,d=3)
      cat("\t\n") 
    }  
  }else{
    if(round(x$pdiscr2[2],3)==0 && round(x$pdiscr2[x$kia+3],3)==0 && (!x$WL)){
      cat("\n\n\t", "Item Effects and Covariance matrix Items ", "\n")
      cat("\n\t", "--- Population Mean Item ---", "\n")
      cat("\t\t\t","EAP","\t","SD","\n")
      cat("\t mu_a","\t\t") 
      printF(x$pdiscr2[1], w=6,d=3)
      cat("\t")	
      printF(x$sepdiscr2[1], w=6,d=3)
      cat("\n\t","mu_phi","\t")
      printF(x$pdiscr2[x$kia+2], w=6,d=3)
      cat("\t")
      printF(x$sepdiscr2[x$kia+2], w=6,d=3)
      
      if(x$simv){
        if(x$predictoria) {
          cat("\n\t\t","Item Difficulty Predictor Effects","\n")
          for(ii in 2:(x$kia+1)){
            if(ii == 2){
              cat("\t\t\t","EAP", "\t", "SD","\t","Sim","\n")
            }
            if(ii == 2) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",ii-2,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t")
            if(ii > 2) {
              printF(x$data$Bia[ii-2], w=6,d=3)
              cat("\t\n")				
            }else{
              cat("\n")
            }	
          }	  
        }
        if(x$predictorit) {
          cat("\t\t","Time Intensity Predictor Effects","\n")
          sst <- 0	
          for(ii in (x$kia+3):(2+x$kit+x$kia)){
            sst<- sst+1	
            if(ii == (x$kia+3)) {
              cat("\t\t\t","EAP", "\t", "SD","\t","Sim","\n")
            }
            if(ii == (x$kia+3)) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",(x$kia-2)+sst,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t")
            if(ii > (x$kia+3)) {
              printF(x$data$Bit[sst-1], w=6,d=3)
              cat("\t\n")				
            }else{
              cat("\n")
            }				
          }
        }
      }else{
        if(x$predictoria) {
          cat("\n\t\t","Item Difficulty Predictor Effects","\n")
          for(ii in 2:(x$kia+1)){
            if(ii == 2){
              cat("\t\t\t","EAP", "\t", "SD","\n")
            }
            if(ii == 2) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",ii-2,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t\n") 
          }
        }
        if(x$predictorit) {
          cat("\t\t","Time Intensity Predictor Effects","\n")
          sst <- 0	
          for(ii in (x$kia+3):(2+x$kit+x$kia)){
            sst<- sst+1	
            if(ii == (x$kia+3)) {
              cat("\t\t\t","EAP", "\t", "SD","\n")
            }
            if(ii == (x$kia+3)) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",(x$kia-2)+sst,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t\n") 
          }
        }
      }	
    }else{
      cat("\n\n\t", "Item Effects and Covariance matrix Items ", "\n")
      cat("\n\t", "--- Population Mean Item ---", "\n")
      cat("\t\t\t","EAP","\t","SD","\n")
      cat("\t mu_a","\t\t") 
      printF(x$pdiscr2[1], w=6,d=3)
      cat("\t")	
      printF(x$sepdiscr2[1], w=6,d=3)
      cat("\n\t","mu_phi","\t")
      printF(x$pdiscr2[x$kia+2], w=6,d=3)
      cat("\t")
      printF(x$sepdiscr2[x$kia+2], w=6,d=3)
      if(x$simv){
        if(x$predictoria) {
          cat("\n\t\t","Item Difficulty Predictor Effects","\n")
          for(ii in 2:(x$kia+1)){
            if(ii == 2){
              cat("\t\t\t","EAP", "\t", "SD","\t","Sim","\n")
            }
            if(ii == 2) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",ii-2,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t")
            if(ii > 2) {
              printF(x$data$Bia[ii-2], w=6,d=3)
              cat("\t\n")				
            }else{
              cat("\n")
            }
          }
        }
        if(x$predictorit) {
          sst <- 0	
          cat("\t\t","Time Intensity Predictor Effects","\n")
          for(ii in (x$kia+3):(2+x$kit+x$kia)){
            sst <- sst+1
            if(ii == (x$kia+3)) {
              cat("\t\t\t","EAP", "\t", "SD","\t","Sim","\n")
            }
            if(ii == (x$kia+3)) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",(x$kia-2)+sst,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t")
            if(ii > (x$kia+3)) {
              printF(x$data$Bit[sst-1], w=6,d=3)
              cat("\t\n")				
            }else{
              cat("\n")
            } 
          }
        }
      }else{
        if(x$predictoria) {
          cat("\n\t\t","Item Difficulty Predictor Effects","\n")
          for(ii in 2:(x$kia+1)){
            if(ii == 2){
              cat("\t\t\t","EAP", "\t", "SD","\n")
            }
            if(ii == 2) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",ii-2,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t\n") 
          }
        }
        if(x$predictorit) {
          sst <- 0	
          cat("\t\t","Time Intensity Predictor Effects","\n")
          for(ii in (x$kia+3):(2+x$kit+x$kia)){
            sst <- sst+1
            if(ii == (x$kia+3)) {
              cat("\t\t\t","EAP", "\t", "SD","\n")
            }
            if(ii == (x$kia+3)) {
              cat("\t Intercept","\t") 
            }else{	
              cat("\t X",(x$kia-2)+sst,"\t\t") 
            }
            printF(x$pdiscr2[ii], w=6,d=3)
            cat("\t") 
            printF(x$sepdiscr2[ii], w=6,d=3)
            cat("\t\n") 
          }
        }
      }	
    }
    if(x$WL){
      cat("\n\n\t", "--- Covariance matrix Items (a,b,error variance,lambda)---", "\n")
      cat("\t","SigmaI","\t\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
    }else{
      cat("\n\n\t", "--- Covariance matrix Items (a,b,phi,lambda)---", "\n")
      cat("\t","SigmaI","\t\t\t", "SD SigmaI","\t\t\t","SigmaI (Correlation)","\n")
    }
    for(ii in 1:4){
      cat("\t") 
      printF(x$pSdiscr2[ii,], w=6,d=3)
      cat("\t")
      if(max(x$sepSdiscr2)<1e02){
        printF(x$sepSdiscr2[ii,], w=6,d=3)
      }else{
        printF(x$sepSdiscr2[ii,], w=8,d=3)
      }
      cat("\t") 
      printF(x$SigmaIcor[ii,], w=6,d=3)
      cat("\t\n") 
    }
  }
  
  
  if(x$nopredictorp){
    cat("\n\n\t", "Mean and Covariance matrix Persons (ability,speed)", "\n")
    cat("\n\t", "--- Population Mean Person (Ability - Speed)---", "\n")
    cat("\t\t\t","muP", "\t", "SD","\n")
    for(ii in 1:2){
      if(ii == 1) {
        cat("\t Ability \t") 
      }else{	
        cat("\t Speed \t\t") 
      }
      printF(x$ppers2[ii], w=6,d=3)
      cat("\t") 
      printF(x$seppers2[ii], w=6,d=3)
      cat("\t\n") 
    }  
  } else {
    if(x$simv){
      cat("\n\n\t", "Person Effects and Covariance matrix Persons (ability,speed)", "\n")
      cat("\n\t", "--- Person Effects (Ability - Speed)---", "\n")
      
      if(x$predictora) {
        cat("\t\t","Ability Predictor Effects","\n")
        for(ii in 1:x$ka){
          if(ii == 1) {
            cat("\t\t","EAP", "\t", "SD","\t", "Sim","\n")
          }
          cat("\t X",ii,"\t") 
  
          printF(x$ppers2[ii], w=6,d=3)
          cat("\t") 
          printF(x$seppers2[ii], w=6,d=3)
          cat("\t") 
          if(!is.null(x$data$Ba[ii])) {
            printF(x$data$Ba[ii], w=6,d=3)
            cat("\t\n") 
          }
          else {
            cat("\n") 
          } 
        } 
      }
      
      if(x$predictort) {
        cat("\t\t","Speed Predictor Effects","\n")
        for(ii in (x$ka+1):(x$kt+x$ka)){
          if(ii == (x$ka+1)) {
            cat("\t\t","EAP", "\t", "SD","\t", "Sim","\n")
          }
            
          cat("\t X",ii,"\t") 
          
          
          #cat("\t X",ii,"\t") 
          printF(x$ppers2[ii], w=6,d=3)
          cat("\t") 
          printF(x$seppers2[ii], w=6,d=3)
          cat("\t")
          if(!is.null(x$data$Bt[ii-x$ka])) {
            printF(x$data$Bt[ii-x$ka], w=6,d=3)
            cat("\t\n") 
          }
          else {
            cat("\n") 
          }
          
        }
      }
    } else {
      cat("\n\n\t", "Person Effects and Covariance matrix Persons (ability,speed)", "\n")
      cat("\n\t", "--- Person Effects (Ability - Speed)---", "\n")
      
      if(x$predictora) {
        cat("\t\t","Ability Predictor Effects","\n")
        for(ii in 1:x$ka){
          if(ii == 1) {
            cat("\t\t","EAP", "\t", "SD","\n")
          }
          cat("\t X",ii,"\t") 
          printF(x$ppers2[ii], w=6,d=3)
          cat("\t") 
          printF(x$seppers2[ii], w=6,d=3)
          cat("\t\n") 
        }  
      }
      
      if(x$predictort) {
        cat("\t\t","Speed Predictor Effects","\n")
        for(ii in (x$ka+1):(x$kt+x$ka)){
          if(ii == (x$ka+1)) {
            cat("\t\t","EAP", "\t", "SD","\n")
          }
          cat("\t X",ii,"\t") 
          printF(x$ppers2[ii], w=6,d=3)
          cat("\t") 
          printF(x$seppers2[ii], w=6,d=3)
          cat("\t\n") 
        }
      }
    }
  }
  
  cat("\n\t", "SigmaP", "\t", "SD SigmaP", "\t", "SigmaP (Correlation)", "\n")
  for (ii in 1:2) {
    cat("\t")
    printF(x$pSpers2[ii, ], w = 6, d = 3)
    cat("\t")
    printF(x$sepSpers2[ii, ], w = 6, d = 3)
    cat("\t")
    printF(x$SigmaPcor[ii, ], w = 6, d = 3)
    cat("\t\n")
  }
  
  if (x$gammamodel) {
    cat("\n\n\t", "--- Shape Parameter Gamma ---", "\n")
    cat("\t", "EAP", "\t", "SD", "\n\t")
    printF(x$estnug, w = 6, d = 3)
    cat("\t")
    printF(x$seestnug[jj], w = 6, d = 3)
    cat("\t")
  }
  
  if ("lZP" %in% names(x)) {
    
    cat("\n\n\n\t", "*** \t\t\t\t\t\t ***", "\n")
    cat("\t", "*** Person Fit Analysis (Log-Normal Speed) ***", "\n")
    cat("\t", "*** \t\t\t\t\t\t ***", "\n")
    
    
    # Percentage Aberrant lZP (chi-square statstic)
    cat("\n\t", "Percentage Outliers Persons (5% level)", "\n")
    cat("\n\t", "lZ", "\n")
    cat("\t", round(100 * length(which(x$lZP < 0.05))/x$N, 2), "%", "\t\n")
    cat("\t 95% Posterior Probability: ", round(100 * sum(x$EAPCP1 > 0.95)/x$N, 2), "%", "\t\n")
    
    # Percentage Persons (Estimation Group, Aberrant Group)
    
    # cat('\n\t', 'Percentage Persons (Estimation Group, Aberrant Group) (5% level)','\n') cat('\n\t', 'Estimation Group', '\t',
    # 'Aberrant','\n') cat('\t',100*round(apply(out$Mingroup,2,mean)[1],3),'\t\t\t',100*round(apply(out$Mingroup,2,mean)[2],3),'\t\n')
    
    
    cat("\n\n\t", "*** Item Fit Analysis ***", "\n")
    
    cat("\n\t", "Misfitting Items (5% level)", "\n")
    set <- which(x$lZI < 0.05)
    if (length(set >= 1)) {
      cat("\n\t", "lI", "\n")
      cat("\t", set, "\t\n")
    } else {
      cat("\t", "No Misfitting Items", "\t\n")
    }
    
    cat("\n\t", "*** Residual Analysis ***", "\n")
    set <- which(x$EAPresid > 0.95)
    if (length(set >= 1)) {
      cat("\n\t", "Percentage Extreme Residuals (.95 Posterior Probability)", "\n")
      cat("\t", round(100 * length(set)/x$N, 4), "%", "\t\n")
      
      ## Identify Extremes
      set <- which(x$EAPresid > 0.95)
      set <- cbind((set%%x$N), (floor(set/x$N) + 1), exp(x$RT[set]))
      dd <- order(set[, 1])
      colnames(set) <- c("Person", "Item", " RT ")
      rownames(set) <- rep(c(""), nrow(set))
      
      cat("\n\t", "Extreme Residuals", "\n")
      cat("\t", "Person", " Item", "\t", " RT ", "\n")
      for (jj in dd) {
        cat("\t")
        printF(set[jj, 1], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 2], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 3], w = 8, d = 4)
        cat("\n")
      }
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKS < 0.05)))/x$K, "%", "of items has non-lognormally distributed residuals", "\t\n")
      set <- which(x$EAPKS < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKS[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    } else {
      cat("\t", "No Extreme Residuals", "\t\n")
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKS < 0.05)))/x$K, "%", "\t\n")
      set <- which(x$EAPKS < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKS[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    }
    
    cat("\n\n\n\t", "*** \t\t\t\t\t\t ***", "\n")
    cat("\t", "*** Person Fit Analysis (IRT Model For Ability) ***", "\n")
    cat("\t", "*** \t\t\t\t\t\t ***", "\n")
    
    # Percentage Aberrant lZP (chi-square statstic) Percentage Aberrant l0 (log-likelihood statstic)
    cat("\n\t", "Percentage Outliers Persons (5% level)", "\n")
    cat("\n\t", "Log-likelihood Statistic", "\n")
    # cat('\t',round(100*length(which(out$lZPA < .05))/x$N,2),'%','\t\n')
    cat("\t", round(100 * length(which(x$PFlp < 0.05))/x$N, 2), "%", "\t\n")
    cat("\t 95% Posterior Probability: ", round(100 * sum(x$EAPCP2 > 0.95)/x$N, 2), "%", "\t\n")
    cat("\t 95% Posterior Probability (Ability and Speed): ", round(100 * sum(x$EAPCP3 > 0.95)/x$N, 2), "%", "\t\n")
    
    
    cat("\n\n\t", "*** Item Fit Analysis ***", "\n")
    
    cat("\n\t", "Misfitting Items (5% level)", "\n")
    # set <- which(out$lZIA < .05)
    set <- which(x$IFlp < 0.05)
    if (length(set >= 1)) {
      cat("\n\t", "lI", "\n")
      cat("\t", set, "\t\n")
    } else {
      cat("\t", "No Misfitting Items", "\t\n")
    }
    
    cat("\n\t", "*** Residual Analysis ***", "\n")
    set <- which(x$EAPresidA > 0.95)
    if (length(set >= 1)) {
      cat("\n\t", "Percentage Extreme Residuals (.95 Posterior Probability)", "\n")
      cat("\t", round(100 * length(set)/x$N, 4), "%", "\t\n")
      
      ## Identify Extremes
      set <- which(x$EAPresidA > 0.95)
      set <- cbind((set%%x$N), (floor(set/x$N) + 1), x$Y[set], x$Mtheta[(set%%x$N)])
      dd <- order(set[, 1])
      colnames(set) <- c("Person", "Item", " Y ", " EAP theta")
      rownames(set) <- rep(c(""), nrow(set))
      
      cat("\n\t", "Extreme Residuals", "\n")
      cat("\t", "Person", " Item", "\t", " Response ", " EAP Theta", "\n")
      for (jj in dd) {
        cat("\t")
        printF(set[jj, 1], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 2], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 3], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 4], w = 8, d = 4)
        cat("\n")
      }
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKSA < 0.05)))/x$K, "%", "of items has non-normally distributed latent residuals", "\t\n")
      set <- which(x$EAPKSA < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKSA[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    } else {
      cat("\t", "No Extreme Residuals", "\t\n")
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKSA < 0.05)))/x$K, "%", "\t\n")
      set <- which(x$EAPKSA < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKSA[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    }
  }  ## close personfit report
  cat("\n\n")
}

#' @method print summary.LNRT
#' @export
print.summary.LNRT <- function(x, ...)
{
  ## Person Fit outcomes
  
  if (x$gammamodel) {
    cat("\n", "Gamma RT Modeling, 2013, J.P. Fox")
  } else {
    cat("\n", "Log-Normal RT Modeling, 2013, J.P. Fox")
  }
  cat("\n", "Summary of results")
  
  
  if (x$gammamodel) {
    if (x$simv) {
      if (ncol(x$Mnug) == 1) {
        cat("\n\n\t", "Time Intensity parameter", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\n")
      }
      if (ncol(x$Mnug) > 1) {
        cat("\n\n\t", "Time Intensity parameter", "\t", "Shape Parameter Gamma", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\n")
      }
    } else {
      
      if (ncol(x$Mnug) == 1) {
        cat("\n\n\t", "Time Intensity parameter", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\n")
      }
      
      if (ncol(x$Mnug) > 1) {
        cat("\n\n\t", "Time Intensity parameter", "\t", "Shape Parameter Gamma", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t\t", "item", "\t", "EAP", "\t", "SD", "\n")
      }
      
    }
  } else {
    if (x$simv) {
      if (x$WL) {
        cat("\n\n\t", "Time Discrimination (Measurement Error)", " ", "Time Intensity", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\n")
      } else {
        cat("\n\n\t", "Time Discrimination", "\t\t", "Time Intensity", "\t\t", "Measurement Error Variance", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "Sim", "\t", "EAP", "\t", 
            "SD", "\t", "Sim", "\n")
      }
    } else {
      if (x$WL) {
        cat("\n\n\t", "Time Discrimination (Measurement Error)", "\t", "Time Intensity", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t\t\t", "item", "\t", "EAP", "\t", "SD", "\t", "\n")
      } else {
        cat("\n\n\t", "Time Discrimination", "\t", "Time Intensity", "\t", "Measurement Error Variance", "\n")
        cat("\t", "item", "\t", "EAP", "\t", "SD", "\t", "item", "\t", "EAP", "\t", "SD", "\t", "EAP", "\t", "SD", "\n")
      }
    }
  }
  
  if (x$gammamodel) {
    
    for (ii in 1:x$K) {
      cat("\n\t", ii, "\t")
      printF(x$tintens[ii], w = 6, d = 3)  # estimated bsp  
      cat("\t")
      printF(x$setintens[ii], w = 6, d = 3)
      if (x$simv) {
        cat("\t")
        printF(x$data$bsp[ii], w = 6, d = 3)  # bsp true value
        # cat('\n\t', ii,'\t')
        if (ncol(x$Mnug) > 1) {
          cat("\t")
          cat("\t", ii, "\t")
          printF(x$estnug[ii], w = 6, d = 3)  # Shape Parameter Gamma EAP
          cat("\t")
          printF(x$seestnug[ii], w = 6, d = 3)
        }
        
      } else {
        if (ncol(x$Mnug) == 1) {
          cat("\t")
        }
        
        if (ncol(x$Mnug) > 1) {
          cat("\t\t", ii, "\t")
          printF(x$estnug[ii], w = 6, d = 3)
          cat("\t")
          printF(x$seestnug[ii], w = 6, d = 3)
        }
      }
    }
  } else {
    for (ii in 1:x$K) {
      cat("\t")
      if (x$simv) {
        cat("\n\t", ii, "\t")
        if (x$WL) {
          printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
          cat("\t")
          printF(1/sqrt(x$data$sigma2[ii]), w = 6, d = 3)  #SIM      
        } else {
          printF(x$tdiscr[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$setdiscr[ii], w = 6, d = 3)  # SD
          cat("\t")
          printF(x$data$ab[ii, 3], w = 6, d = 3)  #SIM   
        }
        cat("\t", ii, "\t")
        printF(x$tintens[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$setintens[ii], w = 6, d = 3)  # SD
        cat("\t")
        printF(x$data$ab[ii, 4], w = 6, d = 3)  # SIM      
        cat("\t")
        if (x$WL) {
          cat(" ")
        } else {
          printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
          cat("\t")
          printF(x$data$sigma2[ii], w = 6, d = 3)  #SIM      
        }
      } else {
        cat("\n\t", ii, "\t")
        if (x$WL) {
          printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$seestsigma2[ii], w = 6, d = 3)  # SD     
          cat("\t\t")
        } else {
          printF(x$tdiscr[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$setdiscr[ii], w = 6, d = 3)  # SD
        }
        cat("\t", ii, "\t")
        printF(x$tintens[ii], w = 6, d = 3)  # EAP
        cat("\t")
        printF(x$setintens[ii], w = 6, d = 3)  # SD
        cat("\t")
        if (x$WL) {
          cat("")
        } else {
          printF(x$estsigma2[ii], w = 6, d = 3)  # EAP
          cat("\t")
          printF(x$seestsigma2[ii], w = 6, d = 3)  #SD     
        }
      }
    }
    
  }
  
  
  cat("\n\n\t", "Mean and Covariance matrix Items (phi,lambda)", "\n")
  cat("\n\t", "--- Population Mean Item ---", "\n")
  cat("\t", "mu_phi", " ", "SD", "\t", "mu_lam", " ", "SD", "\n\t")
  for (jj in c(1, 2)) {
    printF(x$pdiscr[jj], w = 6, d = 3)
    cat("\t ")
    printF(x$sepdiscr[jj], w = 6, d = 3)
    cat("\t")
  }
  cat("\n\n\t", "--- Covariance matrix Items ---", "\n")
  cat("\t", "phi", "\t", "SD", "\t", "Cov", "\t", "SD", "\t", "lambda", " ", "SD", "\n\t")
  for (jj in c(1, 2, 4)) {
    printF(x$pSdiscr[jj], w = 6, d = 3)
    cat("\t")
    printF(x$sepSdiscr[jj], w = 6, d = 3)
    cat("\t")
  }
  
  
  cat("\n\n\t", "Mean and Covariance matrix Persons", "\n")
  cat("\n\t", "--- Population Mean Person ---", "\n")
  cat("\t", "mu_P", "\t ", "SD", "\n\t")
  for (jj in c(1)) {
    printF(x$ppers[jj], w = 6, d = 2)
    cat("\t")
    printF(x$seppers[jj], w = 6, d = 2)
    cat("\t")
  }
  cat("\n\n\t", "--- Covariance matrix Person ---", "\n")
  cat("\t", "Sigma_P", "SD", "\n\t")
  for (jj in c(1)) {
    printF(x$pSpers[jj], w = 6, d = 3)
    cat("\t")
    printF(x$sepSpers[jj], w = 6, d = 3)
    cat("\t")
  }
  
  if ((x$gammamodel) && (ncol(x$Mnug) == 1)) {
    cat("\n\n\t", "--- Shape Parameter Gamma ---", "\n")
    cat("\t", "EAP", "\t", "SD", "\n\t")
    printF(x$estnug, w = 6, d = 3)
    cat("\t")
    printF(x$seestnug[jj], w = 6, d = 3)
    cat("\t")
  }
  
  
  if ("lZP" %in% names(x)) {
    
    cat("\n\n\n\t", "*** Person Fit Analysis ***", "\n")
    
    # Percentage Aberrant lZP3
    cat("\n\t", "Percentage Outliers Persons (5% level)", "\n")
    cat("\n\t", "lZ", "\n")
    cat("\t", round(100 * sum(x$lZP < 0.05)/x$N, 2), "%", "\t\n")
    cat("\t 95% Posterior Probability: ", round(100 * sum(x$EAPCP > 0.95)/x$N, 2), "%", "\t\n")
    
    
    # Percentage Persons (Estimation Group, Aberrant Group)
    
    # cat('\n\t', 'Percentage Persons (Estimation Group, Aberrant Group) (5% level)','\n') cat('\n\t', 'Estimation Group', '\t',
    # 'Aberrant','\n') cat('\t',100*round(apply(x$Mingroup,2,mean)[1],3),'\t\t\t',100*round(apply(x$Mingroup,2,mean)[2],3),'\t\n')
    
    
    cat("\n\n\t", "*** Item Fit Analysis ***", "\n")
    
    cat("\n\t", "Misfitting Items (5% level)", "\n")
    set <- which(x$lZI < 0.05)
    if (length(set >= 1)) {
      cat("\n\t", "lI", "\n")
      cat("\t", set, "\t\n")
    } else {
      cat("\t", "No Misfitting Items", "\t\n")
    }
    
    cat("\n\t", "*** Residual Analysis ***", "\n")
    set <- which(x$EAPresid > 0.95)
    if (length(set >= 1)) {
      cat("\n\t", "Percentage Extreme Residuals (.95 Posterior Probability)", "\n")
      cat("\t", round(100 * length(set)/x$N, 4), "%", "\t\n")
      
      ## Identify Extremes
      set <- which(x$EAPresid > 0.95)
      set <- cbind((set%%x$N), (floor(set/x$N) + 1), exp(x$RT[set]))
      dd <- order(set[, 1])
      colnames(set) <- c("Person", "Item", " RT ")
      rownames(set) <- rep(c(""), nrow(set))
      
      cat("\n\t", "Extreme Residuals", "\n")
      cat("\t", "Person", " Item", "\t", " RT ", "\n")
      for (jj in dd) {
        cat("\t")
        printF(set[jj, 1], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 2], w = 6, d = 0)
        cat("\t")
        printF(set[jj, 3], w = 8, d = 4)
        cat("\n")
      }
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKS < 0.05)))/x$K, "%", "of items has non-lognormally distributed residuals", "\t\n")
      set <- which(x$EAPKS < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKS[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    } else {
      cat("\t", "No Extreme Residuals", "\t\n")
      cat("\n\t", "Kolmogorov Smirnov Test (5% level)", "\n")
      cat("\t", round(100 * length(which(x$EAPKS < 0.05)))/x$K, "%", "\t\n")
      set <- which(x$EAPKS < 0.05)
      if (length(set) > 0) {
        cat("\t", "Item", "\t P-value ", "\n")
        for (jj in 1:length(set)) {
          cat("\t")
          printF(set[jj], w = 6, d = 0)
          cat("\t")
          printF(x$EAPKS[set[jj]], w = 6, d = 3)
          cat("\n")
        }
      }
    }
  }
  cat("\n\n")
   
}