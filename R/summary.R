#' @method summary LNIRT
#' @export
summary.LNIRT <- function(object, ...)
{
  if (class(object$data) == "simLNIRT") {
    simv <- TRUE
  } else {
    simv <- FALSE
  }
  if ("Mnug" %in% names(object)) {
    gammamodel <- TRUE
  } else {
    gammamodel <- FALSE
  }
  
  ## predictors speed and ability
  if(length(object$XPA)> 0){
    predictora <- TRUE
    ka <- ncol(object$XPA)
  } else {
    predictora <- FALSE
    ka <- NULL
  }
  if(length(object$XPT)> 0){
    predictort <- TRUE
    kt <- ncol(object$XPT)
  } else{
    predictort <- FALSE
    kt <- NULL
  }
  
  if(!predictora && !predictort){
    nopredictorp <- TRUE
  } else {
    nopredictorp <- FALSE
  }
  
  ## Item predictors difficulty and time intensity
  if(length(object$XIA)> 0){
    predictoria <- TRUE
    kia <- ncol(object$XIA)
  } else {
    predictoria <- FALSE
    kia <- NULL
  }
  if(length(object$XIT)> 0){
    predictorit <- TRUE
    kit <- ncol(object$XIT)
  } else{
    predictorit <- FALSE
    kit <- NULL
  }
  if(!predictoria && !predictorit){
    nopredictori <- TRUE
  } else {
    nopredictori <- FALSE
  }
  
  N <- length(object$Mtheta[, 1])
  K <- ncol(object$MAB[, , 1])
  XG <- length(object$MAB[, 1, 1])
  bi <- round(0.1 * XG, 0)  #burnin
  
  idiscr <- apply(object$MAB[bi:XG, , 1], 2, mean)
  idiff <- apply(object$MAB[bi:XG, , 2], 2, mean)
  tdiscr <- apply(object$MAB[bi:XG, , 3], 2, mean)
  tintens <- apply(object$MAB[bi:XG, , 4], 2, mean)
  
  seidiscr <- round(sqrt(apply(object$MAB[bi:XG, , 1], 2, var)), 3)
  seidiff <- round(sqrt(apply(object$MAB[bi:XG, , 2], 2, var)), 3)
  setdiscr <- round(sqrt(apply(object$MAB[bi:XG, , 3], 2, var)), 3)
  setintens <- round(sqrt(apply(object$MAB[bi:XG, , 4], 2, var)), 3)
  
  if (object$guess) {
    iguess <- round(apply(object$Mguess, 2, mean), 3)
    seiguess <- round(sqrt(apply(object$Mguess, 2, var)), 3)
  }
  else {
    iguess <- NULL
    seiguess <- NULL
  }
  
  pdiscr2 <- round(apply(object$MmuI[bi:XG, ], 2, mean), 3)
  sepdiscr2 <- round(sqrt(apply(object$MmuI[bi:XG, ], 2, var)), 3)
  
  pSdiscr2 <- matrix(c(round(apply(object$MSI[bi:XG, 1, ], 2, mean), 3), round(apply(object$MSI[bi:XG, 2, ], 2, mean), 3), round(apply(object$MSI[bi:XG, 
                                                                                                                                         3, ], 2, mean), 3), round(apply(object$MSI[bi:XG, 4, ], 2, mean), 3)), ncol = 4, nrow = 4, byrow = TRUE)
  sepSdiscr2 <- matrix(c(round(sqrt(apply(object$MSI[bi:XG, 1, ], 2, var)), 3), round(sqrt(apply(object$MSI[bi:XG, 2, ], 2, var)), 3), round(sqrt(apply(object$MSI[bi:XG, 
                                                                                                                                                          3, ], 2, var)), 3), round(sqrt(apply(object$MSI[bi:XG, 4, ], 2, var)), 3)), ncol = 4, nrow = 4, byrow = TRUE)
  sds <- sqrt(diag(pSdiscr2))
  SigmaIcor <- round(pSdiscr2/(sds %*% t(sds)), 3)
  
  ppers2 <- round(apply(object$MmuP[bi:XG, ], 2, mean), 3)
  seppers2 <- round(sqrt(apply(object$MmuP[bi:XG, ], 2, var)), 3)
  
  pSpers2 <- matrix(c(round(mean(object$MSP[bi:XG, 1, 1]), 3), round(mean(object$MSP[bi:XG, 2, 1]), 3), round(mean(object$MSP[bi:XG, 1, 2]), 3), round(mean(object$MSP[bi:XG, 
                                                                                                                                                           2, 2]), 3)), ncol = 2, nrow = 2)
  sepSpers2 <- matrix(c(round(sqrt(var(object$MSP[bi:XG, 1, 1])), 3), round(sqrt(var(object$MSP[bi:XG, 2, 1])), 3), round(sqrt(var(object$MSP[bi:XG, 
                                                                                                                                     1, 2])), 3), round(sqrt(var(object$MSP[bi:XG, 2, 2])), 3)), ncol = 2, nrow = 2)
  sds <- sqrt(diag(pSpers2))
  SigmaPcor <- round(pSpers2/(sds %*% t(sds)), 3)
  
  estsigma2 <- round(apply(object$Msigma2[bi:XG, ], 2, mean), 3)
  seestsigma2 <- round(sqrt(apply(object$Msigma2[bi:XG, ], 2, var)), 3)
  
  if (gammamodel) {
    estnug <- mean(object$Mnug[bi:XG])
    seestnug <- round(sqrt(var(object$Mnug[bi:XG])), 3)
  }
  else {
    estnug = NULL
    seestnug = NULL
  }
  
  out <- list(Mtheta = object$Mtheta, MTSD = object$MTSD, MAB = object$MAB, MmuP = object$MmuP, MSP = object$MSP, MmuI = object$MmuI, MSI = object$MSI, Mguess = object$Mguess, Msigma2 = object$Msigma2, 
              RT = object$RT, Y = object$Y, simv = simv, gammamodel = gammamodel, WL = object$WL, td = object$td, guess = object$guess, par1 = object$par1, N = N, K = K, XG = XG, bi = bi, idiscr = idiscr, idiff = idiff, tdiscr = tdiscr, tintens = tintens,
              seidiscr = seidiscr, seidiff = seidiff, setdiscr = setdiscr, setintens = setintens, iguess = iguess, seiguess = seiguess, pdiscr2 = pdiscr2, sepdiscr2 = sepdiscr2, pSdiscr2 = pSdiscr2, sepSdiscr2 = sepSdiscr2,
              SigmaIcor = SigmaIcor, ppers2 = ppers2, seppers2 = seppers2, pSpers2 = pSpers2, sepSpers2 = sepSpers2, SigmaPcor = SigmaPcor, estsigma2 = estsigma2, seestsigma2 = seestsigma2, estnug = estnug, seestnug = seestnug, 
              data = object$data, 
              nopredictorp = nopredictorp, nopredictori = nopredictori, predictora = predictora, predictort = predictort,
              predictoria = predictoria, predictorit = predictorit, ka = ka, kt = kt, kia = kia, kit = kit)
  if ("lZP" %in% names(object)) {
    tmp <- list(lZP = object$lZP, lZPT = object$lZPT, lZPA = object$lZPA, lZI = object$lZI, EAPresid = object$EAPresid, EAPresidA = object$EAPresidA, EAPKS = object$EAPKS, EAPKSA = object$EAPKSA, PFl = object$PFl, 
                PFlp = object$PFlp, IFl = object$IFl, IFlp = object$IFlp, EAPl0 = object$EAPl0, EAPCP1 = object$EAPCP1, EAPCP2 = object$EAPCP2, EAPCP3 = object$EAPCP3)
    out <- append(out, tmp)
  }
  
  class(out) <- "summary.LNIRT"
  return(out)
  
  
}

#' @method summary LNRT
#' @export
summary.LNRT <- function(object, ...)
{
  if (class(object$data) == "simLNIRT") {
    simv <- TRUE
  } else {
    simv <- FALSE
  }
  if ("Mnug" %in% names(object)) {
    gammamodel <- TRUE
  } else {
    gammamodel <- FALSE
  }
  
  ## predictors speed
  if(length(object$XPT)> 0){
    predictort <- TRUE
    nopredictorp <- FALSE
    kt <- ncol(object$XPT)
  } else{
    predictort <- FALSE
    nopredictorp <- TRUE
    kt <- NULL
  }
  
  ## Item predictors time intensity
  if(length(object$XIT)> 0){
    predictorit <- TRUE
    nopredictori <- FALSE
    kit <- ncol(object$XIT)
  } else{
    predictorit <- FALSE
    nopredictori <- TRUE
    kit <- NULL
  }

  N <- length(object$Mtheta)
  K <- ncol(object$MAB[, , 1])
  XG <- length(object$MAB[, 1, 1])
  bi <- round(0.1 * XG, 0)  #burnin
  
  ## item parameter estimates
  
  tdiscr <- apply(object$MAB[bi:XG, , 1], 2, mean)
  tintens <- apply(object$MAB[bi:XG, , 2], 2, mean)
  
  setdiscr <- round(sqrt(apply(object$MAB[bi:XG, , 1], 2, var)), 3)
  setintens <- round(sqrt(apply(object$MAB[bi:XG, , 2], 2, var)), 3)
  
  ## item population parameter estimates
  
  pdiscr <- round(apply(object$MmuI[bi:XG, ], 2, mean), 3)
  sepdiscr <- round(sqrt(apply(object$MmuI[bi:XG, ], 2, var)), 3)
  
  pSdiscr <- c(round(apply(object$MSI[bi:XG, , 1], 2, mean), 3), round(apply(object$MSI[bi:XG, , 2], 2, mean), 3))
  sepSdiscr <- c(round(sqrt(apply(object$MSI[bi:XG, , 1], 2, var)), 3), round(sqrt(apply(object$MSI[bi:XG, , 2], 2, var)), 3))
  
  ## person population parameter estimates
  
  if(!is.null(dim(object$MmuP[bi:XG, ])))
    ppers <- round(apply(object$MmuP[bi:XG, ], 2, mean), 3)
  else 
    ppers <- round(mean(object$MmuP[bi:XG, ]), 3)
  if (ncol(out$MmuP) == 1)
    seppers <- round(sqrt(var(object$MmuP[bi:XG, ])), 3)
  else 
    seppers <- round(sqrt(apply(object$MmuP[bi:XG, ], 2, var)), 3) 
  
  pSpers <- round(mean(object$MSP[bi:XG, 1, 1]), 3)
  sepSpers <- round(sqrt(var(object$MSP[bi:XG, 1, 1])), 3)
  
  
  ## Shape parameter Gamma single nu : dim(Mnug) : XG*1 multiple nu : dim(Mnug) : XG*K
  
  if (gammamodel) {
    if (ncol(object$Mnug) == 1) {
      estnug <- mean(object$Mnug[bi:XG])
      seestnug <- round(sqrt(var(object$Mnug[bi:XG])), 3)
    }
    if (ncol(object$Mnug) > 1) {
      estnug <- round(apply(object$Mnug[bi:XG, ], 2, mean), 3)
      seestnug <- round(sqrt(apply(object$Mnug[bi:XG, ], 2, var)), 3)
    }
    
    estsigma2 <- NULL
    seestsigma2 <- NULL
    
  } else {
    ## Measurement error parameter estimates
    
    estsigma2 <- round(apply(object$Msigma2[bi:XG, ], 2, mean), 3)
    seestsigma2 <- round(sqrt(apply(object$Msigma2[bi:XG, ], 2, var)), 3)
    
    estnug <- NULL
    seestnug <- NULL
  }
  
  
  out <- list(Mtheta = object$Mtheta, MTSD = object$MTSD, MAB = object$MAB, MmuP = object$MmuP, MSP = object$MSP, MmuI = object$MmuI, MSI = object$MSI, Msigma2 = object$Msigma2, 
              theta = object$theta, sigma2 = object$sigma2, RT = object$RT, simv = simv, gammamodel = gammamodel, WL = object$WL, Discrimination = object$Discrimination, N = N, K = K, XG = XG, bi = bi, 
              tdiscr = tdiscr, tintens = tintens, setdiscr = setdiscr, setintens = setintens, pdiscr = pdiscr, sepdiscr = sepdiscr, pSdiscr = pSdiscr, sepSdiscr = sepSdiscr,
              ppers = ppers, seppers = seppers, pSpers = pSpers, sepSpers = sepSpers, estsigma2 = estsigma2, seestsigma2 = seestsigma2, estnug = estnug, seestnug = seestnug, data = object$data,
              nopredictorp = nopredictorp, nopredictori = nopredictori, predictort = predictort, predictorit = predictorit, kt = kt, kit = kit)
  if ("lZP" %in% names(object)) {
    tmp <- list(lZP = object$lZP, lZPT = object$lZPT, lZI = object$lZI, EAPresid = object$EAPresid, EAPKS = object$EAPKS, EAPCP = object$EAPCP)
    out <- append(out, tmp)
  }
  
  class(out) <- "summary.LNRT"
  return(out)
  
  
}