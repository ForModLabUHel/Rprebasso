####init Biomass
initBiomasses <- function(pCro,initVarX){
  initVarX<-as.matrix(initVarX) #change vector to matrix when maxlayer=1
  siteType <- initVarX[8,1]
  ##set parameters
  par_betab <- pCro[13,initVarX[1,]]
  par_x <- pCro[19,initVarX[1,]]
  par_beta0 <- pCro[12,initVarX[1,]]
  par_betas <- pCro[14,initVarX[1,]]
  par_mf <- pCro[8,initVarX[1,]]
  par_mr <- pCro[9,initVarX[1,]]
  par_mw <- pCro[10,initVarX[1,]]
  par_alfar <- pCro[20+pmin(siteType,5),initVarX[1,]]
  par_c <- pCro[7,initVarX[1,]]
  par_rhof <- pCro[15,initVarX[1,]]
  par_rhor <- par_alfar * par_rhof
  par_rhow <- pCro[2,initVarX[1,]]
  par_S_branchMod <- pCro[27,initVarX[1,]]
  gammaC <- 0. #initVarX[8,]
  Tbd <- 10 #####to include in the parameters
  
  ###set variables
  A <- initVarX[7,]
  ba <- initVarX[5,]; d <- initVarX[4,]
  N <- ba/(pi*((d/2/100)^2))
  h = initVarX[3,]; hc <- initVarX[6,]
  B = ba/N
  Lc <- h - hc
  betab =  par_betab * Lc^(par_x-1)
  beta0 = par_beta0
  beta1 = (beta0 + betab + par_betas) 
  beta2 = 1. - betab - par_betas 		
  betaC = (beta1 + gammaC * beta2) / par_betas
  wf_STKG <- max(0.,par_rhof * A * N)
  W_froot = max(0.,par_rhor * A * N)  ##to check  ##newX
  W_wsap = max(0.,par_rhow * A * N * (beta1 * h + beta2 * hc)) ##newX
  W_c = max(0.,par_rhow * A * N * hc) #sapwood stem below Crown
  W_s = max(0.,par_rhow * A * N * par_betas * Lc) #sapwood stem within crown
  W_branch =  max(0.,par_rhow * A * N * betab * Lc) #branches biomass
  W_croot = max(0.,par_rhow * beta0 * A * h * N) #W_stem * (beta0 - 1.)	#coarse root biomass
  Wsh = pmax((A+B+sqrt(A*B)) * hc * par_rhow * N/2.9 - W_c,0) #initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
  #initialize Wdb dead branches biomass
  Wdb = pmax(0,ifelse(par_S_branchMod == 1.,Tbd * W_branch * ((0.0337+0.000009749*N)*exp(-0.00456*d^2)+0.00723),
         Tbd * W_branch *((-0.00513+0.000012*N)*exp((0.00000732-0.000000764*N)*d^2)+0.00467)))
  W_stem = max(0.,W_c + W_s + Wsh)
  V = W_stem / par_rhow
  biomasses <- rbind(wf_STKG,W_froot,W_wsap,W_c,W_s,W_branch,W_croot,Wsh,Wdb,W_stem,V)
  return(biomasses)
}
  

#### function to calculate initial sapwood area at crown base (A)
compA <- function(inputs){
  p_ksi = inputs[1]
  p_rhof = inputs[2]
  p_z <- inputs[3]
  Lc = inputs[4]
  A <- max(0.,p_ksi/p_rhof * Lc^p_z)
  return(A)
}


# Function to Compute Hc based on ksi parameter
ksiHcMod <- function(initVar){
  h <- initVar[3] - 1.3
  b <- pi * initVar[4]^2 / 40000
  p_ksi <- pCROBAS[38,initVar[1]]
  p_rhof <- pCROBAS[15,initVar[1]]
  p_z <- pCROBAS[11,initVar[1]]
  Lc <- h* ((p_rhof * b)/(p_ksi * h^p_z))^(1/(p_z-1))
  Hc <- max(0.,(h-Lc))
  return(Hc)
}

###function to replace HC NAs in initial variable initVar
findHcNAs <- function(initVar,pHcMod,HcModV){
  if(is.vector(initVar)){
    if(is.na(initVar[6])){
      if(HcModV==1){
        initVar[6] <- ksiHcMod(initVar)
      }else if(HcModV==2){
        inModHc <- c(pHcMod[,initVar[1]],initVar[3],
                     initVar[4],initVar[2],initVar[5],initVar[5])
        initVar[6] <- model.Hc(inModHc)
      }
    }
  } else if(any(is.na(initVar[6,]))){
    initVar[1,][which(initVar[1,]==0)] <- 1 ###deals with 0 species ID
    HcNAs <- which(is.na(initVar[6,]))
    BAtot <- sum(initVar[5,],na.rm = T)
    if(length(HcNAs)==1){
      if(HcModV==1){
        initVar[6,HcNAs] <- ksiHcMod(initVar[,HcNAs])
      }else if(HcModV==2){
        inModHc <- c(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                     initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],BAtot)
        initVar[6,HcNAs] <- model.Hc(inModHc)
      }
    }else{
      if(HcModV==1){
        initVar[6,HcNAs] <- apply(initVar,2,ksiHcMod)
      }else if(HcModV==2){
        inModHc <- rbind(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                         initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],BAtot)
        initVar[6,HcNAs] <- apply(inModHc,2,model.Hc)
      }
    }
  }
  return(initVar)
}


##Height of the crown base model
model.Hc <- function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
  lnHc_sim <- pValues[1]+pValues[2]*log(H)+pValues[3]*D/H+
    pValues[4]*log(age)+ pValues[5]*log(BA_sp)+
    pValues[6]*(BA_sp/BA_tot)
  Hc_sim <- exp(lnHc_sim)
  return(pmax(Hc_sim,0.,na.rm = T)) 
} 
varNames  <- c('siteID','gammaC','sitetype','species','ETS' ,'P0','age', 'DeadWoodVolume', 'Respi_tot','GPP/1000',
               'H','D', 'BA','Hc_base','Cw','Ac','N','npp','leff','keff','lproj','ET_preles','weight',
               'Wbranch',"WfineRoots",'Litter_fol','Litter_fr','Litter_branch','Litter_wood','V',
               'Wstem','W_croot','wf_STKG', 'wf_treeKG','B_tree','Light',"Vharvested","Wharvested","soilC",
               "aSW","summerSW","Vmort","gross growth", "GPPspecies","Rh species", "NEP sp"," W_wsap","W_c","W_s","Wsh","Wdb","dHc",
               "Wbh","Wcrh","dH")

  getVarNam <- function(){
    return(varNames)
}


  aTmean <- function(TAir,nYears){
    Tmean = colMeans(matrix(TAir,365,nYears))
    return(Tmean)
  }

  aTampl <- function(TAir,nYears){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    TbyYear <- matrix(TAir,365,nYears)
    Tampl = apply(TbyYear, 2, function(x) max(aggregate(x/2,by=list(monthsDays),FUN=mean)) - min(aggregate(x/2,by=list(monthsDays),FUN=mean))  )
    return(Tampl)
  }

  aPrecip <- function(Precip,nYears){
    aP = colSums(matrix(Precip,365,nYears))
    return(aP)
  }
