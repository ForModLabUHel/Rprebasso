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
  wf_STKG <- pmax(0.,par_rhof * A * N)
  W_froot = pmax(0.,par_rhor * A * N)  ##to check  ##newX
  W_wsap = pmax(0.,par_rhow * A * N * (beta1 * h + beta2 * hc)) ##newX
  W_c = pmax(0.,par_rhow * A * N * hc) #sapwood stem below Crown
  W_s = pmax(0.,par_rhow * A * N * par_betas * Lc) #sapwood stem within crown
  W_branch =  pmax(0.,par_rhow * A * N * betab * Lc) #branches biomass
  Wsh = pmax((A+B+sqrt(A*B)) * hc * par_rhow * N/2.9 - W_c,0) #initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
  #initialize Wdb dead branches biomass
  Wdb = pmax(0,ifelse(par_S_branchMod == 1.,Tbd * W_branch * ((0.0337+0.000009749*N)*exp(-0.00456*d^2)+0.00723),
         Tbd * W_branch *((-0.00513+0.000012*N)*exp((0.00000732-0.000000764*N)*d^2)+0.00467)))
  W_stem = pmax(0.,W_c + W_s + Wsh)

  W_crs = par_rhow * beta0 * A * h * N
  W_crh = Wsh * beta0
  W_croot = W_crh + W_crs
  #W_croot = pmax(0.,(par_rhow * Lc * beta0 * A / par_betas * N + (W_c + Wsh) * beta0)) #coarse root biomass
  V = W_stem / par_rhow
  biomasses <- rbind(wf_STKG,W_froot,W_wsap,W_c,W_s,W_branch,W_croot,Wsh,Wdb,W_stem,V,W_crh)

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

varNames  <- c('siteID','gammaC','sitetype','species','ETS' ,'P0','age', 'DeadWoodVolume', 'Respi_tot','GPP/1000',
               'H','D', 'BA','Hc_base','Cw','Ac','N','npp','leff','keff','lproj','ET_preles','weight',
               'Wbranch',"WfineRoots",'Litter_fol','Litter_fr','Litter_fWoody','Litter_cWoody','V',
               'Wstem','W_croot','wf_STKG', 'wf_treeKG','B_tree','Light',"Vharvested","Wharvested","soilC",
               "aSW","dH","Vmort","gross growth", "GPPspecies","Rh species", "NEP sp"," W_wsap","W_c","W_s","Wsh","Wdb","dHc",
               "Wbh","Wcrh")

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

  
  
  
  monthlyGPP <- function(GPPx,func="sum"){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    GPPbyYear <- matrix(GPPx,365,length(GPPx)/365)
    GPPm = unlist(apply(GPPbyYear, 2, function(x) 
      aggregate(x,by=list(monthsDays),FUN=func)$x))
    GPPm <- as.vector(GPPm)
    return(GPPm)
  }
  
  monthlyWeather <- function(varx,func){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    # varByYear <- matrix(varx,365,length(varx)/365)
    varM = unlist(apply(varx, 1, function(x) 
      aggregate(x,by=list(monthsDays),FUN=func)$x))
    varM <- as.vector(varM)
    return(varM)
  }
  
  monthlyFluxes <- function(modOut){
    cueGV <- ###0.5 carbon use efficiency (NPP/GPP) of ground vegetation
    # if(is.matrix(mGPP)){
    ##aggregating total GPP (tree layers + GV) at monthly time step
    mGPPtot <- t(apply(modOut$dailyPRELES[,,1],1,monthlyGPP,"sum")) 

    cueTrees <- modOut$multiOut[,,18,,1]/modOut$multiOut[,,44,,1]
    cueTrees[which(is.na(cueTrees))] <- 0.
    totGPP <- apply(modOut$multiOut[,,44,,1],1:2,sum) + modOut$GVout[,,3]
    cueAll <- pGPP <- array(NA,dim=c(modOut$nSites,modOut$maxYears,(modOut$maxNlayers+1)))
    cueAll[,,1] <- cueGV
    cueAll[,,2:(modOut$maxNlayers+1)] <- cueTrees
    pGPP[,,1] <- modOut$GVout[,,3]/totGPP
    
    for(i in 1:modOut$maxNlayers){
      pGPP[,,(i+1)] <- modOut$multiOut[,,44,i,1] / totGPP
    }
    mGPP <- mNPP <- array(NA,dim=c(modOut$nSites,modOut$maxYears*12,(modOut$maxNlayers+1)))
    for(i in 1:dim(pGPP)[3]){
      for(ij in 1:modOut$nSites){
        mGPP[ij,,i] <- rep(pGPP[ij,,i],each=12) * mGPPtot[ij,]
        mNPP[ij,,i] <- mGPP[ij,,i] * rep(cueAll[ij,,i],each=12)
      }
    }
    
    
    ###litterfall calculations lit is splitted to August and September
    Lf <- aperm(apply(modOut$multiOut[,,26,,1],c(1,3),mLit),c(2,1,3))
    Lfr <- aperm(apply(modOut$multiOut[,,27,,1],c(1,3),mLit),c(2,1,3))
    Lnw <- Lf + Lfr
    Lfw <- aperm(apply(modOut$multiOut[,,28,,1],c(1,3),mLit),c(2,1,3))
    Lw <- aperm(apply(modOut$multiOut[,,29,,1],c(1,3),mLit),c(2,1,3))
    
    nYears <- dim(modOut$multiOut)[2]
    nMonths <- nYears * 12
    nSites <- modOut$nSites
    nLayers <- dim(modOut$multiOut)[4]
    litter <- array(0.,dim=c(nSites,nMonths,nLayers,3))
    litter[,,,1] <- Lnw
    litter[,,,2] <- Lfw
    litter[,,,3] <- Lw
    
    species <- modOut$multiOut[,1,4,,1]
    nSp <- max(species)
    litterSize <- modOut$litterSize[,1:nSp]
    soilC <- array(0.,dim=c(nSites,(nMonths+1),5,3,nLayers))
    soilC[,1,,,] <- modOut$soilC[,1,,,]  ###initialize with soilC from simulations
    nClimID <- modOut$nClimID
    climIDs <- modOut$siteInfo[,2]
    weatherYasso <- array(0,dim=c(nClimID,nMonths,3))
    weatherYasso[,,1] <- t(apply(modOut$weather[,,,2],1,monthlyWeather,"mean")) ###monthly tempearture
    weatherYasso[,,2] <- t(apply(modOut$weather[,,,4],1,monthlyWeather,"sum")) *12 ###monthly precipitation ###*12 rescales to annual mm
    
    soilCtrees <- .Fortran("runYasso",litter=as.array(litter),
                           litterSize=as.array(litterSize),
                           nMonths=as.integer(nMonths), 
                           nLayers=as.integer(nLayers), 
                           nSites=as.integer(nSites),
                           nSp=as.integer(nSp),
                           species=as.matrix(species),
                           nClimID=as.integer(nClimID),
                           climIDs=as.integer(climIDs),
                           pAWEN=as.matrix(parsAWEN),
                           pYasso=as.double(pYAS),
                           climate=as.matrix(weatherYasso),
                           soilC=as.array(soilC)) 
    
    ###calculate soil C for gv
    fAPAR <- modOut$fAPAR
    fAPAR[which(is.na(modOut$fAPAR),arr.ind = T)] <- 0.
    AWENgv <- array(NA,dim=c(dim(modOut$fAPAR),4))
    p0 = modOut$multiOut[,,6,1,1]
    for(ij in 1:nYears){
      AWENgv[,ij,] <- t(sapply(1:nrow(fAPAR), function(i) .Fortran("fAPARgv",fAPAR[i,ij],
                                                                   modOut$ETSy[i,ij],modOut$siteInfo[i,2],
                                                                   0,0,p0[i,ij],rep(0,4))[[7]]))
    }
    
    mAWEN <- array(0,dim=c(nSites,nMonths,5))
    mAWEN[,,1:4] <- aperm(apply(AWENgv,c(1,3),mLit),c(2,1,3))
    mGVit <- t(apply(modOut$GVout[,,2],1,mLit))
    ###calculate steady state soil C per GV
    # ststGV <- matrix(NA,nSites,5)
    soilGV <- array(0,dim=c(nSites,(nMonths+1),5))
    
    
    soilCgv <- .Fortran("runYassoAWENin",
                        mAWEN=as.array(mAWEN),
                        nMonths=as.integer(nMonths),  
                        nSites=as.integer(nSites), 
                        litSize=0.,
                        nClimID=as.integer(nClimID),
                        climIDs=as.integer(climIDs),
                        pYasso=as.double(pYAS),
                        climate=as.matrix(weatherYasso),
                        soilCgv=as.array(soilGV))
    ####add gvsoilc to first layer foliage soilC
    # check in normal runs where ground vegetation soilC is calculated
    soilCtot <- apply(soilCtrees$soilC,1:2,sum) + apply(soilCgv$soilCgv,1:2,sum)
    mRhTot <- soilCtot[,1:nMonths]/10 - soilCtot[,2:(nMonths+1)]/10 +
      apply(litter,1:2,sum)/10 + mGVit/10
    mNPPtot <- apply(mNPP,1:2,sum)
    mRaTot <- mGPPtot - mNPPtot
    mNEPtot <- mNPPtot - mRhTot
    return(list(mGPP=mGPPtot,mNPP=mNPPtot,mRa=mRaTot,mRh=mRhTot,
                mNEP=mNEPtot,soilC=soilCtot))
  } 
  
  
  
  
  
  mLit <- function(aLit){
    nYears <- length(aLit)
    lit <- matrix(0,12,nYears)
    lit[8:9,] <- matrix(aLit/2,2,nYears,byrow = T)
    lit <- as.vector(lit)
    return(lit)
  }
  