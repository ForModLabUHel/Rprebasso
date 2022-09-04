InitMultiSite <- function(nYearsMS,
                          pCROBAS = pCROB,
                          pHcMod = pHcM,
                          pPRELES = pPREL,
                          etmodel = 0,
                          pYASSO =pYAS,
                          pAWEN = parsAWEN,
                          siteInfo = NA,
                          multiInitVar = NA,
                          multiThin = NA,
                          multiNthin = NA,
                          multiInitClearCut = NA,
                          fixBAinitClarcut = 1.,  ###if 1 when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relBa at last year 
                          initCLcutRatio = NA,  ###BA ratio per each species/layer (default is the ba ratio at the begginning of the simulations)
                          areas = NA,
                          PAR,
                          TAir,
                          VPD,
                          Precip,
                          CO2,
                          multiP0=NA,
                          soilC = NA,
                          weatherYasso = NA,
                          litterSize = litterSizeDef,
                          soilCtot = NA,
                          defaultThin = 1.,
                          ClCut = 1.,
                          energyCut = 0.,
                          inDclct = NA,
                          inAclct = NA,
                          yassoRun = 0,
                          smoothP0 = 1,
                          smoothETS = 1,
                          smoothYear=5,
                          HcModV=2,  ####version of model to compute Hc 1 uses the version of based on ksi parameter 2 uses the empirical model
                          tapioPars=pTapio,
                          thdPer = NA,
                          limPer = NA,
                          ftTapioPar = ftTapio,
                          tTapioPar = tTapio,
                          GVrun = 1,
                          thinInt = -999.,
                          mortMod = 1 #flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
){  
  
  nSites <- length(nYearsMS)
  if(length(thinInt)==1) thinInt <- rep(thinInt,nSites)
  if(all(is.na(thdPer))) thdPer <- rep(0.5,nSites)
  if(all(is.na(limPer))) limPer <- rep(0.5,nSites)
  if(all(is.na(areas))) areas <- rep(1.,nSites) ###each site is 1 ha (used to scale regional harvest)
  if(all(is.na(siteInfo))){
    siteInfo = matrix(c(1,1,3,160,0,0,20,3,3,413.,0.45,0.118),nSites,12,byrow = T) ###default values for nspecies and site type = 3
    siteInfo[,1] <- 1:nSites
  }
  colnames(siteInfo) <- c("siteID", "climID", "siteType", "SWinit", "CWinit", 
                          "SOGinit", "Sinit", "nLayers", "nSpecies", "soildepth", 
                          "effective field capacity", "permanent wilting point") 
  
  nLayers <- siteInfo[,8]
  # nSp <- siteInfo[,9]
  climIDs <- siteInfo[,2]
  # if(all(is.na(multiInitVar)) & all(is.na(nSp)) nSp <- rep(3,nSites)
  allSp = ncol(pCROBAS)
  varNam <- getVarNam()
  nVar <- length(varNam)
  
  nClimID <- length(unique(climIDs))
  if(!all((1:nClimID) %in% climIDs) | length(climIDs) != nSites) warning("check consistency between weather inputs and climIDs")
  if(nClimID == 1){
    nClimID = 2
    climIDs[1] <- 2
    siteInfo[1,2] <- 2
    PAR = matrix(PAR,2,length(PAR),byrow = T)
    TAir = matrix(TAir,2,length(TAir),byrow = T)
    VPD = matrix(VPD,2,length(VPD),byrow = T)
    Precip = matrix(Precip,2,length(Precip),byrow = T)
    CO2 <- matrix(CO2,2,length(CO2),byrow = T)
  }
  
  maxYears <- max(nYearsMS)
  maxNlayers <- max(nLayers)
  layerNam <- paste("layer",1:maxNlayers)
  multiOut <- array(0, dim=c(nSites,(maxYears),nVar,maxNlayers,2),
                    dimnames = list(site=NULL,year=NULL,variable=varNam,layer=layerNam,
                                    status=c("stand","thinned")))
  multiEnergyWood <- array(0, dim=c(nSites,(maxYears),maxNlayers,2),
                           dimnames = list(site=NULL,year=NULL,layer=layerNam,
                                           variable=c("volume","biomass")))
  initClearcut = initSeedling.def
  if (all(is.na(multiInitClearCut))) multiInitClearCut <- matrix(initClearcut,nSites,5,byrow = T)
  # multiInitClearCut <- cbind(multiInitClearCut,0.0008025897)
  ###process yasso inputs if missing
  if(all(is.na(soilC))) soilC <- array(0,dim=c(nSites,maxYears,5,3,maxNlayers))
  if(all(is.na(soilCtot))) soilCtot <- matrix(0,nSites,maxYears)
  
  ##process weather inputs for YASSO
  if(all(is.na(weatherYasso))){
    weatherYasso <- array(0,dim=c(nClimID,maxYears,3))
    if(nClimID>1){
      weatherYasso[,,1] <- t(apply(TAir[,1:(maxYears*365)],1,aTmean,maxYears))
      weatherYasso[,,3] <- t(apply(TAir[,1:(maxYears*365)],1,aTampl,maxYears))
      weatherYasso[,,2] <- t(apply(Precip[,1:(maxYears*365)],1,aPrecip,maxYears))
    } else{
      weatherYasso[1,,1] <- aTmean(TAir[1,1:(maxYears*365)],maxYears) 
      weatherYasso[1,,3] <- aTampl(TAir[1,1:(maxYears*365)],maxYears) 
      weatherYasso[1,,2] <- aPrecip(Precip[1,1:(maxYears*365)],maxYears) 
    }
  }
  
  if (length(defaultThin) == 1) defaultThin=as.double(rep(defaultThin,nSites))
  if (length(ClCut) == 1) ClCut=as.double(rep(ClCut,nSites))
  if (length(energyCut) == 1) energyCut=as.double(rep(energyCut,nSites))
  if (length(inDclct) == 1) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == 1) inAclct=matrix(inAclct,nSites,allSp)
  if (length(inDclct) == nSites) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == nSites) inAclct=matrix(inAclct,nSites,allSp)
  if (length(yassoRun) == 1) yassoRun=as.double(rep(yassoRun,nSites))
  # if (length(PREBASversion) == 1) PREBASversion=as.double(rep(PREBASversion,nSites))
  # 
  ###process ETS
  multiETS <- matrix(NA,nClimID,maxYears)
  for(climID in 1:nClimID){
    nYearsX <- max(nYearsMS[which(climIDs==climID)])
    Temp <- TAir[climID,1:(365*nYearsX)]-5
    ETS <- pmax(0,Temp,na.rm=T)
    ETS <- matrix(ETS,365,nYearsX); ETS <- colSums(ETS)
    multiETS[climID,(1:nYearsX)] <- ETS
    
    # xx <- min(10,nYearsX)
    # Ainit = 6 - 0.005*mean(ETS[1:xx]) + 2.25 ## need to add 2*sitetype
    # sitesClimID <- which(climIDs==climID)
    # multiInitClearCut[sitesClimID,5] <- replace(multiInitClearCut[sitesClimID,5],
    #                                             which(is.na(multiInitClearCut[sitesClimID,5])),round(Ainit))
  }
  xx <- min(10,nYearsX)
  Ainits <- multiInitClearCut[,5]
  for(xd in 1:nSites){
    if(is.na(Ainits[xd])) {
      Ainits[xd] = round(6 + 2* siteInfo[xd,3] - 0.005*mean(multiETS[siteInfo[xd,2],1:xx]) + 2.25)
      multiInitClearCut[xd,5] = 999.
    }
  } 
  ETSthres <- 1000; ETSmean <- rowMeans(multiETS)
  if(smoothETS==1. & maxYears > 1){
    for(i in 2:maxYears) multiETS[,i] <- multiETS[,(i-1)] + (multiETS[,i]-multiETS[,(i-1)])/min(i,smoothYear)
  } 
  multiETS[which(is.na(multiETS))] <- 0.
  ####process clearcut
  for(i in 1: nSites){
    if(ClCut[i]==1 & all(is.na(inDclct[i,]))) inDclct[i,] <-
        c(ClCutD_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          NA,NA,NA,NA,NA,NA,NA)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)')
    if(ClCut[i]==1 & all(is.na(inAclct[i,]))) inAclct[i,] <-
        c(ClCutA_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          80,50,13,30,50,13,120)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)')
    if(any(!is.na(inDclct[i,]))) inDclct[i,is.na(inDclct[i,])] <- max(inDclct[i,],na.rm=T)
    if(all(is.na(inDclct[i,]))) inDclct[i,] <- 9999999.99
    if(any(!is.na(inAclct[i,]))) inAclct[i,is.na(inAclct[i,])] <- max(inAclct[i,],na.rm=T)
    if(all(is.na(inAclct[i,]))) inAclct[i,] <- 9999999.99
  }
  
  maxThin <- max(multiNthin)
  ###thinning if missing.  To improve
  if(all(is.na(multiThin))){
    multiNthin <- rep(0,nSites)
    maxThin <- 2
    multiThin <- array(0, dim=c(nSites,maxThin,10))
    multiThin[,,9:10] <- -999
  }
  multiThin[is.na(multiThin)] <- -999
  
  ###PROCESS weather inputs for prebas
  multiweather <- array(-999,dim=c(nClimID,maxYears,365,5))
  
  ##extract weather inputs
  for(i in 1:nClimID){
    nYearsX <- max(nYearsMS[which(climIDs==i)])
    weatherPreles <- array(c(PAR[i,1:(365*nYearsX)],TAir[i,1:(365*nYearsX)],
                             VPD[i,1:(365*nYearsX)],Precip[i,1:(365*nYearsX)],
                             CO2[i,1:(365*nYearsX)]),dim=c(365,nYearsX,5))
    
    weatherPreles <- aperm(weatherPreles, c(2,1,3))
    
    multiweather[i,(1:nYearsX),,] <- weatherPreles
  }
  
  ### compute P0
  ###if P0 is not provided use preles to compute P0
  if(all(is.na(multiP0))){
    multiP0 <- array(NA,dim=c(nClimID,maxYears,2))
    for(climID in 1:nClimID){
      nYearsX <- max(nYearsMS[which(climIDs==climID)])
      P0 <- PRELES(DOY=rep(1:365,nYearsX),PAR=PAR[climID,1:(365*nYearsX)],
                   TAir=TAir[climID,1:(365*nYearsX)],VPD=VPD[climID,1:(365*nYearsX)],
                   Precip=Precip[climID,1:(365*nYearsX)],CO2=rep(380,(365*nYearsX)),
                   fAPAR=rep(1,(365*nYearsX)),LOGFLAG=0,p=pPRELES)$GPP
      P0 <- matrix(P0,365,nYearsX)
      multiP0[climID,(1:nYearsX),1] <- colSums(P0)
    }
    if(smoothP0==1 & maxYears > 1){
      multiP0[,1,2] <- multiP0[,1,1]
      for(i in 2:maxYears) multiP0[,i,2] <- multiP0[,(i-1),2] + (multiP0[,i,1]-multiP0[,(i-1),2])/min(i,smoothYear)
      # multiP0[,,2] <- matrix(rowMeans(multiP0[,,1]),nClimID,maxYears,byrow = F)
    } else{
      multiP0[,,2] <- multiP0[,,1]
    }
  }
  multiP0[which(is.na(multiP0))] <- 0.
  
  if (all(is.na(multiInitVar))){
    multiInitVar <- array(NA,dim=c(nSites,7,maxNlayers))
    multiInitVar[,1,] <- rep(1:maxNlayers,each=nSites)
    multiInitVar[,3,] <- initClearcut[1]; multiInitVar[,4,] <- initClearcut[2]
    multiInitVar[,5,] <- initClearcut[3]/maxNlayers; multiInitVar[,6,] <- initClearcut[4]
    multiInitVar[,2,] <- matrix(Ainits,nSites,maxNlayers)
    for(ikj in 1:maxNlayers){
      p_ksi <- pCROBAS[38,multiInitVar[,1,ikj]]
      p_rhof <- pCROBAS[15,multiInitVar[,1,ikj]]
      p_z <- pCROBAS[11,multiInitVar[,1,ikj]]
      Lc <- multiInitVar[,3,ikj] - multiInitVar[,6,ikj]
      A <- as.numeric(p_ksi/p_rhof * Lc^p_z)
      multiInitVar[,7,ikj] <- A     
    } 
    multiInitVar[which(is.na(multiInitVar))] <- 0.
  }else{
    ####if Height of the crown base is not available use model
    if(maxNlayers==1){
      multiInitVar <- array(aaply(multiInitVar,1,findHcNAs,pHcMod),dim=c(nSites,7,1))
    }else{
      multiInitVar <- aaply(multiInitVar,1,findHcNAs,pHcMod)
    }
    
    ###age cannot be lower than 1 year
    multiInitVar[,2,][which(multiInitVar[,2,]<1)] <- 1
    
    ####compute A
    if(all(is.na(multiInitVar[,7,]))|all(multiInitVar[,7,]==0)){
      for(ikj in 1:maxNlayers){
        not0 <- which(multiInitVar[,3,ikj]>0)
        p_ksi <- pCROBAS[38,multiInitVar[not0,1,ikj]]
        p_rhof <- pCROBAS[15,multiInitVar[not0,1,ikj]]
        p_z <- pCROBAS[11,multiInitVar[not0,1,ikj]]
        Lc <- multiInitVar[not0,3,ikj] - multiInitVar[not0,6,ikj]
        A <- as.numeric(p_ksi/p_rhof * Lc^p_z)
        multiInitVar[not0,7,ikj] <- A     
      } 
    }
    LcCheck <- multiInitVar[,3,] - multiInitVar[,6,]
    if(any(LcCheck<0.)) return("check, some Lc is negative")
    # p_ksi = matrix(pCROBAS[38,multiInitVar[,1,]],nSites,maxNlayers)
    #  p_rhof <- matrix(pCROBAS[15,multiInitVar[,1,]],nSites,maxNlayers)
    #  p_z <- matrix(pCROBAS[11,multiInitVar[,1,]],nSites,maxNlayers)
    #  Lc <- multiInitVar[,3,] - multiInitVar[,6,]
    #  A <- p_ksi/p_rhof * Lc^p_z
    #  multiInitVar[,7,] <- A      # p_ksi=pCROBAS[38,multiInitVar[,1,]]
    # p_rhof <- pCROBAS[15,multiInitVar[,1,]]
    # p_z <- pCROBAS[11,multiInitVar[,1,]]
    # Lc <- multiInitVar[,3,] - multiInitVar[,6,]
    # A <- p_ksi/p_rhof * Lc^p_z
    # multiInitVar[,7,] <- A
    # N = multiInitVar[,5,]/(pi*((multiInitVar[,4,]/2/100)**2))
    # B = multiInitVar[,5,]/N
    # Lc = multiInitVar[,3,] - multiInitVar[,6,]
    # rc = Lc / (multiInitVar[,3,]-1.3) 
    # multiInitVar[,7,] = rc * B
    # multiInitVar[which(is.na(multiInitVar))] <- 0.
    # ops <- which(multiInitVar[,6,]<1.3 & multiInitVar[,3,]>0.,arr.ind = T)
    # if(length(ops)>0.){
    # p_ksi=pCROBAS[38,multiInitVar[,1,][ops]]
    # p_rhof <- pCROBAS[15,multiInitVar[,1,][ops]]
    # p_z <- pCROBAS[11,multiInitVar[,1,][ops]]
    # Lc <- multiInitVar[,3,][ops] - multiInitVar[,6,][ops]
    # A <- p_ksi/p_rhof * Lc^p_z
    # multiInitVar[,7,][ops] <- A
    # }
  }
  
  multiInitVar[,6:7,][which(is.na(multiInitVar[,6:7,]))] <- 0
  if(length(fixBAinitClarcut)==1) fixBAinitClarcut=rep(fixBAinitClarcut,nSites)
  
  if(all(is.na(initCLcutRatio))){
    if(maxNlayers==1){
      initCLcutRatio <- rep(1,nSites)  
    }else{
      initCLcutRatio <- multiInitVar[,5,]/rowSums(multiInitVar[,5,]) 
    }
  }
  initCLcutRatio[which(is.na(initCLcutRatio))] <- 0.
  
  # if(all(is.na(litterSize))){
  #   litterSize <- matrix(0,3,allSp)
  #   litterSize[2,] <- 2
  #   litterSize[1,] <- c(30,30,10)
  #   # siteInfo <- siteInfo[,-c(4,5)]
  # }
  
  ###!!!###initiaize biomasses
  initVarX <- abind(multiInitVar,matrix(siteInfo[,3],nSites,maxNlayers),along=2)
  biomasses <- array(apply(initVarX,1,initBiomasses,pCro=pCROBAS),dim=c(12,maxNlayers,nSites))
  biomasses <- aperm(biomasses,c(3,1,2))
  biomasses[which(is.na(biomasses))] <- 0
  
  multiOut[,1,c(33,25,47:49,24,32,50,51,31,30,54),,1] <- biomasses
  # multiInitVar <- multiInitVar[,1:7,1:maxNlayers]
  
  for(i in 1:maxNlayers){
    sitxx <- which(multiInitVar[,3,i]==0)
    multiOut[sitxx,,,i,] <- 0.
  }
  
  dimnames(multiInitVar) <-  list(site=NULL,
                                  variable=c("SpeciesID","age","H","D","BA","Hc","Ac"),layer=layerNam)
  
  multiSiteInit <- list(
    multiOut = multiOut,
    multiEnergyWood = multiEnergyWood,
    nSites = nSites,
    nClimID = nClimID,
    nLayers = nLayers,
    maxYears = maxYears,
    maxThin = maxThin,
    nYears = nYearsMS,
    areas = areas,
    thinning = multiThin,
    pCROBAS = pCROBAS,
    allSp = allSp,
    siteInfo = siteInfo,
    maxNlayers = maxNlayers,
    nThinning = multiNthin,
    fAPAR = matrix(0.7,nSites,maxYears),
    initClearcut = multiInitClearCut,
    fixBAinitClarcut=fixBAinitClarcut,
    initCLcutRatio = initCLcutRatio,
    ETSy = multiETS,
    P0y = multiP0,
    multiInitVar = multiInitVar,
    weather = multiweather, 
    DOY = 1:365,
    pPRELES = pPRELES,
    etmodel = etmodel,
    soilC = soilC,
    pYASSO = pYASSO,
    pAWEN = pAWEN,
    weatherYasso = weatherYasso,
    litterSize = litterSize,
    soilCtot = soilCtot,
    defaultThin = defaultThin,
    ClCut = ClCut,
    energyCut = energyCut,
    inDclct = inDclct,
    inAclct = inAclct,
    dailyPRELES = array(-999,dim=c(nSites,(maxYears*365),3)),
    yassoRun = yassoRun,
    smoothP0 = smoothP0,
    smoothETS = smoothETS,
    tapioPars=tapioPars,
    thdPer = thdPer,
    limPer = limPer,
    ftTapioPar = ftTapioPar,
    tTapioPar = tTapioPar,
    GVrun=as.integer(GVrun),
    GVout=array(0.,dim = c(nSites,maxYears,5)),
    thinInt = thinInt,
    mortMod = mortMod
  )
  return(multiSiteInit)
}

multiPrebas <- function(multiSiteInit,
                        fertThin = 0,
                        nYearsFert = 20,
                        oldLayer=0){
  ###initialize siteType
  multiSiteInit$multiOut[,,3,,1] <- array(multiSiteInit$siteInfo[,3],
                                          dim=c(multiSiteInit$nSites,
                                                multiSiteInit$maxYears,
                                                multiSiteInit$maxNlayers))
  
  for(ijj in 1:multiSiteInit$maxNlayers){
    siteXs <- which(multiSiteInit$multiInitVar[,1,ijj] %in% 1:ncol(multiSiteInit$pCROBAS))
    multiSiteInit$multiOut[siteXs,,3,ijj,2] =
      matrix(multiSiteInit$pCROBAS[cbind((20+pmin(multiSiteInit$siteInfo[,3],5))[siteXs],
                                         multiSiteInit$multiInitVar[siteXs,1,ijj])],
             length(siteXs),multiSiteInit$maxYears)
  }
  
  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }
  prebas <- .Fortran("multiPrebas",
                     multiOut = as.array(multiSiteInit$multiOut),
                     nSites = as.integer(multiSiteInit$nSites),
                     nClimID = as.integer(multiSiteInit$nClimID),
                     nLayers = as.integer(multiSiteInit$nLayers),######
                     maxYears = as.integer(multiSiteInit$maxYears),
                     maxThin = as.integer(multiSiteInit$maxThin),
                     nYears = as.integer(multiSiteInit$nYears),
                     thinning=as.array(multiSiteInit$thinning),
                     pCROBAS = as.matrix(multiSiteInit$pCROBAS),    ####
                     allSp = as.integer(multiSiteInit$allSp),       ####
                     siteInfo = as.matrix(multiSiteInit$siteInfo[,c(1:7,10:12)]),  ####
                     maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                     nThinning=as.integer(multiSiteInit$nThinning),
                     fAPAR=as.matrix(multiSiteInit$fAPAR),
                     initClearcut=as.matrix(multiSiteInit$initClearcut),
                     fixBAinitClearcut = as.double(multiSiteInit$fixBAinitClarcut),
                     initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                     ETSy=as.matrix(multiSiteInit$ETSy),
                     P0y=as.array(multiSiteInit$P0y),
                     multiInitVar=as.array(multiSiteInit$multiInitVar),
                     weather=as.array(multiSiteInit$weather),
                     DOY= as.integer(multiSiteInit$DOY),
                     pPRELES=as.double(multiSiteInit$pPRELES),
                     etmodel=as.integer(multiSiteInit$etmodel),
                     soilC = as.array(multiSiteInit$soilC),
                     pYASSO=as.double(multiSiteInit$pYASSO),
                     pAWEN = as.matrix(multiSiteInit$pAWEN),
                     weatherYasso = as.array(multiSiteInit$weatherYasso),
                     litterSize = as.array(multiSiteInit$litterSize),
                     soilCtot = as.matrix(multiSiteInit$soilCtot),
                     defaultThin=as.double(multiSiteInit$defaultThin),
                     ClCut=as.double(multiSiteInit$ClCut),
                     energyCut=as.double(multiSiteInit$energyCut),
                     inDclct=as.matrix(multiSiteInit$inDclct),
                     inAclct=as.matrix(multiSiteInit$inAclct),
                     dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                     yassoRun=as.double(multiSiteInit$yassoRun),
                     multiEnergyWood = as.array(multiSiteInit$multiEnergyWood),
                     tapioPars = as.array(multiSiteInit$tapioPars),
                     thdPer=as.double(multiSiteInit$thdPer),
                     limPer=as.double(multiSiteInit$limPer),
                     ftTapioPar = as.array(multiSiteInit$ftTapioPar),
                     tTapioPar = as.array(multiSiteInit$tTapioPar),
                     GVout = as.array(multiSiteInit$GVout),
                     GVrun = as.integer(multiSiteInit$GVrun),
                     thinInt=as.double(multiSiteInit$thinInt),
                     fertThin = as.integer(fertThin),
                     flagFert = as.integer(0),
                     nYearsFert = as.integer(nYearsFert),
                     oldLayer=as.integer(oldLayer),
                     mortMod=as.integer(multiSiteInit$mortMod)
  )
  dimnames(prebas$multiOut) <- dimnames(multiSiteInit$multiOut)
  dimnames(prebas$multiInitVar) <- dimnames(multiSiteInit$multiInitVar)
  names(prebas$siteInfo) <- names(multiSiteInit$siteInfo)
  
  class(prebas) <- "multiPrebas"
  return(prebas)
}


regionPrebas <- function(multiSiteInit,
                         HarvLim = NA,
                         minDharv = 999.,
                         cutAreas = NA,  ### is a matrix: area of cuttings rows are years of simulations
                         ###columns: clcutArea target(1), simulated clCut area(2) (set to 0. will be filled by prebas output);
                         ####precom-thin target(3), sim(4); area firstThin targ(5), sim(6)
                         compHarv = 0.,###flag for compensating harvest if harvest do not reach the desired levels
                         ####compHarv=0 -> no compensation, compHarv=1 compensate harvest with clearcut
                         ### compHarv=2 compensate harvest with thinnings
                         thinFact = 0.25, ####if compHarv = 2 -> thinFact is the percentage of thinning to compansate harvest
                         #######compHarv[1]
                         ageHarvPrior = 0., ####flag used in the IBC-carbon runs of
                         ####the mitigation Scenario and biodiversity protection 
                         ####scenario (protect). If higher then 0. the scenarios is activated and
                         #####the sites are ordered according to the siteType and
                         ###priority is given to the sites where age is lower then ageHarvPrior
                         siteOrder=NA,
                         fertThin = 0.,
                         nYearsFert = 20,
                         oldLayer=0, ####oldLayer == 1 will leave 5-10% basal area at clearcut in the old layer
                         startSimYear=1
){
  
  if(length(HarvLim)==2) HarvLim <- matrix(HarvLim,multiSiteInit$maxYears,2,byrow = T)
  if(all(is.na(HarvLim))) HarvLim <- matrix(0.,multiSiteInit$maxYears,2)
  if(all(is.na(cutAreas))) cutAreas <- matrix(-999.,(multiSiteInit$maxYears),6)
  compHarv <- c(compHarv,thinFact)
  if(ageHarvPrior > 0.){
    sitesCl1 <- which(multiSiteInit$siteInfo[,3]<3.5)
    sitesCl2 <- which(multiSiteInit$siteInfo[,3]>3.5)
    siteOrder1 <- replicate(multiSiteInit$maxYears,sample(sitesCl1))
    siteOrder2 <- replicate(multiSiteInit$maxYears,sample(sitesCl2))
    siteOrder <- rbind(siteOrder1,siteOrder2)
  }else if(all(is.na(siteOrder))){
    siteOrder <- replicate(multiSiteInit$maxYears,sample(1:multiSiteInit$nSites))
  }  
  # reorder first year of siteOreder according to age of the stands, 
  # because in Fortran first year has a bug probably 
  # in the PACK function
if(ageHarvPrior>0){
  domSp <- apply(multiSiteInit$multiInitVar[,5,],1,which.max)
  agesX <- multiSiteInit$multiInitVar[,2,][cbind(1:multiSiteInit$nSites,domSp)]
  newOrdX <- c(which(agesX[siteOrder[,1]] <= ageHarvPrior),
               which(agesX[siteOrder[,1]] > ageHarvPrior))
  siteOrder[,1] <- siteOrder[newOrdX,1]
}  
  
  ###initialize siteType
  multiSiteInit$multiOut[,,3,,1] <- array(multiSiteInit$siteInfo[,3],
                                          dim=c(multiSiteInit$nSites,
                                                multiSiteInit$maxYears,
                                                multiSiteInit$maxNlayers))
  for(ijj in 1:multiSiteInit$maxNlayers){
    siteXs <- which(multiSiteInit$multiInitVar[,1,ijj] %in% 1:ncol(multiSiteInit$pCROBAS))
    multiSiteInit$multiOut[siteXs,,3,ijj,2] =
      matrix(multiSiteInit$pCROBAS[cbind((20+pmin(multiSiteInit$siteInfo[,3],5))[siteXs],
                                         multiSiteInit$multiInitVar[siteXs,1,ijj])],
             length(siteXs),multiSiteInit$maxYears)
  }
  
  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }
  prebas <- .Fortran("regionPrebas",
                     siteOrder = as.matrix(siteOrder),
                     HarvLim = as.matrix(HarvLim),
                     minDharv = as.double(minDharv),
                     multiOut = as.array(multiSiteInit$multiOut),
                     nSites = as.integer(multiSiteInit$nSites),
                     areas = as.double(multiSiteInit$areas),
                     nClimID = as.integer(multiSiteInit$nClimID),
                     nLayers = as.integer(multiSiteInit$nLayers),######
                     maxYears = as.integer(multiSiteInit$maxYears),
                     maxThin = as.integer(multiSiteInit$maxThin),
                     nYears = as.integer(multiSiteInit$nYears),
                     thinning=as.array(multiSiteInit$thinning),
                     pCROBAS = as.matrix(multiSiteInit$pCROBAS),    ####
                     allSp = as.integer(multiSiteInit$allSp),       ####
                     siteInfo = as.matrix(multiSiteInit$siteInfo[,c(1:7,10:12)]),  ####
                     maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                     nThinning=as.integer(multiSiteInit$nThinning),
                     fAPAR=as.matrix(multiSiteInit$fAPAR),
                     initClearcut=as.matrix(multiSiteInit$initClearcut),
                     fixBAinitClarcut = as.double(multiSiteInit$fixBAinitClarcut),
                     initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                     ETSy=as.matrix(multiSiteInit$ETSy),
                     P0y=as.array(multiSiteInit$P0y),
                     multiInitVar=as.array(multiSiteInit$multiInitVar),
                     weather=as.array(multiSiteInit$weather),
                     DOY= as.integer(multiSiteInit$DOY),
                     pPRELES=as.double(multiSiteInit$pPRELES),
                     etmodel=as.integer(multiSiteInit$etmodel),
                     soilC = as.array(multiSiteInit$soilC),
                     pYASSO=as.double(multiSiteInit$pYASSO),
                     pAWEN = as.matrix(multiSiteInit$pAWEN),
                     weatherYasso = as.array(multiSiteInit$weatherYasso),
                     litterSize = as.array(multiSiteInit$litterSize),
                     soilCtot = as.matrix(multiSiteInit$soilCtot),
                     defaultThin=as.double(multiSiteInit$defaultThin),
                     ClCut=as.double(multiSiteInit$ClCut),
                     energyCut=as.double(multiSiteInit$energyCut),
                     inDclct=as.matrix(multiSiteInit$inDclct),
                     inAclct=as.matrix(multiSiteInit$inAclct),
                     dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                     yassoRun=as.double(multiSiteInit$yassoRun),
                     multiEnergyWood = as.array(multiSiteInit$multiEnergyWood),
                     tapioPars = as.array(multiSiteInit$tapioPars),
                     thdPer=as.double(multiSiteInit$thdPer),
                     limPer=as.double(multiSiteInit$limPer),
                     ftTapioPar = as.array(multiSiteInit$ftTapioPar),
                     tTapioPar = as.array(multiSiteInit$tTapioPar),
                     GVout = as.array(multiSiteInit$GVout),
                     GVrun = as.integer(multiSiteInit$GVrun),
                     cutAreas=as.matrix(cutAreas),
                     compHarv=as.double(compHarv),
                     thinInt=as.double(multiSiteInit$thinInt),
                     ageHarvPrior = as.double(ageHarvPrior),
                     fertThin = as.integer(fertThin),
                     flagFert = as.integer(rep(0,multiSiteInit$nSites)),
                     nYearsFert = as.integer(nYearsFert),
                     oldLayer=as.integer(oldLayer),
                     mortMod=as.integer(multiSiteInit$mortMod),
                     startSimYear = as.integer(startSimYear)
  )
  class(prebas) <- "regionPrebas"
  if(prebas$maxNlayers>1){
    rescalVbyArea <- prebas$multiOut[,,37,,1] * prebas$areas
    prebas$totHarv <- apply(rescalVbyArea,2,sum)
  }else{
    prebas$totHarv <- colSums(prebas$multiOut[,,37,1,1]*prebas$areas)
  }
  
  dimnames(prebas$multiOut) <- dimnames(multiSiteInit$multiOut)
  dimnames(prebas$multiInitVar) <- dimnames(multiSiteInit$multiInitVar)
  names(prebas$siteInfo) <- names(multiSiteInit$siteInfo)
  return(prebas)
}




reStartRegionPrebas <- function(multiSiteInit,
                         HarvLim,
                         minDharv = 999.,
                         cutAreas,  ### is a matrix: area of cuttings rows are years of simulations
                         ###columns: clcutArea target(1), simulated clCut area(2) (set to 0. will be filled by prebas output);
                         ####precom-thin target(3), sim(4); area firstThin targ(5), sim(6)
                         compHarv,###flag for compensating harvest if harvest do not reach the desired levels
                         ####compHarv=0 -> no compensation, compHarv=1 compensate harvest with clearcut
                         ### compHarv=2 compensate harvest with thinnings
                         thinFact=0.25, ####if compHarv = 2 -> thinFact is the percentage of thinning to compansate harvest
                         #######compHarv[1]
                         ageHarvPrior = 0, ####flag used in the IBC-carbon runs of
                         ####the mitigation Scenario and biodiversity protection 
                         ####scenario (protect). If higher then 0. the scenarios is activated and
                         #####the sites are ordered according to the siteType and
                         ###priority is given to the sites where age is lower then ageHarvPrior
                         siteOrder=NA,
                         fertThin = 0.,
                         nYearsFert = 20,
                         oldLayer=0, ####oldLayer == 1 will leave 5-10% basal area at clearcut in the old layer
                         startSimYear
){
  
  if(length(HarvLim)==2) HarvLim <- matrix(HarvLim,multiSiteInit$maxYears,2,byrow = T)
  if(all(is.na(HarvLim))) HarvLim <- matrix(0.,multiSiteInit$maxYears,2)
  if(all(is.na(cutAreas))) cutAreas <- matrix(-999.,(multiSiteInit$maxYears),6)
  compHarv <- c(compHarv,thinFact)
  if(ageHarvPrior > 0.){
    sitesCl1 <- which(multiSiteInit$siteInfo[,3]<3.5)
    sitesCl2 <- which(multiSiteInit$siteInfo[,3]>3.5)
    siteOrder1 <- replicate(multiSiteInit$maxYears,sample(sitesCl1))
    siteOrder2 <- replicate(multiSiteInit$maxYears,sample(sitesCl2))
    siteOrder <- rbind(siteOrder1,siteOrder2)
  }else if(all(is.na(siteOrder))){
    siteOrder <- replicate(multiSiteInit$maxYears,sample(1:multiSiteInit$nSites))
  }  
  # reorder first year of siteOreder according to age of the stands, 
  # because in Fortran first year has a bug probably 
  # in the PACK function
  if(ageHarvPrior>0){
    domSp <- apply(multiSiteInit$multiInitVar[,5,],1,which.max)
    agesX <- multiSiteInit$multiInitVar[,2,][cbind(1:multiSiteInit$nSites,domSp)]
    newOrdX <- c(which(agesX[siteOrder[,1]] <= ageHarvPrior),
                 which(agesX[siteOrder[,1]] > ageHarvPrior))
    siteOrder[,1] <- siteOrder[newOrdX,1]
  }  
  
  ###initialize siteType
  multiSiteInit$multiOut[,,3,,1] <- array(multiSiteInit$siteInfo[,3],
                                          dim=c(multiSiteInit$nSites,
                                                multiSiteInit$maxYears,
                                                multiSiteInit$maxNlayers))
  for(ijj in 1:multiSiteInit$maxNlayers){
    siteXs <- which(multiSiteInit$multiInitVar[,1,ijj] %in% 1:ncol(multiSiteInit$pCROBAS))
    multiSiteInit$multiOut[siteXs,,3,ijj,2] =
      matrix(multiSiteInit$pCROBAS[cbind((20+pmin(multiSiteInit$siteInfo[,3],5))[siteXs],
                                         multiSiteInit$multiInitVar[siteXs,1,ijj])],
             length(siteXs),multiSiteInit$maxYears)
  }
  
  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }
  prebas <- .Fortran("regionPrebas",
                     siteOrder = as.matrix(siteOrder),
                     HarvLim = as.matrix(HarvLim),
                     minDharv = as.double(minDharv),
                     multiOut = as.array(multiSiteInit$multiOut),
                     nSites = as.integer(multiSiteInit$nSites),
                     areas = as.double(multiSiteInit$areas),
                     nClimID = as.integer(multiSiteInit$nClimID),
                     nLayers = as.integer(multiSiteInit$nLayers),######
                     maxYears = as.integer(multiSiteInit$maxYears),
                     maxThin = as.integer(multiSiteInit$maxThin),
                     nYears = as.integer(multiSiteInit$nYears),
                     thinning=as.array(multiSiteInit$thinning),
                     pCROBAS = as.matrix(multiSiteInit$pCROBAS),    ####
                     allSp = as.integer(multiSiteInit$allSp),       ####
                     siteInfo = as.matrix(multiSiteInit$siteInfo),  ####
                     maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                     nThinning=as.integer(multiSiteInit$nThinning),
                     fAPAR=as.matrix(multiSiteInit$fAPAR),
                     initClearcut=as.matrix(multiSiteInit$initClearcut),
                     fixBAinitClarcut = as.double(multiSiteInit$fixBAinitClarcut),
                     initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                     ETSy=as.matrix(multiSiteInit$ETSy),
                     P0y=as.array(multiSiteInit$P0y),
                     multiInitVar=as.array(multiSiteInit$multiInitVar),
                     weather=as.array(multiSiteInit$weather),
                     DOY= as.integer(multiSiteInit$DOY),
                     pPRELES=as.double(multiSiteInit$pPRELES),
                     etmodel=as.integer(multiSiteInit$etmodel),
                     soilC = as.array(multiSiteInit$soilC),
                     pYASSO=as.double(multiSiteInit$pYASSO),
                     pAWEN = as.matrix(multiSiteInit$pAWEN),
                     weatherYasso = as.array(multiSiteInit$weatherYasso),
                     litterSize = as.array(multiSiteInit$litterSize),
                     soilCtot = as.matrix(multiSiteInit$soilCtot),
                     defaultThin=as.double(multiSiteInit$defaultThin),
                     ClCut=as.double(multiSiteInit$ClCut),
                     energyCut=as.double(multiSiteInit$energyCut),
                     inDclct=as.matrix(multiSiteInit$inDclct),
                     inAclct=as.matrix(multiSiteInit$inAclct),
                     dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                     yassoRun=as.double(multiSiteInit$yassoRun),
                     multiEnergyWood = as.array(multiSiteInit$multiEnergyWood),
                     tapioPars = as.array(multiSiteInit$tapioPars),
                     thdPer=as.double(multiSiteInit$thdPer),
                     limPer=as.double(multiSiteInit$limPer),
                     ftTapioPar = as.array(multiSiteInit$ftTapioPar),
                     tTapioPar = as.array(multiSiteInit$tTapioPar),
                     GVout = as.array(multiSiteInit$GVout),
                     GVrun = as.integer(multiSiteInit$GVrun),
                     cutAreas=as.matrix(cutAreas),
                     compHarv=as.double(compHarv),
                     thinInt=as.double(multiSiteInit$thinInt),
                     ageHarvPrior = as.double(ageHarvPrior),
                     fertThin = as.integer(fertThin),
                     flagFert = as.integer(rep(0,multiSiteInit$nSites)),
                     nYearsFert = as.integer(nYearsFert),
                     oldLayer=as.integer(oldLayer),
                     mortMod=as.integer(multiSiteInit$mortMod),
                     startSimYear = as.integer(startSimYear)
  )
  class(prebas) <- "regionPrebas"
  if(prebas$maxNlayers>1){
    rescalVbyArea <- prebas$multiOut[,,37,,1] * prebas$areas
    prebas$totHarv <- apply(rescalVbyArea,2,sum)
  }else{
    prebas$totHarv <- colSums(prebas$multiOut[,,37,1,1]*prebas$areas)
  }
  
  dimnames(prebas$multiOut) <- dimnames(multiSiteInit$multiOut)
  dimnames(prebas$multiInitVar) <- dimnames(multiSiteInit$multiInitVar)
  names(prebas$siteInfo) <- names(multiSiteInit$siteInfo)
  return(prebas)
}

