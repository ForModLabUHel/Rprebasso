#' Title
#'
#' @param nYearsMS 
#' @param pCROBAS 
#' @param pHcMod 
#' @param pPRELES 
#' @param etmodel 
#' @param pYASSO 
#' @param pAWEN 
#' @param siteInfo 
#' @param multiInitVar 
#' @param multiThin 
#' @param multiNthin 
#' @param multiInitClearCut 
#' @param fixBAinitClarcut 
#' @param initCLcutRatio 
#' @param areas 
#' @param PAR 
#' @param TAir 
#' @param VPD 
#' @param Precip 
#' @param CO2 
#' @param multiP0 
#' @param soilC 
#' @param weatherYasso 
#' @param litterSize 
#' @param soilCtot 
#' @param defaultThin 
#' @param ClCut 
#' @param energyCut 
#' @param inDclct 
#' @param inAclct 
#' @param yassoRun 
#' @param smoothP0 
#' @param smoothETS 
#' @param smoothYear 
#' @param HcModV  version of model to compute Hc 1 uses the version of based on ksi parameter 2 uses the empirical model; default value (HcModV_def) is 1
#' @param tapioPars 
#' @param thdPer 
#' @param limPer 
#' @param ftTapioPar 
#' @param tTapioPar 
#' @param GVrun 
#' @param thinInt 
#' @param mortMod 
#' @param ECMmod 
#' @param pECMmod 
#' @param layerPRELES 
#' @param LUEtrees 
#' @param LUEgv 
#' @param alpharNcalc #alphar calculations based on Nitrogen availability. deafault value is FALSE (no nitrogen impact). =1calculates N uptake
#' @param p0currClim # vector of average annual P0 for the climIDs at current climate. if NA the first five years of the simulations will be used to calculate it.
#' @param TcurrClim # vector of average annual temperature for the climIDs at current climate. if NA the first five years of the simulations will be used to calculate it.
#' @param PcurrClim # vector of average annual precipitation for the climIDs current climate. if NA the first five years of the simulations will be used to calculate it.
#'
#' @return Initialize PREBAS and return an object list that can be inputted to multiPrebas and regionPrebas functions to run PREBAS 
#' @export
#'
#' @examples
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
                          multiInitClearCut = NA, ###A matrix (rows are sites columns are the variables) with initial stand variables after clearcut: H, D, BA, Hc, Ainit. Ainit is the year when the stand reaches measurable size. If NA the default values from initSeedling.def are used Ainit and is automatically computed using air temperature.
                          fixBAinitClarcut = 1.,  ###if 1, when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relative basal area at last year 
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
                          HcModV=HcModV_def,  ####version of model to compute Hc 1 uses the version of based on ksi parameter 2 uses the empirical model; default value (HcModV_def) is 1
                          tapioPars=pTapio,
                          thdPer = NA,
                          limPer = NA,
                          ftTapioPar = ftTapio,
                          tTapioPar = tTapio,
                          GVrun = 1,
                          thinInt = -999.,
                          mortMod = 1, #flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
                          ECMmod=0, #flag for ECM modelling MAkela et al.2022
                          pECMmod = parsECMmod,
                          layerPRELES = 0,
                          LUEtrees = pLUEtrees,
                          LUEgv = pLUEgv,
                          alpharNcalc=FALSE,
                          p0currClim = NA,
                          TcurrClim = NA,
                          PcurrClim = NA
                          ){  
  
  nSites <- length(nYearsMS)
  if(length(mortMod)==1) mortMod <- rep(mortMod,2)
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
      Ainits[xd] = max(round(6 + 2* siteInfo[xd,3] - 0.005*mean(multiETS[siteInfo[xd,2],1:xx]) + 2.25+2),2)
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
          NA,NA,NA,NA,NA,NA,NA,NA)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)','quil')
    if(ClCut[i]==1 & all(is.na(inAclct[i,]))) inAclct[i,] <-
        c(ClCutA_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          80,50,13,30,50,13,120,100)  ###"fasy","pipi","eugl","rops","popu",'eugrur','piab(DE)','quil')
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
  
  if(all(is.na(multiInitVar))){
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
      multiInitVar <- array(aaply(multiInitVar,1,findHcNAs,pHcMod,pCROBAS,HcModV),dim=c(nSites,7,1))
    }else{
      multiInitVar <- aaply(multiInitVar,1,findHcNAs,pHcMod,pCROBAS,HcModV)
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
  

  # ###set PRELES parameters
  # if(layerPRELES==0){
  #   if(!is.vector(pPRELES)) stop("check pPRELES parameters, it should be a vector")
  #   pPRELES <- matrix(pPRELES, nrow =length(pPRELES), ncol=ncol(pCROBAS))
  # }
  # if(layerPRELES==1){
  #   if(ncol(pPRELES) != ncol(pCROBAS)) stop("check consistency in species column between pPRELES and pCROBAS")
  # }
  # 
  # ###use the most common species to calculate P0
  # if(maxNlayers > 1){
  #   domSp <- multiInitVar[,1,][matrix(c(1:nSites,apply(multiInitVar[,5,],1,which.max)),ncol=2)]
  # }else{
  #   domSp <- multiInitVar[,1,1]
  # }
  # domSp <- as.numeric(names(which.max(table(domSp))))
  # pPRELESx <- pPRELES[,domSp]
  # 
  ### compute P0
  ###if P0 is not provided use preles to compute P0
  if(all(is.na(multiP0))){
    multiP0 <- array(NA,dim=c(nClimID,maxYears,2))
    for(climID in 1:nClimID){
      nYearsX <- max(nYearsMS[which(climIDs==climID)])
      P0 <- PRELES(DOY=rep(1:365,nYearsX),PAR=PAR[climID,1:(365*nYearsX)],
                   TAir=TAir[climID,1:(365*nYearsX)],VPD=VPD[climID,1:(365*nYearsX)],
                   Precip=Precip[climID,1:(365*nYearsX)],CO2=CO2[climID,1:(365*nYearsX)],
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
  
  ###initialize siteType
  multiOut[,,3,,1] <- array(siteInfo[,3],dim=c(nSites,maxYears,maxNlayers))
  ###initialize alfar
  for(ijj in 1:maxNlayers){
    siteXs <- which(multiInitVar[,1,ijj] %in% 1:ncol(pCROBAS))
    multiOut[siteXs,,3,ijj,2] =
      matrix(pCROBAS[cbind((20+pmin(siteInfo[,3],5))[siteXs],
              multiInitVar[siteXs,1,ijj])],
             length(siteXs),maxYears)
  }
  
if(alpharNcalc){
  ###initialize alfar
  if(all(is.na(p0currClim))) p0currClim <- rowMeans(multiP0[,1:min(maxYears,10),1])
  p0ratio <- multiP0[,,1]/p0currClim
  if(all(is.na(TcurrClim))) TcurrClim <- apply(weatherYasso[,1:min(10,maxYears),1],1,mean)
  if(all(is.na(PcurrClim))) PcurrClim <- apply(weatherYasso[,1:min(10,maxYears),2],1,mean)
  fT0 <- fTfun(TcurrClim,PcurrClim)
  fT <- fTfun(weatherYasso[,,1],weatherYasso[,,2])
  fTratio <- fT/fT0 
  alpharNfact <- p0ratio * fTratio
  
  for(ijj in 1:nClimID){
      # siteXs <- which(siteInfo[,2] == ijj)
      siteXs <- which(siteInfo[,2]==ijj)
      if(length(siteXs)==1 & maxNlayers==1) multiOut[siteXs,,3,,2] <- multiOut[siteXs,,3,,2] * alpharNfact[ijj,]
      if(length(siteXs)==1 & maxNlayers>1) multiOut[siteXs,,3,,2] <- sweep(multiOut[siteXs,,3,,2],1,alpharNfact[ijj,],FUN="*") 
      if(length(siteXs)>1) multiOut[siteXs,,3,,2] <- sweep(multiOut[siteXs,,3,,2],2,alpharNfact[ijj,],FUN="*") 
  }
  ####alphar is smoothed using a running average of 10 years
  if(maxNlayers==1) multiOut[,1,3,,2] <-  apply(multiOut[,1:10,3,1,2],1,mean)
  if(maxNlayers>1) multiOut[,1,3,,2] <-  apply(multiOut[,1:10,3,,2],c(1,3),mean)
  for(ijj in 2:maxYears){
    multiOut[,ijj,3,,2] <- multiOut[,(ijj-1),3,,2] + (multiOut[,ijj,3,,2] - multiOut[,(ijj-1),3,,2])/10
  }
} 
  
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
    mortMod = mortMod,
    ECMmod = ECMmod,
    pECMmod = pECMmod,
    layerPRELES = layerPRELES,
    LUEtrees = LUEtrees,
    LUEgv = LUEgv,
    alpharNcalc=alpharNcalc
  )
  return(multiSiteInit)
}

#' Title
#'
#' @param multiSiteInit output from the InitMultiSite function: initialized prebas inputs
#' @param fertThin flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 
#' @param nYearsFert number of years after thinnings for which the fertilization is effective. default values is 20 years
#' @param oldLayer flag for retention trees after clearcut (randomly 5-10 percent basal area is left after clearcut)
#'
#' @return
#' #'  soilC Initial soil carbon compartments for each layer. Array with dimentions = c(nSites, nYears,5,3,nLayers). 1st dim = # of sites in the simulations; 2nd dim = # of simulation years; The 3rd dimention (5) corresponds to the AWENH pools; the 4th dimention (3) corresponds to the tree organs (foliage, branch and stem). \cr
#'   \cr
#'  soilCtot matrix of stand total annual soilcarbon per year. rows are # of sites; columns corresponds to # of years \cr
#'   \cr
#'  multiOut  An array with annual model outputs. 1st dimension corresponds to # of sites; 2nd dimension corresponds to the number of years of the simulation (nYears); 3nd dimension corresponds to the output variables (see list below); 4th dimension corresponds to the number of layers in the stand (nLayers); 5th dimensions reports the state of the stand (1) and (2) the variables of the harvested trees (2). \cr
#' Output variables: \cr
#' 1."siteID" \cr
#' 2."gammaC" internal parameter \cr
#' 3."sitetype" site fertility class \cr
#' 4."species" \cr
#' 5."ETS" effective temperature sums \cr
#' 6."P0" Potential annual gross primary production (gC m-2 y-1) \cr
#' 7."age" Age of the layer (years) \cr
#' 8."DeadWoodVolume" Dead wood volume (m3 ha-1) \cr
#' 9."Respi_tot" Autotrophic respiration (gC m-2 y-1) \cr
#' 10."GPP/1000" Total tree GPP  (kgC m-2 y-1) \cr
#' 11."H" Layer average height (m) \cr
#' 12."D" Layer average diameter at breast height (cm) \cr
#' 13."BA" Layer basal area (m-2 ha-1) \cr
#' 14."Hc_base" Layer Base of crown height (m) \cr
#' 15."Cw" Crown width (m) \cr
#' 16."A" Sapwood area of average tree at crown base (m2) \cr
#' 17."N" Layer density \cr
#' 18."npp" net primary production (gC m-2 y-1) \cr
#' 19."leff" Effective leaf area \cr
#' 20."keff" Effective light extintion coefficient \cr
#' 21."lproj" Projected leaf area \cr
#' 22."ET_preles" Annual evapotranspiration (mm y-1) \cr
#' 23."weight" Layer weight on photosynthesis \cr
#' 24."Wbranch" Living Branch biomass (kgC ha-1) \cr
#' 25."WfineRoots" Fine roots biomass (kgC ha-1) \cr
#' 26."Litter_fol" Foliage litter (kgC ha-1) \cr
#' 27."Litter_fr" Fine root litter (kgC ha-1) \cr
#' 28."Litter_fWoody" fine woody litter (kgC ha-1)\cr 
#' 29."Litter_cWood" coarse woody litter (kgC ha-1) \cr
#' 30."V" Layer volume (m3 ha-1) \cr
#' 31."Wstem" Stem Biomass (kgC ha-1) \cr
#' 32."W_croot" Course root Biomass (kgC ha-1) \cr
#' 33."wf_STKG" Foliage biomass (kgC ha-1) \cr
#' 34."wf_treeKG" Foliage biomass of the average tree (kgC ha-1) \cr
#' 35."B_tree" Basal area of average tree (m2) \cr
#' 36."Light" The proportion of light that has not been intercepted by canopy (meanlight)  \cr
#' 37."VroundWood" harvested round wood volume (m3 ha-1) \cr
#' 38."WroundWood" haversted round wood biomass (kgC ha-1) \cr
#' 39."soilC" totaal soil carbon (kgC ha-1) \cr
#' 40."aSW" average available soil water (mm) \cr
#' 41."dH" height growth (m) \cr
#' 42."Vmort" volume of dead trees (m3 ha-1) \cr
#' 43."grossGrowth" gross growth (m3 ha-1 y-1) \cr
#' 44."GPPtrees" Gross primary production per tree layer (gC m-2 y-1) \cr
#' 45."Rh" heterotrophic respiration (gC m-2 y-1) \cr
#' 46."NEP" Net ecosystem exchange (gC m-2 y-1), note 1st layer include fluxes from ground vegetation \cr
#' 47." W_wsap" sapwood biomass (kgC ha-1) \cr
#' 48."W_c" sapwood stem below Crown (kgC ha-1) \cr
#' 49."W_s" sapwood stem within crown (kgC ha-1) \cr
#' 50."Wsh" biomass of stem heartwood  (kgC ha-1) \cr
#' 51."Wdb" biomass of dead branches on living trees (kgC ha-1) \cr
#' 52."dHc" Height of the crown base change (m) \cr
#' 53."Wbh" biomass of branches heartwood (kgC ha-1) \cr
#' 54."Wcrh" biomass of coarse root heartwood (kgC ha-1)\cr
#' \cr
#'  multiEnergyWood  An array with annual energywood harvested. \cr
#'  1st dimension corresponds to the number of site of the simulation (nSites); \cr
#'  2nd dimension corresponds to the number of years of the simulation (nYears); \cr
#'  3rd dimension corresponds to the layers;\cr
#'  4th dimension has 2 elements that correspond to volume and biomasses respectively\cr
#' 
#' \cr
#'  dailyPRELES  An array with PRELES output: \cr
#'  1st dimension corresponds to the number of climIDs (site weather conditions defined in siteInfo, note that multiple stands can share the same climate) ; \cr
#'  2nd dimension corresponds to the days ; \cr
#'  3rd dimension corresponds to the PRELES outputs: GPP, ET and SW. \cr
#' \cr
#'  GVout  A array with the ground vegetation model output: \cr
#'  1st dimension corresponds to the number of site of the simulation (nSites); \cr
#'  2nd dimension corresponds to the years ; \cr
#'  3rd dimension corresponds to the GV outputs: fAPAR_gv, litGV, photoGV, Wgv,GVnpp. \cr
#' @export
#'
#' @examples
multiPrebas <- function(multiSiteInit,
                        fertThin = 0,
                        nYearsFert = 20,
                        oldLayer=0){
  
  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }

  ####avoid species = 0  replace with species 1 when layer is empty
  multiSiteInit$multiInitVar[,1,][which(multiSiteInit$multiInitVar[,1,]==0)] <- 1
  multiSiteInit$multiOut[,,4,,1][which(multiSiteInit$multiOut[,,4,,1]==0)] = 1

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
                     pPRELES=as.matrix(multiSiteInit$pPRELES),
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
                     mortMod=as.double(multiSiteInit$mortMod),
                     ECMmod=as.integer(multiSiteInit$ECMmod),
                     pECMmod=as.double(multiSiteInit$pECMmod),
                     layerPRELES = as.integer(multiSiteInit$layerPRELES),
                     LUEtrees = as.double(multiSiteInit$LUEtrees),
                     LUEgv = as.double(multiSiteInit$LUEgv)
  )
  dimnames(prebas$multiOut) <- dimnames(multiSiteInit$multiOut)
  dimnames(prebas$multiInitVar) <- dimnames(multiSiteInit$multiInitVar)
  names(prebas$siteInfo) <- names(multiSiteInit$siteInfo)
  prebas$alpharNcalc = multiSiteInit$alpharNcalc
  
  class(prebas) <- "multiPrebas"
  return(prebas)
}


#' Title
#'
#' @param multiSiteInit 
#' @param HarvLim 
#' @param minDharv 
#' @param cutAreas 
#' @param compHarv 
#' @param thinFact 
#' @param ageHarvPrior 
#' @param siteOrder 
#' @param fertThin 
#' @param nYearsFert 
#' @param oldLayer 
#' @param startSimYear 
#'
#' @return
#' @export
#'
#' @examples
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
  # if(length(startSimYear)==1) startSimYear <- rep(startSimYear,multiSiteInit$nSites)
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
  

  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }

  ####avoid species = 0  replace with species 1 when layer is empty
  multiSiteInit$multiInitVar[,1,][which(multiSiteInit$multiInitVar[,1,]==0)] <- 1
  multiSiteInit$multiOut[,,4,,1][which(multiSiteInit$multiOut[,,4,,1]==0)] = 1

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
                     pPRELES=as.matrix(multiSiteInit$pPRELES),
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
                     mortMod=as.double(multiSiteInit$mortMod),
                     startSimYear = as.integer(startSimYear),
                     ECMmod=as.integer(multiSiteInit$ECMmod),
                     pECMmod=as.double(multiSiteInit$pECMmod),
                     layerPRELES = as.integer(multiSiteInit$layerPRELES),
                     LUEtrees = as.double(multiSiteInit$LUEtrees),
                     LUEgv = as.double(multiSiteInit$LUEgv)
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
  prebas$alpharNcalc = multiSiteInit$alpharNcalc
  return(prebas)
}




#' Title
#'
#' @param multiSiteInit 
#' @param HarvLim 
#' @param minDharv 
#' @param cutAreas 
#' @param compHarv 
#' @param thinFact 
#' @param ageHarvPrior 
#' @param siteOrder 
#' @param fertThin 
#' @param nYearsFert 
#' @param oldLayer 
#' @param startSimYear 
#'
#' @return
#' @export
#'
#' @examples
reStartRegionPrebas <- function(multiSiteInit,
                         HarvLim = NA,
                         minDharv = 999.,
                         cutAreas = NA,  ### is a matrix: area of cuttings rows are years of simulations
                         ###columns: clcutArea target(1), simulated clCut area(2) (set to 0. will be filled by prebas output);
                         ####precom-thin target(3), sim(4); area firstThin targ(5), sim(6)
                         compHarv=0,###flag for compensating harvest if harvest do not reach the desired levels
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
  
  # if(length(startSimYear)==1) startSimYear <- rep(startSimYear,multiSiteInit$nSites)
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
  
  # if(FALSE){
  if(oldLayer==1){
    multiSiteInit <- addOldLayer(multiSiteInit)
  }
  ####avoid species = 0  replace with species 1 when layer is empty
  multiSiteInit$multiInitVar[,1,][which(multiSiteInit$multiInitVar[,1,]==0)] <- 1
  multiSiteInit$multiOut[,,4,,1][which(multiSiteInit$multiOut[,,4,,1]==0)] = 1
  
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
                     pPRELES=as.matrix(multiSiteInit$pPRELES),
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
                     mortMod=as.double(multiSiteInit$mortMod),
                     startSimYear = as.integer(startSimYear),
                     ECMmod=as.integer(multiSiteInit$ECMmod),
                     pECMmod=as.double(multiSiteInit$pECMmod),
                     layerPRELES = as.integer(multiSiteInit$layerPRELES),
                     LUEtrees = as.double(multiSiteInit$LUEtrees),
                     LUEgv = as.double(multiSiteInit$LUEgv)
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


