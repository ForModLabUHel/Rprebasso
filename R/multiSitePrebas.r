InitMultiSite <- function(nYearsMS,
                          pCROBAS = pCROB,
                          pPRELES = pPREL,
                          PREBASversion = 0,
                          etmodel = 0,
                          pYASSO =pYAS,
                          pAWEN = parsAWEN,
                          siteInfo = NA,
                          multiInitVar = NA,
                          multiThin = NA,
                          multiNthin = NA,
                          multiInitClearCut = NA,
                          fixBAinitClarcut = 1.,
                          initCLcutRatio = NA,
                          PAR,
                          TAir,
                          VPD,
                          Precip,
                          CO2,
                          multiP0=NA,
                          soilC = NA,
                          weatherYasso = NA,
                          litterSize = NA,
                          soilCtot = NA,
                          defaultThin = 1.,
                          ClCut = 1.,
                          inDclct = NA,
                          inAclct = NA,
                          yassoRun = 0,
                          lukeRuns){

  nSites <- length(nYearsMS)
  if(all(is.na(siteInfo))){
    siteInfo = matrix(c(1,1,3,160,0,0,20,3,3),nSites,9,byrow = T) ###default values for nspecies and site type = 3
    siteInfo[,1] <- 1:nSites
  }
  nLayers <- siteInfo[,8]
  if(length(fixBAinitClarcut)==1) fixBAinitClarcut=rep(fixBAinitClarcut,nSites)
  if(all(is.na(initCLcutRatio))){
    initCLcutRatio <- matrix(0.,nSites,max(nLayers))
    for(iz in 1:nSites) initCLcutRatio[iz,1:nLayers[iz]] <- rep(1/nLayers[iz],nLayers[iz])
  }
  # nSp <- siteInfo[,9]
  climIDs <- siteInfo[,2]
  # if(all(is.na(multiInitVar)) & all(is.na(nSp)) nSp <- rep(3,nSites)
  allSp = ncol(pCROBAS)
  varNam <- getVarNam()
  nVar <- length(varNam)

  nClimID <- length(unique(climIDs))
  if(!all((1:nClimID) %in% climIDs) | length(climIDs) != nSites) return("check consistency between weather inputs and climIDs")

  maxYears <- max(nYearsMS)
  maxNlayers <- max(nLayers)
  layerNam <- paste("layer",1:maxNlayers)
  multiOut <- array(0, dim=c(nSites,(maxYears),nVar,maxNlayers,2),
                    dimnames = list(NULL,NULL,varNam,layerNam,
                                    c("stand","thinned")))
  initClearcut = c(1.5,0.5,0.0431969,0.,NA)
  if (all(is.na(multiInitClearCut))) multiInitClearCut <- matrix(initClearcut,nSites,5,byrow = T)

  ###process yasso inputs if missing
  if(is.na(soilC)) soilC <- array(0,dim=c(nSites,maxYears,5,3,maxNlayers))
  if(is.na(soilCtot)) soilCtot <- matrix(0,nSites,maxYears)

  ##process weather inputs for YASSO
  if(all(is.na(weatherYasso))){
    weatherYasso <- array(0,dim=c(nClimID,maxYears,3))
    weatherYasso[,,1] <- t(apply(TAir[,1:(maxYears*365)],1,aTmean,maxYears))
    weatherYasso[,,3] <- t(apply(TAir[,1:(maxYears*365)],1,aTampl,maxYears))
    weatherYasso[,,2] <- t(apply(Precip[,1:(maxYears*365)],1,aPrecip,maxYears))
  }

  if (length(defaultThin) == 1) defaultThin=as.double(rep(defaultThin,nSites))
  if (length(ClCut) == 1) ClCut=as.double(rep(ClCut,nSites))
  if (length(inDclct) == 1) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == 1) inAclct=matrix(inAclct,nSites,allSp)
  if (length(inDclct) == nSites) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == nSites) inAclct=matrix(inAclct,nSites,allSp)
  if (length(yassoRun) == 1) yassoRun=as.double(rep(yassoRun,nSites))
  if (length(PREBASversion) == 1) PREBASversion=as.double(rep(PREBASversion,nSites))

  ###process ETS
  multiETS <- matrix(NA,nClimID,maxYears)
  for(climID in 1:nClimID){
    nYearsX <- max(nYearsMS[which(climIDs==climID)])
    Temp <- TAir[climID,1:(365*nYearsX)]-5
    ETS <- pmax(0,Temp,na.rm=T)
    ETS <- matrix(ETS,365,nYearsX); ETS <- colSums(ETS)
    multiETS[climID,(1:nYearsX)] <- ETS

    xx <- min(10,nYearsX)
    Ainit = 6 + 2*3.5 - 0.005*mean(ETS[1:xx]) + 2.25 ## this is not dependent to site type? check with Annikki
    sitesClimID <- which(climIDs==climID)
    multiInitClearCut[sitesClimID,5] <- replace(multiInitClearCut[sitesClimID,5],
                                                which(is.na(multiInitClearCut[sitesClimID,5])),round(Ainit))
  }
  ETSthres <- 1000; ETSmean <- rowMeans(multiETS)

  ####process clearcut
  for(i in 1: nSites){
    if(ClCut[i]==1 & all(is.na(inDclct[i,]))) inDclct[i,] <-
        c(ClCutD_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]))
    if(ClCut[i]==1 & all(is.na(inAclct[i,]))) inAclct[i,] <-
        c(ClCutA_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]))
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
    multiThin <- array(0, dim=c(nSites,maxThin,8))
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

    multiweather[i,(1:nYearsMS[i]),,] <- weatherPreles
  }

  ### compute P0
  ###if P0 is not provided use preles to compute P0
  if(all(is.na(multiP0))){
    multiP0 <- matrix(NA,nClimID,maxYears)
    for(climID in 1:nClimID){
      nYearsX <- max(nYearsMS[which(climIDs==climID)])
      P0 <- PRELES(DOY=rep(1:365,nYearsX),PAR=PAR[climID,1:(365*nYearsX)],
                   TAir=TAir[climID,1:(365*nYearsX)],VPD=VPD[climID,1:(365*nYearsX)],
                   Precip=Precip[climID,1:(365*nYearsX)],CO2=rep(380,(365*nYearsX)),
                   fAPAR=rep(1,(365*nYearsX)),LOGFLAG=0,p=pPRELES)$GPP
      P0 <- matrix(P0,365,nYearsX)
      multiP0[climID,(1:nYearsX)] <- colSums(P0)
    }}

  if (all(is.na(multiInitVar))){
    multiInitVar <- array(NA,dim=c(nSites,6,maxNlayers))
    multiInitVar[,1,] <- rep(1:maxNlayers,each=nSites)
    multiInitVar[,3,] <- initClearcut[1]; multiInitVar[,4,] <- initClearcut[2]
    multiInitVar[,5,] <- initClearcut[3]/maxNlayers; multiInitVar[,6,] <- initClearcut[4]
    multiInitVar[,2,] <- matrix(multiInitClearCut[,5],nSites,maxNlayers)
  }
  if(all(is.na(litterSize))){
    litterSize <- matrix(0,3,allSp)
    litterSize[2,] <- 2
    litterSize[1,] <- c(30,30,10)
    # siteInfo <- siteInfo[,-c(4,5)]
  }

  multiSiteInit <- list(
    multiOut = multiOut,
    nSites = nSites,
    nClimID = nClimID,
    nLayers = nLayers,
    maxYears = maxYears,
    maxThin = maxThin,
    nYears = nYearsMS,
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
    inDclct = inDclct,
    inAclct = inAclct,
    dailyPRELES = array(-999,dim=c(nSites,(maxYears*365),3)),
    yassoRun = yassoRun,
    lukeRuns = lukeRuns,
    PREBASversion = PREBASversion)
  return(multiSiteInit)
}

multiPrebas <- function(multiSiteInit){
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
                     siteInfo = as.matrix(multiSiteInit$siteInfo[,1:7]),  ####
                     maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                     nThinning=as.integer(multiSiteInit$nThinning),
                     fAPAR=as.matrix(multiSiteInit$fAPAR),
                     initClearcut=as.matrix(multiSiteInit$initClearcut),
                     fixBAinitClearcut = as.double(multiSiteInit$fixBAinitClarcut),
                     initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                     ETSy=as.matrix(multiSiteInit$ETSy),
                     P0y=as.matrix(multiSiteInit$P0y),
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
                     inDclct=as.matrix(multiSiteInit$inDclct),
                     inAclct=as.matrix(multiSiteInit$inAclct),
                     dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                     yassoRun=as.double(multiSiteInit$yassoRun),
                     PREBASversion=as.double(multiSiteInit$PREBASversion),
                     lukeRuns=as.double(multiSiteInit$lukeRuns))
  class(prebas) <- "multiPrebas"
  return(prebas)
}


regionPrebas <- function(multiSiteInit,
                        HarvLim = NA,
                        minDharv = 15){

  if(length(HarvLim)==1) HarvLim <- rep(HarvLim,multiSiteInit$maxYears)
  if(any(is.na(HarvLim))) HarvLim[which(is.na(HarvLim))] <- 0.

  siteOrder <- matrix(1:multiSiteInit$nSites,multiSiteInit$nSites,multiSiteInit$maxYears)
  siteOrder <- apply(siteOrder,2,sample,multiSiteInit$nSites)

  prebas <- .Fortran("regionPrebas",
                   siteOrder = as.matrix(siteOrder),
                   HarvLim = as.double(HarvLim),
                   minDharv = as.double(minDharv),
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
                   siteInfo = as.matrix(multiSiteInit$siteInfo[,1:7]),  ####
                   maxNlayers = as.integer(multiSiteInit$maxNlayers), ####
                   nThinning=as.integer(multiSiteInit$nThinning),
                   fAPAR=as.matrix(multiSiteInit$fAPAR),
                   initClearcut=as.matrix(multiSiteInit$initClearcut),
                   fixBAinitClarcut = as.double(multiSiteInit$fixBAinitClarcut),
                   initCLcutRatio = as.matrix(multiSiteInit$initCLcutRatio),
                   ETSy=as.matrix(multiSiteInit$ETSy),
                   P0y=as.matrix(multiSiteInit$P0y),
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
                   inDclct=as.matrix(multiSiteInit$inDclct),
                   inAclct=as.matrix(multiSiteInit$inAclct),
                   dailyPRELES = as.array(multiSiteInit$dailyPRELES),
                   yassoRun=as.double(multiSiteInit$yassoRun),
                   PREBASversion=as.double(multiSiteInit$PREBASversion),
                   lukeRuns=as.double(multiSiteInit$lukeRuns))
class(prebas) <- "regionPrebas"
if(prebas$maxNlayers>1){
    prebas$totHarv <- apply(prebas$multiOut[,,37,,1],2,sum)
  }else{
    prebas$totHarv <- prebas$multiOut[,,37,,1]
  }
return(prebas)
}



