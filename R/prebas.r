
prebas <- function(nYears,
                   pCROBAS = pCROB,
                   pPRELES = pPREL,
                   PREBASversion = 0,
                   etmodel = 0,
                   pYASSO = pYAS,
                   pAWEN = parsAWEN,
                   siteInfo = NA,
                   thinning=NA,
                   initClearcut = c(1.5,0.5,0.0431969,0.,0.),
                   fixBAinitClarcut = 1.,
                   initCLcutRatio = NA,
                   PAR,TAir,VPD,Precip,CO2,
                   P0=NA,
                   initVar = NA,
                   soilC = NA,
                   weatherYasso = NA,
                   litterSize = NA,
                   soilCtot = NA,
                   defaultThin = 1.,
                   ClCut = 1.,
                   inDclct = NA,
                   inAclct = NA,
                   yassoRun = 0){

  ###process weather###
  if(length(PAR) >= (nYears*365)){
    PAR = PAR[1:(nYears*365)]
    TAir = TAir[1:(nYears*365)]
    VPD = VPD[1:(nYears*365)]
    Precip = Precip[1:(nYears*365)]
    CO2 = CO2[1:(nYears*365)]
  } else{
    stop("daily weather inputs < nYears*365")
  }
  ###

  ###proc thinnings##
  if(all(is.na(thinning))) thinning=matrix(0,1,8)
  thinning[is.na(thinning)] <- -999
  nThinning = max(1,nrow(thinning))
  thinning <- thinning[order(thinning[,2],thinning[,1],thinning[,3]),]
  ###
  if(all(is.na(initVar))) {
    nLayers <- 3 }else {
    nLayers <- ifelse(is.null(ncol(initVar)),1,ncol(initVar))
  }
  nSp = ncol(pCROBAS)
  if(anyNA(siteInfo)) siteInfo = c(1,1,3,160,0,0,20) ###default values for nspecies and site type = 3

  if(all(is.na(initCLcutRatio))){
    initCLcutRatio <- rep(1/nLayers,nLayers)
  }

  varNam <- getVarNam()
  nVar <- length(varNam)

  layerNam <- paste("layer",1:nLayers)
  output <- array(0, dim=c((nYears),nVar,nLayers,2),
        dimnames = list(NULL,varNam,layerNam,c("stand","thinned")))
  fAPAR <- rep(0.7,nYears)

  ###compute ETS year
  Temp <- TAir[1:(365*nYears)]-5
  ETS <- pmax(0,Temp,na.rm=T)
  ETS <- matrix(ETS,365,nYears); ETS <- colSums(ETS)

  ###if P0 is not provided use preles to compute P0
  if(is.na(P0)){
    P0 <- PRELES(DOY=rep(1:365,nYears),
                 PAR=PAR,TAir=TAir,VPD=VPD,Precip=Precip,CO2=CO2,
                 fAPAR=rep(1,length(PAR)),LOGFLAG=0,p=pPRELES)$GPP
    P0 <- matrix(P0,365,nYears);P0 <- colSums(P0)
  }

  ETSthres <- 1000; ETSmean <- mean(ETS)

  ####process clearcut
  if(any(!is.na(c(inDclct,inAclct)))){
    if(is.na(inDclct)) inDclct <- 9999999.99
    if(is.na(inAclct)) inAclct <- 9999999.99
  }
  # if(ClCut==1 & all(is.na(initVar)) & is.na(inDclct)) inDclct <-
  if(ClCut==1 & all(is.na(inDclct))) inDclct <-
    c(ClCutD_Pine(ETSmean,ETSthres,siteInfo[3]),
      ClCutD_Spruce(ETSmean,ETSthres,siteInfo[3]),
      ClCutD_Birch(ETSmean,ETSthres,siteInfo[3]))
  # if(ClCut==1 & all(is.na(initVar)) & is.na(inAclct)) inAclct <-
  if(ClCut==1 & all(is.na(inAclct))) inAclct <-
    c(ClCutA_Pine(ETSmean,ETSthres,siteInfo[3]),
      ClCutA_Spruce(ETSmean,ETSthres,siteInfo[3]),
      ClCutA_Birch(ETSmean,ETSthres,siteInfo[3]))
  if(any(is.na(inDclct))) inDclct[is.na(inDclct)] <- 9999999.99
  if(length(inDclct)==1) inDclct<- rep(inDclct,nSp)
  if(any(is.na(inAclct))) inAclct[is.na(inAclct)] <- 9999999.99
  if(length(inAclct)==1) inAclct<- rep(inAclct,nSp)

###if any initial value is given the model is initialized from plantation
  if (all(is.na(initVar))){
    initVar <- matrix(NA,6,nLayers)
    initVar[1,] <- 1:nLayers
    initVar[3,] <- initClearcut[1]; initVar[4,] <- initClearcut[2]
    initVar[5,] <- initClearcut[3]/nLayers; initVar[6,] <- initClearcut[4]
  }

  xx <- min(10,nYears)
  Ainit = 6 + 2*3.5 - 0.005*(sum(ETS[1:xx])/xx) + 2.25
  initVar[2,which(is.na(initVar[2,]))] <- initClearcut[5] <- round(Ainit)

  ####process weather PRELES (!!to check 365/366 days per year)
  weatherPreles <- array(c(PAR,TAir,VPD,Precip,CO2),dim=c(365,nYears,5))
  weatherPreles <- aperm(weatherPreles, c(2,1,3))

  ###initialise soil inputs
  if(all(is.na(soilCtot))) soilCtot = numeric(nYears)
  if(all(is.na(soilC))) soilC = array(0,dim = c(nYears,5,3,nLayers))
  if(all(is.na(litterSize))){
    litterSize = matrix(0,3,nLayers)
    litterSize[2,] <- 2
    for (i in 1:nLayers) litterSize[1,i] <- ifelse(initVar[1,i]==3,10,30)
  }

##process weather inputs for YASSO
  if(all(is.na(weatherYasso))){
    weatherYasso = matrix(0,nYears,3)
    weatherYasso[,1] = aTmean(TAir,nYears)
    weatherYasso[,3] = aTampl(TAir,nYears)
    weatherYasso[,2] = aPrecip(Precip,nYears)
  }

  PREBASversion <- paste("prebas_v",PREBASversion,sep='')

  prebas <- .Fortran(PREBASversion,
                     nYears=as.integer(nYears),
                     nLayers=as.integer(nLayers),
                     nSp=as.integer(nSp),
                     siteInfo = as.numeric(siteInfo),
                     pCROBAS = as.matrix(pCROBAS),
                     initVar=as.matrix(initVar),
                     thinning=as.matrix(thinning),
                     output=as.array(output),
                     nThinning=as.integer(nThinning),
                     maxYearSite=as.integer(nYears),
                     fAPAR=as.numeric(fAPAR),
                     initClearcut=as.numeric(initClearcut),
                     fixBAinitClarcut=as.numeric(fixBAinitClarcut),
                     initCLcutRatio = as.double(initCLcutRatio),
                     ETS = as.numeric(ETS),
                     P0 = as.numeric(P0),
                     weather=as.array(weatherPreles),
                     DOY= as.integer(1:365),
                     pPRELES=as.numeric(pPRELES),
                     etmodel = as.integer(etmodel),
                     soilC = as.array(soilC),
                     pYASSO=as.numeric(pYASSO),
                     pAWEN = as.matrix(pAWEN),
                     weatherYasso = as.matrix(weatherYasso),
                     litterSize = as.matrix(litterSize),
                     soilCtot=as.numeric(soilCtot),
                     defaultThin=as.double(defaultThin),
                     ClCut=as.double(ClCut),
                     inDclct=as.double(inDclct),
                     inAclct=as.double(inAclct),
                     dailyPRELES = matrix(-999,(nYears*365),3),
                     yassoRun=as.double(yassoRun))
  class(prebas) <- "prebas"
  return(prebas)
}

