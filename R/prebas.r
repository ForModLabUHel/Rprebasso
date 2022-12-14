prebas <- function(nYears,
                   pCROBAS = pCROB,
                   pHcMod = pHcM,
                   pPRELES = pPREL,
                   pYASSO = pYAS,
                   pAWEN = parsAWEN,
                   # PREBASversion = 0,
                   etmodel = 0,
                   siteInfo = NA,
                   thinning=NA,
                   initClearcut = initSeedling.def,
                   fixBAinitClarcut = 1.,
                   initCLcutRatio = NA,
                   PAR,TAir,VPD,Precip,CO2,
                   P0=NA,
                   initVar = NA,
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
                   tapioPars=pTapio,
                   thdPer=0.5,
                   limPer=0.5,
                   ftTapioPar = ftTapio,
                   tTapioPar = tTapio,
                   GVrun = 1, ###flag for Ground vegetation model 1-> runs the GV model
                   thinInt=-999.,
                   fertThin=0.,
                   nYearsFert=20,
                   protect=0,
                   mortMod=1,
                   ECMmod=0 #flag for ECM modelling MAkela et al.2022
                   ){
  
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
  if(all(is.na(thinning))){
    thinning=matrix(0,1,10)
    thinning[,9:10] <- -999
  } 
  thinning[is.na(thinning)] <- -999
  nThinning = max(1,nrow(thinning))
  thinning <- thinning[order(thinning[,2],thinning[,1],thinning[,3]),]
  ###
  if(all(is.na(initVar))) {
    nLayers <- 3 }else {
      nLayers <- ifelse(is.null(ncol(initVar)),1,ncol(initVar))
    }
  nSp = ncol(pCROBAS)
  if(anyNA(siteInfo)) siteInfo = c(1,1,3,160,0,0,20,413.,0.45,0.118) ###default values for nspecies and site type = 3
                                  
  if(all(is.na(initCLcutRatio))){
    initCLcutRatio <- rep(1/nLayers,nLayers)
  }
  
  varNam <- getVarNam()
  nVar <- length(varNam)
  
  layerNam <- paste("layer",1:nLayers)
  output <- array(0, dim=c((nYears),nVar,nLayers,2),
                  dimnames = list(year=NULL,variable=varNam,layer=layerNam,status=c("stand","thinned")))
  energyWood <- array(0, dim=c((nYears),nLayers,2),
                      dimnames = list(year=NULL,layer=layerNam,variable=c("volume","biomass")))
  fAPAR <- rep(0.7,nYears)
  
  ###compute ETS year
  Temp <- TAir[1:(365*nYears)]-5
  ETS <- pmax(0,Temp,na.rm=T)
  ETS <- matrix(ETS,365,nYears); ETS <- colSums(ETS)
  if(smoothETS==1. & nYears>1) for(i in 2:nYears) ETS[i] <- ETS[(i-1)] + (ETS[i]-ETS[(i-1)])/min(i,smoothYear)
  
  ###if P0 is not provided use preles to compute P0
  if(is.na(P0)){
    P0 <- PRELES(DOY=rep(1:365,nYears),
                 PAR=PAR,TAir=TAir,VPD=VPD,Precip=Precip,CO2=CO2,
                 fAPAR=rep(1,length(PAR)),LOGFLAG=0,p=pPRELES)$GPP
    P0 <- matrix(P0,365,nYears);P0 <- colSums(P0)
  }
  P0 <- matrix(P0,nYears,2)
  if(smoothP0==1.& nYears>1){
    P0[1,2] <- P0[1,1]
    for(i in 2:nYears) P0[i,2] <- P0[(i-1),2] + (P0[i,1]-P0[(i-1),2])/min(i,smoothYear)
  } 
  
  ETSthres <- 1000; ETSmean <- mean(ETS)
  
  ####process clearcut
  if(any(!is.na(c(inDclct,inAclct)))){
    if(is.na(inDclct)) inDclct <- 9999999.99
    if(is.na(inAclct)) inAclct <- 9999999.99
  }
  # if(ClCut==1 & all(is.na(initVar)) & is.na(inDclct)) inDclct <-
  if(ClCut==1 & all(is.na(inDclct))) inDclct <-
    c(ClCutD_Pine(ETSmean,ETSthres,siteInfo[3]), ####pine in Finland
      ClCutD_Spruce(ETSmean,ETSthres,siteInfo[3]), ####spruce in Finland
      ClCutD_Birch(ETSmean,ETSthres,siteInfo[3]), ####birch in Finland
      NA,NA,NA,NA)  ###"fasy","pipi","eugl","rops"
  # if(ClCut==1 & all(is.na(initVar)) & is.na(inAclct)) inAclct <-
  if(ClCut==1 & all(is.na(inAclct))) inAclct <-
    c(ClCutA_Pine(ETSmean,ETSthres,siteInfo[3]),   ####pine in Finland
      ClCutA_Spruce(ETSmean,ETSthres,siteInfo[3]), ####spruce in Finland 
      ClCutA_Birch(ETSmean,ETSthres,siteInfo[3]),  ####birch in Finland
      80,50,13,30)  ###"fasy","pipi","eugl","rops"  
  if(any(is.na(inDclct))) inDclct[is.na(inDclct)] <- 9999999.99
  if(length(inDclct)==1) inDclct<- rep(inDclct,nSp)
  if(any(is.na(inAclct))) inAclct[is.na(inAclct)] <- 9999999.99
  if(length(inAclct)==1) inAclct<- rep(inAclct,nSp)
  
  ###if any initial value is given the model is initialized from plantation
  if (all(is.na(initVar))){
    initVar <- matrix(NA,7,nLayers)
    initVar[1,] <- 1:nLayers
    initVar[3,] <- initClearcut[1]; initVar[4,] <- initClearcut[2]
    initVar[5,] <- initClearcut[3]/nLayers; initVar[6,] <- initClearcut[4]
  }
  
  
  ####if Height of the crown base is not available use model
  initVar <- findHcNAs(initVar,pHcMod)
  # initialize A
  for(ikj in 1:nLayers){
    p_ksi=pCROBAS[38,initVar[1,ikj]]
    p_rhof <- pCROBAS[15,initVar[1,ikj]]
    p_z <- pCROBAS[11,initVar[1,ikj]]
    Lc <- initVar[3,ikj] - initVar[6,ikj]
    A <- p_ksi/p_rhof * Lc^p_z
    initVar[7,ikj] <- A     
  } 
  
  # print(initVar)
  xx <- min(10,nYears)
  Ainit = 6 + 2*siteInfo[3] - 0.005*(sum(ETS[1:xx])/xx) + 2.25
  if(is.na(initClearcut[5])) initClearcut[5] <- Ainit
  if(length(initVar[2,which(is.na(initVar[2,]))])>0){
     initVar[2,which(is.na(initVar[2,]))] <- as.numeric(round(Ainit))
    }
  # print(initVar)
  
  ####process weather PRELES (!!to check 365/366 days per year)
  weatherPreles <- array(c(PAR,TAir,VPD,Precip,CO2),dim=c(365,nYears,5))
  weatherPreles <- aperm(weatherPreles, c(2,1,3))
  
  ###initialise soil inputs
  if(all(is.na(soilCtot))) soilCtot = numeric(nYears)
  if(all(is.na(soilC))) soilC = array(0,dim = c(nYears,5,3,nLayers))
  # if(all(is.na(litterSize))){
  #   litterSize = matrix(0,3,nLayers)
  #   litterSize[2,] <- 2
  #   for (i in 1:nLayers) litterSize[1,i] <- ifelse(initVar[1,i]==3,10,30)
  # }
  
  ##process weather inputs for YASSO
  if(all(is.na(weatherYasso))){
    weatherYasso = matrix(0,nYears,3)
    weatherYasso[,1] = aTmean(TAir,nYears)
    weatherYasso[,3] = aTampl(TAir,nYears)
    weatherYasso[,2] = aPrecip(Precip,nYears)
  }
  
  ###init biomasses
  if(nLayers==1)initVarX <- matrix(c(initVar,siteInfo[3]),8,1)
  if(nLayers>1) initVarX <- rbind(initVar,siteInfo[3])
  biomasses <- initBiomasses(pCROBAS,as.matrix(initVarX))
  biomasses[which(is.na(biomasses))] <- 0.
  output[1,c(33,25,47:49,24,32,50,51,31,30,54),,1] <- biomasses
  # print(biomasses)
  initVar <- as.matrix(initVar[1:7,])
  # PREBASversion <- paste("prebas_v",PREBASversion,sep='')
  
  ###initialize siteType and alfar parameter
  output[,3,,1] <- siteInfo[3]
  for(ijj in 1:nLayers) output[,3,ijj,2] = pCROBAS[(20+min(siteInfo[3],5)),initVar[1,ijj]]
  
  prebas <- .Fortran("prebas",
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
                     P0 = as.matrix(P0),
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
                     energyCut=as.double(energyCut),
                     inDclct=as.double(inDclct),
                     inAclct=as.double(inAclct),
                     dailyPRELES = matrix(-999,(nYears*365),3),
                     yassoRun=as.double(yassoRun),
                     energyWood = as.array(energyWood),
                     tapioPars = as.array(tapioPars),
                     thdPer = as.double(thdPer),
                     limPer = as.double(limPer),
                     ftTapioPar = as.array(ftTapioPar),
                     tTapioPar = as.array(tTapioPar),
                     GVout = matrix(0,nYears,5),
                     GVrun = as.integer(GVrun),
                     thinInt = as.double(thinInt),
                     fertThin = as.integer(fertThin),
                     flagFert = as.integer(0),
                     nYearsFert = as.integer(nYearsFert),
                     protect = as.integer(protect),
                     mortMod = as.integer(mortMod),
                     ECMmod = as.integer(ECMmod)
                     )
  class(prebas) <- "prebas"
  return(prebas)
}

