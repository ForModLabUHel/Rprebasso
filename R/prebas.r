#' @title prebas model function
#' @description function to run the PREBAS model for a single site!
#'
#' @param nYears number of years of the simulations
#' @param pCROBAS (49 x nSpecies(11) matrix) Matrix of parameter sets of CROBAS growth model, each column corresponds to a species. Default values pCROBAS = pCROB are the parameter sets for Scots pine (Pinus sylvestris), Norway spruce (Picea abies), Silver birch (Betula pendula), European beech (Fagus sylvatica), Maritime pine (Pinus pinaster), Blue gum (Eucalyptus globulus), Black locust (Robinia pseudoacacia), Populus(in Romania), Eucalyptus grandis x Eucalyptus urophylla (in Paraguay), and Norway spruce(in Germany), Quercus Ilex. Default is pCROB, print(pCROB) to see the parameter values and names.
#' @param pHcMod  (7 x nSpecies(11) matrix) Matrix of parameter sets Height of the crown base model.
#' @param pPRELES A vector of PRELES parameter. Default is the boreal forest version pPREL.  
#' @param pYASSO A vector of YASSO parameters. Default is global pYAS
#' @param pAWEN (12 x nSpecies(11) matrix) Matrix of parameter sets for litterfal decomposition in YASSO pools. Default is parsAWEN
#' @param etmodel Evapotranspiration model for PRELES. Default etmodel = 0. Possible values -1, 0, 1, 2
#' @param siteInfo vector of site info SiteID, climID, siteType, SWinit (initial soil water), CWinit (initial crown water), SOGinit (initial snow on ground), Sinit (initial temperature acclimation state), soildepth, effective field capacity, permanent wilthing point. Default = c(1,1,3,160,0,0,20,413.,0.45,0.118), i.e. siteType = 3. 
#' @param thinning A matrix with thinnig inputs. Rows correspond to a thinning event. Column 1 year from the start of the simulation; column 2 is siteID; column 3 layer where thinnings are carried out; column 4 to 7 stand variables (H, D, B, Hc); column 8 parameter that indicates if the stand variables (column 4:7) are provided as fraction of the actual model outputs (value=1 means that fraction is used); column 9 is the stand density after thinning if its value is not -999; colum 10 is Sapwood area of average tree at crown base (m2) if its value is not -999 (see examples).
#' @param initClearcut A numeric vector with initial stand variables after clearcut: H, D, BA, Hc, Ainit. Ainit is the year when the stand reaches measurable size. If NA the default values from initSeedling.def are used Ainit and is automatically computed using air temperature.
#' @param fixBAinitClarcut If 1, when clearcuts occur the species initial biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relative basal area at last year 
#' @param initCLcutRatio vector (of the length of the number of layers) with basal area fraction (to the total basal area)  at replanting it fixBAinitClarcut == 1
#' @param PAR A numeric vector of daily sums of photosynthetically active radiation, mmol/m2
#' @param TAir A numeric vector of daily mean temperature, degrees C.
#' @param VPD A numeric vector of daily mean vapour pressure deficits, kPa.
#' @param Precip A numeric vector of daily rainfall, mm
#' @param CO2 A numeric vector of air CO2, ppm
#' @param P0 A numeric vector with the annual potential photosynthesis (gC m-2 y-1). If P0 is not provided PRELES is used to compute P0 using fAPAR = 1.
#' @param initVar initial state of the forest. matrix with VarsIn,nLayers dimentions: VarsIn = initial state input variables (speciesID, age(if NA is calculated by initialAgeSeedl function),h,dbh,ba(by layer),hc,Ac(NA))
#' @param soilC Initial soil carbon compartments for each layer. Array with dimentions = c(nYears,5,3,nLayers). The second dimention (5) corresponds to the AWENH pools; the third dimention (3) corresponds to the tree organs (foliage, branch and stem).
#' @param weatherYasso Annual weather inputs for Yasso, matrix with nYears and 3 columns. If NA it is internally calculated using the daily weather inputs
#' @param litterSize Marix with litter inputs for YASSO. Rows are tree organs, columns correspond to the species of pCROBAS.
#' @param soilCtot total annual soilcarbon per year
#' @param defaultThin If defaultThin = 1 (default) Finnish standard managment practices are applied (ref). 
#' @param ClCut If ClCut = 1 clearcuts are applied. If inDclct = NA and inAclct = NA Finnish standard clearcut practices are applied (ref).
#' @param energyCut if==1 energy wood is harvested at thinning and clearcut
#' @param inDclct Vector of Diameter (cm) threshold for clearcut. Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
#' @param inAclct Vector of Age (year) threshold for clearcut.  Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
#' @param yassoRun flag for YASSO calculations. If yassoRun=1 the YASSO model is run to compute the carbon balance of the soil.
#' @param smoothP0 if 1 P0 is smoothed using an average window of smoothYear (see below) # years 
#' @param smoothETS if 1 ETS is smoothed using an average window of smoothYear (see below) # years 
#' @param smoothYear # of years for smoothing P0 and ETS
#' @param tapioPars parameters for the tapio rules applied at harvests in Finalnd. tapioPars(sitetype, conif/decid, south/center/north, thinning parameters), and parameters for modifying thinnig limits and thresholds
#' @param thdPer value varying between 0 and 1 (default=0.5) defining the intensities of the thinnings based of Finnish default rules. 
#' @param limPer value varying between 0 and 1 (default=0.5) defining the BASAL area threshold used for thinning implementation according to Finnish rules.
#' @param ftTapioPar   first thinning parameters. 
#' @param tTapioPar  Tending  thinning parameters. 
#' @param GVrun flag for Ground vegetation model 1-> runs the GV model
#' @param thinInt parameter that determines the thinning intensity; from below (thinInt>1) or above (thinInt<1);
#' @param fertThin flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 
#' @param nYearsFert number of years after thinnings for which the fertilization is effective 
#' @param protect flag for retention trees after clearcut (randomly 5-10 percent basal area is left after clearcut)
#' @param mortMod flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
#' @param ECMmod flag for activation of the ectomycoryzal model (ECMmod=1). see Makela et al. 2022
#' @param pECMmod parameters of ECM model 
#' @param layerPRELES flag indicating if preles is going to be run by layer with species specific parameters or using 1 parameter set for the whole forest
#' @param LUEtrees light use efficiency parameters for tree species
#' @param LUEgv light use efficiency parameter for ground vegetation
#' @param aplharNcalc #alphar calculations based on Nitrogen availability. deafault value is FALSE (no nitrogen impact). =1calculates N uptake
#'
#' @return
#'  soilC Initial soil carbon compartments for each layer. Array with dimentions = c(nYears,5,3,nLayers). The second dimention (5) corresponds to the AWENH pools; the third dimention (3) corresponds to the tree organs (foliage, branch and stem). \cr
#'   \cr
#'  soilCtot stand total annual soilcarbon per year. \cr
#'   \cr
#'  output  An array with annual model outputs. 1st dimension corresponds to the number of years of the simulation (nYears); 2nd dimension corresponds to the output variables (see list below); 3rd dimension corresponds to the number of layers in the stand (nLayers); 4th dimensions reports the state of the stand (1) and (2) the variables of the harvested trees (2). \cr
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
#'  energyWood  An array with annual energywood harvested. \cr
#'  1st dimension corresponds to the number of years of the simulation (nYears); \cr
#'  2nd dimension corresponds to the layers;\cr
#'  3rd dimension has 2 elements that correspond to volume and biomasses respectively\cr
#' 
#' \cr
#'  dailyPRELES  A matrix with PRELES output: \cr
#'  1st dimension corresponds to the days ; \cr
#'  2nd dimension corresponds to the PRELES outputs: GPP, ET and SW. \cr
#' \cr
#'  GVout  A matrix with the ground vegetation model output: \cr
#'  1st dimension corresponds to the years ; \cr
#'  2nd dimension corresponds to the GV outputs: fAPAR_gv, litGV, photoGV, Wgv,GVnpp. \cr
#' 
#' @export
#'
#' @examples
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
                   ECMmod=0, #flag for ECM modelling MAkela et al.2022
                   pECMmod=parsECMmod,
                   layerPRELES = 0,
                   LUEtrees = pLUEtrees,
                   LUEgv = pLUEgv
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
  
  ###set PRELES parameters
  # if(layerPRELES==0){
  #   if(!is.vector(pPRELES)) stop("check pPRELES parameters, it should be a vector")
  #   pPRELES <- matrix(pPRELES, nrow =length(pPRELES), ncol=ncol(pCROBAS))
  # }
    
  ###if P0 is not provided use preles to compute P0
  # domSp <- initVar[1,which.max(initVar[5,])]
  # pPRELESx <- pPRELES[,domSp]
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
  Ainit = max((6 + 2*siteInfo[3] - 0.005*(sum(ETS[1:xx])/xx) + 2.25 +2),2)
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
  
  if(aplharNcalc){
    ###initialize alfar
    p0currClim <- mean(P0[1:min(maxYears,5),1])
    p0ratio <- P0[,1]/p0currClim
    T0 <- mean(weatherYasso[1:min(5,maxYears),1])
    precip0 <- mean(weatherYasso[1:min(5,maxYears),2])
    fT0 <- fTfun(T0,precip0)
    fT <- fTfun(weatherYasso[,1],weatherYasso[,2])
    fTratio <- fT/fT0 
    alpharNfact <- p0ratio * fTratio
    if(maxNlayers==1) output[,3,,2] <- output[,3,,2] * alpharNfact
    if(maxNlayers>1) output[,3,,2] <- sweep(output[,3,,2],1,alpharNfact,FUN="*") 
    }
  } 
  
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
                     ECMmod = as.integer(ECMmod),
                     pECMmod = as.numeric(pECMmod),
                     layerPRELES = as.integer(layerPRELES),
                     LUEtrees = as.double(LUEtrees),
                     LUEgv = as.double(LUEgv)
                     )
  class(prebas) <- "prebas"
  return(prebas)
}

