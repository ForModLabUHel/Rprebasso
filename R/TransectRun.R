#' @title Transect runs for pine, spruce, or mixed forests
#' @description The simulations of pine or spruce or mixed forests at 7 selected points over a North-South transect. Start from seedlings!
#'
#' @param SiteType  A value between 1 to 5.
#' @param initVar   initial state of the forest. Array with nSites,VarsIn,nLayers dimentions: nSites=number of sites (7); VarsIn = initial state input variables (speciesID, age(if NA is calculated by initialAgeSeedl function),h,dbh,ba(by layer),hc,Ac(NA))
#' @param species   If the initial state of the forest is not provided (initVar) Should be 'Pine', or 'Spruce', or 'Birch' or 'Mixed'
#' @param nYears    Number of years to run the model.Default value is 100.
#' @param pPRELES   A vector of PRELES parameter. Default is the boreal forest version pPREL. 
#' @param pCROBAS   (47 x nSpecies matrix) Matrix of parameter sets, each column corresponds to a species. Default values pCROBAS = pCROB are the parameter sets for Scots pine (Pinus sylvestris), Norway spruce (Picea abies), Silver birch (Betula pendula), European beech (Fagus sylvatica), Maritime pine (Pinus pinaster), Blue gum (Eucalyptus globulus), Black locust (Robinia pseudoacacia), Populus(in Romania), Eucalyptus grandis x Eucalyptus urophylla (in Paraguay), and Norway spruce(in Germany). Default is pCROB, print(pCROB) to see the parameter values and names.
#' @param AgeEffect Carbon allocation strategy for seedlings. When True, the early growth will decrease due to age effect,by setting pCROBAS[45,]<-0.3.
#' @param defaultThin If defaultThin = 1 (default) Finnish standard management practices are applied (ref).
#' @param multiThin 	A array with thinnig inputs. Three dimensions: nSites x maxThin x 9. The first dimension is the number of sites. The second dimension, maxThin, is the number of thinnings. For the third demention, element 1 is year from the start of the simulation; element 2 is siteID; element 3 layer where thinnings are carried out; element 4 to 7 stand variables (H, D, B, Hc); element 8 parameter that indicates if the stand variables (column 4:7) are provided as fraction of the actual model outputs (value=1 means that fraction is used); element 9 is the stand density after thinning if its value is not -999.
#' @param ClCut     Vector of Diameter (cm) threshold for clearcut. Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
#' @param energyCut Energy cutting strategy will be applied if set to 1. Default value is 0.
#' @param mortMod flag for the mortality model selection (1= Reineke, 2= random (Siilipehto, 2020), 3= both models)
#' @param ECMmod flag for the ECM modelling activation 1 -> model ECM according to Makela et al. 2022, 0 -> no ECM modelling
#' @param multiInitClearCut A Matrix: matrix(initClearcut,NoOfSites,5,byrow = T), where initClearcut includes those 5 variables H,dbh,BA,HC,AC, same with 'initSeedling.def'
#' @param multiNthin A matrix with thinning inputs. Rows correspond to a thinning event. Column 1 year from the start of the simulation; column 2 is siteID; column 3 layer where thinnings are carried out; column 4 to 7 stand variables (H, D, B, Hc); column 8 parameter that indicates if the stand variables (column 4:7) are provided as fraction of the actual model outputs (value=1 means that fraction is used); column 9 is the stand density after thinning if its value is not -999; colum 10 is Sapwood area of average tree at crown base (m2) if its value is not -999 (see examples).
#' @param GVrun 
#' @param pHcMod 
#' @param etmodel 
#' @param pYASSO 
#' @param pAWEN 
#' @param fixAinit = NA,
#' @param fixBAinitClearcut 
#' @param initCLcutRatio 
#' @param multiP0 
#' @param soilC 
#' @param weatherInput ##list of weather inputs for PRELES: each variables is a matrix with nrow=nSites,columns=number of days in the simulations
#' @param weatherYasso 
#' @param litterSize 
#' @param soilCtot 
#' @param inDclct 
#' @param inAclct 
#' @param inHclct 
#' @param yassoRun 
#' @param smoothP0 
#' @param smoothETS 
#' @param smoothYear 
#' @param HcModV flag for the Hc model: 1 use the pipe model defined in the HcPipeMod function, different from 1 uses empirical models; default value (HcModV_def) is 1
#' @param tapioPars 
#' @param thdPer 
#' @param limPer 
#' @param ftTapioPar 
#' @param tTapioPar 
#' @param thinInt 
#' @param layerPRELES 
#' @param LUEtrees light use efficiency parameters for tree species
#' @param LUEgv light use efficiency parameter for ground vegetation
#' @param alpharNcalc #alphar calculations based on Nitrogen availability. deafault value is FALSE (no nitrogen impact). =1calculates N uptake
#' @param p0currClim # vector of average annual P0 for the climIDs at current climate. if NA the first five years of the simulations will be used to calculate it.
#' @param TcurrClim # vector of average annual temperature for the climIDs at current climate. if NA the first five years of the simulations will be used to calculate it.
#' @param PcurrClim # vector of average annual precipitation for the climIDs current climate. if NA the first five years of the simulations will be used to calculate it.
#' @param ingrowth # flag to simulate ingrowth
#' @param soilPar # input a matrix (dim=nSites,3 ) with soil depth, FC, WP, for each site if NA uses the default values
#' @param modVersion # model version to use in the simulations it can be multiSite or region
#' @param nYearsFert number of years after thinnings for which the fertilization is effective. default values is 20 years
#' @param yearFert simulation year when fertilization occurred
#' @param deltaSiteTypeFert fertilization impact
#' @param fertThin flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 
#' @param oldLayer flag for retention trees after clearcut (randomly 5-10 percent basal area is left after clearcut)
#' @param latitude latitude of the site
#' @param TsumSBBs initial temperature sums for bark beetle risk for the two years before the first year if not available it will be calculated using the first year
#' @param SMIt0 site vector of initial SoilMoirture index
#' @param TminTmax array(climaIDs,ndays,2) with daily Tmin Tmax values for each climID, Tmin and Tmax will be used to calculate the Nesterov Index that will be used in the fire risk calculations  
#' @param soilC_steadyState flag for soilC at steady state calculations. if true the soilC at st st is calculated with the average litterfall of the simulations and soilC balance is computed for each year
#' @param disturbanceON flag for activating disturbance modules. can be one of "wind", "fire",  "bb" or a combination of the three, ex. c("fire", "bb") 
#' @param CO2model CO2 model for PRELES. Default CO2model = 1 (Launaniemi) ; CO2model = 2 (Kolari) 
#' @param lightnings used in fire disturbance module. is the frequency of lightning-caused ignition events (ha-1 d-1) used in the fire module it is a matrix of dimensions nSites,ndays 
#' @param popden used in fire disturbance module. It is the population density (individuals km-2). it is a matrix of dimensions nSites,ndays 
#' @param a_nd used in fire disturbance module. a(ND) is a parameter expressing the propensity of people to produce ignition events (ignitions individual-1 d-1). site specific parameter. vector of lenght nSites
#' @param NIout flag to return the nesterov index
#' @param FDIout flag to return the fire danger index instead of SW daily preles, set to 1 to return the FDI
#' 
#' @importFrom plyr aaply
#'
#' @return The output from multiPrebas()
#' 
#' @export
#'
#' @examples 
#' # Qucik examples
#' Pine.1 <- TransectRun(SiteType = 1, species = "Pine")
#' plot(Pine.1$multiOut[1, , 30, 1, 1], type = "l")
#' 
#' #' An example of initail sates
#' InitStands<-matrix(c(1,NA,1.5,0.5,0.0431969,0.,0.),nrow = 7,ncol=7,byrow = T)
#' Pine.i <-TransectRun(SiteType = 1, initVar = InitStands)
#' plot(Pine.i$multiOut[1, , 30, 1, 1], type = "l")
#' 
#' #An example of thinning
#' thinnings<-array(0,c(7,2,10))
#' thinnings[,,1]<- c(rep(20,7),rep(40,7)) #year from the start of the simulation
#' thinnings[,,2]<- 1:7 #siteID
#' thinnings[,,3]<- 1 #layer where thinnings are carried out
#' thinnings[,,4]<- 1.1 #H
#' thinnings[,,5]<- 1.1#D
#' thinnings[,,6]<- 0.5#B
#' thinnings[,,7]<- 1#Hc
#' thinnings[,,8]<- 1 #if the stand variables (column 4:7) are provided as fraction of the actual model outputs
#' thinnings[,,9:10] <- -999 #9 the stand density after thinning if its value is not -999.
#' Pine.thin <-TransectRun(SiteType = 1, initVar = InitStands,multiThin = thinnings,multiNthin = rep(2,7))
#' plot(Pine.thin$multiOut[1, , 30, 1, 1], type = "l")
#' 
#' #An example of new settings for seedlings after a clear cut
#' initClearcut<-c(1.5100000,0.8000000,0.0631969,0.2000000,NA) # H,dbh,BA,HC,AC
#' afterClearCut <- matrix(initClearcut,7,5,byrow = T)
#' Pine.newseed <-TransectRun(SiteType = 1, initVar = InitStands,multiThin = thinnings,multiNthin = rep(2,7),multiInitClearCut = afterClearCut)
#' plot(Pine.newseed$multiOut[1, , 30, 1, 1], type = "l")
#' 


TransectRun <- function(SiteType = NA, initVar = NA, species = NA, nYears = 100, pCROBAS = pCROB,
                        pPRELES= NA,AgeEffect = F, defaultThin = 1,multiThin = NA, multiNthin = NA, ClCut = 1, energyCut = 0,
                        GVrun=1,
                        pHcMod = pHcM,
                        etmodel = 0,
                        pYASSO =pYAS,
                        pAWEN = parsAWEN,
                        multiInitClearCut = NA,
                        fixAinit = 0, ###fix initial year age (vector of length # sites) 0,initial age is calculated by the model, other wise use years 
                        fixBAinitClearcut = 1.,  ###if 1 when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relBa at last year 
                        initCLcutRatio = NA,  ###BA ratio per each species/layer (default is the ba ratio at the begginning of the simulations)
                        multiP0=NA,
                        soilC = NA,
                        weatherInput=NULL,
                        weatherYasso = NA,
                        litterSize = litterSizeDef,
                        soilCtot = NA,
                        inDclct = NA,
                        inAclct = NA,
                        inHclct = NA,
                        yassoRun = 1,
                        smoothP0 = 1,
                        smoothETS = 1,
                        smoothYear=5,
                        HcModV=HcModV_def,  ####version of model to compute Hc 1 uses the version of based on ksi parameter 2 uses the empirical model; default value (HcModV_def) is 1
                        tapioPars=pTapio,
                        thdPer = NA,
                        limPer = NA,
                        ftTapioPar = ftTapio,
                        tTapioPar = tTapio,
                        thinInt = -999.,
                        mortMod = 1, #flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
                        ECMmod = 0,
                        layerPRELES = 0,
                        LUEtrees = pLUEtrees,
                        LUEgv = pLUEgv,
                        alpharNcalc=FALSE,
                        p0currClim = NA,
                        TcurrClim = NA,
                        PcurrClim = NA,
                        ingrowth = FALSE,
                        soilPar = NA, #### input a matrix with soil depth, FC, WP, for each site if NA uses the default values
                        siteInfoDist = NA,
                        modVersion = "multiSite",
                        nYearsFert = 20,
                        yearFert=NULL,
                        deltaSiteTypeFert = 1,
                        fertThin=0.,
                        oldLayer=0,
                        latitude = c(60.295,60.959,61.377,62.647,64.441,66.143,68.203),
                        TsumSBBs = matrix(-999.,7,4),
                        SMIt0 = rep(-999,7),
                        TminTmax = NA,
                        soilC_steadyState=FALSE,
                        disturbanceON = NA,
                        CO2model = 2,
                        lightnings = NA,
                        popden = NA,
                        a_nd = NA,
                        NIout = F,
                        FDIout = 0
) {
  
  if(!CO2model %in% 1:2) stop(paste0("set CO2model 1 or 2"))
  if(all(is.na(pPRELES))){
    pPRELES <- pPREL
    pPRELES[12:13] <- pCO2model[CO2model,]
  }
  
  if(all(!is.na(soilC))){
    if(soilC_steadyState) warning("soilC at steady state was not computed because initial soilC was inputed")
    soilC_steadyState = FALSE
    yassoRun = 1
  }
  if(soilC_steadyState) yassoRun = 1
  if(nrow(pCROBAS)!=53) stop("check that pCROBAS has 53 parameters, see pCROB to compare")
  if(!modVersion %in% c("multiSite","region")) stop("modVersion must be region or multiSite")
  
  nSites <- 7
  siteInfo <- matrix(c(NA, NA, NA, 160, 0, 0, 20, 3, 3, 413, 0.45, 0.118,3), nSites, 13, byrow = T)
  if (all(is.na(SiteType))) {
    SiteType <- 3
    warning("siteType 3 was assigned to all sites since SiteType was not provided")
  }
  siteInfo[, 3] <- SiteType
  siteInfo[, 2] <- siteInfo[, 1] <- 1:nSites
  
  if (is.na(species) & all(is.na(initVar))) {
    species <- "Pine"
    warning("Species was assigned to 'Pine' since species was not provided")
  }
  
  if (!all(is.na(initVar))) {
    siteInfo[, 8:9] <- 1
    if(length(dim(initVar))>2) siteInfo[, 8:9] <- dim(initVar)[3]
  }
  
  if (species == "Pine" & all(is.na(initVar))) {
    initVar <- aperm(array(c(1, NA, initSeedling.def), dim = c(7, 7, 1)), c(2, 1, 3))
    siteInfo[, 8:9] <- 1 # even-aged pure forests
  }
  
  if (species == "Spruce" & all(is.na(initVar))) {
    initVar <- aperm(array(c(2, NA, initSeedling.def), dim = c(7, 7, 1)), c(2, 1, 3))
    siteInfo[, 8:9] <- 1 # even-aged pure forests
  }
  
  if (species == "Birch" & all(is.na(initVar))) {
    initVar <- aperm(array(c(3, NA, initSeedling.def), dim = c(7, 7, 1)), c(2, 1, 3))
    siteInfo[, 8:9] <- 1 # even-aged pure forests
  }
  
  if(!all(is.na(soilPar))) siteInfo[,10:12] <- soilPar 
  
  if (species == "Mixed" & all(is.na(initVar))) {
    initVar <- array(NA, dim = c(7, 7, 3))
    initVar[, , 1] <- matrix(c(1, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, , 2] <- matrix(c(2, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, , 3] <- matrix(c(3, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, 5, 1:2] <- initVar[, 5, 1:2] * 0.45 ## 45% basal area is pine and 45% spruce
    initVar[, 5, 3] <- initVar[, 5, 3] * 0.1 ## 10% basal area is birch
  }
  
  if (AgeEffect == T) {
    pCROBAS[45, ] <- 0.3
  }
  if(is.null(weatherInput)){
    path.weatherTran <-
      system.file("transect", "weatherTran.rda", package = "Rprebasso")
    load(file = path.weatherTran)
    nRep <- ceiling(nYears / (dim(PARtran)[2] / 365))
    PAR = do.call(cbind, replicate(nRep, PARtran, simplify = FALSE))
    TAir = do.call(cbind, replicate(nRep, TAirtran, simplify = FALSE))
    VPD = do.call(cbind, replicate(nRep, VPDtran, simplify = FALSE))
    Precip = do.call(cbind, replicate(nRep, Preciptran, simplify = FALSE))
    CO2 = do.call(cbind, replicate(nRep, CO2tran, simplify = FALSE))
  }else{
    PAR = weatherInput$PAR
    TAir = weatherInput$TAir
    VPD = weatherInput$VPD
    Precip = weatherInput$Precip
    CO2 = weatherInput$CO2
  }
  
  initPrebas <- InitMultiSite(
    nYearsMS = rep(nYears, nSites),
    siteInfo = siteInfo,
    pCROBAS = pCROBAS,
    pPRELES = pPRELES,
    multiInitVar = as.array(initVar),
    PAR = PAR,
    TAir = TAir,
    VPD = VPD,
    Precip = Precip,
    CO2 = CO2,
    # litterSize = litterSize,
    # pAWEN = parsAWEN,
    defaultThin = defaultThin,
    multiThin = multiThin,
    multiNthin = multiNthin,
    ClCut = ClCut,
    yassoRun = yassoRun,
    energyCut = energyCut,
    GVrun=GVrun,
    pHcMod = pHcM,
    etmodel = etmodel,
    pYASSO =pYASSO,
    pAWEN = pAWEN,
    multiInitClearCut = multiInitClearCut,
    fixAinit = fixAinit,
    fixBAinitClearcut = fixBAinitClearcut,  ###if 1 when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relBa at last year 
    initCLcutRatio = initCLcutRatio,  ###BA ratio per each species/layer (default is the ba ratio at the begginning of the simulations)
    multiP0=multiP0,
    soilC = soilC,
    weatherYasso = weatherYasso,
    litterSize = litterSize,
    soilCtot = soilCtot,
    inDclct = inDclct,
    inAclct = inAclct,
    inHclct = inHclct,
    smoothP0 = smoothP0,
    smoothETS = smoothETS,
    smoothYear=smoothYear,
    HcModV=HcModV,  ####version of model to compute Hc 1 uses the version of based on ksi parameter 2 uses the empirical model
    tapioPars=tapioPars,
    thdPer = thdPer,
    limPer = limPer,
    ftTapioPar = ftTapioPar,
    tTapioPar = tTapioPar,
    thinInt = thinInt,
    mortMod = mortMod,#flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
    ECMmod=ECMmod,
    layerPRELES=layerPRELES,
    LUEtrees = LUEtrees,
    LUEgv = LUEgv,
    alpharNcalc=alpharNcalc,
    p0currClim = p0currClim,
    TcurrClim = TcurrClim,
    PcurrClim = PcurrClim,
    ingrowth = ingrowth,
    siteInfoDist = siteInfoDist,
    latitude = latitude,
    TsumSBBs = TsumSBBs,
    SMIt0 = SMIt0,
    TminTmax = TminTmax,
    disturbanceON = disturbanceON,
    CO2model=CO2model,
    lightnings = lightnings,
    popden = popden,
    a_nd = a_nd,
    NIout = NIout,
    FDIout = FDIout
    )

  initPrebas$multiInitVar[, 2, ] <- initialAgeSeedl(initPrebas$siteInfo[, 3], rowMeans(initPrebas$ETS)) # Initial age

  
  if(modVersion=="region"){
    TransectOut <- regionPrebas(initPrebas,
                                fertThin = fertThin,
                                nYearsFert = nYearsFert,
                                yearFert=yearFert,
                                deltaSiteTypeFert = deltaSiteTypeFert,
                                oldLayer=oldLayer)
  } 
  if(modVersion=="multiSite"){
    TransectOut <- multiPrebas(initPrebas,
                               fertThin = fertThin,
                               nYearsFert = nYearsFert,
                               yearFert=yearFert,
                               deltaSiteTypeFert = deltaSiteTypeFert,
                               oldLayer=oldLayer)
  }  
  
  if(soilC_steadyState){
    stst_soilC <- stXX_GV(TransectOut,GVrun = 1)
    TransectOut <-yassoPREBASin(prebOut=TransectOut,initSoilC=stst_soilC)
    TransectOut$stst_soilC <- stst_soilC
  }
  return(TransectOut)
}
