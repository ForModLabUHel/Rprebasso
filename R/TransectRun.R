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
#' 
#' @importFrom plyr aaply
#'
#' @return The output from multiPrebas()
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
                        pPRELES= pPREL,AgeEffect = F, defaultThin = 1,multiThin = NA, multiNthin = NA, ClCut = 1, energyCut = 0,
                        GVrun=1,
                        pHcMod = pHcM,
                        etmodel = 0,
                        pYASSO =pYAS,
                        pAWEN = parsAWEN,
                        multiInitClearCut = NA,
                        fixBAinitClarcut = 1.,  ###if 1 when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relBa at last year 
                        initCLcutRatio = NA,  ###BA ratio per each species/layer (default is the ba ratio at the begginning of the simulations)
                        multiP0=NA,
                        soilC = NA,
                        weatherYasso = NA,
                        litterSize = litterSizeDef,
                        soilCtot = NA,
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
                        thinInt = -999.,
                        mortMod = 1 #flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
) {
  nSites <- 7
  siteInfo <- matrix(c(NA, NA, NA, 160, 0, 0, 20, 3, 3, 413, 0.45, 0.118), nSites, 12, byrow = T)
  if (is.na(SiteType)) {
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

  if (species == "Mixed" & all(is.na(initVar))) {
    initVar <- array(NA, dim = c(7, 7, 3))
    initVar[, , 1] <- matrix(c(1, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, , 2] <- matrix(c(2, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, , 3] <- matrix(c(3, NA, initSeedling.def), 7, 7, byrow = T)
    initVar[, 5, 1:2] <- initVar[, 5, 1:2] * 0.45 ## 45% basal area is pine and 45% spruce
    initVar[, 5, 3] <- initVar[, 5, 3] * 0.1 ## 10% basal area is birch
  }

  path.par <-
    system.file("transect", "PARtran.csv", package = "Rprebasso")
  path.tair <-
    system.file("transect", "TAirtran.csv", package = "Rprebasso")
  path.precip <-
    system.file("transect", "Preciptran.csv", package = "Rprebasso")
  path.vpd <-
    system.file("transect", "VPDtran.csv", package = "Rprebasso")
  path.co2 <-
    system.file("transect", "CO2tran.csv", package = "Rprebasso")

  PARtran <- as.matrix(read.csv(
    file = path.par,
    header = T,
    row.names = 1
  ))

  Preciptran <- as.matrix(
    read.csv(
      file = path.precip,
      header = T,
      row.names = 1
    )
  )
  VPDtran <- as.matrix(
    read.csv(
      file = path.vpd,
      header = T,
      row.names = 1
    )
  )
  TAirtran <- as.matrix(
    read.csv(
      file = path.tair,
      header = T,
      row.names = 1
    )
  )
  CO2tran <- as.matrix(
    read.csv(
      file = path.co2,
      header = T,
      row.names = 1
    )
  )

  nRep <- ceiling(nYears / (dim(PARtran)[2] / 365))

  if (AgeEffect == T) {
    pCROBAS[45, ] <- 0.3
  }

  initPrebas <- InitMultiSite(
    nYearsMS = rep(nYears, nSites),
    siteInfo = siteInfo,
    pCROBAS = pCROBAS,
    pPRELES = pPRELES,
    multiInitVar = as.array(initVar),
    PAR = do.call(cbind, replicate(nRep, PARtran, simplify = FALSE)),
    TAir = do.call(cbind, replicate(nRep, TAirtran, simplify = FALSE)),
    VPD = do.call(cbind, replicate(nRep, VPDtran, simplify = FALSE)),
    Precip = do.call(cbind, replicate(nRep, Preciptran, simplify = FALSE)),
    CO2 = do.call(cbind, replicate(nRep, CO2tran, simplify = FALSE)),
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
    fixBAinitClarcut = fixBAinitClarcut,  ###if 1 when clearcut occur the species inital biomass is fixed at replanting using the values in initCLcutRatio else at replanting the replanting follows species relBa at last year 
    initCLcutRatio = initCLcutRatio,  ###BA ratio per each species/layer (default is the ba ratio at the begginning of the simulations)
    multiP0=multiP0,
    soilC = soilC,
    weatherYasso = weatherYasso,
    litterSize = litterSize,
    soilCtot = soilCtot,
    inDclct = inDclct,
    inAclct = inAclct,
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
    mortMod = mortMod #flag for mortality model selection 1= reineke model; 2: random mort mod based on Siilipehto et al.2020; 3 = both models
  )
  initPrebas$multiInitVar[, 2, ] <- initialAgeSeedl(initPrebas$siteInfo[, 3], rowMeans(initPrebas$ETS)) # Initial age
  TransectOut <- multiPrebas(initPrebas)
  return(TransectOut)
}
