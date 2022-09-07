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
#' @param ClCut     Vector of Diameter (cm) threshold for clearcut. Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
#' @param energyCut Energy cutting strategy will be applied if set to 1. Default value is 0.
#' @param mortMod flag for the mortality model selection (1= Reineke, 2= random (Siilipehto, 2020), 3= both models)
#' @param ECMmod flag for the ECM modelling activation 1 -> model ECM according to Makela et al. 2022, 0 -> no ECM modelling
#'
#' @importFrom plyr aaply
#'
#' @return The output from multiPrebas()
#' @export
#'
#' @examples Pine.1 <- TransectRun(SiteType = 1, species = "Pine")
#' plot(Pine.1$multiOut[1, , 30, 1, 1], type = "l")

#' 
TransectRun <- function(SiteType = NA, initVar = NA, species = NA, nYears = 100, pCROBAS = pCROB,
                        pPRELES= pPREL,AgeEffect = F, defaultThin = 1, ClCut = 1, energyCut = 0,
                        mortMod=1,GVrun=1) {
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
    defaultThin = defaultThin,
    multiInitVar = as.array(initVar),
    PAR = do.call(cbind, replicate(nRep, PARtran, simplify = FALSE)),
    TAir = do.call(cbind, replicate(nRep, TAirtran, simplify = FALSE)),
    VPD = do.call(cbind, replicate(nRep, VPDtran, simplify = FALSE)),
    Precip = do.call(cbind, replicate(nRep, Preciptran, simplify = FALSE)),
    CO2 = do.call(cbind, replicate(nRep, CO2tran, simplify = FALSE)),
    # litterSize = litterSize,
    # pAWEN = parsAWEN,
    ClCut = ClCut,
    yassoRun = 1,
    energyCut = energyCut,
    mortMod=mortMod,
    GVrun=GVrun
  )
  initPrebas$multiInitVar[, 2, ] <- initialAgeSeedl(initPrebas$siteInfo[, 3], rowMeans(initPrebas$ETS)) # Initial age
  TransectOut <- multiPrebas(initPrebas)
  return(TransectOut)
}
