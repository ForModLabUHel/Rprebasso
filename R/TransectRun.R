#' @title Transect runs for pine, spruce, or mixed forests
#' @description The simulations of pine or spruce or mixed forests at 7 selected points over a North-South transect
#' @param SiteType  A value between 1 to 5.
#' @param species   Should be 'Pine', or 'Spruce', or 'Mixed'
#' @param nYears    Number of years to run the model.Default value is 100.
#' @param pCROBAS   (47 x nSpecies matrix) Matrix of parameter sets, each column corresponds to a species. Default values pCROBAS = pCROB are the parameter sets for Scots pine (Pinus sylvestris), Norway spruce (Picea abies), Silver birch (Betula pendula), European beech (Fagus sylvatica), Maritime pine (Pinus pinaster), Blue gum (Eucalyptus globulus), Black locust (Robinia pseudoacacia), Populus(in Romania), Eucalyptus grandis x Eucalyptus urophylla (in Paraguay), and Norway spruce(in Germany). Default is pCROB, print(pCROB) to see the parameter values and names. 
#' @param AgeEffect Carbon allocation strategy for seedlings. When True, the early growth will decrease due to age effect.
#' @param newInitSeedling Change the initial state of seedlings. Default value is False.
#' @param H0_new    Initial height of seedlings. Taking effect only when newInitSeedling == T.
#' @param Hc0_new   Initial crown base height of seedlings. Taking effect only when newInitSeedling == T.
#' @param N0_new    Initial stand density of seedlings. Taking effect only when newInitSeedling == T.
#' @param defaultThin If defaultThin = 1 (default) Finnish standard management practices are applied (ref).
#'
#' @return The output from multiPrebas()
#' @export
#'
#' @examples Pine.1<-TransectRun(SiteType = 1,species = 'Pine')
#'Pine.1.ns<-TransectRun(SiteType = 1,species = 'Pine',newInitSeedling = T)
#'Pine.1.ns.ae<-TransectRun(SiteType = 1,species = 'Pine',newInitSeedling = T,AgeEffect = T)
#'           
#'plot(Pine.1$multiOut[1,,30,1,1],type = 'l')
#'lines(Pine.1.ns$multiOut[1,,30,1,1],col='red')
#'lines(Pine.1.ns.ae$multiOut[1,,30,1,1],col='blue')
#' 
TransectRun <- function(SiteType = NA, species = NA, nYears = 100, pCROBAS = pCROB,
                         AgeEffect = F, newInitSeedling = F, H0_new=1.35, Hc0_new=0.2, N0_new=2000,
                        defaultThin = 1) {
  nSites <- 7
  siteInfo <- matrix(c(NA, NA, SiteType, 160, 0, 0, 20, 3, 3, 413, 0.45, 0.118), nSites, 12, byrow = T)
  siteInfo[, 2] <- siteInfo[, 1] <- 1:nSites

  if (species == "Pine") {
    path <-
      system.file("transect", "initVarP.csv", package = "Rprebasso")
    initVar <-
      read.csv(
        file = path,
        header = T,
        row.names = 1
      )
    initVar <-
      aperm(array(as.matrix(initVar), dim = c(7, 1, 7)), c(3, 1, 2))
    if(newInitSeedling ==T){
      initVar[,3,1]<-H0_new
      initVar[,6,1]<-Hc0_new
      initVar[,5,1]<-3.14*0.525*0.525/40000*N0_new
      }
    siteInfo[, 8:9] <- 1 # even-aged pure forests
  }

  if (species == "Spruce") {
    path <-
      system.file("transect", "initVarSp.csv", package = "Rprebasso")
    initVar <-
      read.csv(
        file = path,
        header = T,
        row.names = 1
      )
    initVar <-
      aperm(array(as.matrix(initVar), dim = c(7, 1, 7)), c(3, 1, 2))
    if(newInitSeedling ==T){
      initVar[,3,1]<-H0_new
      initVar[,6,1]<-Hc0_new
      initVar[,5,1]<-3.14*0.5*0.5/40000*N0_new
    }
    siteInfo[, 8:9] <- 1 # even-aged pure forests
  }

  if (species == "Mixed") {
    path <-
      system.file("transect", "initVarAll.csv", package = "Rprebasso")
    initVar <-
      read.csv(
        file = path,
        header = T,
        row.names = 1
      )
    initVar <-
      aperm(array(as.matrix(initVar), dim = c(7, 3, 7)), c(3, 1, 2))
    if(newInitSeedling ==T){
      initVar[,3,1]<-H0_new
      initVar[,6,1]<-Hc0_new
      initVar[,5,1]<-3.14*0.525*0.525/40000*N0_new
      initVar[,3,2]<-H0_new
      initVar[,6,2]<-Hc0_new
      initVar[,5,2]<-3.14*0.5*0.5/40000*N0_new
      initVar[,3,2]<-H0_new
      initVar[,6,2]<-Hc0_new
      initVar[,5,2]<-3.14*0.5*0.5/40000*N0_new
    }
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
  
  if(AgeEffect==T){pCROBAS[45,]<-0.3}
  
  initPrebas <- InitMultiSite(
    nYearsMS = rep(nYears, nSites),
    siteInfo = siteInfo,
    pCROBAS = pCROBAS,
    defaultThin = defaultThin,
    multiInitVar = as.array(initVar),
    PAR = do.call(cbind, replicate(nRep, PARtran, simplify = FALSE)),
    TAir = do.call(cbind, replicate(nRep, TAirtran, simplify = FALSE)),
    VPD = do.call(cbind, replicate(nRep, VPDtran, simplify = FALSE)),
    Precip = do.call(cbind, replicate(nRep, Preciptran, simplify = FALSE)),
    CO2 = do.call(cbind, replicate(nRep, CO2tran, simplify = FALSE)),
    # litterSize = litterSize,
    # pAWEN = parsAWEN,
    # ClCut = 1.,
    yassoRun = 1
  )
  initPrebas$multiInitVar[,2,] <- round(6 + 2*SiteType - 0.005*rowMeans(initPrebas$ETSy) + 2.25) ##Initial age
  TransectOut <- multiPrebas(initPrebas)
  return(TransectOut)
}
