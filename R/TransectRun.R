TransectRun <- function(SiteType = NA,species = NA, nYears = 100, pCROBAS=pCROB,
                        defaultThin = 1) {
  
  nSites <- 7
  siteInfo <- matrix(c(NA,NA,SiteType,160,0,0,20,3,3,413,0.45,0.118),nSites,12,byrow = T)
  siteInfo[,2] <- siteInfo[,1] <- 1:nSites 
  
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
    siteInfo[,8:9] <- 1  # even-aged pure forests
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
    siteInfo[,8:9] <- 1  # even-aged pure forests
  }
  
  if (species == 'Mixed'){
    path <-
      system.file("transect", "initVarAll.csv", package = "Rprebasso")
    initVar <-
      read.csv(
        file = path, 
        header = T,
        row.names = 1
      )
    initVar <-
      aperm(array(as.matrix(initVar),dim=c(7,3,7)),c(3,1,2))
  
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
  
  PARtran <-as.matrix(    read.csv(
      file = path.par, 
      header = T,
      row.names = 1
    ))

  Preciptran <-as.matrix(
    read.csv(
      file = path.precip, 
      header = T,
      row.names = 1
    ))
  VPDtran <-as.matrix(
    read.csv(
      file = path.vpd, 
      header = T,
      row.names = 1
    ))
  TAirtran <-as.matrix(
    read.csv(
      file = path.tair, 
      header = T,
      row.names = 1
    ))
  CO2tran <-as.matrix(
    read.csv(
      file = path.co2, 
      header = T,
      row.names = 1
    ))
  
  nRep <- ceiling(nYears/(dim(PARtran)[2]/365))
  
  initPrebas <- InitMultiSite(  nYearsMS = rep(nYears,nSites),
                                siteInfo=siteInfo,                                     
                                pCROBAS = pCROBAS,
                                defaultThin=defaultThin,
                                multiInitVar = as.array(initVar),
                                PAR = do.call(cbind, replicate(nRep, PARtran, simplify=FALSE)),
                                TAir=do.call(cbind, replicate(nRep, TAirtran, simplify=FALSE)),
                                VPD=do.call(cbind, replicate(nRep, VPDtran, simplify=FALSE)),
                                Precip=do.call(cbind, replicate(nRep, Preciptran, simplify=FALSE)),
                                CO2=do.call(cbind, replicate(nRep, CO2tran, simplify=FALSE)),
                                # litterSize = litterSize,
                                # pAWEN = parsAWEN,
                                # ClCut = 1.,
                                yassoRun = 1)
  
  TransectOut<- multiPrebas(initPrebas)
  return(TransectOut)
  }
