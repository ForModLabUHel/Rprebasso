####Sweeden BAU populus###
sw_bau_pop <- function(initPrebas,pop_sites,
                       speciesID=8,
                       ftDens_before=1800, ####minimum density before precommercial thinning (first thinning)
                       ftH_before=5, ####minimum H before precommercial thinning (first thinning)
                       ftDens_target=1600, ####target density after precommercial thinning (first thinning)
                       tDens_before=1800, ####minimum density before commercial thinning
                       tH_before=14, ####minimum H before commercial thinning
                       tDens_target=700, ####target density after commercial thinning
                       age_Clcut = 70){
  
  initPrebas$ftTapioPar[1,speciesID,1,1:3] <- 99999.
  initPrebas$tTapioPar[,speciesID,1,1:3] <- 99999.
  
  initPrebas$ftTapioPar[1,speciesID,1,1] <- ftDens_before 
  initPrebas$ftTapioPar[1,speciesID,1,2] <- ftH_before
  initPrebas$ftTapioPar[1,speciesID,1,3] <- ftDens_target
  
  initPrebas$tTapioPar[1,speciesID,1,1] <- tDens_before
  initPrebas$tTapioPar[1,speciesID,1,2] <- tH_before
  initPrebas$tTapioPar[1,speciesID,1,3] <- tDens_target
  
  initPrebas$defaultThin[pop_sites] <- speciesID #note that this id will be used to choose the thinning parameters in B_prebas code (alternative_chooseThin)
  initPrebas$clct_pars[pop_sites,,] <- 9999
  initPrebas$clct_pars[pop_sites,,2] <- age_Clcut
  return(initPrebas)
}


####Sweeden BAU fague sylvestris###
thinning_def_fagus <- matrix(0,13,11)
# names(thinning_def_fagus) <- c("year","species","layer","H","D","B","Hc","frac_flag","N","Ac","pHarvTrees")

#commercial thinnings 1
thinning_def_fagus[1,] <- c(20,1,1,1.05,1.02,0.8,1,1,-999,-999,1)
thinning_def_fagus[2,] <- c(30,1,1,1.05,1.02,0.8,1,1,-999,-999,1)
thinning_def_fagus[3,] <- c(40,1,1,1.05,1.02,0.8,1,1,-999,-999,1)

#commercial thinnings 2
thinning_def_fagus[4,] <- c(50,1,1,0.98,0.98,0.9,1,1,-999,-999,1)
thinning_def_fagus[5,] <- c(60,1,1,0.98,0.98,0.9,1,1,-999,-999,1)
thinning_def_fagus[6,] <- c(70,1,1,0.98,0.98,0.9,1,1,-999,-999,1)
thinning_def_fagus[7,] <- c(80,1,1,0.98,0.98,0.9,1,1,-999,-999,1)
thinning_def_fagus[8,] <- c(90,1,1,0.98,0.98,0.9,1,1,-999,-999,1)
thinning_def_fagus[9,] <- c(100,1,1,0.98,0.98,0.9,1,1,-999,-999,1)

####clearcut
thinning_def_fagus[10,] <- c(110,1,1,1,1,0.75,1,1,-999,-999,1)
thinning_def_fagus[11,] <- c(120,1,1,1,1,0.7,1,1,-999,-999,1)
thinning_def_fagus[12,] <- c(130,1,1,1,1,0.,1,1,-999,-999,1)  #final cut

# thinning_def_fagus[13,] <- c(135,1,1,1.5,0.5,0.035,0.1,0,-999,-999,1) #replanting

sw_bau_fagsy <- function(initPrebas,siteXs,
                         ClCut_fag=130, #130 ####if ClCut_fag is NA the ingrowth layers are not removed at clearcut
                         thin_def=thinning_def_fagus ####default thinning matrix for fagus
){
  
  initPrebas$defaultThin[siteXs] <- 0
  if(is.na(ClCut_fag)){
    initPrebas$ClCut[siteXs] <- 0
    ClCut_fag <- thin_def[nrow(thin_def),1]
  }else{
    ####if ClCut_fag is active remove all trees also from ingrowth
    initPrebas$ClCut[siteXs] <- 1
    initPrebas$clct_pars[siteXs,,2] <- ClCut_fag
    initPrebas$clct_pars[siteXs,,c(1,3)] <- 9999
    thin_def <- thin_def[-nrow(thin_def),]
  }
  
  
  nlayers <- dim(initPrebas$multiInitVar)[3]
  maxYear <- initPrebas$maxYears
  
  layerX <- 1
  thinningX <- list()
  for(ij in siteXs){
    
    thin_def2 <- thin_def
    thin_def2[,1] <- thin_def[,1] + max(thin_def[,1])
    thinning_def_X <- rbind(thin_def,thin_def2)
    yearSims <- thinning_def_X[,1] - min(initPrebas$multiInitVar[ij,2,layerX],ClCut_fag)
    lines_to_move <- which(yearSims<1)
    yearSims_lasts <- max(yearSims) + thinning_def_X[lines_to_move,1]
    thinning_def_X[,1] <- yearSims
    sel <- which(thinning_def_X[,1]>0 & thinning_def_X[,1]<maxYear+1)
    thinning_def_X <- thinning_def_X[sel,]
    
    for(i in 1:initPrebas$nLayers[ij]){
      # if(initPrebas$multiInitVar[ij,5,i]>0){ #check if active site
      thinning_layer <- thinning_def_X
      thinning_layer[,2] <- initPrebas$multiInitVar[ij,1,i]
      thinning_layer[,3] <- i
      if(i==1){
        thinning_new=thinning_layer
      }else{
        thinning_new <- rbind(thinning_new,thinning_layer)
      }
      # }
    }
    thinningAll <- rbind(thinning_new,initPrebas$thinning[ij,1:initPrebas$nThinning[ij],])
    thinningAll <- thinningAll[order(thinningAll[,1], thinningAll[,3]), ]
    
    thinningX[[as.character(ij)]] <- thinningAll
  }
  
  maxNthin <- max(sapply(thinningX, nrow))
  maxNthin <- max(initPrebas$maxThin,maxNthin)
  
  initPrebas$maxThin <- maxNthin
  
  thinning <- array(0,dim=c(initPrebas$nSites,maxNthin,11))
  thinning[,,9:10] <- -999
  thinning[,,11] <- 1
  
  for(i in 1:initPrebas$nSites){
    thinning[i,1:initPrebas$nThinning[i],] <- initPrebas$thinning[i,1:initPrebas$nThinning[i],]
  }
  for(i in siteXs){
    initPrebas$nThinning[i] <- nrow(thinningX[[as.character(i)]])
    thinning[i,1:initPrebas$nThinning[i],] <- thinningX[[as.character(i)]]
  }
  initPrebas$thinning <- thinning
  
  return(initPrebas)
}


####Sweeden BAU piCo###
sw_bau_piCo <- function(initPrebas,siteXs,
                        speciesID=5,
                        # ftDens_before=1800, ####minimum density before precommercial thinning (first thinning)
                        # ftH_before=5, ####minimum H before precommercial thinning (first thinning)
                        # ftDens_target=1600, ####target density after precommercial thinning (first thinning)
                        tDens_before=900, ####minimum density before commercial thinning
                        tAge_before=30, ####minimum H before commercial thinning
                        tDens_target=800, ####target density after commercial thinning
                        age_Clcut = 55){
  
  initPrebas$ftTapioPar[1,speciesID,1,1:3] <- 99999.
  initPrebas$tTapioPar[,speciesID,1,1:3] <- 99999.
  
  initPrebas$tTapioPar[1,speciesID,1,1] <- tDens_before
  initPrebas$tTapioPar[1,speciesID,2,2] <- tAge_before
  initPrebas$tTapioPar[1,speciesID,1,3] <- tDens_target
  
  initPrebas$defaultThin[siteXs] <- speciesID #note that this id will be used to choose the thinning parameters in B_prebas code (alternative_chooseThin)
  initPrebas$clct_pars[siteXs,,] <- 9999
  initPrebas$clct_pars[siteXs,,2] <- age_Clcut
  return(initPrebas)
}



####Sweeden BAU alnus sp###
thinning_def_alnus <- matrix(0,3,11)
# names(thinning_def_fagus) <- c("year","species","layer","H","D","B","Hc","frac_flag","N","Ac","pHarvTrees")

#commercial thinnings 1
thinning_def_alnus[1,] <- c(25,1,1,1.0,1.0,0.75,1,1,-999,-999,1)
thinning_def_alnus[2,] <- c(35,1,1,1.0,1.0,0.75,1,1,-999,-999,1)

####clearcut
thinning_def_alnus[3,] <- c(60,1,1,1,1,0.,1,1,-999,-999,1)  #final cut

# thinning_def_alnus[4,] <- c(65,1,1,1.5,0.5,0.035,0.1,0,-999,-999,1) #replanting


sw_bau_AlnSp <- function(initPrebas,siteXs,
                         ClCut_age=60, #130 ####if ClCut_age is NA the ingrowth layers are not removed at clearcut
                         thin_def=thinning_def_alnus ####default thinning matrix for fagus
){
  
  initPrebas$defaultThin[siteXs] <- 0
  if(is.na(ClCut_age)){
    initPrebas$ClCut[siteXs] <- 0
    ClCut_age <- thin_def[3,1]
  }else{
    ####if ClCut_age is active remove all trees also from ingrowth
    initPrebas$ClCut[siteXs] <- 1
    initPrebas$clct_pars[siteXs,,2] <- ClCut_age
    initPrebas$clct_pars[siteXs,,c(1,3)] <- 9999
    thin_def <- thin_def[-3,]
  }
  
  
  nlayers <- dim(initPrebas$multiInitVar)[3]
  maxYear <- initPrebas$maxYears
  
  layerX <- 1
  thinningX <- list()
  for(ij in siteXs){
    
    thinning_def2 <- thin_def
    thinning_def2[,1] <- thin_def[,1] + ClCut_age
    thinning_defX <- rbind(thin_def,thinning_def2)
    yearSims <- thinning_defX[,1] - min(initPrebas$multiInitVar[ij,2,layerX],ClCut_age)
    lines_to_move <- which(yearSims<1)
    yearSims_lasts <- max(yearSims) + thinning_defX[lines_to_move,1]
    thinning_defX[,1] <- yearSims
    sel <- which(thinning_defX[,1]>0 & thinning_defX[,1]<maxYear+1)
    thinning_defX <- thinning_defX[sel,]
    
    for(i in 1:initPrebas$nLayers[ij]){
      # if(initPrebas$multiInitVar[ij,5,i]>0){ #check if active site
      thinning_layer <- thinning_defX
      thinning_layer[,2] <- initPrebas$multiInitVar[ij,1,i]
      thinning_layer[,3] <- i
      if(i==1){
        thinning_new=thinning_layer
      }else{
        thinning_new <- rbind(thinning_new,thinning_layer)
      }
      # }
    }
    thinningAll <- rbind(thinning_new,initPrebas$thinning[ij,1:initPrebas$nThinning[ij],])
    thinningAll <- thinningAll[order(thinningAll[,1], thinningAll[,3]), ]
    
    thinningX[[as.character(ij)]] <- thinningAll
  }
  
  maxNthin <- max(sapply(thinningX, nrow))
  maxNthin <- max(initPrebas$maxThin,maxNthin)
  
  initPrebas$maxThin <- maxNthin
  
  thinning <- array(0,dim=c(initPrebas$nSites,maxNthin,11))
  thinning[,,9:10] <- -999
  thinning[,,11] <- 1
  
  for(i in 1:initPrebas$nSites){
    thinning[i,1:initPrebas$nThinning[i],] <- initPrebas$thinning[i,1:initPrebas$nThinning[i],]
  }
  for(i in siteXs){
    initPrebas$nThinning[i] <- nrow(thinningX[[as.character(i)]])
    thinning[i,1:initPrebas$nThinning[i],] <- thinningX[[as.character(i)]]
  }
  initPrebas$thinning <- thinning
  
  return(initPrebas)
}




####Sweeden BAU quercus robur###
thinning_def_quercus <- matrix(0,5,11)
# names(thinning_def_quercus) <- c("year","species","layer","H","D","B","Hc","frac_flag","N","Ac","pHarvTrees")

#commercial thinnings
thinning_def_quercus[1,] <- c(30,1,1,1.0,1.0,1.,1,1,400,-999,1)
thinning_def_quercus[2,] <- c(60,1,1,1.0,1.0,1.,1,1,200,-999,1)

####clearcut
thinning_def_quercus[3,] <- c(110,1,1,1,1,0.2,1,1,-999,-999,1)
thinning_def_quercus[4,] <- c(115,1,1,1,1,0.5,1,1,-999,-999,1)
thinning_def_quercus[5,] <- c(120,1,1,1,1,0.,1,1,-999,-999,1)  #final cut

# thinning_def_fagus[13,] <- c(135,1,1,1.5,0.5,0.035,0.1,0,-999,-999,1) #replanting

sw_bau_QueRob <- function(initPrebas,siteXs,
                          ClCut_year=120, #130 ####if ClCut_year is NA the ingrowth layers are not removed at clearcut
                          thin_def=thinning_def_quercus ####default thinning matrix for fagus
){
  
  initPrebas$defaultThin[siteXs] <- 0
  if(is.na(ClCut_year)){
    initPrebas$ClCut[siteXs] <- 0
    ClCut_year <- thin_def[nrow(thin_def),1]
  }else{
    ####if ClCut_year is active remove all trees also from ingrowth
    initPrebas$ClCut[siteXs] <- 1
    initPrebas$clct_pars[siteXs,,2] <- ClCut_year
    initPrebas$clct_pars[siteXs,,c(1,3)] <- 9999
    thin_def <- thin_def[-nrow(thin_def),]
  }
  
  
  nlayers <- dim(initPrebas$multiInitVar)[3]
  maxYear <- initPrebas$maxYears
  
  layerX <- 1
  thinningX <- list()
  for(ij in siteXs){
    
    thin_def2 <- thin_def
    thin_def2[,1] <- thin_def[,1] + max(thin_def[,1])
    thinning_def_X <- rbind(thin_def,thin_def2)
    yearSims <- thinning_def_X[,1] - min(initPrebas$multiInitVar[ij,2,layerX],ClCut_year)
    lines_to_move <- which(yearSims<1)
    yearSims_lasts <- max(yearSims) + thinning_def_X[lines_to_move,1]
    thinning_def_X[,1] <- yearSims
    sel <- which(thinning_def_X[,1]>0 & thinning_def_X[,1]<maxYear+1)
    thinning_def_X <- thinning_def_X[sel,]
    
    for(i in 1:initPrebas$nLayers[ij]){
      # if(initPrebas$multiInitVar[ij,5,i]>0){ #check if active site
      thinning_layer <- thinning_def_X
      thinning_layer[,2] <- initPrebas$multiInitVar[ij,1,i]
      thinning_layer[,3] <- i
      if(i==1){
        thinning_new=thinning_layer
      }else{
        thinning_new <- rbind(thinning_new,thinning_layer)
      }
      # }
    }
    thinningAll <- rbind(thinning_new,initPrebas$thinning[ij,1:initPrebas$nThinning[ij],])
    thinningAll <- thinningAll[order(thinningAll[,1], thinningAll[,3]), ]
    
    thinningX[[as.character(ij)]] <- thinningAll
  }
  
  maxNthin <- max(sapply(thinningX, nrow))
  maxNthin <- max(initPrebas$maxThin,maxNthin)
  
  initPrebas$maxThin <- maxNthin
  
  thinning <- array(0,dim=c(initPrebas$nSites,maxNthin,11))
  thinning[,,9:10] <- -999
  thinning[,,11] <- 1
  
  for(i in 1:initPrebas$nSites){
    thinning[i,1:initPrebas$nThinning[i],] <- initPrebas$thinning[i,1:initPrebas$nThinning[i],]
  }
  for(i in siteXs){
    initPrebas$nThinning[i] <- nrow(thinningX[[as.character(i)]])
    thinning[i,1:initPrebas$nThinning[i],] <- thinningX[[as.character(i)]]
  }
  initPrebas$thinning <- thinning
  
  return(initPrebas)
}




####bau processing using a thinning matrix inputs
####this function was used for estonia BAU
bau_in_thinningMatrix <- function(initPrebas,siteXs,
                         ClCut_age,
                         ClCut_D = 9999,
                         ClCut_H = 9999,
                         nTree_seedlings,
                         year_seedling,
                         yearThin=NA,
                         baThin=NA,
                         hThin = 1,
                         dbhThin = 1, 
                         hcThin = 1,
                         fracThin = 1,
                         dens_after_Thin = -999,
                         acThin = -999,
                         pHarvTreeThin = 1,
                         sitesToThin = NA
){
  
  initPrebas$defaultThin[siteXs] <- 0
  initPrebas$initClearcut[siteXs,3] <- pi*(initPrebas$initClearcut[siteXs,2]/200)^2 * nTree_seedlings
  initPrebas$initClearcut[siteXs,5] <- year_seedling
  initPrebas$thinning[siteXs,,4:8] <- 0
  initPrebas$nThinning[siteXs] <- 0
  initPrebas$maxThin <- max(initPrebas$nThinning)
  initPrebas$thinning <- initPrebas$thinning[,1:initPrebas$maxThin,]
      
  if(is.na(ClCut_age)){
    initPrebas$ClCut[siteXs] <- 0
  }else{
    ####if ClCut_age is active remove all trees also from ingrowth
    initPrebas$ClCut[siteXs] <- 1
    initPrebas$clct_pars[siteXs,,1] <- ClCut_D
    initPrebas$clct_pars[siteXs,,2] <- ClCut_age
    initPrebas$clct_pars[siteXs,,3] <- ClCut_H
  }

  if(all(!is.na(yearThin))){
    nThinMat <- length(yearThin)
    
    # names(thinning) <- c("year","species","layer","H","D","B","Hc","frac_flag","N","Ac","pHarvTrees")
    thin_def <- matrix(0,(nThinMat+1),11) ###create thinning matrix
    thin_def[1:nThinMat,1] <- yearThin
    # thin_def[1:nThinMat,2] <- 1 #species
    # thin_def[1:nThinMat,3] <- 1 #layer
    thin_def[1:nThinMat,4] <- hThin
    thin_def[1:nThinMat,5] <- dbhThin
    thin_def[1:nThinMat,6] <- baThin
    thin_def[1:nThinMat,7] <- hcThin
    thin_def[1:nThinMat,8] <- fracThin
    thin_def[1:nThinMat,9] <- dens_after_Thin
    thin_def[1:nThinMat,10] <- acThin
    thin_def[1:nThinMat,11] <- pHarvTreeThin
    ####clearcut
    thin_def[(nThinMat+1),] <- c(ClCut_age,1,1,1,1,0.,1,1,-999,-999,1)  #final cut
    
    if(is.na(ClCut_age)){
      # initPrebas$ClCut[siteXs] <- 0
      ClCut_age <- thin_def[nrow(thin_def),1]
    }else{
      ####if ClCut_age is active remove all trees also from ingrowth
      # initPrebas$ClCut[siteXs] <- 1
      # initPrebas$clct_pars[siteXs,,1] <- ClCut_D
      # initPrebas$clct_pars[siteXs,,2] <- ClCut_age
      # initPrebas$clct_pars[siteXs,,3] <- ClCut_H
      thin_def <- thin_def[-nrow(thin_def),]
    }
    if(!is.matrix(thin_def)) thin_def <- matrix(thin_def,nrow=1)
    
    nlayers <- dim(initPrebas$multiInitVar)[3]
    maxYear <- initPrebas$maxYears
    
    layerX <- 1
    thinningX <- list()
    if(all(is.na(sitesToThin))){
      sitesThinned <- siteXs
    }else if(length(sitesToThin)>1){
      sitesThinned <- sitesToThin
    }else if(sitesToThin < 1){
      n <- round(length(siteXs) * sitesToThin) 
      sitesThinned <- sample(siteXs, n)
    }
    
    for(ij in sitesThinned){
      
      thinning_def2 <- thinning_def3 <- thin_def#matrix(thin_def,nrow = 1)
      thinning_def2[,1] <- thin_def[,1] + ClCut_age
      thinning_def3[,1] <- thin_def[,1] + 2*ClCut_age
      thinning_defX <- rbind(thin_def,thinning_def2,thinning_def3)
      yearSims <- thinning_defX[,1] - min(initPrebas$multiInitVar[ij,2,layerX],ClCut_age)
      lines_to_move <- which(yearSims<1)
      yearSims_lasts <- max(yearSims) + thinning_defX[lines_to_move,1]
      thinning_defX[,1] <- yearSims
      sel <- which(thinning_defX[,1]>0 & thinning_defX[,1]<maxYear+1)
      thinning_defX <- thinning_defX[sel,]
      
      for(i in 1:initPrebas$nLayers[ij]){
        # if(initPrebas$multiInitVar[ij,5,i]>0){ #check if active site
        thinning_layer <- matrix(thinning_defX,ncol=11)
        thinning_layer[,2] <- initPrebas$multiInitVar[ij,1,i]
        thinning_layer[,3] <- i
        if(i==1){
          thinning_new=thinning_layer
        }else{
          thinning_new <- rbind(thinning_new,thinning_layer)
        }
        # }
      }
      thinningAll <- thinning_new#rbind(thinning_new,initPrebas$thinning[ij,1:initPrebas$nThinning[ij],])
      thinningAll <- thinningAll[order(thinningAll[,1], thinningAll[,3]), ]
      
      thinningX[[as.character(ij)]] <- thinningAll
    }
    
    nthin <- sapply(thinningX, nrow)
    maxNthin <- max(sapply(thinningX, nrow))
    maxNthin <- max(initPrebas$maxThin,maxNthin)
    
    initPrebas$maxThin <- maxNthin
    
    thinning <- array(0,dim=c(initPrebas$nSites,maxNthin,11))
    thinning[,,9:10] <- -999
    thinning[,,11] <- 1
    
    for(i in 1:initPrebas$nSites){
      if(initPrebas$nThinning[i]>0) thinning[i,1:initPrebas$nThinning[i],] <- initPrebas$thinning[i,1:initPrebas$nThinning[i],]
    }
    initPrebas$nThinning[sitesThinned] <- nthin
    for(i in sitesThinned){
      initPrebas$nThinning[i] <- nrow(thinningX[[as.character(i)]])
      if(initPrebas$nThinning[i]>0) thinning[i,1:initPrebas$nThinning[i],] <- thinningX[[as.character(i)]]
    }
    initPrebas$thinning <- thinning
  }
  
  return(initPrebas)
}

#### Estonian BAU parameters
#### parameters for the bau_in_thinningMatrix function 
# picab_pinsy used for PicAb_PinSy_CC parameters
est_bau_pars_def <- list()
est_bau_pars_def$picab_pinsy$ClCut_age = 90
est_bau_pars_def$picab_pinsy$nTree_seedlings = 2200
est_bau_pars_def$picab_pinsy$year_seedling = 5
est_bau_pars_def$picab_pinsy$baThin = c(0.75,0.75)
est_bau_pars_def$picab_pinsy$yearThin = c(35,45)

# picab_decid used for: 
# PicAb_PopTr_CC, PicAb_BetSp_CC, PicAb_BetSp_CC (!!!aklnus removed at second thinning !!to be implemented)
# PicAb_hf_CC
est_bau_pars_def$picab_decid$ClCut_age = 80
est_bau_pars_def$picab_decid$nTree_seedlings = 2200
est_bau_pars_def$picab_decid$year_seedling = 4
est_bau_pars_def$picab_decid$yearThin = c(30,40)
est_bau_pars_def$picab_decid$baThin = c(0.75,0.75)

#PicAb_lf_CC, PinSy_hf_CC
est_bau_pars_def$conif$ClCut_age = 90
est_bau_pars_def$conif$nTree_seedlings = 2200
est_bau_pars_def$conif$year_seedling = 5
est_bau_pars_def$conif$yearThin = c(30,40)
est_bau_pars_def$conif$baThin = c(0.75,0.75)

#PinSy_mf_CC
est_bau_pars_def$PinSy_mf_CC$ClCut_age = 100
est_bau_pars_def$PinSy_mf_CC$nTree_seedlings = 2200
est_bau_pars_def$PinSy_mf_CC$year_seedling = 5
est_bau_pars_def$PinSy_mf_CC$yearThin = c(30,40)
est_bau_pars_def$PinSy_mf_CC$baThin = c(0.75,0.75)

#PinSy_lf_CC
est_bau_pars_def$PinSy_lf_CC$ClCut_age = 115
est_bau_pars_def$PinSy_lf_CC$nTree_seedlings = 2200
est_bau_pars_def$PinSy_lf_CC$year_seedling = 5
est_bau_pars_def$PinSy_lf_CC$yearThin = c(30,40)
est_bau_pars_def$PinSy_lf_CC$baThin = c(0.75,0.75)

#BetPe_PopTr_CC, BetPe_AlnSp_CC, BetSp_CC_High (to set min dbh)
#AlnSp_PopTr_CC, AlnSp_BetSp_CC, AlnSp_CC
est_bau_pars_def$BetPe_Al$ClCut_age = 60
est_bau_pars_def$BetPe_Al$nTree_seedlings = 2200
est_bau_pars_def$BetPe_Al$year_seedling = 5
est_bau_pars_def$BetPe_Al$yearThin = c(25,35)
est_bau_pars_def$BetPe_Al$baThin = c(0.75,0.75)

#BetSp_CC_Low (!!!to set min dbh)
est_bau_pars_def$BetSp_CC_Low$ClCut_age = 70
est_bau_pars_def$BetSp_CC_Low$nTree_seedlings = 2200
est_bau_pars_def$BetSp_CC_Low$year_seedling = 5
est_bau_pars_def$BetSp_CC_Low$yearThin = c(25,35)
est_bau_pars_def$BetSp_CC_Low$baThin = c(0.75,0.75)

#PopTr_CC, PopTr_BetPe_CC, PopTr_AlnGl_CC
est_bau_pars_def$PopTr_sp$ClCut_age = 43
est_bau_pars_def$PopTr_sp$nTree_seedlings = 2200
est_bau_pars_def$PopTr_sp$year_seedling = 3
est_bau_pars_def$PopTr_sp$yearThin = c(20)
est_bau_pars_def$PopTr_sp$baThin = c(0.75)

#PopTr_PicAb_CC
est_bau_pars_def$PopTr_PicAb_CC$ClCut_age = 75
est_bau_pars_def$PopTr_PicAb_CC$nTree_seedlings = 2200
est_bau_pars_def$PopTr_PicAb_CC$year_seedling = 4
est_bau_pars_def$PopTr_PicAb_CC$yearThin = c(30,40)
est_bau_pars_def$PopTr_PicAb_CC$baThin = c(0.75,0.75)

#PopTr_PicAb_CC
est_bau_pars_def$PopTr_PicAb_CC$ClCut_age = 75
est_bau_pars_def$PopTr_PicAb_CC$nTree_seedlings = 2200
est_bau_pars_def$PopTr_PicAb_CC$year_seedling = 4
est_bau_pars_def$PopTr_PicAb_CC$yearThin = c(30,40)
est_bau_pars_def$PopTr_PicAb_CC$baThin = c(0.75,0.75)

#AlnSp_PicAb_CC
est_bau_pars_def$AlnSp_PicAb_CC$ClCut_age = 70
est_bau_pars_def$AlnSp_PicAb_CC$nTree_seedlings = 2200
est_bau_pars_def$AlnSp_PicAb_CC$year_seedling = 5
est_bau_pars_def$AlnSp_PicAb_CC$yearThin = c(25,35)
est_bau_pars_def$AlnSp_PicAb_CC$baThin = c(0.75,0.75)

#QueSp_CC
est_bau_pars_def$QueSp_CC$ClCut_age = 105
est_bau_pars_def$QueSp_CC$nTree_seedlings = 2000
est_bau_pars_def$QueSp_CC$year_seedling = 5
est_bau_pars_def$QueSp_CC$yearThin = c(30,60)
est_bau_pars_def$QueSp_CC$baThin = c(1,1)
est_bau_pars_def$QueSp_CC$dens_after_Thin = c(400,200)

#QueSp_SW
est_bau_pars_def$QueSp_SW$ClCut_age = 110
est_bau_pars_def$QueSp_SW$nTree_seedlings = 1500
est_bau_pars_def$QueSp_SW$year_seedling = 5
est_bau_pars_def$QueSp_SW$yearThin = c(30,60,100,105)
est_bau_pars_def$QueSp_SW$baThin = c(1,1,0.2,0.5)
est_bau_pars_def$QueSp_SW$dens_after_Thin = c(400,200,-999,-999)


#denmark parameters
den_bau_pars_def <- list()

##PicAb_CC, PinSy_CC
den_bau_pars_def$conif$ClCut_age = 65
den_bau_pars_def$conif$ClCut_D = 45
den_bau_pars_def$conif$nTree_seedlings = 2700
den_bau_pars_def$conif$year_seedling = 3
den_bau_pars_def$conif$yearThin = c(35,45)
den_bau_pars_def$conif$baThin = c(0.75,0.75)
# den_bau_pars_def$conif$dens_after_Thin = c(-999,-999)

##FagSy_SW_High
den_bau_pars_def$FagSy_SW_High$ClCut_age = 115
# den_bau_pars_def$FagSy_SW_High$ClCut_D = 65
den_bau_pars_def$FagSy_SW_High$nTree_seedlings = 3500
den_bau_pars_def$FagSy_SW_High$year_seedling = 3
den_bau_pars_def$FagSy_SW_High$yearThin = c(15,25,35,45,55,65,75,85,95,105)
den_bau_pars_def$FagSy_SW_High$baThin = c(0.80,0.80,0.80,0.9,0.9,0.9,0.9,0.9,0.75,0.7)
den_bau_pars_def$FagSy_SW_High$hThin = c(1.02,1.02,1.02,0.98,0.98,0.98,0.98,0.98,0.75,0.7)
den_bau_pars_def$FagSy_SW_High$dbhThin = c(1.02,1.02,1.02,0.98,0.98,0.98,0.98,0.98,0.75,0.7)
# den_bau_pars_def$conif$dens_after_Thin = c(-999,-999)

##FagSy_SW_Low
den_bau_pars_def$FagSy_SW_Low$ClCut_age = 130
# den_bau_pars_def$FagSy_SW_Low$ClCut_D = 65
den_bau_pars_def$FagSy_SW_Low$nTree_seedlings = 3500
den_bau_pars_def$FagSy_SW_Low$year_seedling = 3
den_bau_pars_def$FagSy_SW_Low$yearThin = c(15,25,35,45,55,65,75,85,95,105,110,120)
den_bau_pars_def$FagSy_SW_Low$baThin = c(0.80,0.80,0.80,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.75,0.7)
den_bau_pars_def$FagSy_SW_Low$hThin = c(1.02,1.02,1.02,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.75,0.7)
den_bau_pars_def$FagSy_SW_Low$dbhThin = c(1.02,1.02,1.02,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.75,0.7)

##PopTr_CC
###same as sweden

##QueSp_SW
###same as sweden

# BetSp_CC
# set default thinning =1 and clearcut = 1 for the betula sites

# LarDe_CC
den_bau_pars_def$LarDe_CC$speciesID=3
den_bau_pars_def$LarDe_CC$ftDens_before=2000 ####minimum density before precommercial thinning (first thinning)
den_bau_pars_def$LarDe_CC$ftH_before=5 ####minimum H before precommercial thinning (first thinning)
den_bau_pars_def$LarDe_CC$ftDens_target=1300 ####target density after precommercial thinning (first thinning)
den_bau_pars_def$LarDe_CC$tDens_before=800 ####minimum density before commercial thinning
den_bau_pars_def$LarDe_CC$tH_before=13 ####minimum H before commercial thinning
den_bau_pars_def$LarDe_CC$tDens_target=700 ####target density after commercial thinning
den_bau_pars_def$LarDe_CC$age_Clcut = 70
  
##FraEx_CC
den_bau_pars_def$FraEx_CC$ClCut_age = 110
den_bau_pars_def$FraEx_CC$nTree_seedlings = 1500
den_bau_pars_def$FraEx_CC$year_seedling = 3
den_bau_pars_def$FraEx_CC$yearThin = c(31)
den_bau_pars_def$FraEx_CC$baThin = c(1)
den_bau_pars_def$FraEx_CC$dens_after_Thin = 700

## DouFi_SW 
den_bau_pars_def$DouFi_SW$ClCut_age = 120
den_bau_pars_def$DouFi_SW$nTree_seedlings = 2700
den_bau_pars_def$DouFi_SW$year_seedling = 3
den_bau_pars_def$DouFi_SW$yearThin = c(35,45,55,90,105)
den_bau_pars_def$DouFi_SW$baThin = c(0.70,0.70,0.70,0.5,0.5)
den_bau_pars_def$DouFi_SW$hThin = c(0.98,0.98,0.98,1,1)
den_bau_pars_def$DouFi_SW$dbhThin = c(0.98,0.98,0.98,1,1)

# AceSp_CC
den_bau_pars_def$AceSp_CC$ClCut_age = 85
den_bau_pars_def$AceSp_CC$nTree_seedlings = 2500
den_bau_pars_def$AceSp_CC$year_seedling = 3
den_bau_pars_def$AceSp_CC$yearThin = c(30,40)
den_bau_pars_def$AceSp_CC$baThin = c(0.8,0.8)
den_bau_pars_def$AceSp_CC$hThin = c(0.98,0.98)
den_bau_pars_def$AceSp_CC$dbhThin = c(0.98,0.98)

# AlnSp_CC
den_bau_pars_def$AlnSp_CC$ClCut_age = 60
den_bau_pars_def$AlnSp_CC$nTree_seedlings = 2500
den_bau_pars_def$AlnSp_CC$year_seedling = 3
den_bau_pars_def$AlnSp_CC$yearThin = c(25,35)
den_bau_pars_def$AlnSp_CC$baThin = c(0.75,0.75)


#### Irish BAU parameters
#### parameters for the bau_in_thinningMatrix function 
irl_bau_pars_def <- list()

# PicSi_good = PicSi_CC_good, PicSi_CC_CT_good
irl_bau_pars_def$PicSi_good$ClCut_age = 40
irl_bau_pars_def$PicSi_good$nTree_seedlings = 2500
irl_bau_pars_def$PicSi_good$year_seedling = 3
irl_bau_pars_def$PicSi_good$baThin = c(0.75)
irl_bau_pars_def$PicSi_good$yearThin = c(18)
irl_bau_pars_def$PicSi_good$sitesToThin =0.3

# PicSi_low = PicSi_CC_low, PicSi_CC_CT_low
irl_bau_pars_def$PicSi_low$ClCut_age = 50
irl_bau_pars_def$PicSi_low$nTree_seedlings = 2500
irl_bau_pars_def$PicSi_low$year_seedling = 4
irl_bau_pars_def$PicSi_low$baThin = c(0.75)
irl_bau_pars_def$PicSi_low$yearThin = c(25)
irl_bau_pars_def$PicSi_low$sitesToThin =0.3
irl_bau_pars_def$PicSi_low$hThin = c(1.02)
irl_bau_pars_def$PicSi_low$dbhThin = c(1.02)

# PicSidom_CC
irl_bau_pars_def$PicSidom_CC$ClCut_age = 50
irl_bau_pars_def$PicSidom_CC$nTree_seedlings = 2500
irl_bau_pars_def$PicSidom_CC$year_seedling = 3
irl_bau_pars_def$PicSidom_CC$baThin = c(0.75)
irl_bau_pars_def$PicSidom_CC$yearThin = c(20)
irl_bau_pars_def$PicSidom_CC$sitesToThin =0.65

# PicSidom_CT_CC
irl_bau_pars_def$PicSidom_CT_CC$ClCut_age = 50
irl_bau_pars_def$PicSidom_CT_CC$nTree_seedlings = 2500
irl_bau_pars_def$PicSidom_CT_CC$year_seedling = 3
irl_bau_pars_def$PicSidom_CT_CC$baThin = c(0.75)
irl_bau_pars_def$PicSidom_CT_CC$yearThin = c(20)
irl_bau_pars_def$PicSidom_CT_CC$sitesToThin =0.55
irl_bau_pars_def$PicSidom_CT_CC$hThin = c(1.02)
irl_bau_pars_def$PicSidom_CT_CC$dbhThin = c(1.02)

# PinCo_CC_good
irl_bau_pars_def$PinCo_CC_good$ClCut_age = 40
irl_bau_pars_def$PinCo_CC_good$nTree_seedlings = 2500
irl_bau_pars_def$PinCo_CC_good$year_seedling=3

# PinCo_CC_low
irl_bau_pars_def$PinCo_CC_low$ClCut_age = 50
irl_bau_pars_def$PinCo_CC_low$nTree_seedlings = 1800
irl_bau_pars_def$PinCo_CC_low$year_seedling=4

# BL_softwood_limited
irl_bau_pars_def$BL_softwood_limited$ClCut_age = 60
irl_bau_pars_def$BL_softwood_limited$nTree_seedlings = 2500
irl_bau_pars_def$BL_softwood_limited$year_seedling=3

# BL_hardwood_limited
irl_bau_pars_def$BL_hardwood_limited$ClCut_age = 90
irl_bau_pars_def$BL_hardwood_limited$nTree_seedlings = 2500
irl_bau_pars_def$BL_hardwood_limited$year_seedling = 4
irl_bau_pars_def$BL_hardwood_limited$baThin = c(0.7)
irl_bau_pars_def$BL_hardwood_limited$yearThin = c(25)
irl_bau_pars_def$BL_hardwood_limited$hThin = c(1.02)
irl_bau_pars_def$BL_hardwood_limited$dbhThin = c(1.02)

#### Latvia BAU parameters
#### parameters for the bau_in_thinningMatrix function 
lat_bau_pars_def <- list()

forType <- "PicAb_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 85
lat_bau_pars_def[[forType]]$nTree_seedlings = 2000
lat_bau_pars_def[[forType]]$year_seedling = 3
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "PinSy_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 111
lat_bau_pars_def[[forType]]$nTree_seedlings = 3000
lat_bau_pars_def[[forType]]$year_seedling = 4
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "PinSy_BetPe_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 100
lat_bau_pars_def[[forType]]$nTree_seedlings = 3000
lat_bau_pars_def[[forType]]$year_seedling = 4
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "PicAb_dec" #"PicAb_PopTr_CC","PicAb_BetSp_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 80
lat_bau_pars_def[[forType]]$nTree_seedlings = 2500
lat_bau_pars_def[[forType]]$year_seedling = 4
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "PicAb_PinSy_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 90
lat_bau_pars_def[[forType]]$nTree_seedlings = 2500
lat_bau_pars_def[[forType]]$year_seedling = 4
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(35,45)

forType <- "decid" #"BetSp_CC", "AlnGl_CC", "AlnIn_CC","AlnSp_BetSp_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 61
lat_bau_pars_def[[forType]]$nTree_seedlings = 2000
lat_bau_pars_def[[forType]]$year_seedling = 3
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(25,35)

forType <- "AlnSp_PicAb_CC" 
lat_bau_pars_def[[forType]]$ClCut_age = 70
lat_bau_pars_def[[forType]]$nTree_seedlings = 2000
lat_bau_pars_def[[forType]]$year_seedling = 2
lat_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lat_bau_pars_def[[forType]]$yearThin = c(25,35)

forType <- "PopTr"# "PopTr_CC", "PopTr_BetPe_CC"
lat_bau_pars_def[[forType]]$ClCut_age = 45
lat_bau_pars_def[[forType]]$nTree_seedlings = 2000
lat_bau_pars_def[[forType]]$year_seedling = 1
lat_bau_pars_def[[forType]]$baThin = c(0.75)
lat_bau_pars_def[[forType]]$yearThin = c(20)

forType <- "FraEx_CC" 
lat_bau_pars_def[[forType]]$ClCut_age = 85
lat_bau_pars_def[[forType]]$nTree_seedlings = 1500
lat_bau_pars_def[[forType]]$year_seedling = 3
lat_bau_pars_def[[forType]]$baThin = c(1)
lat_bau_pars_def[[forType]]$yearThin = c(31)
lat_bau_pars_def[[forType]]$dens_after_Thin <- 700

forType <- "QueSp_SW"
lat_bau_pars_def[[forType]]$ClCut_age = 121
lat_bau_pars_def[[forType]]$nTree_seedlings = 1500
lat_bau_pars_def[[forType]]$year_seedling = 5
lat_bau_pars_def[[forType]]$yearThin = c(30,60,101,106)
lat_bau_pars_def[[forType]]$baThin = c(1,1,0.2,0.5)
lat_bau_pars_def[[forType]]$dens_after_Thin = c(400,200,-999,-999)

forType <- "TilSp_SW" 
lat_bau_pars_def[[forType]]$ClCut_age = 95
lat_bau_pars_def[[forType]]$nTree_seedlings = 2000
lat_bau_pars_def[[forType]]$year_seedling = 3
lat_bau_pars_def[[forType]]$baThin = c(1,0.3,0.5)
lat_bau_pars_def[[forType]]$yearThin = c(31,81,86)
lat_bau_pars_def[[forType]]$dens_after_Thin <- c(700,-999,-999)


#### Lithuania BAU parameters
#### parameters for the bau_in_thinningMatrix function 
lit_bau_pars_def <- list()

forType <- "PicAb_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 96
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75)
lit_bau_pars_def[[forType]]$yearThin = c(16,30,40)

forType <- "PicAb_GSC"
lit_bau_pars_def[[forType]]$ClCut_age = 125
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75,0.50,0.5)
lit_bau_pars_def[[forType]]$yearThin = c(16,30,40,95,110)

forType <- "PicAb_BeSp_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 80
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lit_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "PinSy_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 111
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75)
lit_bau_pars_def[[forType]]$yearThin = c(16,30,40)

forType <- "PinSy_SW"
lit_bau_pars_def[[forType]]$ClCut_age = 125
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75,0.3)
lit_bau_pars_def[[forType]]$yearThin = c(16,30,40,110)

forType <- "PinSy_BetPe_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 100
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lit_bau_pars_def[[forType]]$yearThin = c(30,40)

forType <- "BetSp_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 76
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75)
lit_bau_pars_def[[forType]]$yearThin = c(31)

forType <- "BetSp_AlnGl_SW"
lit_bau_pars_def[[forType]]$ClCut_age = 75
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.2)
lit_bau_pars_def[[forType]]$yearThin = c(31,65)

forType <- "AlnGl_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 75
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75)
lit_bau_pars_def[[forType]]$yearThin = c(26)

forType <- "AlnIn_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 41
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75)
lit_bau_pars_def[[forType]]$yearThin = c(26)

forType <- "AlnSp" #"AlnSp_PinSy_CC", "AlnSp_BetSp_CC", "BetPe_AlnSp_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 60
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.75)
lit_bau_pars_def[[forType]]$yearThin = c(25,35)

forType <- "PopTr" #"PopTr_CC","PopTr_BetPe_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 51
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75)
lit_bau_pars_def[[forType]]$yearThin = c(26)

forType <- "PopTr_AlnIn_GSC"
lit_bau_pars_def[[forType]]$ClCut_age = 65
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(0.75,0.3)
lit_bau_pars_def[[forType]]$yearThin = c(26,50)

forType <- "FraEx_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 110
lit_bau_pars_def[[forType]]$nTree_seedlings = 2500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(1)
lit_bau_pars_def[[forType]]$yearThin = c(31)
lit_bau_pars_def[[forType]]$dens_after_Thin = 700

forType <- "QueSp_CC"
lit_bau_pars_def[[forType]]$ClCut_age = 131
lit_bau_pars_def[[forType]]$nTree_seedlings = 1500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(1,1)
lit_bau_pars_def[[forType]]$yearThin = c(31,60)
lit_bau_pars_def[[forType]]$dens_after_Thin = c(400,200)

forType <- "QueSp_SW"
lit_bau_pars_def[[forType]]$ClCut_age = 131
lit_bau_pars_def[[forType]]$nTree_seedlings = 1500
lit_bau_pars_def[[forType]]$year_seedling = 3
lit_bau_pars_def[[forType]]$baThin = c(1,1,0.2,0.5)
lit_bau_pars_def[[forType]]$yearThin = c(31,60,121,126)
lit_bau_pars_def[[forType]]$dens_after_Thin = c(400,200,-999,-999)

#### Germany BAU parameters
#### parameters for the bau_in_thinningMatrix function 
de_bau_pars_def <- list()

forType <- "General_CC"
de_bau_pars_def[[forType]]$ClCut_age = 95
de_bau_pars_def[[forType]]$nTree_seedlings = 2000
de_bau_pars_def[[forType]]$year_seedling = 3
de_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75)
de_bau_pars_def[[forType]]$yearThin = c(35,45,55)
de_bau_pars_def[[forType]]$hThin = c(0.98,0.98,0.98)
de_bau_pars_def[[forType]]$dbhThin = c(0.98,0.98,0.98)

forType <- "BLSCON_GSC"
de_bau_pars_def[[forType]]$ClCut_age = 125
de_bau_pars_def[[forType]]$nTree_seedlings = 2500
de_bau_pars_def[[forType]]$year_seedling = 3
de_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75,0.5,0.5)
de_bau_pars_def[[forType]]$yearThin = c(35,55,75,95,110)
de_bau_pars_def[[forType]]$hThin = c(0.98,0.98,0.98,1,1)
de_bau_pars_def[[forType]]$dbhThin = c(0.98,0.98,0.98,1,1)

forType <- "PicAb_CC"
de_bau_pars_def[[forType]]$ClCut_age = 90
de_bau_pars_def[[forType]]$nTree_seedlings = 2500
de_bau_pars_def[[forType]]$year_seedling = 3
de_bau_pars_def[[forType]]$baThin = c(0.8,0.8,0.8,0.8)
de_bau_pars_def[[forType]]$yearThin = c(20,35,50,65)
de_bau_pars_def[[forType]]$hThin = c(1.02,1.02,0.98,0.98)
de_bau_pars_def[[forType]]$dbhThin = c(1.02,1.02,0.98,0.98)

forType <- "QueSp_FagSy_SS"
de_bau_pars_def[[forType]]$ClCut_age = 125
de_bau_pars_def[[forType]]$nTree_seedlings = 2500
de_bau_pars_def[[forType]]$year_seedling = 3
de_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75,0.2,0.5)
de_bau_pars_def[[forType]]$yearThin = c(35,55,75,105,115)
de_bau_pars_def[[forType]]$hThin = c(0.98,0.98,0.98,1,1)
de_bau_pars_def[[forType]]$dbhThin = c(0.98,0.98,0.98,1,1)

forType <- "QueSp_SS"
de_bau_pars_def[[forType]]$ClCut_age = 130
de_bau_pars_def[[forType]]$nTree_seedlings = 2500
de_bau_pars_def[[forType]]$year_seedling = 3
de_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.75,0.2,0.5)
de_bau_pars_def[[forType]]$yearThin = c(40,80,110,120,130)
de_bau_pars_def[[forType]]$hThin = c(0.98,0.98,0.98,1,1)
de_bau_pars_def[[forType]]$dbhThin = c(0.98,0.98,0.98,1,1)

#### Poland BAU parameters
#### parameters for the bau_in_thinningMatrix function 
pl_bau_pars_def <- list()

forType <- "FagSy_SW"
pl_bau_pars_def[[forType]]$ClCut_age = 130
pl_bau_pars_def[[forType]]$nTree_seedlings = 7000
pl_bau_pars_def[[forType]]$year_seedling = 3
pl_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.7,0.5)
pl_bau_pars_def[[forType]]$yearThin = c(30,45,100,110)
pl_bau_pars_def[[forType]]$hThin = c(1,1,1,1)
pl_bau_pars_def[[forType]]$dbhThin = c(1,1,1,1)

forType <- "QueSp_SW"
pl_bau_pars_def[[forType]]$ClCut_age = 150
pl_bau_pars_def[[forType]]$nTree_seedlings = 7000
pl_bau_pars_def[[forType]]$year_seedling = 3
pl_bau_pars_def[[forType]]$baThin = c(0.75,0.75,0.6,0.6)
pl_bau_pars_def[[forType]]$yearThin = c(30,45,120,130)
pl_bau_pars_def[[forType]]$hThin = c(0.98,0.98,1,1)
pl_bau_pars_def[[forType]]$dbhThin = c(0.98,0.98,1,1)

forType <- "CON_CC"
pl_bau_pars_def[[forType]]$ClCut_age = 80
pl_bau_pars_def[[forType]]$nTree_seedlings = 4000
pl_bau_pars_def[[forType]]$year_seedling = 3
pl_bau_pars_def[[forType]]$baThin = c(0.7,0.80,0.85,0.9,0.95)
pl_bau_pars_def[[forType]]$yearThin = c(25,35,45,55,65)
pl_bau_pars_def[[forType]]$hThin = c(1,1,1,1,1)
pl_bau_pars_def[[forType]]$dbhThin = c(1,1,1,1,1)


#' management function updater (ForestNavigator)
#'
#' @param initPrebas Rprebasso initialization object for multisite created by the InitMultiSite function
#' @param forest_type_management_tab Table with forest management (for_man) by forest type (forest_type) and siteID (site)
#' @param country Country of simulations. Choose between: Sweden
#' @param management type of management. choose between: BAU
#'
#' @returns
#' @export
#'
#' @examples
forest_management_update <- function(initPrebas, 
                                     forest_type_management_tab, 
                                     country, 
                                     management,
                                     est_bau_pars=est_bau_pars_def,
                                     den_bau_pars=den_bau_pars_def,
                                     irl_bau_pars=irl_bau_pars_def,
                                     lat_bau_pars=lat_bau_pars_def,
                                     lit_bau_pars=lit_bau_pars_def,
                                     de_bau_pars=de_bau_pars_def,
                                     pl_bau_pars=pl_bau_pars_def){
  available_countries <- c("Sweden","Finland","Estonia","Denmark","Ireland","Latvia","Lithuania","Germany","Poland","United Kingdom")
  available_managements <- c("bau", "noman")
  if(!country %in% available_countries) stop(cat("This country: ", country,
                                                 " is not between the available countries: ", available_countries,fill = TRUE))
  if(!management %in% available_managements) stop(cat("This management: ", management, 
                                                      " is not between the available managements: ", available_managements,fill = TRUE))
  if(country == "Sweden" & management=="bau"){
    if(dim(initPrebas$ftTapioPar)[2] < 12){
      dims <- dim(initPrebas$ftTapioPar)
      ftTapioPar <- array(999,dim = c(5,12,3,7))
      ftTapioPar[1:dims[1],1:dims[2],1:dims[3],1:dims[4]] <- initPrebas$ftTapioPar
      initPrebas$ftTapioPar <- ftTapioPar
    }
    if(dim(initPrebas$tTapioPar)[2] < 12){
      dims <- dim(initPrebas$tTapioPar)
      tTapioPar <- array(999,dim = c(5,12,3,7))
      tTapioPar[1:dims[1],1:dims[2],1:dims[3],1:dims[4]] <- initPrebas$tTapioPar
      initPrebas$tTapioPar <- tTapioPar
    }
    ##find the sites with alternative management##
    pop_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PopTr_CC")])
    alnus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnSp_CC")])
    quercus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    pinco_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinCo_CC")])
    fagus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FagSy_SW")])
    ##----##
    
    ## update the initialization##
    if(length(pop_sites)>0) initPrebas <- sw_bau_pop(initPrebas,pop_sites)
    if(length(alnus_sites)>0) initPrebas <- sw_bau_AlnSp(initPrebas,alnus_sites)
    if(length(quercus_sites)>0) initPrebas <- sw_bau_QueRob(initPrebas,quercus_sites)
    if(length(pinco_sites)>0) initPrebas <- sw_bau_piCo(initPrebas,pinco_sites)
    if(length(fagus_sites)>0) initPrebas <- sw_bau_fagsy(initPrebas,fagus_sites)
    ##----##
  }
  
  if(country == "Estonia" & management=="bau"){
    
    ##find the sites with alternative management##
    picab_pinsy_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicAb_PinSy_CC")])
    picab_decid_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                            c("PicAb_PopTr_CC", "PicAb_BetSp_CC", "PicAb_BetSp_CC", "PicAb_hf_CC"))]) 
    conif_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                            c("PicAb_lf_CC", "PinSy_hf_CC"))])
    PinSy_mf_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_mf_CC")])
    PinSy_lf_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_lf_CC")])
    BetPe_Al_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                  c("BetPe_PopTr_CC", "BetPe_AlnSp_CC", "BetSp_CC_High", "AlnSp_PopTr_CC", "AlnSp_BetSp_CC", "AlnSp_CC"))]) 
    BetSp_CC_Low_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BetSp_CC_Low")])
    PopTr_sp_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                  c("PopTr_CC", "PopTr_BetPe_CC", "PopTr_AlnGl_CC"))]) 
    PopTr_PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PopTr_PicAb_CC")])
    AlnSp_PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnSp_PicAb_CC")])
    QueSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_CC")])
    QueSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    ##----##
    
    ## update the initialization##
    if(length(picab_pinsy_sites)>0) initPrebas <- bau_in_thinningMatrix(
                    initPrebas=initPrebas,
                    siteXs=picab_pinsy_sites,
                    ClCut_age= est_bau_pars$picab_pinsy$ClCut_age,
                    nTree_seedlings=est_bau_pars$picab_pinsy$nTree_seedlings,
                    year_seedling=est_bau_pars$picab_pinsy$year_seedling,
                    yearThin=est_bau_pars$picab_pinsy$yearThin,
                    baThin=est_bau_pars$picab_pinsy$baThin,
                    hThin = 1,
                    dbhThin = 1,
                    hcThin = 1,
                    fracThin = 1,
                    dens_after_Thin = -999,
                    acThin = -999,
                    pHarvTreeThin = 1)
    if(length(picab_decid_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=picab_decid_sites,
      ClCut_age= est_bau_pars$picab_decid$ClCut_age,
      nTree_seedlings=est_bau_pars$picab_decid$nTree_seedlings,
      year_seedling=est_bau_pars$picab_decid$year_seedling,
      yearThin=est_bau_pars$picab_decid$yearThin,
      baThin=est_bau_pars$picab_decid$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    if(length(conif_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=conif_sites,
      ClCut_age= est_bau_pars$conif$ClCut_age,
      nTree_seedlings=est_bau_pars$conif$nTree_seedlings,
      year_seedling=est_bau_pars$conif$year_seedling,
      yearThin=est_bau_pars$conif$yearThin,
      baThin=est_bau_pars$conif$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    if(length(PinSy_mf_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_mf_CC_sites,
      ClCut_age= est_bau_pars$PinSy_mf_CC$ClCut_age,
      nTree_seedlings=est_bau_pars$PinSy_mf_CC$nTree_seedlings,
      year_seedling=est_bau_pars$PinSy_mf_CC$year_seedling,
      yearThin=est_bau_pars$PinSy_mf_CC$yearThin,
      baThin=est_bau_pars$PinSy_mf_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(PinSy_lf_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_lf_CC_sites,
      ClCut_age= est_bau_pars$PinSy_lf_CC$ClCut_age,
      nTree_seedlings=est_bau_pars$PinSy_lf_CC$nTree_seedlings,
      year_seedling=est_bau_pars$PinSy_lf_CC$year_seedling,
      yearThin=est_bau_pars$PinSy_lf_CC$yearThin,
      baThin=est_bau_pars$PinSy_lf_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(BetPe_Al_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BetPe_Al_sites,
      ClCut_age= est_bau_pars$BetPe_Al$ClCut_age,
      nTree_seedlings=est_bau_pars$BetPe_Al$nTree_seedlings,
      year_seedling=est_bau_pars$BetPe_Al$year_seedling,
      yearThin=est_bau_pars$BetPe_Al$yearThin,
      baThin=est_bau_pars$BetPe_Al$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(BetSp_CC_Low_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BetSp_CC_Low_sites,
      ClCut_age= est_bau_pars$BetSp_CC_Low$ClCut_age,
      nTree_seedlings=est_bau_pars$BetSp_CC_Low$nTree_seedlings,
      year_seedling=est_bau_pars$BetSp_CC_Low$year_seedling,
      yearThin=est_bau_pars$BetSp_CC_Low$yearThin,
      baThin=est_bau_pars$BetSp_CC_Low$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(PopTr_sp_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_sp_sites,
      ClCut_age= est_bau_pars$PopTr_sp$ClCut_age,
      nTree_seedlings=est_bau_pars$PopTr_sp$nTree_seedlings,
      year_seedling=est_bau_pars$PopTr_sp$year_seedling,
      yearThin=est_bau_pars$PopTr_sp$yearThin,
      baThin=est_bau_pars$PopTr_sp$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(PopTr_PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_PicAb_CC_sites,
      ClCut_age= est_bau_pars$PopTr_PicAb_CC$ClCut_age,
      nTree_seedlings=est_bau_pars$PopTr_PicAb_CC$nTree_seedlings,
      year_seedling=est_bau_pars$PopTr_PicAb_CC$year_seedling,
      yearThin=est_bau_pars$PopTr_PicAb_CC$yearThin,
      baThin=est_bau_pars$PopTr_PicAb_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(AlnSp_PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnSp_PicAb_CC_sites,
      ClCut_age= est_bau_pars$AlnSp_PicAb_CC$ClCut_age,
      nTree_seedlings=est_bau_pars$AlnSp_PicAb_CC$nTree_seedlings,
      year_seedling=est_bau_pars$AlnSp_PicAb_CC$year_seedling,
      yearThin=est_bau_pars$AlnSp_PicAb_CC$yearThin,
      baThin=est_bau_pars$AlnSp_PicAb_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(QueSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_CC_sites,
      ClCut_age= est_bau_pars$QueSp_CC$ClCut_age,
      nTree_seedlings=est_bau_pars$QueSp_CC$nTree_seedlings,
      year_seedling=est_bau_pars$QueSp_CC$year_seedling,
      yearThin=est_bau_pars$QueSp_CC$yearThin,
      baThin=est_bau_pars$QueSp_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = est_bau_pars$QueSp_CC$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(QueSp_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_SW_sites,
      ClCut_age= est_bau_pars$QueSp_SW$ClCut_age,
      nTree_seedlings=est_bau_pars$QueSp_SW$nTree_seedlings,
      year_seedling=est_bau_pars$QueSp_SW$year_seedling,
      yearThin=est_bau_pars$QueSp_SW$yearThin,
      baThin=est_bau_pars$QueSp_SW$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = est_bau_pars$QueSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)

    ##----##
  }

  if(country == "Denmark" & management=="bau"){
    
    ##find the sites with alternative management##
    ##PicAb_CC, PinSy_CC
    conif_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                                                c("PicAb_CC", "PinSy_CC"))])
    FagSy_SW_High_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FagSy_SW_High")])
    FagSy_SW_Low_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FagSy_SW_Low")])
    PopTr_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PopTr_CC")])
    QueSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    BetSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BetSp_CC")])
    LarDe_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "LarDe_CC")])
    FraEx_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FraEx_CC")])
    AceSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AceSp_CC")])
    AlnSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnSp_CC")])
    DouFi_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "DouFi_SW")])
    ##----##
    
    ## update the initialization##
    if(length(conif_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=conif_sites,
      ClCut_age= den_bau_pars$conif$ClCut_age,
      nTree_seedlings=den_bau_pars$conif$nTree_seedlings,
      year_seedling=den_bau_pars$conif$year_seedling,
      yearThin=den_bau_pars$conif$yearThin,
      baThin=den_bau_pars$conif$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(FagSy_SW_High_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas = initPrebas,
      siteXs = FagSy_SW_High_sites,
      ClCut_age = den_bau_pars$FagSy_SW_High$ClCut_age,
      nTree_seedlings = den_bau_pars$FagSy_SW_High$nTree_seedlings,
      year_seedling=den_bau_pars$FagSy_SW_High$year_seedling,
      yearThin=den_bau_pars$FagSy_SW_High$yearThin,
      baThin=den_bau_pars$FagSy_SW_High$baThin,
      hThin = den_bau_pars$FagSy_SW_High$hThin,
      dbhThin = den_bau_pars$FagSy_SW_High$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(FagSy_SW_Low_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas = initPrebas,
      siteXs = FagSy_SW_Low_sites,
      ClCut_age = den_bau_pars$FagSy_SW_Low$ClCut_age,
      nTree_seedlings = den_bau_pars$FagSy_SW_Low$nTree_seedlings,
      year_seedling=den_bau_pars$FagSy_SW_Low$year_seedling,
      yearThin=den_bau_pars$FagSy_SW_Low$yearThin,
      baThin=den_bau_pars$FagSy_SW_Low$baThin,
      hThin = den_bau_pars$FagSy_SW_Low$hThin,
      dbhThin = den_bau_pars$FagSy_SW_Low$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    ##PopTr_CC
    if(length(PopTr_CC_sites)>0) initPrebas <- sw_bau_pop(initPrebas,PopTr_CC_sites)
    
    ##QueSp_SW
    if(length(QueSp_SW_sites)>0) initPrebas <- sw_bau_QueRob(initPrebas,QueSp_SW_sites)
    
    # BetSp_CC
    if(length(BetSp_CC_sites)>0){
      initPrebas$ClCut[BetSp_CC_sites] <- 1
      initPrebas$defaultThin[BetSp_CC_sites] <- 1
    }
    
    if(length(LarDe_CC_sites)>0) initPrebas <- sw_bau_pop(initPrebas = initPrebas,
            pop_sites = LarDe_CC_sites,
            speciesID=den_bau_pars$LarDe_CC$speciesID,
            ftDens_before=den_bau_pars$LarDe_CC$ftDens_before, ####minimum density before precommercial thinning (first thinning)
            ftH_before=den_bau_pars$LarDe_CC$ftH_before, ####minimum H before precommercial thinning (first thinning)
            ftDens_target=den_bau_pars$LarDe_CC$ftDens_target, ####target density after precommercial thinning (first thinning)
            tDens_before=den_bau_pars$LarDe_CC$tDens_before, ####minimum density before commercial thinning
            tH_before=den_bau_pars$LarDe_CC$tH_before, ####minimum H before commercial thinning
            tDens_target=den_bau_pars$LarDe_CC$tDens_target, ####target density after commercial thinning
            age_Clcut = den_bau_pars$LarDe_CC$age_Clcut)
    
    if(length(FraEx_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=FraEx_CC_sites,
      ClCut_age= den_bau_pars$FraEx_CC$ClCut_age,
      nTree_seedlings=den_bau_pars$FraEx_CC$nTree_seedlings,
      year_seedling=den_bau_pars$FraEx_CC$year_seedling,
      yearThin=den_bau_pars$FraEx_CC$yearThin,
      baThin=den_bau_pars$FraEx_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = den_bau_pars$FraEx_CC$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(AceSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AceSp_CC_sites,
      ClCut_age= den_bau_pars$AceSp_CC$ClCut_age,
      nTree_seedlings=den_bau_pars$AceSp_CC$nTree_seedlings,
      year_seedling=den_bau_pars$AceSp_CC$year_seedling,
      yearThin=den_bau_pars$AceSp_CC$yearThin,
      baThin=den_bau_pars$AceSp_CC$baThin,
      hThin = den_bau_pars$AceSp_CC$hThin,
      dbhThin = den_bau_pars$AceSp_CC$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    if(length(AlnSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnSp_CC_sites,
      ClCut_age= den_bau_pars$AlnSp_CC$ClCut_age,
      nTree_seedlings=den_bau_pars$AlnSp_CC$nTree_seedlings,
      year_seedling=den_bau_pars$AlnSp_CC$year_seedling,
      yearThin=den_bau_pars$AlnSp_CC$yearThin,
      baThin=den_bau_pars$AlnSp_CC$baThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    
    if(length(DouFi_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=DouFi_SW_sites,
      ClCut_age= den_bau_pars$DouFi_SW$ClCut_age,
      nTree_seedlings=den_bau_pars$DouFi_SW$nTree_seedlings,
      year_seedling=den_bau_pars$DouFi_SW$year_seedling,
      yearThin=den_bau_pars$DouFi_SW$yearThin,
      baThin=den_bau_pars$DouFi_SW$baThin,
      hThin = den_bau_pars$DouFi_SW$hThin,
      dbhThin = den_bau_pars$DouFi_SW$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    ##----##
  }

  if(country == "Ireland" & management=="bau"){
    
    ##find the sites with alternative management##
    PicSi_good_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                           c("PicSi_CC_good", "PicSi_CC_CT_good","PicSi_CT_CC_good"))])
    PicSi_low_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                           c("PicSi_CC_low", "PicSi_CC_CT_low"))])
    PicSidom_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicSidom_CC")])
    PicSidom_CT_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicSidom_CT_CC")])
    PinCo_CC_good_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinCo_CC_good")])
    PinCo_CC_low_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinCo_CC_low")])
    BL_softwood_limited_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BL_softwood_limited")])
    BL_hardwood_limited_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BL_hardwood_limited")])
    ##----##
    
    ## update the initialization##
    # PicSi_good
    if(length(PicSi_good_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicSi_good_sites,
      ClCut_age= irl_bau_pars$PicSi_good$ClCut_age,
      nTree_seedlings=irl_bau_pars$PicSi_good$nTree_seedlings,
      year_seedling=irl_bau_pars$PicSi_good$year_seedling,
      yearThin=irl_bau_pars$PicSi_good$yearThin,
      baThin=irl_bau_pars$PicSi_good$baThin,
      sitesToThin = irl_bau_pars$PicSi_good$sitesToThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PicSi_low
    if(length(PicSi_low_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicSi_low_sites,
      ClCut_age= irl_bau_pars$PicSi_low$ClCut_age,
      nTree_seedlings=irl_bau_pars$PicSi_low$nTree_seedlings,
      year_seedling=irl_bau_pars$PicSi_low$year_seedling,
      baThin=irl_bau_pars$PicSi_low$baThin,
      yearThin=irl_bau_pars$PicSi_low$yearThin,
      sitesToThin = irl_bau_pars$PicSi_low$sitesToThin,
      hThin = irl_bau_pars$PicSi_low$hThin,
      dbhThin = irl_bau_pars$PicSi_low$hThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    # PicSidom_CC
    if(length(PicSidom_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicSidom_CC_sites,
      ClCut_age= irl_bau_pars$PicSidom_CC$ClCut_age,
      nTree_seedlings=irl_bau_pars$PicSidom_CC$nTree_seedlings,
      year_seedling=irl_bau_pars$PicSidom_CC$year_seedling,
      baThin=irl_bau_pars$PicSidom_CC$baThin,
      yearThin=irl_bau_pars$PicSidom_CC$yearThin,
      sitesToThin = irl_bau_pars$PicSidom_CC$sitesToThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PicSidom_CT_CC
    if(length(PicSidom_CT_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicSidom_CT_CC_sites,
      ClCut_age= irl_bau_pars$PicSidom_CT_CC$ClCut_age,
      nTree_seedlings=irl_bau_pars$PicSidom_CT_CC$nTree_seedlings,
      year_seedling=irl_bau_pars$PicSidom_CT_CC$year_seedling,
      baThin=irl_bau_pars$PicSidom_CT_CC$baThin,
      yearThin=irl_bau_pars$PicSidom_CT_CC$yearThin,
      sitesToThin = irl_bau_pars$PicSidom_CT_CC$sitesToThin,
      hThin = irl_bau_pars$PicSidom_CT_CC$hThin,
      dbhThin = irl_bau_pars$PicSidom_CT_CC$hThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PinCo_CC_good
    if(length(PinCo_CC_good_sites)>0){initPrebas <- bau_in_thinningMatrix(
        initPrebas=initPrebas,
        siteXs=PinCo_CC_good_sites,
        ClCut_age= irl_bau_pars$PinCo_CC_good$ClCut_age,
        nTree_seedlings=irl_bau_pars$PinCo_CC_good$nTree_seedlings,
        year_seedling=irl_bau_pars$PinCo_CC_good$year_seedling)
    }
    
    # PinCo_CC_low
    if(length(PinCo_CC_low_sites)>0){initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinCo_CC_low_sites,
      ClCut_age= irl_bau_pars$PinCo_CC_low$ClCut_age,
      nTree_seedlings=irl_bau_pars$PinCo_CC_low$nTree_seedlings,
      year_seedling=irl_bau_pars$PinCo_CC_low$year_seedling)
    }

    # BL_softwood_limited
    if(length(BL_softwood_limited_sites)>0){initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BL_softwood_limited_sites,
      ClCut_age= irl_bau_pars$BL_softwood_limited$ClCut_age,
      nTree_seedlings=irl_bau_pars$BL_softwood_limited$nTree_seedlings,
      year_seedling=irl_bau_pars$BL_softwood_limited$year_seedling)
    }

    # BL_hardwood_limited
    if(length(BL_hardwood_limited_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BL_hardwood_limited_sites,
      ClCut_age= irl_bau_pars$BL_hardwood_limited$ClCut_age,
      nTree_seedlings=irl_bau_pars$BL_hardwood_limited$nTree_seedlings,
      year_seedling=irl_bau_pars$BL_hardwood_limited$year_seedling,
      baThin=irl_bau_pars$BL_hardwood_limited$baThin,
      yearThin=irl_bau_pars$BL_hardwood_limited$yearThin,
      hThin = irl_bau_pars$BL_hardwood_limited$hThin,
      dbhThin = irl_bau_pars$BL_hardwood_limited$hThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    ##----##
  }

  if(country == "Latvia" & management=="bau"){
    
    ##find the sites with alternative management##
    PicAb_dec_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                              c("PicAb_PopTr_CC","PicAb_BetSp_CC"))])
    decid_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                              c("BetSp_CC", "AlnGl_CC", "AlnIn_CC","AlnSp_BetSp_CC"))])
    PopTr_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                              c("PopTr_CC", "PopTr_BetPe_CC"))])
    PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicAb_CC")])
    PinSy_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_CC")])
    PinSy_BetPe_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_BetPe_CC")])
    PicAb_PinSy_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicAb_PinSy_CC")])
    AlnSp_PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnSp_PicAb_CC")])
    FraEx_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FraEx_CC")])
    QueSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    TilSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "TilSp_SW")])
    ##----##
    
    ## update the initialization##
    # PicAb_CC
    if(length(PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_CC_sites,
      ClCut_age= lat_bau_pars$PicAb_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$PicAb_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$PicAb_CC$year_seedling,
      baThin=lat_bau_pars$PicAb_CC$baThin,
      yearThin=lat_bau_pars$PicAb_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    # PicAb_CC
    if(length(PinSy_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_CC_sites,
      ClCut_age= lat_bau_pars$PinSy_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$PinSy_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$PinSy_CC$year_seedling,
      baThin=lat_bau_pars$PinSy_CC$baThin,
      yearThin=lat_bau_pars$PinSy_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    # PinSy_BetPe_CC
    if(length(PinSy_BetPe_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_BetPe_CC_sites,
      ClCut_age= lat_bau_pars$PinSy_BetPe_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$PinSy_BetPe_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$PinSy_BetPe_CC$year_seedling,
      baThin=lat_bau_pars$PinSy_BetPe_CC$baThin,
      yearThin=lat_bau_pars$PinSy_BetPe_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PicAb_dec
    if(length(PicAb_dec_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_dec_sites,
      ClCut_age= lat_bau_pars$PicAb_dec$ClCut_age,
      nTree_seedlings=lat_bau_pars$PicAb_dec$nTree_seedlings,
      year_seedling=lat_bau_pars$PicAb_dec$year_seedling,
      baThin=lat_bau_pars$PicAb_dec$baThin,
      yearThin=lat_bau_pars$PicAb_dec$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PicAb_PinSy_CC
    if(length(PicAb_PinSy_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_PinSy_CC_sites,
      ClCut_age= lat_bau_pars$PicAb_PinSy_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$PicAb_PinSy_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$PicAb_PinSy_CC$year_seedling,
      baThin=lat_bau_pars$PicAb_PinSy_CC$baThin,
      yearThin=lat_bau_pars$PicAb_PinSy_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # decid
    if(length(decid_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=decid_sites,
      ClCut_age= lat_bau_pars$decid$ClCut_age,
      nTree_seedlings=lat_bau_pars$decid$nTree_seedlings,
      year_seedling=lat_bau_pars$decid$year_seedling,
      baThin=lat_bau_pars$decid$baThin,
      yearThin=lat_bau_pars$decid$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # AlnSp_PicAb_CC
    if(length(AlnSp_PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnSp_PicAb_CC_sites,
      ClCut_age= lat_bau_pars$AlnSp_PicAb_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$AlnSp_PicAb_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$AlnSp_PicAb_CC$year_seedling,
      baThin=lat_bau_pars$AlnSp_PicAb_CC$baThin,
      yearThin=lat_bau_pars$AlnSp_PicAb_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # PopTr
    if(length(PopTr_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_sites,
      ClCut_age= lat_bau_pars$PopTr$ClCut_age,
      nTree_seedlings=lat_bau_pars$PopTr$nTree_seedlings,
      year_seedling=lat_bau_pars$PopTr$year_seedling,
      baThin=lat_bau_pars$PopTr$baThin,
      yearThin=lat_bau_pars$PopTr$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # FraEx_CC
    if(length(FraEx_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=FraEx_CC_sites,
      ClCut_age= lat_bau_pars$FraEx_CC$ClCut_age,
      nTree_seedlings=lat_bau_pars$FraEx_CC$nTree_seedlings,
      year_seedling=lat_bau_pars$FraEx_CC$year_seedling,
      baThin=lat_bau_pars$FraEx_CC$baThin,
      yearThin=lat_bau_pars$FraEx_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lat_bau_pars$FraEx_CC$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # QueSp_SW
    if(length(QueSp_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_SW_sites,
      ClCut_age= lat_bau_pars$QueSp_SW$ClCut_age,
      nTree_seedlings=lat_bau_pars$QueSp_SW$nTree_seedlings,
      year_seedling=lat_bau_pars$QueSp_SW$year_seedling,
      baThin=lat_bau_pars$QueSp_SW$baThin,
      yearThin=lat_bau_pars$QueSp_SW$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lat_bau_pars$QueSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # TilSp_SW
    if(length(TilSp_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=TilSp_SW_sites,
      ClCut_age= lat_bau_pars$TilSp_SW$ClCut_age,
      nTree_seedlings=lat_bau_pars$TilSp_SW$nTree_seedlings,
      year_seedling=lat_bau_pars$TilSp_SW$year_seedling,
      baThin=lat_bau_pars$TilSp_SW$baThin,
      yearThin=lat_bau_pars$TilSp_SW$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lat_bau_pars$TilSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
    
    
    ##----##
  }
  
  
  if(country == "Lithuania" & management=="bau"){
    
    ##find the sites with alternative management##
    PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicAb_CC")])
    PicAb_GSC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PicAb_GSC")])
    
    PicAb_BeSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% c("PicAb_BeSp_CC","PicAb_BetSp_CC"))])
    PinSy_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_CC")])
    PinSy_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_SW")])
    PinSy_BetPe_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinSy_BetPe_CC")])
    BetSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BetSp_CC")])
    BetSp_AlnGl_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "BetSp_AlnGl_SW")])
    AlnGl_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnGl_CC")])
    AlnIn_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnIn_CC")])
    AlnSp_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                          c("AlnSp_PinSy_CC", "AlnSp_BetSp_CC", "BetPe_AlnSp_CC"))])
    PopTr_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                          c("PopTr_CC","PopTr_BetPe_CC"))])
    PopTr_AlnIn_GSC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PopTr_AlnIn_GSC")])
    FraEx_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FraEx_CC")])
    QueSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_CC")])
    QueSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    ##----##
    
    ## update the initialization##
    # PicAb_CC
    if(length(PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_CC_sites,
      ClCut_age= lit_bau_pars$PicAb_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PicAb_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$PicAb_CC$year_seedling,
      baThin=lit_bau_pars$PicAb_CC$baThin,
      yearThin=lit_bau_pars$PicAb_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PicAb_GSC
    if(length(PicAb_GSC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_GSC_sites,
      ClCut_age= lit_bau_pars$PicAb_GSC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PicAb_GSC$nTree_seedlings,
      year_seedling=lit_bau_pars$PicAb_GSC$year_seedling,
      baThin=lit_bau_pars$PicAb_GSC$baThin,
      yearThin=lit_bau_pars$PicAb_GSC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PicAb_BeSp_CC
    if(length(PicAb_BeSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_BeSp_CC_sites,
      ClCut_age= lit_bau_pars$PicAb_BeSp_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PicAb_BeSp_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$PicAb_BeSp_CC$year_seedling,
      baThin=lit_bau_pars$PicAb_BeSp_CC$baThin,
      yearThin=lit_bau_pars$PicAb_BeSp_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PinSy_CC
    if(length(PinSy_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_CC_sites,
      ClCut_age= lit_bau_pars$PinSy_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PinSy_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$PinSy_CC$year_seedling,
      baThin=lit_bau_pars$PinSy_CC$baThin,
      yearThin=lit_bau_pars$PinSy_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PinSy_SW
    if(length(PinSy_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_SW_sites,
      ClCut_age= lit_bau_pars$PinSy_SW$ClCut_age,
      nTree_seedlings=lit_bau_pars$PinSy_SW$nTree_seedlings,
      year_seedling=lit_bau_pars$PinSy_SW$year_seedling,
      baThin=lit_bau_pars$PinSy_SW$baThin,
      yearThin=lit_bau_pars$PinSy_SW$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PinSy_BetPe_CC
    if(length(PinSy_BetPe_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PinSy_BetPe_CC_sites,
      ClCut_age= lit_bau_pars$PinSy_BetPe_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PinSy_BetPe_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$PinSy_BetPe_CC$year_seedling,
      baThin=lit_bau_pars$PinSy_BetPe_CC$baThin,
      yearThin=lit_bau_pars$PinSy_BetPe_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # BetSp_CC
    if(length(BetSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BetSp_CC_sites,
      ClCut_age= lit_bau_pars$BetSp_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$BetSp_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$BetSp_CC$year_seedling,
      baThin=lit_bau_pars$BetSp_CC$baThin,
      yearThin=lit_bau_pars$BetSp_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # BetSp_AlnGl_SW
    if(length(BetSp_AlnGl_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BetSp_AlnGl_SW_sites,
      ClCut_age= lit_bau_pars$BetSp_AlnGl_SW$ClCut_age,
      nTree_seedlings=lit_bau_pars$BetSp_AlnGl_SW$nTree_seedlings,
      year_seedling=lit_bau_pars$BetSp_AlnGl_SW$year_seedling,
      baThin=lit_bau_pars$BetSp_AlnGl_SW$baThin,
      yearThin=lit_bau_pars$BetSp_AlnGl_SW$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # AlnGl_CC
    if(length(AlnGl_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnGl_CC_sites,
      ClCut_age= lit_bau_pars$AlnGl_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$AlnGl_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$AlnGl_CC$year_seedling,
      baThin=lit_bau_pars$AlnGl_CC$baThin,
      yearThin=lit_bau_pars$AlnGl_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # AlnIn_CC
    if(length(AlnIn_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnIn_CC_sites,
      ClCut_age= lit_bau_pars$AlnIn_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$AlnIn_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$AlnIn_CC$year_seedling,
      baThin=lit_bau_pars$AlnIn_CC$baThin,
      yearThin=lit_bau_pars$AlnIn_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # AlnSp
    if(length(AlnSp_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=AlnSp_sites,
      ClCut_age= lit_bau_pars$AlnSp$ClCut_age,
      nTree_seedlings=lit_bau_pars$AlnSp$nTree_seedlings,
      year_seedling=lit_bau_pars$AlnSp$year_seedling,
      baThin=lit_bau_pars$AlnSp$baThin,
      yearThin=lit_bau_pars$AlnSp$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PopTr
    if(length(PopTr_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_sites,
      ClCut_age= lit_bau_pars$PopTr$ClCut_age,
      nTree_seedlings=lit_bau_pars$PopTr$nTree_seedlings,
      year_seedling=lit_bau_pars$PopTr$year_seedling,
      baThin=lit_bau_pars$PopTr$baThin,
      yearThin=lit_bau_pars$PopTr$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # PopTr_AlnIn_GSC
    if(length(PopTr_AlnIn_GSC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_AlnIn_GSC_sites,
      ClCut_age= lit_bau_pars$PopTr_AlnIn_GSC$ClCut_age,
      nTree_seedlings=lit_bau_pars$PopTr_AlnIn_GSC$nTree_seedlings,
      year_seedling=lit_bau_pars$PopTr_AlnIn_GSC$year_seedling,
      baThin=lit_bau_pars$PopTr_AlnIn_GSC$baThin,
      yearThin=lit_bau_pars$PopTr_AlnIn_GSC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # FraEx_CC
    if(length(FraEx_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=FraEx_CC_sites,
      ClCut_age= lit_bau_pars$FraEx_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$FraEx_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$FraEx_CC$year_seedling,
      baThin=lit_bau_pars$FraEx_CC$baThin,
      yearThin=lit_bau_pars$FraEx_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lit_bau_pars$QueSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
    # QueSp_CC
    if(length(QueSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_CC_sites,
      ClCut_age= lit_bau_pars$QueSp_CC$ClCut_age,
      nTree_seedlings=lit_bau_pars$QueSp_CC$nTree_seedlings,
      year_seedling=lit_bau_pars$QueSp_CC$year_seedling,
      baThin=lit_bau_pars$QueSp_CC$baThin,
      yearThin=lit_bau_pars$QueSp_CC$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lit_bau_pars$QueSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
    # QueSp_SW
    if(length(QueSp_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_SW_sites,
      ClCut_age= lit_bau_pars$QueSp_SW$ClCut_age,
      nTree_seedlings=lit_bau_pars$QueSp_SW$nTree_seedlings,
      year_seedling=lit_bau_pars$QueSp_SW$year_seedling,
      baThin=lit_bau_pars$QueSp_SW$baThin,
      yearThin=lit_bau_pars$QueSp_SW$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = lit_bau_pars$QueSp_SW$dens_after_Thin,
      acThin = -999,
      pHarvTreeThin = 1)
  }
  
  if(country == "United Kingdom" & management=="bau"){
    
    ##find the sites with alternative management##
    PopTr_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("FRA_PopTr","FRA_PopX"))])
    BL_softwood_limited_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("FRA_BLS09","FRA_BLS01","IRL_BLS02","FRA_BLS02","FRA_BLSCON02","FRA_CONBLS02"))])
    BL_hardwood_limited_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("FRA_MIX06","FRA_BLSCON03","FRA_CasSa","FRA_BetPe","FRA_MIX03","IRL_BLS01",
                                  "FRA_MIX04", "FRA_MIX11", "IRL_MIX02"))])
    PinCo_CC_good_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("FRA_PinNi","IRL_MIX01"))])
    PinCo_CC_low_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("FRA_CONBLS05","IRL_CON02"))])
    PicSi_good_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                c("IRL_PicSi"))])
    FagSy_SW_High_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                      c("FRA_FagSy","FRA_FagSy-AbiAl"))])
    conif_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                    c("IRL_CON01", "FRA_PinSy","FRA_CON03","IRL_CONBLS02"))])
    QueSp_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FRA_QueRo-PinSy")])
    ##----##
    
    ## update the initialization##
    # PopTr
    if(length(PopTr_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PopTr_sites,
      ClCut_age= lit_bau_pars$PopTr$ClCut_age,
      nTree_seedlings=lit_bau_pars$PopTr$nTree_seedlings,
      year_seedling=lit_bau_pars$PopTr$year_seedling,
      baThin=lit_bau_pars$PopTr$baThin,
      yearThin=lit_bau_pars$PopTr$yearThin,
      hThin = 1,
      dbhThin = 1,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    # BL_softwood_limited
    if(length(BL_softwood_limited_sites)>0){initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BL_softwood_limited_sites,
      ClCut_age= irl_bau_pars$BL_softwood_limited$ClCut_age,
      nTree_seedlings=irl_bau_pars$BL_softwood_limited$nTree_seedlings,
      year_seedling=irl_bau_pars$BL_softwood_limited$year_seedling)
    }
    # BL_hardwood_limited
    if(length(BL_hardwood_limited_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BL_hardwood_limited_sites,
      ClCut_age= irl_bau_pars$BL_hardwood_limited$ClCut_age,
      nTree_seedlings=irl_bau_pars$BL_hardwood_limited$nTree_seedlings,
      year_seedling=irl_bau_pars$BL_hardwood_limited$year_seedling,
      baThin=irl_bau_pars$BL_hardwood_limited$baThin,
      yearThin=irl_bau_pars$BL_hardwood_limited$yearThin,
      hThin = irl_bau_pars$BL_hardwood_limited$hThin,
      dbhThin = irl_bau_pars$BL_hardwood_limited$hThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
  # PinCo_CC_good
  if(length(PinCo_CC_good_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas=initPrebas,
    siteXs=PinCo_CC_good_sites,
    ClCut_age= irl_bau_pars$PinCo_CC_good$ClCut_age,
    nTree_seedlings=irl_bau_pars$PinCo_CC_good$nTree_seedlings,
    year_seedling=irl_bau_pars$PinCo_CC_good$year_seedling)
  
  # PinCo_CC_low
  if(length(PinCo_CC_low_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas=initPrebas,
    siteXs=PinCo_CC_low_sites,
    ClCut_age= irl_bau_pars$PinCo_CC_low$ClCut_age,
    nTree_seedlings=irl_bau_pars$PinCo_CC_low$nTree_seedlings,
    year_seedling=irl_bau_pars$PinCo_CC_low$year_seedling)
  
  # PicSi_good
  if(length(PicSi_good_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas=initPrebas,
    siteXs=PicSi_good_sites,
    ClCut_age= irl_bau_pars$PicSi_good$ClCut_age,
    nTree_seedlings=irl_bau_pars$PicSi_good$nTree_seedlings,
    year_seedling=irl_bau_pars$PicSi_good$year_seedling,
    yearThin=irl_bau_pars$PicSi_good$yearThin,
    baThin=irl_bau_pars$PicSi_good$baThin,
    sitesToThin = irl_bau_pars$PicSi_good$sitesToThin,
    hThin = 1,
    dbhThin = 1,
    hcThin = 1,
    fracThin = 1,
    dens_after_Thin = -999,
    acThin = -999,
    pHarvTreeThin = 1)
  # FagSy
  if(length(FagSy_SW_High_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas = initPrebas,
    siteXs = FagSy_SW_High_sites,
    ClCut_age = den_bau_pars$FagSy_SW_High$ClCut_age,
    nTree_seedlings = den_bau_pars$FagSy_SW_High$nTree_seedlings,
    year_seedling=den_bau_pars$FagSy_SW_High$year_seedling,
    yearThin=den_bau_pars$FagSy_SW_High$yearThin,
    baThin=den_bau_pars$FagSy_SW_High$baThin,
    hThin = den_bau_pars$FagSy_SW_High$hThin,
    dbhThin = den_bau_pars$FagSy_SW_High$dbhThin,
    hcThin = 1,
    fracThin = 1,
    dens_after_Thin = -999,
    acThin = -999,
    pHarvTreeThin = 1)
# conif
  if(length(conif_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas=initPrebas,
    siteXs=conif_sites,
    ClCut_age= est_bau_pars$conif$ClCut_age,
    nTree_seedlings=est_bau_pars$conif$nTree_seedlings,
    year_seedling=est_bau_pars$conif$year_seedling,
    yearThin=est_bau_pars$conif$yearThin,
    baThin=est_bau_pars$conif$baThin,
    hThin = 1,
    dbhThin = 1,
    hcThin = 1,
    fracThin = 1,
    dens_after_Thin = -999,
    acThin = -999,
    pHarvTreeThin = 1)
  #QueSp_CC
  if(length(QueSp_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
    initPrebas=initPrebas,
    siteXs=QueSp_CC_sites,
    ClCut_age= est_bau_pars$QueSp_CC$ClCut_age,
    nTree_seedlings=est_bau_pars$QueSp_CC$nTree_seedlings,
    year_seedling=est_bau_pars$QueSp_CC$year_seedling,
    yearThin=est_bau_pars$QueSp_CC$yearThin,
    baThin=est_bau_pars$QueSp_CC$baThin,
    hThin = 1,
    dbhThin = 1,
    hcThin = 1,
    fracThin = 1,
    dens_after_Thin = est_bau_pars$QueSp_CC$dens_after_Thin,
    acThin = -999,
    pHarvTreeThin = 1)
 }
  if(country == "Germany" & management=="bau"){
    
    ##find the sites with alternative management##
    General_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                      c("General_CC"))])
    BLSCON_GSC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                      c("BLSCON_GSC"))])
    PicAb_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                      c("PicAb_CC"))])
    QueSp_FagSy_SS_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_FagSy_SS")])
    QueSp_SS_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SS")])
    ##----##
    
    ## update the initialization##
    # General_CC
    if(length(General_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=General_CC_sites,
      ClCut_age= de_bau_pars$General_CC$ClCut_age,
      nTree_seedlings=de_bau_pars$General_CC$nTree_seedlings,
      year_seedling=de_bau_pars$General_CC$year_seedling,
      baThin=de_bau_pars$General_CC$baThin,
      yearThin=de_bau_pars$General_CC$yearThin,
      hThin = de_bau_pars$General_CC$hThin,
      dbhThin = de_bau_pars$General_CC$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # BLSCON_GSC
    if(length(BLSCON_GSC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=BLSCON_GSC_sites,
      ClCut_age= de_bau_pars$BLSCON_GSC$ClCut_age,
      nTree_seedlings=de_bau_pars$BLSCON_GSC$nTree_seedlings,
      year_seedling=de_bau_pars$BLSCON_GSC$year_seedling,
      baThin=de_bau_pars$BLSCON_GSC$baThin,
      yearThin=de_bau_pars$BLSCON_GSC$yearThin,
      hThin = de_bau_pars$BLSCON_GSC$hThin,
      dbhThin = de_bau_pars$BLSCON_GSC$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    

    # PicAb_CC
    if(length(PicAb_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=PicAb_CC_sites,
      ClCut_age= de_bau_pars$PicAb_CC$ClCut_age,
      nTree_seedlings=de_bau_pars$PicAb_CC$nTree_seedlings,
      year_seedling=de_bau_pars$PicAb_CC$year_seedling,
      baThin=de_bau_pars$PicAb_CC$baThin,
      yearThin=de_bau_pars$PicAb_CC$yearThin,
      hThin = de_bau_pars$PicAb_CC$hThin,
      dbhThin = de_bau_pars$PicAb_CC$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    # QueSp_FagSy_SS
    if(length(QueSp_FagSy_SS_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_FagSy_SS_sites,
      ClCut_age= de_bau_pars$QueSp_FagSy_SS$ClCut_age,
      nTree_seedlings=de_bau_pars$QueSp_FagSy_SS$nTree_seedlings,
      year_seedling=de_bau_pars$QueSp_FagSy_SS$year_seedling,
      baThin=de_bau_pars$QueSp_FagSy_SS$baThin,
      yearThin=de_bau_pars$QueSp_FagSy_SS$yearThin,
      hThin = de_bau_pars$QueSp_FagSy_SS$hThin,
      dbhThin = de_bau_pars$QueSp_FagSy_SS$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)

    # QueSp_SS
    if(length(QueSp_SS_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_SS_sites,
      ClCut_age= de_bau_pars$QueSp_SS$ClCut_age,
      nTree_seedlings=de_bau_pars$QueSp_SS$nTree_seedlings,
      year_seedling=de_bau_pars$QueSp_SS$year_seedling,
      baThin=de_bau_pars$QueSp_SS$baThin,
      yearThin=de_bau_pars$QueSp_SS$yearThin,
      hThin = de_bau_pars$QueSp_SS$hThin,
      dbhThin = de_bau_pars$QueSp_SS$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
  }
  
  if(country == "Poland" & management=="bau"){
    
    ##find the sites with alternative management##
    FagSy_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                                                     c("FagSy_SW"))])
    QueSp_SW_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                                                     c("QueSp_SW"))])
    CON_CC_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man %in% 
                                                                   c("CON_CC"))])
    ##----##
    
    ## update the initialization##
    # FagSy_SW
    if(length(FagSy_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=FagSy_SW_sites,
      ClCut_age= pl_bau_pars$FagSy_SW$ClCut_age,
      nTree_seedlings=pl_bau_pars$FagSy_SW$nTree_seedlings,
      year_seedling=pl_bau_pars$FagSy_SW$year_seedling,
      baThin=pl_bau_pars$FagSy_SW$baThin,
      yearThin=pl_bau_pars$FagSy_SW$yearThin,
      hThin = pl_bau_pars$FagSy_SW$hThin,
      dbhThin = pl_bau_pars$FagSy_SW$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    # QueSp_SW
    if(length(QueSp_SW_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=QueSp_SW_sites,
      ClCut_age= pl_bau_pars$QueSp_SW$ClCut_age,
      nTree_seedlings=pl_bau_pars$QueSp_SW$nTree_seedlings,
      year_seedling=pl_bau_pars$QueSp_SW$year_seedling,
      baThin=pl_bau_pars$QueSp_SW$baThin,
      yearThin=pl_bau_pars$QueSp_SW$yearThin,
      hThin = pl_bau_pars$QueSp_SW$hThin,
      dbhThin = pl_bau_pars$QueSp_SW$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
    
    # CON_CC
    if(length(CON_CC_sites)>0) initPrebas <- bau_in_thinningMatrix(
      initPrebas=initPrebas,
      siteXs=CON_CC_sites,
      ClCut_age= pl_bau_pars$CON_CC$ClCut_age,
      nTree_seedlings=pl_bau_pars$CON_CC$nTree_seedlings,
      year_seedling=pl_bau_pars$CON_CC$year_seedling,
      baThin=pl_bau_pars$CON_CC$baThin,
      yearThin=pl_bau_pars$CON_CC$yearThin,
      hThin = pl_bau_pars$CON_CC$hThin,
      dbhThin = pl_bau_pars$CON_CC$dbhThin,
      hcThin = 1,
      fracThin = 1,
      dens_after_Thin = -999,
      acThin = -999,
      pHarvTreeThin = 1)
    
  }
  
  return(initPrebas)    
}

