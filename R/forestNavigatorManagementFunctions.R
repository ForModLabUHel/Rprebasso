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
                                  ClCut_age, #130 ####if ClCut_age is NA the ingrowth layers are not removed at clearcut
                                  nTree_seedlings,
                                  year_seedling,
                                  yearThin,
                                  baThin,
                                  hThin = 1,
                                  dbhThin = 1, 
                                  hcThin = 1,
                                  fracThin = 1,
                                  dens_after_Thin = -999,
                                  acThin = -999,
                                  pHarvTreeThin = 1
){
  
  
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
  
  initPrebas$initClearcut[siteXs,3] <- pi*(initPrebas$initClearcut[siteXs,2]/200)^2 * nTree_seedlings
  initPrebas$initClearcut[siteXs,5] <- year_seedling
  initPrebas$defaultThin[siteXs] <- 0
  if(is.na(ClCut_age)){
    initPrebas$ClCut[siteXs] <- 0
    ClCut_age <- thin_def[nrow(thin_def),1]
  }else{
    ####if ClCut_age is active remove all trees also from ingrowth
    initPrebas$ClCut[siteXs] <- 1
    initPrebas$clct_pars[siteXs,,2] <- ClCut_age
    initPrebas$clct_pars[siteXs,,c(1,3)] <- 9999
    thin_def <- thin_def[-nrow(thin_def),]
  }
  if(!is.matrix(thin_def)) thin_def <- matrix(thin_def,nrow=1)
  
  nlayers <- dim(initPrebas$multiInitVar)[3]
  maxYear <- initPrebas$maxYears
  
  layerX <- 1
  thinningX <- list()
  for(ij in siteXs){
    
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
                                     est_bau_pars=est_bau_pars_def){
  available_countries <- c("Sweden","Finland","Estonia")
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
    
    # if(dim(initPrebas$ftTapioPar)[2] < 12){
    #   dims <- dim(initPrebas$ftTapioPar)
    #   ftTapioPar <- array(999,dim = c(5,12,3,7))
    #   ftTapioPar[1:dims[1],1:dims[2],1:dims[3],1:dims[4]] <- initPrebas$ftTapioPar
    #   initPrebas$ftTapioPar <- ftTapioPar
    # }
    # if(dim(initPrebas$tTapioPar)[2] < 12){
    #   dims <- dim(initPrebas$tTapioPar)
    #   tTapioPar <- array(999,dim = c(5,12,3,7))
    #   tTapioPar[1:dims[1],1:dims[2],1:dims[3],1:dims[4]] <- initPrebas$tTapioPar
    #   initPrebas$tTapioPar <- tTapioPar
    # }
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
  
  return(initPrebas)    
}

