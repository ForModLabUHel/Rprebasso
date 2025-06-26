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
    yearSims <- thinning_def_X[,1] - initPrebas$multiInitVar[ij,2,layerX]
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
    yearSims <- thinning_defX[,1] - initPrebas$multiInitVar[ij,2,layerX]
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
    yearSims <- thinning_def_X[,1] - initPrebas$multiInitVar[ij,2,layerX]
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
forest_management_update <- function(initPrebas, forest_type_management_tab, country, management){
  available_countries <- c("sweden","finland")
  available_managements <- c("bau", "noman")
  if(!country %in% available_countries) stop(cat("This country: ", country,
                                                 " is not between the available countries: ", available_countries,fill = TRUE))
  if(!management %in% available_managements) stop(cat("This management: ", management, 
                                                      " is not between the available managements: ", available_managements,fill = TRUE))
  if(country == "sweden" & management=="bau"){
    ##find the sites with alternative management##
    pop_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PopTr_CC")])
    alnus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "AlnSp_CC")])
    quercus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "QueSp_SW")])
    pinco_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "PinCo_CC")])
    fagus_sites <- sort(forest_type_management_tab$site[which(forest_type_management_tab$for_man == "FagSy_SW")])
    ##----##
    
    ## update the initialization##
    initPrebas <- sw_bau_pop(initPrebas,pop_sites)
    initPrebas <- sw_bau_AlnSp(initPrebas,alnus_sites)
    initPrebas <- sw_bau_QueRob(initPrebas,quercus_sites)
    initPrebas <- sw_bau_piCo(initPrebas,pinco_sites)
    initPrebas <- sw_bau_fagsy(initPrebas,fagus_sites)
    ##----##
  }
  return(initPrebas)    
}
