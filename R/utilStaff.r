###function to select subset sites in a initPrebas object.
###useful to test the model locally . you can select few sites and reinitilize the model just for those sites to get lighter computational load
#' Title
#'
#' @param initPrebas multisite Initialization object
#' @param sitex sites to be selected
#'
#' @returns a multisite prebas initialization object with just the sitex sites
#' @export
#'
#' @examples
initPrebas_subsetter <- function(initPrebas,sitex){
  notsubset <- c("nSites", "nClimID", "maxYears",
                 "maxThin","pCROBAS", "allSp","maxNlayers",
                 "ETSy", "P0y", "weather", "DOY", "pPRELES", "etmodel", 
                 "pYASSO",  "pAWEN",   "weatherYasso",
                 "litterSize","smoothP0","smoothETS", "tapioPars", 
                 "ftTapioPar","tTapioPar", "GVrun",
                 "mortMod", "ECMmod",  "pECMmod",
                 "layerPRELES", "LUEtrees","LUEgv",  
                 "alpharNcalc", "siteInfoDist", 
                 "dist_flag", "CO2model" )
  
  # Subset all elements except 'notX'
  subset_initPrebas <- lapply(names(initPrebas), function(name) {
    x <- initPrebas[[name]]
    if (!name %in% notsubset) {
      if (is.vector(x)) {
        x = x[sitex]
      } else if (length(dim(x))==2) {
        x = x[sitex,, drop = FALSE]
      } else if (length(dim(x))==3) {
        x = x[sitex,,, drop = FALSE]
      } else if (length(dim(x))==4) {
        x = x[sitex,,,, drop = FALSE]
      } else if (length(dim(x))==5) {
        x = x[sitex,,,,, drop = FALSE]
      }
    }
    return(x)
  })
  
  # Preserve names
  names(subset_initPrebas) <- names(initPrebas)
  
  subset_initPrebas$nSites <- dim(subset_initPrebas$multiOut)[1]
  # process weather inputs
  subset_initPrebas$nClimID <- length(unique(subset_initPrebas$siteInfo[,2]))
  climIDsunique <- sort(unique(subset_initPrebas$siteInfo[,2]))
  if(length(climIDsunique)==1){
    subset_initPrebas$siteInfo[,2] <- match(subset_initPrebas$siteInfo[,2], climIDsunique)
    subset_initPrebas$siteInfo[2,2] <- 2
    climIDsunique <- c(climIDsunique,climIDsunique)
    subset_initPrebas$ETSy <- subset_initPrebas$ETSy[climIDsunique,]
    subset_initPrebas$P0y <- subset_initPrebas$P0y[climIDsunique,,]
    subset_initPrebas$weather <- subset_initPrebas$weather[climIDsunique,,,]
    subset_initPrebas$weatherYasso <- subset_initPrebas$weatherYasso[climIDsunique,,]
  }else{
    subset_initPrebas$siteInfo[,2] <- match(subset_initPrebas$siteInfo[,2], climIDsunique)
    subset_initPrebas$ETSy <- subset_initPrebas$ETSy[climIDsunique,]
    subset_initPrebas$P0y <- subset_initPrebas$P0y[climIDsunique,,]
    subset_initPrebas$weather <- subset_initPrebas$weather[climIDsunique,,,]
    subset_initPrebas$weatherYasso <- subset_initPrebas$weatherYasso[climIDsunique,,]
  }
  
  subset_initPrebas$maxYears <- max(subset_initPrebas$nYears)
  subset_initPrebas$maxThin <- max(subset_initPrebas$nThinning)
  subset_initPrebas$thinning <- subset_initPrebas$thinning[,1:subset_initPrebas$maxThin,]
  subset_initPrebas$maxNlayers <- max(subset_initPrebas$nLayers)
  return(subset_initPrebas)
}


###regression model for the drained peatland forested sites (paper reference)
peat_regression_model <- function(BA,Tseason,siteType,peat_reg_pars=peat_regression_pars,maxsiteType = 5){
  siteType <- min(siteType,maxsiteType)
  p_st <- peat_reg_pars$p_st[siteType]
  p_ba <- peat_reg_pars$p_ba
  p_Tseason <- peat_reg_pars$p_Tseason
  rh <- p_st + p_ba * BA + p_Tseason * Tseason
  return(rh)
}

###multisite version of regression model for the drained peatland forested sites (paper reference)
peat_regression_model_multiSite <- function(modOut,peat_sites,peat_regr_pars=peat_regression_pars,max_siteType = 5){
  
  peat_sites <- sort(peat_sites)
  
  siteTypes <- modOut$siteInfo[peat_sites,3] 
  if(modOut$maxNlayers==1){
    BA <- modOut$multiOut[peat_sites,,13,1,1]
  }else{
    BA <- apply(modOut$multiOut[peat_sites,,13,,1],1:2,sum)
  }
  
  Tseason <- apply(modOut$weather[modOut$siteInfo[peat_sites,"climID"],,121:304,2],1:2,mean)
  
  Rh_peat <- peat_regression_model(BA,Tseason,siteTypes,peat_regr_pars,max_siteType) * 12/44 #converts CO2 equivalents to gC m-2-y
  Rh_mineral <- apply(modOut$multiOut[peat_sites,,45,,1],1:2,sum)
  NEP_mineral <- apply(modOut$multiOut[peat_sites,,46,,1],1:2,sum)
  
  Rh_net <- Rh_peat - Rh_mineral
  NEP_peat <- NEP_mineral - Rh_net
  
  ####repleace values
  modOut$multiOut[peat_sites,,45,,1] <- 0
  modOut$multiOut[peat_sites,,46,,1] <- 0
  modOut$multiOut[peat_sites,,45,1,1] <- Rh_peat
  modOut$multiOut[peat_sites,,46,1,1] <- NEP_peat
  
  return(modOut)
}


# Derive HC from foliage mass this function is useful for databases like MS-NFI
fHc_fol <- function(D,B,H, pCROB){
  
  # B basal area (m2)
  # D diameter (cm)
  # H height (m)
  
  VHc <- c()
  
  for(specid in 1:3){
    
    if(specid == 1){Wf <- 792.307 + 47.304 * B} # Pine
    if(specid == 2) {Wf <- 15000 * (B / (B + 40))} # Spruce
    if(specid == 3){Wf <- 0.5 *(792.307 + 47.304 * B)} # Birch
    
    
    z <- pCROB[11,specid]
    ksi <- pCROB[38,specid]
    N <- B/(pi/4*(D/100)^2) # Number of trees
    
    wf <- Wf/N
    L <- (wf/ksi)^(1/z)
    VHc[specid] <- max(0.2,H-L) # Hc height to crown base
    
  }
  return(VHc)
}




#' Initial age of the seedlings
#'
#' @param SiteType Site type
#' @param ETS Effective temperature sum
#'
#' @return Initial age of the seedlings
#' @export
#'
#' @examples initialAgeSeedl(2, 1200)
initialAgeSeedl <- function(SiteType,ETS){
  round(6 + 2*SiteType - 0.005*ETS + 2.25) ##Initial age
}

#' initBiomasses
#'
#' @param pCro matrix of CROBAS parameters where each column correspond to a different species
#' @param initVarX initial state of the forests.
#'
#' @return the biomasses at initialization 
#' @export
#'
#' @examples
initBiomasses <- function(pCro,initVarX){
  initVarX<-as.matrix(initVarX) #change vector to matrix when maxlayer=1
  siteType <- initVarX[8,1]
  layerXs <- which(initVarX[1,] %in% 1:ncol(pCROB))
  ##set parameters
  par_betab <- pCro[13,initVarX[1,layerXs]]
  par_x <- pCro[19,initVarX[1,layerXs]]
  par_beta0 <- pCro[12,initVarX[1,layerXs]]
  par_betas <- pCro[14,initVarX[1,layerXs]]
  par_mf <- pCro[8,initVarX[1,layerXs]]
  par_mr <- pCro[9,initVarX[1,layerXs]]
  par_mw <- pCro[10,initVarX[1,layerXs]]
  par_c <- pCro[7,initVarX[1,layerXs]]
  par_rhof <- pCro[15,initVarX[1,layerXs]]
  par_rhow <- pCro[2,initVarX[1,layerXs]]
  par_S_branchMod <- pCro[27,initVarX[1,layerXs]]
  gammaC <- 0. #initVarX[8,]
  Tbd <- 10 #####to include in the parameters
  
  ###set variables
  A <- initVarX[7,layerXs]
  ba <- initVarX[5,layerXs]; d <- initVarX[4,layerXs]
  N <- ba/(pi*((d/2/100)^2))
  h = initVarX[3,layerXs]; hc <- initVarX[6,layerXs]
  B = ba/N
  Lc <- h - hc
  
  ### Fine root allocation of early growth  
  age_factor <-  (1. - (1- pCro[44,initVarX[1,layerXs]])/ (1. + exp((-h+ pCro[45,initVarX[1,layerXs]])/ pCro[46,initVarX[1,layerXs]])))/ pCro[44,initVarX[1,layerXs]] 
  par_alfar <- pCro[20+pmin(siteType,5),initVarX[1,layerXs]]* age_factor
  par_rhor <- par_alfar * par_rhof
  beta0 <- par_beta0 * age_factor
  
  betab =  par_betab * Lc^(par_x-1)
  beta1 = (beta0 + betab + par_betas) 
  beta2 = 1. - betab - par_betas 		
  betaC = (beta1 + gammaC * beta2) / par_betas
  wf_STKG <- pmax(0.,par_rhof * A * N)
  W_froot = pmax(0.,par_rhor * A * N)  ##to check  ##newX
  W_wsap = pmax(0.,par_rhow * A * N * (beta1 * h + beta2 * hc)) ##newX
  W_c = pmax(0.,par_rhow * A * N * hc) #sapwood stem below Crown
  W_s = pmax(0.,par_rhow * A * N * par_betas * Lc) #sapwood stem within crown
  W_branch =  pmax(0.,par_rhow * A * N * betab * Lc) #branches biomass
  Wsh = pmax((A+B+sqrt(A*B)) * hc * par_rhow * N/2.9 - W_c,0) #initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
  #initialize Wdb dead branches biomass
  Wdb = pmax(0,ifelse(par_S_branchMod == 1.,Tbd * W_branch * ((0.0337+0.000009749*N)*exp(-0.00456*d^2)+0.00723),
         Tbd * W_branch *((-0.00513+0.000012*N)*exp((0.00000732-0.000000764*N)*d^2)+0.00467)))
  W_stem = pmax(0.,W_c + W_s + Wsh)

  W_crs = par_rhow * beta0 * A * h * N
  W_crh = Wsh * beta0
  W_croot = W_crh + W_crs
  #W_croot = pmax(0.,(par_rhow * Lc * beta0 * A / par_betas * N + (W_c + Wsh) * beta0)) #coarse root biomass
  V = W_stem / par_rhow
  biomassesX <- rbind(wf_STKG,W_froot,W_wsap,W_c,W_s,W_branch,W_croot,Wsh,Wdb,W_stem,V,W_crh)
  biomasses<-matrix(0,nrow = nrow(biomassesX), ncol = ncol(initVarX))
  biomasses[,layerXs]<-biomassesX
  # biomasses[which(is.na(biomasses))] <- 0.
  return(biomasses)
}

#### function to calculate initial sapwood area at crown base (A)
#' compA
#'
#' @param inputs 
#'
#' @return
#' @export
#'
#' @examples
compA <- function(inputs){
  p_ksi = inputs[1]
  p_rhof = inputs[2]
  p_z <- inputs[3]
  Lc = inputs[4]
  A <- max(0.,p_ksi/p_rhof * Lc^p_z)
  return(A)
}

varNames  <- c('siteID','gammaC','sitetype','species','ETS' ,'P0','age', 'DeadWoodVolume', 'Respi_tot','GPPTot/1000',
               'H','D', 'BA','Hc_base','Cw','Ac','N','npp','leff','keff','lproj','ET_preles','weight',
               'Wbranch',"WfineRoots",'Litter_fol','Litter_fr','Litter_fWoody','Litter_cWoody','V',
               'Wstem','W_croot','wf_STKG', 'wf_treeKG','B_tree','Light',"VroundWood","WroundWood","soilC",
               "aSW","dH","Vmort","grossGrowth/bb BA disturbed", "GPPtrees","Rh/SBBpob[layer_1]", "NEP/SMI[layer_1]","W_wsap/fireRisk[layer_1]","W_c","W_s","Wsh","Wdb","dHc",
               "Wbh","Wcrh")

  getVarNam <- function(){
    return(varNames)
}


  aTmean <- function(TAir,nYears){
    Tmean = colMeans(matrix(TAir,365,nYears))
    return(Tmean)
  }

  aTampl <- function(TAir,nYears){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    TbyYear <- matrix(TAir,365,nYears)
    Tampl = apply(TbyYear, 2, function(x) max(aggregate(x/2,by=list(monthsDays),FUN=mean)) - min(aggregate(x/2,by=list(monthsDays),FUN=mean))  )
    return(Tampl)
  }

  aPrecip <- function(Precip,nYears){
    aP = colSums(matrix(Precip,365,nYears))
    return(aP)
  }

  
  
  ### Function for the calculation of monthly mean values from a daily model output variable GPP 
#' monthlyGPP
#'
#' @param GPPx 
#' @param func 
#'
#' @return A vector of the mean daily GPP for each month.
#' @export
#'
#' @examples
  monthlyGPP <- function(GPPx,func="sum"){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),             # Number of days in each month
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    GPPbyYear <- matrix(GPPx,365,length(GPPx)/365)                                           # Divide GPP data into a matrix so that rows are days and columns are years
    GPPm = unlist(apply(GPPbyYear, 2, function(x)                                            # Applies aggregation on the columns of input matrix of GPPbyYear. 
      aggregate(x,by=list(monthsDays),FUN=func)$x))                                          # Aggregation is done by groups in the variable monthsDays, as there is no timestamp on datapoints of GPP
    GPPm <- as.vector(GPPm)
    return(GPPm)
  }
  
#' monthlyWeather
#'
#' @param varx 
#' @param func 
#'
#' @return A vector of the mean daily weather for each month. 
#' @export
#'
#' @examples
  monthlyWeather <- function(varx,func){                                                     # This function works similarly as the monthlyGPP above but with different data set
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    # varByYear <- matrix(varx,365,length(varx)/365)
    if(is.null(dim(varx))){
      varM <- aggregate(varX,by=list(monthsDays),FUN=func)$x
    }else{
      varM = unlist(apply(varx, 1, function(x) 
        aggregate(x,by=list(monthsDays),FUN=func)$x))
    }
    varM <- as.vector(varM)
    return(varM)
  }
  
  ### A post-process function for modifying prebas model output for meeting the Digital Twin Earth project requirements
#' monthlyFluxes
#'
#' @param modOut PREBAS output of a multisite run
#' @param weatherOption 
#'
#' @return
#' @export
#'
#' @examples
  monthlyFluxes <- function(modOut,weatherOption=1){
    if(class(modOut)!="prebas"){
      cueGV <- 0.5 ###carbon use efficiency (NPP/GPP) of ground vegetation
      nYears <- dim(modOut$multiOut)[2]
      nMonths <- nYears * 12
      nSites <- modOut$nSites
      nLayers <- dim(modOut$multiOut)[4]
      
      # if(is.matrix(mGPP)){
      ##aggregating total GPP (tree layers + GV) at monthly time step
      mGPPtot <- t(apply(modOut$dailyPRELES[,,1],1,monthlyGPP,"sum"))                           # Daily GPP is converted into daily GPP using the function monthlyGPP, which is defined above
      # $dailyPRELES is a output matrix from the prebas model
      ### Calculate CUE for different layers
      cueTrees <- modOut$multiOut[,,18,,1]/modOut$multiOut[,,44,,1]                             # CUE = calculate npp / GPPspecies, but include data only from the standing trees
      cueTrees[which(is.na(cueTrees))] <- 0.                                                    # Replace NAs with a value
      cueAll <- pGPP <- array(NA,dim=c(modOut$nSites,modOut$maxYears,(modOut$maxNlayers+1)))    # Create a pair of 3D arrays: x=1:7, y=1:150, z=1:4
      cueAll[,,1] <- cueGV                                                                      # When layer is 1, then CUE is cueGV?
      cueAll[,,2:(modOut$maxNlayers+1)] <- cueTrees                                             # ... and otherwise its cueTrees?
      
      ### Calculate total GPP (totGPP) and ratio of GPP at a layer/total GPP (pGPP)    
      if(nYears==1 & max(nLayers)==1){
        totGPP <- modOut$multiOut[,1,44,1,1] + modOut$GVout[,,3]
      }else{
        totGPP <- apply(modOut$multiOut[,,44,,1],1:2,sum) + modOut$GVout[,,3]                     # Sum of GPP for each year and site in a 2d array
      }
      pGPP[,,1] <- modOut$GVout[,,3]/totGPP                                                     # Ratio between CUE of vegetation and CUE total?
      
      for(i in 1:modOut$maxNlayers){                                                            # Ratio between CUE of vegetation and pine or spruce or mixed
        pGPP[,,(i+1)] <- modOut$multiOut[,,44,i,1] / totGPP
      }
      ### Create arrays for monthly GPP and NPP values
      mGPP <- mNPP <- array(NA,dim=c(modOut$nSites,modOut$maxYears*12,(modOut$maxNlayers+1)))   # Create another 3D array, similar as before.
      for(i in 1:dim(pGPP)[3]){                                                                 # Loop goes through forest layers 1-4 with i
        for(ij in 1:modOut$nSites){                                                             # ... and for every layers it goes through different sites ij
          mGPP[ij,,i] <- rep(pGPP[ij,,i],each=12) * mGPPtot[ij,]                                # Calculate monthly GPP for different layers by multiplying layer's share of GPP with monthly total GPP
          mNPP[ij,,i] <- mGPP[ij,,i] * rep(cueAll[ij,,i],each=12)                               # Calculate monthly NPP for different layers by multiplying monthly total GPP with layers yearly CUE
        }
      }
      
      
      ###litterfall calculations lit is splitted to August and September
      # Woody litter
      # These 5 variables below are different litter categories, have 3d structure: site,litter variable,layer
      if(nYears==1 & max(nLayers)==1){
        Lf <- t(matrix(mLit(modOut$multiOut[,,26,,1],months=8:9),nMonths,nSites))
        Lfr <- t(matrix(mLit(modOut$multiOut[,,27,,1],months=8:9),nMonths,nSites))
        Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
        Lfw <- t(matrix(mLit(modOut$multiOut[,,28,,1],months=8:9),nMonths,nSites))
        Lw <- t(matrix(mLit(modOut$multiOut[,,29,,1],months=8:9),nMonths,nSites))
      }else{
        if(nLayers>1){
          Lf <- aperm(apply(modOut$multiOut[,,26,,1],c(1,3),mLit,months=8:9),c(2,1,3))    # Runs model output variable Litter_fol through the 'mLit' function and transforms the containing array by switching places between rows and columns
          Lfr <- aperm(apply(modOut$multiOut[,,27,,1],c(1,3),mLit,months=8:9),c(2,1,3))   # The mLit function moves litter fall to the August and September time period in the model output
          Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
          Lfw <- aperm(apply(modOut$multiOut[,,28,,1],c(1,3),mLit,months=8:9),c(2,1,3))   # fine woody litter
          Lw <- aperm(apply(modOut$multiOut[,,29,,1],c(1,3),mLit,months=8:9),c(2,1,3))    # Coarse woody litter
        }else{
          Lf <- t(apply(modOut$multiOut[,,26,,1],1,mLit,months=8:9))    # Runs model output variable Litter_fol through the 'mLit' function and transforms the containing array by switching places between rows and columns
          Lfr <- t(apply(modOut$multiOut[,,27,,1],1,mLit,months=8:9))   # The mLit function moves litter fall to the August and September time period in the model output
          Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
          Lfw <- t(apply(modOut$multiOut[,,28,,1],1,mLit,months=8:9))   # fine woody litter
          Lw <- t(apply(modOut$multiOut[,,29,,1],1,mLit,months=8:9))    # Coarse woody litter
        }
      }
      
      ### Create a input array for the Yasso model with the litter fall information
      litter <- array(0.,dim=c(nSites,nMonths,nLayers,3))
      litter[,,,1] <- Lnw   # Non-woody litter
      litter[,,,2] <- Lfw   # Branch litter
      litter[,,,3] <- Lw    # Woody litter
      # litter <- litter*12
      ### Prepare also other initialization information for Yasso
      species <- modOut$multiOut[,1,4,,1]
      nSp <- max(species)
      litterSize <- modOut$litterSize[,1:nSp]
      soilC <- array(0.,dim=c(nSites,(nMonths+1),5,3,nLayers))
      soilC[,1,,,] <- modOut$soilC[,1,,,1:nLayers]  ###initialize with soilC from simulations
      nClimID <- modOut$nClimID
      climIDs <- modOut$siteInfo[,2]
      weatherYasso <- array(0,dim=c(nClimID,nMonths,3))
      if(weatherOption==1){
        weatherYasso[,,1] <- t(apply(modOut$weather[,,,2],1,monthlyWeather,"mean")) ###monthly tempearture
        weatherYasso[,,2] <- t(apply(modOut$weather[,,,4],1,monthlyWeather,"sum")) *12 ###monthly precipitation ###*12 rescales to annual mm
      }else if(weatherOption==2){
        weatherYasso[,,1] <- t(apply(modOut$weather[,,,2],1,monthlyWeather,"mean")) ###monthly tempearture
        weatherYasso[,,2] <- t(apply(modOut$weather[,,,4],1,monthlyWeather,"sum")) *12 ###monthly precipitation ###*12 rescales to annual mm
        weatherYasso[,,3] <- (t(apply(modOut$weather[,,,2],1,monthlyWeather,"max")) - 
                                t(apply(modOut$weather[,,,2],1,monthlyWeather,"min"))) /2
      }else if(weatherOption==3){
        for(i in 1:nClimID){
          for(j in 1:nYears){
            weatherYasso[i,(((j-1)*12)+1):(j*12),1] <- modOut$weatherYasso[i,j,1]
            weatherYasso[i,(((j-1)*12)+1):(j*12),2] <- modOut$weatherYasso[i,j,2]
            weatherYasso[i,(((j-1)*12)+1):(j*12),3] <- modOut$weatherYasso[i,j,3]
          }
        }
      }
      ### Run the Yasso model, which is a function in the src/A_routines.90 file
      soilCtrees <- .Fortran("runYassoMonthly",litter=as.array(litter*12), ###*12 rescale monthly litter to annual units
                             litterSize=as.array(litterSize),
                             nMonths=as.integer(nMonths), 
                             nLayers=as.integer(nLayers), 
                             nSites=as.integer(nSites),
                             nSp=as.integer(nSp),
                             species=as.matrix(species),
                             nClimID=as.integer(nClimID),
                             climIDs=as.integer(climIDs),
                             pAWEN=as.matrix(parsAWEN),
                             pYasso=as.double(pYAS),
                             climate=as.matrix(weatherYasso),
                             soilC=as.array(soilC)) 
      
      ###calculate soil C for gv
      fAPAR <- modOut$fAPAR                                                                                  # fAPAR from the preles model within prebas. 3d array, value for each year, but no layer dimension
      fAPAR[which(is.na(modOut$fAPAR),arr.ind = T)] <- 0.                                                    # Convert NAs to values
      AWENgv <- array(NA,dim=c(dim(modOut$fAPAR),4))                                                         # Add the layer dimension on the fAPAR array
      p0 = as.matrix(modOut$multiOut[,,6,1,1])                                                                          # Annual potential photosynthesis, 3d array: site,year,value
      ETSy = as.matrix(modOut$multiOut[,,5,1,1])
      
      for(ij in 1:nYears){
        AWENgv[,ij,] <- t(sapply(1:nrow(fAPAR), function(i) 
          .Fortran("fAPARgv",fAPAR[i,ij],
                   ETSy[i,ij],modOut$siteInfo[i,2],
                   0,0,p0[i,ij],rep(0,4),0)[[7]]))
      }
      
      mAWEN <- array(0,dim=c(nSites,nMonths,5))
      mAWEN[,,1:4] <- aperm(apply(AWENgv,c(1,3),mLit,months=6:9),c(2,1,3))*12 ###*12 rescale monthly litter to annual units
      if(nYears==1 & max(nLayers)==1){
        mGVit <- t(matrix(mLit(modOut$GVout[,,2],months=6:9),nMonths,nSites))
      }else{
        mGVit <- t(apply(modOut$GVout[,,2],1,mLit,months=6:9))
      }
      
      ###calculate steady state soil C per GV
      # ststGV <- matrix(NA,nSites,5)
      soilGV <- array(0,dim=c(nSites,(nMonths+1),5))
      
      soilCgv <- .Fortran("runYassoAWENinMonthly",
                          mAWEN=as.array(mAWEN),
                          nMonths=as.integer(nMonths),  
                          nSites=as.integer(nSites), 
                          litSize=0.,
                          nClimID=as.integer(nClimID),
                          climIDs=as.integer(climIDs),
                          pYasso=as.double(pYAS),
                          climate=as.matrix(weatherYasso),
                          soilCgv=as.array(soilGV))
      ####add gvsoilc to first layer foliage soilC
      # check in normal runs where ground vegetation soilC is calculated
      soilCtot <- apply(soilCtrees$soilC,1:2,sum) + apply(soilCgv$soilCgv,1:2,sum)
      mRhTot <- soilCtot[,1:nMonths]/10 - soilCtot[,2:(nMonths+1)]/10 +
        apply(litter,1:2,sum)/10 + mGVit/10
      mNPPtot <- apply(mNPP,1:2,sum)
      mRaTot <- mGPPtot - mNPPtot
      mNEPtot <- mNPPtot - mRhTot
      return(list(mGPP=mGPPtot,mNPP=mNPPtot,mRa=mRaTot,mRh=mRhTot,
                  mNEP=mNEPtot,soilC=soilCtot))
    
    }else{
      
      cueGV <- 0.5 ###carbon use efficiency (NPP/GPP) of ground vegetation
      nYears <- dim(modOut$output)[1]
      nMonths <- nYears * 12
      nSites <- 1
      nLayers <- dim(modOut$output)[3]
      
      # if(is.matrix(mGPP)){
      ##aggregating total GPP (tree layers + GV) at monthly time step
      mGPPtot <- monthlyGPP(modOut$dailyPRELES[,1])                           # Daily GPP is converted into daily GPP using the function monthlyGPP, which is defined above
      # $dailyPRELES is a output matrix from the prebas model
      ### Calculate CUE for different layers
      cueTrees <- modOut$output[,18,,1]/modOut$output[,44,,1]                             # CUE = calculate npp / GPPspecies, but include data only from the standing trees
      cueTrees[which(is.na(cueTrees))] <- 0.                                                    # Replace NAs with a value
      cueAll <- pGPP <- array(NA,dim=c(nYears,(nLayers+1)))    # Create a pair of 3D arrays: x=1:7, y=1:150, z=1:4
      cueAll[,1] <- cueGV                                                                      # When layer is 1, then CUE is cueGV?
      cueAll[,2:(nLayers+1)] <- cueTrees                                             # ... and otherwise its cueTrees?
      
      ### Calculate total GPP (totGPP) and ratio of GPP at a layer/total GPP (pGPP)    
      if(nYears==1 & nLayers==1){
        totGPP <- modOut$output[1,44,1,1] + modOut$GVout[,3]
      }else{
        totGPP <- rowSums(modOut$output[,44,,1]) + modOut$GVout[,3]                     # Sum of GPP for each year and site
      }
      pGPP[,1] <- modOut$GVout[,3]/totGPP                                                     # Ratio between GPP of ground vegetation and GPP total
      
      for(i in 1:nLayers){                                                            # Ratio between CUE of vegetation and pine or spruce or mixed
        pGPP[,(i+1)] <- modOut$output[,44,i,1] / totGPP
      }
      ### Create arrays for monthly GPP and NPP values
      mGPP <- mNPP <- array(NA,dim=c(nYears*12,(nLayers+1)))   # Create another matrix, similar as before.
      for(i in 1:dim(pGPP)[2]){                                                                 # Loop goes through forest layers 1-4 with i
        mGPP[,i] <- rep(pGPP[,i],each=12) * mGPPtot                                # Calculate monthly GPP for different layers by multiplying layer's share of GPP with monthly total GPP
        mNPP[,i] <- mGPP[,i] * rep(cueAll[,i],each=12)                               # Calculate monthly NPP for different layers by multiplying monthly total GPP with layers yearly CUE
      }
      
      
      ###litterfall calculations lit is splitted to August and September
      # Woody litter
      # These 5 variables below are different litter categories, have 3d structure: site,litter variable,layer
      if(nYears==1 & nLayers==1){
        Lf <- mLit(modOut$output[,26,1,1],months=8:9)
        Lfr <- mLit(modOut$output[,27,1,1],months=8:9)
        Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
        Lfw <- mLit(modOut$output[,28,1,1],months=8:9)
        Lw <- mLit(modOut$output[,29,1,1],months=8:9)
      }else{
        if(nLayers>1){
          Lf <- apply(modOut$output[,26,,1],2,mLit,months=8:9)    # Runs model output variable Litter_fol through the 'mLit' function and transforms the containing array by switching places between rows and columns
          Lfr <- apply(modOut$output[,27,,1],2,mLit,months=8:9)   # The mLit function moves litter fall to the August and September time period in the model output
          Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
          Lfw <- apply(modOut$output[,28,,1],2,mLit,months=8:9)   # fine woody litter
          Lw <- apply(modOut$output[,29,,1],2,mLit,months=8:9)    # Coarse woody litter
        }else{
          Lf <- mLit(modOut$output[,26,1,1],months=8:9)
          Lfr <- mLit(modOut$output[,27,1,1],months=8:9)
          Lnw <- Lf + Lfr                                                                 # Non-woody litter = The sum of Foliage litter (var. Litter_fol) and Fine root litter (var. Litter_fr) prebas output variables
          Lfw <- mLit(modOut$output[,28,1,1],months=8:9)
          Lw <- mLit(modOut$output[,29,1,1],months=8:9)
        }
      }
      
      ### Create a input array for the Yasso model with the litter fall information
      litter <- array(0.,dim=c(nSites,nMonths,nLayers,3))
      litter[1,,,1] <- Lnw   # Non-woody litter
      litter[1,,,2] <- Lfw   # Branch litter
      litter[1,,,3] <- Lw    # Woody litter
      litter[which(is.na(litter))] <- 0
      # litter <- litter*12
      ### Prepare also other initialization information for Yasso
      species <- modOut$output[1,4,,1]
      nSp <- max(species)
      litterSize <- modOut$litterSize[,1:nSp]
      soilC <- array(0.,dim=c(nSites,(nMonths+1),5,3,nLayers))
      soilC[1,1,,,] <- modOut$soilC[1,,,1:nLayers]  ###initialize with soilC from simulations
      nClimID <- 1#modOut$nClimID
      climIDs <- 1#modOut$siteInfo[,2]
      weatherYasso <- array(0,dim=c(nClimID,nMonths,3))
      if(weatherOption==1){
        weatherYasso[1,,1] <- monthlyWeather(modOut$weather[,,2],"mean") ###monthly tempearture
        weatherYasso[1,,2] <- monthlyWeather(modOut$weather[,,4],"sum") *12 ###monthly precipitation ###*12 rescales to annual mm
      }else if(weatherOption==2){
        weatherYasso[,,1] <- monthlyWeather(modOut$weather[,,2],"mean") ###monthly tempearture
        weatherYasso[,,2] <- monthlyWeather(modOut$weather[,,4],"sum") *12 ###monthly precipitation ###*12 rescales to annual mm
        weatherYasso[,,3] <- (monthlyWeather(modOut$weather[,,2],"max") - 
                                monthlyWeather(modOut$weather[,,2],"min")) /2
      }else if(weatherOption==3){
        for(j in 1:nYears){
          weatherYasso[1,(((j-1)*12)+1):(j*12),1] <- modOut$weatherYasso[j,1]
          weatherYasso[1,(((j-1)*12)+1):(j*12),2] <- modOut$weatherYasso[j,2]
          weatherYasso[1,(((j-1)*12)+1):(j*12),3] <- modOut$weatherYasso[j,3]
        }
      }
      ### Run the Yasso model, which is a function in the src/A_routines.90 file
      soilCtrees <- .Fortran("runYassoMonthly",litter=as.array(litter*12), ###*12 rescale monthly litter to annual units
                             litterSize=as.array(litterSize),
                             nMonths=as.integer(nMonths), 
                             nLayers=as.integer(nLayers), 
                             nSites=as.integer(nSites),
                             nSp=as.integer(nSp),
                             species=as.matrix(species),
                             nClimID=as.integer(nClimID),
                             climIDs=as.integer(climIDs),
                             pAWEN=as.matrix(parsAWEN),
                             pYasso=as.double(pYAS),
                             climate=as.matrix(weatherYasso),
                             soilC=as.array(soilC)) 
      
      ###calculate soil C for gv
      fAPAR <- modOut$fAPAR                                                                                  # fAPAR from the preles model within prebas. 3d array, value for each year, but no layer dimension
      fAPAR[which(is.na(modOut$fAPAR),arr.ind = T)] <- 0.                                                    # Convert NAs to values
      AWENgv <- array(NA,dim=c(nYears,4))                                                         # Add the layer dimension on the fAPAR array
      p0 = as.matrix(modOut$output[,6,1,1])                                                                          # Annual potential photosynthesis, 3d array: site,year,value
      ETSy = as.matrix(modOut$output[,5,1,1])
      
      for(ij in 1:nYears){
        AWENgv[ij,] <- .Fortran("fAPARgv",fAPAR[ij],
                                ETSy[ij],climIDs,
                                0,0,p0[ij],rep(0,4),0)[[7]]
      }
      
      mAWEN <- array(0,dim=c(nSites,nMonths,5))
      mAWEN[1,,1:4] <- apply(AWENgv,2,mLit,months=6:9)*12 ###*12 rescale monthly litter to annual units
      if(nYears==1 & nLayers==1){
        mGVit <- mLit(modOut$GVout[1,2],months=6:9)
      }else{
        mGVit <- mLit(modOut$GVout[,2],months=6:9)
      }
      
      ###calculate steady state soil C per GV
      # ststGV <- matrix(NA,nSites,5)
      soilGV <- array(0,dim=c(nSites,(nMonths+1),5))
      
      soilCgv <- .Fortran("runYassoAWENinMonthly",
                          mAWEN=as.array(mAWEN),
                          nMonths=as.integer(nMonths),  
                          nSites=as.integer(nSites), 
                          litSize=0.,
                          nClimID=as.integer(nClimID),
                          climIDs=as.integer(climIDs),
                          pYasso=as.double(pYAS),
                          climate=as.matrix(weatherYasso),
                          soilCgv=as.array(soilGV))
      ####add gvsoilc to first layer foliage soilC
      # check in normal runs where ground vegetation soilC is calculated
      soilCtot <- apply(soilCtrees$soilC,1:2,sum) + apply(soilCgv$soilCgv,1:2,sum)
      mRhTot <- soilCtot[,1:nMonths]/10 - soilCtot[,2:(nMonths+1)]/10 +
        apply(litter,1:2,sum)/10 + mGVit/10
      mNPPtot <- apply(mNPP,1,sum)
      mRaTot <- mGPPtot - mNPPtot
      mNEPtot <- mNPPtot - mRhTot
      return(list(mGPP=mGPPtot,mNPP=mNPPtot,mRa=mRaTot,mRh=mRhTot[1,],
                  mNEP=mNEPtot[1,],soilC=soilCtot[1,]))      
    }
  } 
  
  
   
  ### This function is used for adding yearling litter fall to the fall period on a timeseries data set, used in the monthlyFluxes function
#' mLit
#'
#' @param aLit Number of years
#' @param months August and September 
#'
#' @return A time series with the annual autumn litter fall.  
#' @export
#'
#' @examples
  mLit <- function(aLit,months=8:9){
    nYears <- length(aLit)                                                        # Gets the length of the array
    nMonths <- length(months)
    lit <- matrix(0,12,nYears)                                                    # Creates a matrix by 12 rows, nYears columns and fills it with zeros
    lit[months,] <- matrix(aLit/nMonths,nMonths,nYears,byrow = T)                                # This appears to create a new matrix within a matrix, on rows 8-9. Is this just used to fill these rows with data?
    lit <- as.vector(lit)                                                         # Matrix is converted intoa vector, which removes the table structure and appears to create a timeseries data set
    return(lit)                                                                   
  }
  
  
  ###Function to calculate basal weighted mean of a PREBAS output
#' baWmean
#'
#' @param modOut PREBAS output of a multisite run
#' @param varX index of the variable for which the basal area weighted mean needs to be calculated
#'
#' @return The weighted mean basal area for a PREBAS output (including both multi site and single site runs)
#' @export
#'
#' @examples
  baWmean <- function(modOut,varX){
    ###calculates basal area weighted mean for single site runs
    if(class(modOut)=="prebas"){
      weightXs <- as.matrix(apply(modOut$output[,13,,1],1,FUN=function(vec)vec/sum(vec,na.rm=T)))
      weightXs <- t(weightXs)
      weightXs[which(is.na(weightXs))] <- 0.
      weigthedMean <- rowSums(modOut$output[,varX,,1]*weightXs)
    }else{
      ###calculates basal area weighted mean for multi site runs
      weightXs <- apply(modOut$multiOut[,,13,,1],1:2,FUN=function(vec)vec/sum(vec,na.rm=T))
      weightXs <- aperm(weightXs,c(2:3,1))
      weightXs[which(is.na(weightXs))] <- 0.
      weigthedMean <- apply(modOut$multiOut[,,varX,,1]*weightXs,1:2,sum,na.rm=T)
    }
    return(weigthedMean)
  }
  
####function to calculate quadratic mean diameter
#' DqFun
#'
#' @param D average diameter
#' @param N stand density
#'
#' @return The quadratic mean diameter for a stand. 
#' @export
#'
#' @examples
  DqFun <- function(D,N){ 
    sqrt(sum(D^2*N)/sum(N))
  }
  
####Calculate quadratic mean diameter using a multisite PREBAS run
#' DqFunMS
#'
#' @param modOut PREBAS output of a multisite run
#'
#' @return The quadratic mean diameter of a multisite PREBAS run 
#' @export
#'
#' @examples
  DqFunMS <- function(modOut){
    D <-modOut[,,12,,1]
    N <-modOut[,,17,,1]
    Dq <- sqrt(apply(D^2*N,1:2,sum)/apply(N,1:2,sum))
    return(Dq)
  }
  
####Calculate Reineke Stand density index
#' reineke1
#'
#' @param Dq quadratic mean diameter of a multisite PREBAS run
#' @param k 
#' @param spID 
#'
#' @return Reineke Stand density index
#' @export
#'
#' @examples
  reineke1 <- function(Dq,k,spID){
    Nmax <- k/Dq^1.605
    return(Nmax)
  }
  
  ####Calculate Reineke Stand density index (PREBAS version) using PREBAS parameters
  reineke <- function(D,N,pCrobas,spID){
    Ntot <- sum(N,na.rm=T)
    par_kRein <- pCrobas[17,spID]
    reinekeLayer = Ntot*(D/25.)**(1.66)
    reinX = reinekeLayer / par_kRein
    return(reinX)
  }
  
  ####MUltisite version for Reineke Stand density index (PREBAS version) using PREBAS parameters
  reinekeMS <- function(modOut,pCrobas){
    D <- modOut[,,12,,1]
    N <- Ntot <- par_kRein <- modOut[,,17,,1]
    spID <- modOut[,,4,,1]
    nLayers <- dim(N)[3]
    nYears <- dim(N)[2]
    Ntot[,,1] <- apply(N,1:2,sum,na.rm=T)
    if(dim(N)[3]>1){
      for(i in 2:nLayers) Ntot[,,i] <- Ntot[,,1]
    }
    zeros <- which(spID==0.,arr.ind=T)
    spID[zeros] <- 7
    for(j in 1:nYears){
      for(i in 1:nLayers){
        par_kRein[,j,i] <-  pCrobas[17,spID[,j,i]]
      }
    }
    par_kRein <- array(pCrobas[17,spID],dim=dim(N))
    par_kRein[zeros] <- 0.
    spID[zeros] <- 0
    reinekeLayer = Ntot*(D/25.)**(1.66)
    reinX = reinekeLayer / par_kRein
    return(reinX)
  }
  
  # 
  # D <- data.all$dbh
  # Ntot <- data.all$ba/(pi * (data.all$dbh/200)^2)
  # reinekeLayer = Ntot*(D/25.)**(1.66)
  # reinX1 = reinekeLayer / pCrobas[17,1]
  # reinX2 = reinekeLayer / pCrobas[17,2]
  # reinX3 = reinekeLayer / pCrobas[17,3]
  ####MUltisite version for Reineke Stand density index (PREBAS version) using PREBAS parameters
  reinekeMSinit <- function(multiInitVar,pCrobas=pCROB){
    D <- multiInitVar[,4,]
    N <- Ntot <- par_kRein <- multiInitVar[,5,]/(pi*(multiInitVar[,4,]/200)^2)
    spID <- multiInitVar[,1,]
    nLayers <- dim(N)[2]
    if(!is.null(nLayers)){
      Ntot[,1] <- apply(N,1,sum,na.rm=T)
      for(i in 2:nLayers) Ntot[,i] <- Ntot[,1]
    }
    zeros <- which(spID==0.,arr.ind=T)
    spID[zeros] <- 7
    if(!is.null(nLayers)){
      for(i in 1:nLayers){
        par_kRein[,i] <-  pCrobas[17,spID[,i]]
      }
    }else{
      par_kRein <-  pCrobas[17,spID]
    }
    par_kRein[zeros] <- 0.
    spID[zeros] <- 0
    reinekeLayer = Ntot*(D/25.)**(1.66)
    reinX = reinekeLayer / par_kRein
    return(reinX)
  }
  
  
  
####steady state solution for Dead wood volume calculations
  ### arguments: 
  ## modRun -> multiRun from PREBAS
  ## siteSel -> selection of sites to include in the steady state calculation
  ## yearX -> selection of years to include in the steady state calculation
  ## yMix -> if T randomize the years
  
  
  initDeadW <- function(modRun,siteSel=NA,yearX=NA, yMix=F){
    nSim <- dim(modRun$multiOut)[2]
    if(all(is.na(siteSel))) siteSel <- 1:dim(modRun$multiOut)[1]
    if(all(is.na(yearX))) yearX <- 1:nSim
    if(yMix) yearX <- sample(yearX)
    modOut <- modRun$multiOut[siteSel,yearX,,,]
    pCrobas <- modRun$pCROBAS
    nSites <- dim(modOut)[1]
    nLayers <- dim(modOut)[4]
    nYears <- dim(modOut)[2]
    # year=2
    # ij=1 #species
    nX <- ceiling((max(150,nYears*3))/nYears)
    modOutX <- selX <- modOut[,,c(12,42,4),,1]
    for(i in 1:(nX-1)){
      modOutX <- abind(modOutX,selX,along=2)
    }
    sims <- modRun$multiOut[siteSel,,c(12,42,4),,1]
    sims[,,2,] <- 0.
    modOutX <- abind(modOutX,sims,along=2)
    nYearsX <- dim(modOutX)[2]
    # modOutX[,(nYears*+1):(nYears*3),2,] <- 0
    deadWV <- array(0,dim=c(nSites,(nYearsX+2),nLayers))
    #if(dN<0.) then
    for(siteX in 1:nSites){
      for(year in 1:nYearsX){
        for(ij in 1:nLayers){
          if(year == 1){
            D <- modOutX[siteX,year,1,ij]
          }else{
            D <- modOutX[siteX,(year-1),1,ij]
          }
          # V <- modOutX[siteX,year,30,ij]
          Vmort <- modOutX[siteX,year,2,ij]
          species <- modOutX[siteX,year,3,ij]
          if(Vmort>0 & D>pCrobas[48,species]){
            # D <- modOutX[siteX,year,12,ij,1]
            perVmort <- exp(-exp(pCrobas[35,species] + pCrobas[36,species]*(1:(nYearsX-year)) + 
                                   pCrobas[37,species]*D + pCrobas[44,species]))
            perVmort[which(perVmort<pCrobas[49,species])] <- 0
            deadWV[siteX,year,ij] = Vmort + deadWV[siteX,year,ij]
            deadWV[siteX,(year+1):nYearsX,ij] = deadWV[siteX,(year+1):nYearsX,ij] +
              Vmort *perVmort
          }
        }
      }
    }
    # deadWVx <- apply(deadWV,1:2,sum)
    deadWVss <- apply(deadWV,2:3,mean)[(nYearsX-nSim+1):(nYearsX),]
    # modRun$multiOut[selX,,8,,1] <- modRun$multiOut[selX,,8,,1] + deadWVss
    return(list(ssDeadW=deadWVss,deadWV=deadWV))
  }
  
  ####Function used in the oldLayer option 
  ####the function is used to add an additional empty layer
  ##for the oldLayer
  ###add layer for old layer in initial state
  addOldLayer <- function(multiSiteInit){
    maxNlayers = multiSiteInit$maxNlayers = 
      multiSiteInit$maxNlayers + 1
    ###add layer for old layer in multiOut array
    dims=dim(multiSiteInit$multiOut);dims[4] <- maxNlayers
    multiOut = array(0,dim=dims)
    multiOut[,,,1:(maxNlayers-1),] = multiSiteInit$multiOut
    multiSiteInit$multiOut <- multiOut
    ###add layer for old layer in initial state
    dims=dim(multiSiteInit$multiInitVar);dims[3] <- maxNlayers
    multiInitVar = array(0,dim=dims)
    multiInitVar[,,1:(maxNlayers-1)] = multiSiteInit$multiInitVar
    multiSiteInit$multiInitVar <- multiInitVar
    ###add layer
    multiSiteInit$nLayers = multiSiteInit$nLayers + 1
    ###add layer to nLayers in siteInfo
    multiSiteInit$siteInfo[,8] = multiSiteInit$siteInfo[,8] + 1
    ###add layer to initClCutRatio
    initCLcutRatio <- matrix(0.,multiSiteInit$nSites,
                             maxNlayers)
    initCLcutRatio[,1:(maxNlayers-1)] = multiSiteInit$initCLcutRatio
    multiSiteInit$initCLcutRatio = initCLcutRatio
    
    dims=dim(multiSiteInit$soilC);dims[5] <- maxNlayers
    soilC = array(0,dim=dims)
    soilC[,,,,1:(maxNlayers-1)] = multiSiteInit$soilC
    multiSiteInit$soilC <- soilC
    
    dims=dim(multiSiteInit$multiEnergyWood);dims[3] <- maxNlayers
    multiEnergyWood = array(0,dim=dims)
    multiEnergyWood[,,1:(maxNlayers-1),] <- multiSiteInit$multiEnergyWood
    multiSiteInit$multiEnergyWood <- multiEnergyWood
    
    return(multiSiteInit)
  }
  
  fTfun <- function(TAir, precip){
    fT <- exp(0.059*TAir - 0.001*TAir^2)* (1-exp(-1.858*precip))
    return(fT)
  }
  
  
  alpharN <- function(alphar0, p0,p00,TAir0,TAir,precip0,precip){
    fT0 <- fTfun(TAir0,precip0)
    fT <- fTfun(TAir,precip)
    alphar <- alphar0 * p0/p00*fT0/fT
  }
  
  
  ###wrapper funtion to call in R the Fortran funtion used to run the ground vegetation model and reporting the outputs for each vegetation type  
  GVbyVtypes <- function(inputs){
    
    output <-rep(0,13)
    test <- .Fortran("fAPARgvByVtypes",
                     inputs = as.double(inputs),
                     output = as.double(output))
    output <- test$output
    return(output)
  }

  ###function to run the GVmodel and report the output by vegetation type: grass and herbs, shrubs and moss&lichens  
  ###modOut is a PREBAS run (muöltisite or region) it needs to be adapted to single site runs
  ####GPP and NPP are calculated based on fAPAR and biomass ratios respectively
  GVmodByVegType <- function(modOut){
    if(class(modOut)!="prebas"){
      
      nYears <- modOut$maxYears
      fAPAR <- modOut$fAPAR
      ets <- modOut$multiOut[,,5,1,1]
      siteTypeX <- modOut$multiOut[,,3,1,1]
      GVout <- array(NA,dim=c(19,modOut$nSites,nYears))
      
      inputs <- abind(fAPAR,ets,siteTypeX,along=3)
      
      GVout[1:13,,] <- apply(inputs,1:2,GVbyVtypes)
      
      oo <- apply(GVout[1:3,,],2:3,sum)
      fAPARratio <- sweep(GVout[1:3,,],2:3,apply(GVout[1:3,,],2:3,sum),FUN="/")
      GVgpp <-  sweep(fAPARratio,2:3,modOut$GVout[,,3],FUN="*") 
      
      Wtot <- GVout[9:11,,]  ##above ground
      Wtot[1:2,,] <- Wtot[1:2,,] + GVout[12:13,,]     ##add belowground 
      Wratio <- sweep(Wtot,2:3,apply(Wtot,2:3,sum),FUN="/")
      GVnpp <-  sweep(Wratio,2:3,modOut$GVout[,,5],FUN="*") 
      
      GVout[14:16,,] <- GVgpp
      GVout[17:19,,] <- GVnpp
      
      namesX <- c(paste0("fAPAR_",c("gra&her", "srhb", "mos&lic")),
                  paste0("liag_",c("gra&her", "srhb", "mos&lic")),
                  paste0("litbg_",c("gra&her", "srhb")),
                  paste0("Wag_",c("gra&her", "srhb", "mos&lic")),
                  paste0("Wbg_",c("gra&her", "srhb")),
                  paste0("GPP_",c("gra&her", "srhb", "mos&lic")),
                  paste0("NPP_",c("gra&her", "srhb", "mos&lic"))
      )
      dimnames(GVout) <- list(GVmod= namesX,sites=NULL,years=NULL)
      return(GVout)
    }else{
      nYears <- modOut$nYears
      fAPAR <- modOut$fAPAR
      ets <- modOut$ETS
      siteTypeX <- modOut$output[,3,1,1]
      GVout <- array(NA,dim=c(19,nYears))
      
      inputs <- abind(fAPAR,ets,siteTypeX,along=2)
      
      GVout[1:13,] <- apply(inputs,1,GVbyVtypes)
      
      oo <- apply(GVout[1:3,],2,sum)
      fAPARratio <- sweep(GVout[1:3,],2,apply(GVout[1:3,],2,sum),FUN="/")
      GVgpp <-  sweep(fAPARratio,2,modOut$GVout[,3],FUN="*") 
      
      Wtot <- GVout[9:11,]  ##above ground
      Wtot[1:2,] <- Wtot[1:2,] + GVout[12:13,]     ##add belowground 
      Wratio <- sweep(Wtot,2,apply(Wtot,2,sum),FUN="/")
      GVnpp <-  sweep(Wratio,2,modOut$GVout[,5],FUN="*") 
      
      GVout[14:16,] <- GVgpp
      GVout[17:19,] <- GVnpp
      
      namesX <- c(paste0("fAPAR_",c("gra&her", "srhb", "mos&lic")),
                  paste0("liag_",c("gra&her", "srhb", "mos&lic")),
                  paste0("litbg_",c("gra&her", "srhb")),
                  paste0("Wag_",c("gra&her", "srhb", "mos&lic")),
                  paste0("Wbg_",c("gra&her", "srhb")),
                  paste0("GPP_",c("gra&her", "srhb", "mos&lic")),
                  paste0("NPP_",c("gra&her", "srhb", "mos&lic"))
      )
      dimnames(GVout) <- list(GVmod= namesX,years=NULL)
      return(GVout)
    }
  }
  
  #' Nesterov Index
  #' A cumulative function of daily Tmax and dew-point temperature Tdew, eq. 5 in TH2010
  #' @param rain daily precipitation (mm)
  #' @param tmin daily minimum temperature (ºC)
  #' @param tmaX daily maximum temperature (ºC)
  #'
  #' @return the Nesterov Index 
  #' @export
  #'
  #' @examples
  NesterovInd <- function(rain, tmin, tmax){
    nDays <- length(rain)
    NI <- rep(0,nDays)
    daysX <- (which(rain[2:nDays]< 3 & (tmin[2:nDays]-4)>=0))+1 #do not consider the first day
    if(length(daysX)>0){
      NI[daysX] = (tmax[daysX]*(tmax[daysX]-tmin[daysX]-4.))+
        (tmax[daysX-1]*(tmax[daysX-1]-tmin[daysX-1]-4))  
    }
    return(NI)
  }
  