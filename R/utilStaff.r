####init Biomass
initBiomasses <- function(pCro,initVarX){
  initVarX<-as.matrix(initVarX) #change vector to matrix when maxlayer=1
  siteType <- initVarX[8,1]
  ##set parameters
  par_betab <- pCro[13,initVarX[1,]]
  par_x <- pCro[19,initVarX[1,]]
  par_beta0 <- pCro[12,initVarX[1,]]
  par_betas <- pCro[14,initVarX[1,]]
  par_mf <- pCro[8,initVarX[1,]]
  par_mr <- pCro[9,initVarX[1,]]
  par_mw <- pCro[10,initVarX[1,]]
  par_alfar <- pCro[20+pmin(siteType,5),initVarX[1,]]
  par_c <- pCro[7,initVarX[1,]]
  par_rhof <- pCro[15,initVarX[1,]]
  par_rhor <- par_alfar * par_rhof
  par_rhow <- pCro[2,initVarX[1,]]
  par_S_branchMod <- pCro[27,initVarX[1,]]
  gammaC <- 0. #initVarX[8,]
  Tbd <- 10 #####to include in the parameters
  
  ###set variables
  A <- initVarX[7,]
  ba <- initVarX[5,]; d <- initVarX[4,]
  N <- ba/(pi*((d/2/100)^2))
  h = initVarX[3,]; hc <- initVarX[6,]
  B = ba/N
  Lc <- h - hc
  betab =  par_betab * Lc^(par_x-1)
  beta0 = par_beta0
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
  biomasses <- rbind(wf_STKG,W_froot,W_wsap,W_c,W_s,W_branch,W_croot,Wsh,Wdb,W_stem,V,W_crh)

  return(biomasses)
}

#### function to calculate initial sapwood area at crown base (A)
compA <- function(inputs){
  p_ksi = inputs[1]
  p_rhof = inputs[2]
  p_z <- inputs[3]
  Lc = inputs[4]
  A <- max(0.,p_ksi/p_rhof * Lc^p_z)
  return(A)
}

varNames  <- c('siteID','gammaC','sitetype','species','ETS' ,'P0','age', 'DeadWoodVolume', 'Respi_tot','GPP/1000',
               'H','D', 'BA','Hc_base','Cw','Ac','N','npp','leff','keff','lproj','ET_preles','weight',
               'Wbranch',"WfineRoots",'Litter_fol','Litter_fr','Litter_fWoody','Litter_cWoody','V',
               'Wstem','W_croot','wf_STKG', 'wf_treeKG','B_tree','Light',"Vharvested","Wharvested","soilC",
               "aSW","dH","Vmort","gross growth", "GPPspecies","Rh species", "NEP sp"," W_wsap","W_c","W_s","Wsh","Wdb","dHc",
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
  monthlyGPP <- function(GPPx,func="sum"){
    monthsDays <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),             # Number of days in each month
                    rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
    GPPbyYear <- matrix(GPPx,365,length(GPPx)/365)                                           # Divide GPP data into a matrix so that rows are days and columns are years
    GPPm = unlist(apply(GPPbyYear, 2, function(x)                                            # Applies aggregation on the columns of input matrix of GPPbyYear. 
      aggregate(x,by=list(monthsDays),FUN=func)$x))                                          # Aggregation is done by groups in the variable monthsDays, as there is no timestamp on datapoints of GPP
    GPPm <- as.vector(GPPm)
    return(GPPm)
  }
  
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
  monthlyFluxes <- function(modOut,weatherOption=1){
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
                    0,0,p0[i,ij],rep(0,4))[[7]]))
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
  } 
  
  
   
  ### This function is used for adding yearling litter fall to the fall period on a timeseries data set, used in the monthlyFluxes function
  mLit <- function(aLit,months=8:9){
    nYears <- length(aLit)                                                        # Gets the length of the array
    nMonths <- length(months)
    lit <- matrix(0,12,nYears)                                                    # Creates a matrix by 12 rows, nYears columns and fills it with zeros
    lit[months,] <- matrix(aLit/nMonths,nMonths,nYears,byrow = T)                                # This appears to create a new matrix within a matrix, on rows 8-9. Is this just used to fill these rows with data?
    lit <- as.vector(lit)                                                         # Matrix is converted intoa vector, which removes the table structure and appears to create a timeseries data set
    return(lit)                                                                   
  }
  
  
  ###Function to calculate basal weighted mean of a PREBAS output
  ### modOut = multisite PREBAS output
  ### varX = index of the variable for which the basal area weighted mean needs to be calculated
  baWmean <- function(modOut,varX){
    ###calculates basal area weighted mean for single site runs
    if(class(modOut)=="prebas"){
      weightXs <- as.matrix(apply(modOut$output[,13,,1],1,FUN=function(vec)vec/sum(vec)))
      weightXs <- t(weightXs)
      weightXs[which(is.na(weightXs))] <- 0.
      weigthedMean <- rowSums(modOut$output[,varX,,1]*weightXs)
    }else{
      ###calculates basal area weighted mean for multi site runs
      weightXs <- apply(modOut$multiOut[,,13,,1],1:2,FUN=function(vec)vec/sum(vec))
      weightXs <- aperm(weightXs,c(2:3,1))
      weightXs[which(is.na(weightXs))] <- 0.
      weigthedMean <- apply(modOut$multiOut[,,varX,,1]*weightXs,1:2,sum)
    }
    return(weigthedMean)
  }
  
####function to calculate quadratic mean diameter
  ####D -> average diameter 
  ####N stand density
  DqFun <- function(D,N) sqrt(sum(D^2*N)/sum(N))
  
####Calculate quadratic mean diameter using a multisite PREBAS run
  ###modOut PREBAS output of a multisite run
  DqFunMS <- function(modOut){
    D <-modOut$multiOut[,,12,,1]
    N <-modOut$multiOut[,,17,,1]
    Dq <- sqrt(sum(D^2*N)/sum(N))
    return(Dq)
  }
  
####Calculate Reineke Stand density index
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
  reinekeMSinit <- function(initPrebas){
    pCrobas <- initPrebas$pCROBAS
    D <- initPrebas$multiInitVar[,4,]
    N <- Ntot <- par_kRein <- initPrebas$multiInitVar[,5,]/(pi*(initPrebas$multiInitVar[,4,]/200)^2)
    spID <- initPrebas$multiInitVar[,1,]
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
  
  