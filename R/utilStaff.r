
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
  W_crs = par_rhow * beta0 * A * H * N
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
               'Wbranch',"WfineRoots",'Litter_fol','Litter_fr','Litter_branch','Litter_wood','V',
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
  
  
