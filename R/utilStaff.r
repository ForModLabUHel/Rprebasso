###default parameters
pCROB <- matrix(NA,40,3,dimnames = list(NULL,c("pisy","piab","beal")))
pCROB[,1] <- c(0.22406007,195.16571750,20.00788979,0.25003240,3.90739583,0.90777990,
 0.29990356,0.40000000,0.50000000,0.03002172,1.79017546,0.28690180,
 0.39976573,0.38006730,180.23481660,0.01052194,855.00952220,0.97559558,
 0.80000000,0.68700000,0.40000000,0.44000000,0.47000000,0.64000000,
 0.84000000,1.00000000,1.00000000,1.00000000,1.40000000,1250.00000000,
 0.00000000, 33.00000000,0.03000000,0.00066900, -2.65300000,0.05500000,
 -0.03000000,0.07000000,20.00788979,0.00000000)

pCROB[,2] <- c(1.585611e-01,1.833353e+02,2.002536e+01,2.503559e-01,9.728822e+00,1.764810e+00,
 2.694873e-01,4.000000e-01,5.000000e-01,3.005746e-02,1.700265e+00,2.567875e-01,
 4.996872e-01,4.660972e-01,2.000728e+02,5.568057e-03,1.037631e+03,4.001300e-01,
 6.000000e-01,8.740000e-01,1.000000e-01,2.800000e-01,3.800000e-01,4.800000e-01,
 5.800000e-01,2.000000e+00,2.000000e+00,1.000000e+00,1.400000e+00,1.250000e+03,
 0.000000e+00,3.700000e+01,3.000000e-02,3.270000e-04,-2.948000e+00,5.900000e-02,
 -3.000000e-02,8.000000e-02,2.002536e+01,0.000000e+00)

pCROB[,3] <- c(0.17815383,216.13110400, 40.56651551,0.31592870,0.96149996,1.43787917,
  0.21899819,0.40000000,0.50000000,0.03011309,1.94740316,0.49585742,
  0.38946665,0.48527250,101.14578740,0.02440534,1064.46260000,0.40215544,
  0.80000000,0.00000000,0.35000000,0.50000000,0.64000000,0.75000000,
  0.94000000,1.00000000,1.00000000,2.00000000,1.40000000,1250.00000000,
  0.00000000,37.00000000,0.03000000,0.00064400,-3.32400000,0.13500000,
  -0.03000000,0.01500000,40.56651551,0.00000000)
pHcM <- matrix(NA,6,3,dimnames = list(NULL,c("pisy","piab","beal")))
pHcM[,1] <- c(-1.67918379,1.16550784,-0.23744806,0.19196957,0.07822739,-0.26698427)
pHcM[,2] <- c(-2.9597995,0.8591922,-0.3549810,0.3362082,0.2469336,0.1316909)
pHcM[,3] <- c(-1.9287089,1.0760549,-0.1079130,0.1922377,0.1363654,-0.3804504)
pPREL = c(413.000000, 0.450000, 0.118000, 3.000000, 0.745700, 10.930000, -3.063000, 17.720000,
            -0.102700, 0.036730, 0.777900, 0.500000, -0.364000, 0.271500, 0.835100, 0.073480,
            0.999600, 0.442800, 1.200000, 0.330000, 4.970496, 0.000000, 0.000000, 160.000000,
            0.000000, 0.000000, 20.000000, -999.000000, -999.000000, -999.000000)
pYAS = c(4.897147e-01, 4.913873e+00, 2.419735e-01, 9.487642e-02, 4.362893e-01, 2.499740e-01,
          9.151269e-01, 9.925823e-01, 8.385374e-02, 1.147678e-02, 6.083150e-04, 4.761282e-04,
          6.603773e-02, 7.713417e-04, 1.040174e-01, 6.488076e-01, -1.548718e-01, -1.956802e-02,
          -9.171713e-01, -4.035943e-04, -1.670727e-04, 9.059805e-02, -2.144096e-04, 4.877247e-02,
          -7.913602e-05, 3.518549e-02, -2.089906e-04, -1.808920e+00, -1.172547e+00, -1.253595e+01,
          4.596472e-03, 1.302583e-03, -4.389227e-01, 1.267467e+00, 2.569142e-01)
parsAWEN <- matrix(NA,12,3,dimnames = list(NULL,c("pisy","piab","beal")))
parsAWEN[,1] <- c(0.518000,0.177300,0.088700,0.216000,0.474660,0.019012,0.078308,
                  0.430248,0.670000,0.022500,0.007500,0.285000)
parsAWEN[,2] <- c(0.482600,0.131700,0.065800,0.319900,0.474660,0.019012,
                  0.078308,0.430248,0.665000,0.017500,0.002500,0.305000)
parsAWEN[,3] <- c(0.407900,0.198000,0.099000,0.295100,0.474660,0.019012,
                  0.078308,0.430248,0.715000,0.015000,0.000000,0.275000)

ClCut_birch <- matrix(NA,2,4)
ClCut_birch[1,] <- c(30.0,60,28.5,60)
ClCut_birch[2,] <- c(28.5,60,27.0,60)
ClCut_pine <- matrix(NA,3,4)
ClCut_pine[1,] <- c(29.0,70,26.0,80)
ClCut_pine[2,] <- c(27.5,80,25.0,90)
ClCut_pine[3,] <- c(24.0,90,23.5,100)
ClCut_spruce <- matrix(NA,2,4)
ClCut_spruce[1,] <- c(30,60,28.0,70)
ClCut_spruce[2,] <- c(28,70,26.5,80)


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
  W_croot = pmax(0.,(Lc * beta0 * A / par_betas * N + (W_c + Wsh) * beta0)) #coarse root biomass
  V = W_stem / par_rhow
  biomasses <- rbind(wf_STKG,W_froot,W_wsap,W_c,W_s,W_branch,W_croot,Wsh,Wdb,W_stem,V)
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


# Function to Compute Hc based on ksi parameter
ksiHcMod <- function(initVar){
  h <- initVar[3] - 1.3
  b <- pi * initVar[4]^2 / 40000
  p_ksi <- pCROBAS[38,initVar[1]]
  p_rhof <- pCROBAS[15,initVar[1]]
  p_z <- pCROBAS[11,initVar[1]]
  Lc <- h* ((p_rhof * b)/(p_ksi * h^p_z))^(1/(p_z-1))
  Hc <- max(0.,(h-Lc))
  return(Hc)
}

###function to replace HC NAs in initial variable initVar
findHcNAs <- function(initVar,pHcMod,HcModV){
  if(is.vector(initVar)){
    if(is.na(initVar[6])){
      if(HcModV==1){
        initVar[6] <- ksiHcMod(initVar)
      }else if(HcModV==2){
        inModHc <- c(pHcMod[,initVar[1]],initVar[3],
                     initVar[4],initVar[2],initVar[5],initVar[5])
        initVar[6] <- model.Hc(inModHc)
      }
    }
  } else if(any(is.na(initVar[6,]))){
    initVar[1,][which(initVar[1,]==0)] <- 1 ###deals with 0 species ID
    HcNAs <- which(is.na(initVar[6,]))
    BAtot <- sum(initVar[5,],na.rm = T)
    if(length(HcNAs)==1){
      if(HcModV==1){
        initVar[6,HcNAs] <- ksiHcMod(initVar[,HcNAs])
      }else if(HcModV==2){
        inModHc <- c(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                     initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],BAtot)
        initVar[6,HcNAs] <- model.Hc(inModHc)
      }
    }else{
      if(HcModV==1){
        initVar[6,HcNAs] <- apply(initVar,2,ksiHcMod)
      }else if(HcModV==2){
        inModHc <- rbind(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                         initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],BAtot)
        initVar[6,HcNAs] <- apply(inModHc,2,model.Hc)
      }
    }
  }
  return(initVar)
}


##Height of the crown base model
model.Hc <- function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
  lnHc_sim <- pValues[1]+pValues[2]*log(H)+pValues[3]*D/H+
    pValues[4]*log(age)+ pValues[5]*log(BA_sp)+
    pValues[6]*(BA_sp/BA_tot)
  Hc_sim <- exp(lnHc_sim)
  return(pmax(Hc_sim,0.,na.rm = T)) 
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
  
  
