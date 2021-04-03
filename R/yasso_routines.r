library(Matrix)


Yasso15Parameters<- c(4.8971473e-01,   4.9138734e+00,   2.4197346e-01,   9.4876416e-02,
                      4.3628932e-01,   2.4997402e-01,   9.1512685e-01,   9.9258227e-01,   8.3853738e-02,   1.1476783e-02,   6.0831497e-04,   4.7612821e-04,   6.6037729e-02,   7.7134168e-04,   1.0401742e-01,   6.4880756e-01,
                      -1.5487177e-01,  -1.9568024e-02,  -9.1717130e-01,  -4.0359430e-04,  -1.6707272e-04,
                      9.0598047e-02,  -2.1440956e-04,   4.8772465e-02,  -7.9136021e-05,   3.5185492e-02,  -2.0899057e-04,  -1.8089202e+00,  -1.1725473e+00,  -1.2535951e+01,
                      4.5964720e-03,   1.3025826e-03,
                      -4.3892271e-01,   1.2674668e+00,   2.5691424e-01)

# This file is the Yasso15 model core code that is an improved and updated version based on
# Yasso07 description Tuomi & Liski 17.3.2008  (Yasso07.pdf)
# and Taru Palosuo's code in December 2011
#
# This version uses the separate temperature/precipitation dependencies for the N and H compartments
# The parameters were estimated using e.g. the additional global scale Zinke data set
#
# Possibility to compute model prediction for steady state conditions has been included this version.
#
# Last edited 24.8.2015
# - Marko J??rvenp????


#  Instructions  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# 1) first run the source code for the function with "source(....r)"
# 2) then you can use the Yasso-function by just calling it Yasso15_R_version(...)
# 3) Input for the function:
# 	1. Yasso15Parameters - Yasso parameters as a vector, length 35
#		1-16 matrix A entries: 4*alpha, 12*p
#		17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION
#		22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2
#		24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2
#		26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2
#		28-30 Precipitation-dependence parameters for AWE, N and H fractions: gamma, gamma_N, gamma_H
#		31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)
#		33-35 Woody parameters: theta_1, theta_2, r
#	  2. SimulationTime - time when the result is requested [a]
#   3. MeanTemperature - mean annual temperature [C]
#   4. TemperatureAmplitude - temperature amplitude i.e. (T_max-T_min)/2, [C]
#   5. Precipitation - annual precipitation [mm]
#   6. InitialCPool - initial C pools of model compartments, length 5, [whatever]
#   7. LitterInput - mean litter input, 5 columns AWENH, must be the same unit as InitialCpool per year
#   8. WoodySize - size of woody litter (for non-woody litter this is 0) [cm]
#   9. Steadystate_pred - set to 1 if ignore 'SimulationTime' and compute solution
#      in steady-state conditions (which sould give equal solution as if time is set large enough)
# 4) The function returns the amount of litter as 5-vector (AWENH compartments) at SimulationTime
#    (or as time is infinity)

# NOTE that this function eats only one type of material at the time. So, non-woody and different woody litter
# materials needs to be calculated separately (and finally count together if desired).

# The output of the function is the vector AWENH compartments at the given time since the simulation start


# Basics  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# additional R libraries (as needed)
#library(Matrix)  # tai Matrix tms


# Function definition   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

Yasso15_R_version <- function(Yasso15Parameters, SimulationTime, MeanTemperature, TemperatureAmplitude,
                              Precipitation, InitialCPool, LitterInput, WoodySize, Steadystate_pred = 0)
{

  # using shorter names for input
  theta <- Yasso15Parameters
  t <- SimulationTime
  climate <- c(MeanTemperature,Precipitation,TemperatureAmplitude)
  init <- InitialCPool
  b <- LitterInput
  d <- WoodySize

  # temperature annual cycle approximation
  te1 <- climate[1]+4*climate[3]/pi*(1/sqrt(2)-1)
  te2 <- climate[1]-4*climate[3]/(sqrt(2)*pi)
  te3 <- climate[1]+4*climate[3]/pi*(1-1/sqrt(2))
  te4 <- climate[1]+4*climate[3]/(sqrt(2)*pi)

  # Average temperature dependence
  te <- c(te1,te2,te3,te4)
  tem <- mean(exp(theta[22]*te+theta[23]*te^2))
  temN <- mean(exp(theta[24]*te+theta[25]*te^2))
  temH <- mean(exp(theta[26]*te+theta[27]*te^2))

  # Precipitation dependence
  tem <- tem*(1.0-exp(theta[28]*climate[2]/1000.0)) # precipitation as m
  temN <- temN*(1.0-exp(theta[29]*climate[2]/1000.0))
  temH <- temH*(1.0-exp(theta[30]*climate[2]/1000.0))

  # Size class dependence -- no effect if d == 0.0
  size_dep <- min(1.0,(1+theta[33]*d+theta[34]*d^2)^(-abs(theta[35])))

  # check rare case where no decomposition happens for some of the compartments
  # (basically, if no rain)
  if(tem <= 1e-16)
  {
    xt <- init + b*t;
    return(xt);
  }

  # Calculating matrix A (will work ok despite the sign of alphas)

  alpha <- abs(c(theta[1], theta[2], theta[3], theta[4], theta[32]))  # Vector of decomposition rates

  # Creating the matrix A_p
  row1 <- c(-1, theta[5], theta[6], theta[7], 0)
  row2 <- c(theta[8], -1, theta[9], theta[10], 0)
  row3 <- c(theta[11], theta[12], -1, theta[13], 0)
  row4 <- c(theta[14], theta[15], theta[16], -1, 0)
  row5 <- c(theta[31], theta[31], theta[31], theta[31], -1)
  A_p <- matrix(c(row1, row2, row3, row4, row5), 5, 5, byrow=T) # A_p is now computed

  # computing the diagonal coefficient matrix k
  k <- diag(c(tem*alpha[1:3]*size_dep,temN*alpha[4]*size_dep,temH*alpha[5])) # no size effect in humus

  A <- A_p%*%k # coefficient matrix A is now computed

  # Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init
  if (Steadystate_pred)
  {
    # Solve DE directly in steady state conditions (time = infinity)
    # using the formula 0 = x'(t) = A*x + b => x = -A^-1*b
    xt <- solve(-A,b)
    xt <- as.matrix(xt)
  } else
  {
    # Solve DE in given time using the analytical formula
    z1 <- A %*% init + b;
    mexpAt <- expm(A*t)
    z2 <- mexpAt %*% z1 - b
    xt <- as.matrix(solve(A,z2))
  }
  return(xt)

} # end of Yasso15 function


# Note for Birch Betula pubenscens and brown leaves is used

foliage.AWEN <- function(Lf, spec) {

  fol.AWEN <- matrix(0,nrow=length(Lf), ncol=4)

  ma <- (1:length(Lf))[spec==1]

  ku <- (1:length(Lf))[spec==2]

  ko <- (1:length(Lf))[spec==3]

  fol.AWEN[,1][ma] <- 0.518*Lf[ma]

  fol.AWEN[,1][ku] <- 0.4826*Lf[ku]

  fol.AWEN[,1][ko] <- 0.4079*Lf[ko]

  fol.AWEN[,2][ma] <- 0.1773*Lf[ma]

  fol.AWEN[,2][ku] <- 0.1317*Lf[ku]

  fol.AWEN[,2][ko] <- 0.198*Lf[ko]

  fol.AWEN[,3][ma] <- 0.0887*Lf[ma]

  fol.AWEN[,3][ku] <- 0.0658*Lf[ku]

  fol.AWEN[,3][ko] <- 0.099*Lf[ko]

  fol.AWEN[,4][ma] <- 0.216*Lf[ma]

  fol.AWEN[,4][ku] <- 0.3199*Lf[ku]

  fol.AWEN[,4][ko] <- 0.2951*Lf[ko]

  return(fol.AWEN)

}

## Branches are here

# It seems that there is only valiues for pine (these are applied for others as well)



branches.AWEN <- function(Lb) {

  fb.AWEN <- matrix(0,nrow=length(Lb), ncol=4)

  a <- c(0.4763,0.4933,0.4289,0.5068,0.4607,0.5047,0.4642,0.5307,0.5256,0.4661,0.5060,

         0.4941,0.4848,0.4158,0.4605,0.4423,0.4811,0.4434,0.5141,0.4312,0.4867,0.3997,0.4758,0.4741,0.4996)

  w <- c(0.0196,0.0105,0.0197,0.0120,0.0107,0.0106,0.0130,0.0126,0.0116,0.0195,0.0180,

         0.0257,0.0219,0.0295,0.0242,0.0198,0.0242,0.0263,0.0188,0.0218,0.0207,0.0234,0.0176,0.0248,0.0188)

  e <- c(0.0870,0.0659,0.1309,0.0506,0.0874,0.0519,0.0840,0.0382,0.0394,0.0996,0.0647,

         0.0905,0.0633,0.1131,0.0874,0.1101,0.0681,0.1108,0.0561,0.1128,0.0452,0.1161,0.0678,0.0698,0.0470)

  n <- c(0.4170,0.4303,0.4205,0.4306,0.4412,0.4328,0.4388,0.4186,0.4234,0.4148,0.4112,

         0.4456,0.4300,0.4416,0.4279,0.4278,0.4266,0.4195,0.4110,0.4341,0.4474,0.4608,0.4388,0.4313,0.4346)

  fb.AWEN[,1] <- mean(a)*Lb

  fb.AWEN[,2] <- mean(w)*Lb

  fb.AWEN[,3] <- mean(e)*Lb

  fb.AWEN[,4] <- mean(n)*Lb

  return(fb.AWEN)

}



stem.AWEN <- function(Lst, spec) {

  st.AWEN <- matrix(0,nrow=length(Lst), ncol=4)

  ma <- (1:length(Lst))[spec==1]

  ku <- (1:length(Lst))[spec==2]

  ko <- (1:length(Lst))[spec==3]

  st.AWEN[,1][ma] <- 0.5*(0.66+0.68)*Lst[ma]

  st.AWEN[,1][ku] <- 0.5*(0.63+0.7)*Lst[ku]

  st.AWEN[,1][ko] <- 0.5*(0.65+0.78)*Lst[ko]

  st.AWEN[,2][ma] <- 0.5*(0.03+0.015)*Lst[ma]

  st.AWEN[,2][ku] <- 0.5*(0.03+0.005)*Lst[ku]

  st.AWEN[,2][ko] <- 0.5*(0.03+0)*Lst[ko]

  st.AWEN[,3][ma] <- 0.5*(0+0.015)*Lst[ma]

  st.AWEN[,3][ku] <- 0.5*(0+0.005)*Lst[ku]

  st.AWEN[,3][ko] <- 0

  st.AWEN[,4][ma] <- 0.5*(0.28+0.29)*Lst[ma]

  st.AWEN[,4][ku] <- 0.5*(0.33+0.28)*Lst[ku]

  st.AWEN[,4][ko] <- 0.5*(0.22+0.33)*Lst[ko]

  return(st.AWEN)

}

#YAsso function that to be used in R function apply
yasso <- function(input,pYasso,yassoOut){
  SimulationTime <- t
  Steadystate_pred <- Steadystate_pred
  MeanTemperature <- input[1]
  TemperatureAmplitude <- input[2]
  Precipitation <- input[3]
  InitialCPool <- input[4:8]
  LitterWp <- c(stem.AWEN(input[9],1),0)
  LitterfWp <- c(branches.AWEN(input[10]),0)
  LitterFp <- c(foliage.AWEN(input[11],1),0)
  LitterWsp <- c(stem.AWEN(input[12],2),0)
  LitterfWsp <- c(branches.AWEN(input[13]),0)
  LitterFsp <- c(foliage.AWEN(input[14],2),0)
  LitterWb <- c(stem.AWEN(input[15],3),0)
  LitterfWb <- c(branches.AWEN(input[16]),0)
  LitterFb <- c(foliage.AWEN(input[17],3),0)
  sizeWp <- input[18]
  sizefWp <- input[19]
  sizeWsp <- input[20]
  sizefWsp <- input[21]
  sizeWb <- input[22]
  sizefWb <- input[23]

  # run yasso for Woody in pine stands
  yassoOut[1,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[1,], LitterWp, sizeWp, Steadystate_pred)
  # run yasso for fine Woody in pine stands
  yassoOut[2,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[2,], LitterfWp, sizefWp, Steadystate_pred)
  # run yasso for foliage in pine stands
  yassoOut[3,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[3,], LitterFp, 0, Steadystate_pred)
  # run yasso for Woody in spruce stands
  yassoOut[4,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[4,], LitterWsp, sizeWsp, Steadystate_pred)
  # run yasso for fine Woody in spruce stands
  yassoOut[5,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[5,], LitterfWsp, sizefWsp, Steadystate_pred)
  # run yasso for foliage in spruce stands
  yassoOut[6,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[6,], LitterFsp, 0, Steadystate_pred)
  # run yasso for Woody in birch stands
  yassoOut[7,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[7,], LitterWb, sizeWb, Steadystate_pred)
  # run yasso for fine Woody in birch stands
  yassoOut[8,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[8,], LitterfWb, sizefWb, Steadystate_pred)
  # run yasso for foliage in birch stands
  yassoOut[9,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[9,], LitterFb, 0, Steadystate_pred)
  return(yassoOut)
}


compAWENH <- function(Lit,parsAWEN,spec,litType){
  #litType if 1 uses Foliage pars, 2 branches, 3 woody
  #spec vector with species
  #parsAWEN matrix with parameters and columns = species
  #Lit litter
  AWENH = numeric(5)
  AWENH[1] = parsAWEN[((litType-1)*4)+1,spec]*Lit
  AWENH[2] = parsAWEN[((litType-1)*4)+2,spec]*Lit
  AWENH[3] = parsAWEN[((litType-1)*4)+3,spec]*Lit
  AWENH[4] = parsAWEN[((litType-1)*4)+4,spec]*Lit
  return(AWENH)
}


#YAsso function that to be used in R function apply
yassoStSt <- function(SimulationTime,Steadystate_pred,
                   MeanTemperature,TemperatureAmplitude,Precipitation,
                   Lit.W,lit.fW,litt.F,InitialCPool,
                   pYasso){

  apply(ciao,stem.AWEN,1)
  LitterWp <- c(stem.AWEN(input[9],1),0)
  LitterfWp <- c(branches.AWEN(input[10]),0)
  LitterFp <- c(foliage.AWEN(input[11],1),0)
  LitterWsp <- c(stem.AWEN(input[12],2),0)
  LitterfWsp <- c(branches.AWEN(input[13]),0)
  LitterFsp <- c(foliage.AWEN(input[14],2),0)
  LitterWb <- c(stem.AWEN(input[15],3),0)
  LitterfWb <- c(branches.AWEN(input[16]),0)
  LitterFb <- c(foliage.AWEN(input[17],3),0)
  sizeWp <- input[18]
  sizefWp <- input[19]
  sizeWsp <- input[20]
  sizefWsp <- input[21]
  sizeWb <- input[22]
  sizefWb <- input[23]

  # run yasso for Woody in pine stands
  yassoOut[1,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[1,], LitterWp, sizeWp, Steadystate_pred)
  # run yasso for fine Woody in pine stands
  yassoOut[2,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[2,], LitterfWp, sizefWp, Steadystate_pred)
  # run yasso for foliage in pine stands
  yassoOut[3,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[3,], LitterFp, 0, Steadystate_pred)
  # run yasso for Woody in spruce stands
  yassoOut[4,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[4,], LitterWsp, sizeWsp, Steadystate_pred)
  # run yasso for fine Woody in spruce stands
  yassoOut[5,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[5,], LitterfWsp, sizefWsp, Steadystate_pred)
  # run yasso for foliage in spruce stands
  yassoOut[6,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[6,], LitterFsp, 0, Steadystate_pred)
  # run yasso for Woody in birch stands
  yassoOut[7,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[7,], LitterWb, sizeWb, Steadystate_pred)
  # run yasso for fine Woody in birch stands
  yassoOut[8,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[8,], LitterfWb, sizefWb, Steadystate_pred)
  # run yasso for foliage in birch stands
  yassoOut[9,] <- Yasso15_R_version(pYasso, SimulationTime, MeanTemperature, TemperatureAmplitude,
                                    Precipitation, yassoOut[9,], LitterFb, 0, Steadystate_pred)

  return(yassoOut)
}




####below are the functions to compute soil steady state carbon from prebas output
StStYasso <- function(litter,parsAWEN,spec,Tmean,Tamp,Precip,litterSize,litType,pYasso,
                      t=1, stst=1,soilCin=rep(0,5)){
  AWEN <- compAWENH(litter,parsAWEN,spec,litType)

  soilC <- Yasso15_R_version(pYasso, SimulationTime=t, Tmean,
                                 Tamp, Precip,soilCin,
                                 AWEN, litterSize,
                                 Steadystate_pred=stst)
  return(soilC)
}

ageEndRot <- function(x){
  if(inherits(x,"prebas")){
    nYearsStst = x$output[(which(x$output[,30,1,1]==0)[1]),7,1,1]
  }
  if(inherits(x,"multiPrebas") | inherits(x,"regionPrebas")){
    endRot <- apply(x$multiOut[,,30,1,1],1,function(aa) which(aa==0)[1])
    index <- matrix(c(1:x$nSites,endRot,rep(7,x$nSites),rep(1,x$nSites),rep(1,x$nSites)),7,5)
    nYearsStst <- x$multiOut[index]
  }
  return(nYearsStst)
}


LitterforYassoStSt <- function(x,rot=1,years=NA){
  ###Function to extract litter from PREBAS output and
  ###compute the average litter input for YASSO
  ## x is a prebas output; rot = 1 the steady state soilC is calculated
  ##on the rotation length;
  ## years = number of years from which compute the average annual litter inputs;
  if(inherits(x,"prebas")){
    if (rot==1) {nYearsStst = ageEndRot(x)} else if(!is.na(years)){
      nYearsStst=years}else{
        nYearsStst=length(x$output[,30,1,1])
      }
    litterSize <- x$litterSize
    nLayers <- x$nLayers
    if(nLayers==1){
      input <-  c(sum(colSums(x$output[1:nYearsStst,26:27,,1])),colSums(x$output[1:nYearsStst,28:29,,1]))/nYearsStst#array(0.,3,nLayers)
      litter <- data.table(input)
      setnames(litter,"layer 1")
      litter[,"litterSize 1":=litterSize[3:1,1]][]
      litter[,litType:=1:3]
      #    class(litter) <- "litterPrebas"
      return(litter)
    } else{
    input <-  rbind(colSums(colSums(x$output[1:nYearsStst,26:27,,1])),colSums(x$output[1:nYearsStst,28:29,,1]))/nYearsStst#array(0.,3,nLayers)
    litter <- data.table(input)
    for(j in 1:nLayers) litter[,paste("litterSize",j):=litterSize[3:1,j]][]
    litter[,litType:=1:3]
#    class(litter) <- "litterPrebas"
    return(litter)}
  }
  # if(inherits(x,"multiPrebas")){
  #   nSites <- x$nSites
  #   if (rot==1) {nYearsStst = ageEndRot(x)} else if(!is.na(years)){
  #     nYearsStst=years}else{
  #       nYearsStst=length(x$output[,30,1,1])
  #     }
  #   litterSize <- x$litterSize
  #   nLayers <- x$nLayers
  #   folLit <- x$multiOut[,,26,,1] + x$multiOut[,,27,,1]
  #   litter <-
  #   input <-  apply(x$multiOut,1,function(ops) rbind(colSums(colSums(ops[1:nYearsStst,26:27,,1])),colSums(ops[1:nYearsStst,28:29,,1])))#/nYearsStst#array(0.,3,nLayers)
  #   litter <- data.table(input)
  #   for(j in 1:nLayers) litter[,paste("litterSize",j):=litterSize[3:1,j]][]
  #   litter[,litType:=1:3]
  #   #    class(litter) <- "litterPrebas"
  #   return(litter)
  # }
}




soilCstst <- function(litter,Tmean,Tamp,Precip,species, ###species is a vector of species code with length = to nLayers
                      pAWEN = parsAWEN,pYasso=pYAS,
                      t=1,stst=1,soilCin=NA){

  if(length(dim(litter))==2){
    litter <- data.table(litter)
    nLayers <- (ncol(litter)-1)/2
    if(is.na(soilCin)) soilCin <- array(0,dim=c(5,3,nLayers))
    soilC = array(0,dim = c(5,3,nLayers))

    setnames(litter, c(paste("layer", 1:nLayers),paste("litterSize", 1:nLayers),"litType"))
    layersNam <- names(litter[,1:nLayers])
    litterSizeNam <- names(litter[,(nLayers+1):(nLayers*2)])

    for(j in 1:nLayers) soilC[,,j] <- matrix(unlist(litter[,
                                    .(list(StStYasso(get(layersNam[j]),
                                    parsAWEN=pAWEN,spec=species[j],Tmean,Tamp,Precip,
                                    get(litterSizeNam[j]),`litType`,pYasso,
                                    t=t,stst=stst,soilCin = soilCin[,,j]))),
                                    by=1:nrow(litter)][,2]),5,3)
    return(soilC)
  }
}



sCststPrebasOut <- function(x){
    litter <- LitterforYassoStSt(x)[]
    Tmean <- mean(x$weatherYasso[,1])
    Tamp <- mean(x$weatherYasso[,3])
    Precip <- mean(x$weatherYasso[,2])
    soilStSt <- soilCstst(litter, Tmean,Tamp,Precip,x$initVar[1,])
    return(soilStSt)
}


####Wrapper function for YASSO runs (fortran version) with PREBAS inputs
####Wrapper function for YASSO runs (fortran version) with PREBAS inputs
yassoPREBASinOldversion <- function(prebOut,initSoilC,pYASSO = pYAS, litterSize = NA, pAWEN=parsAWEN){
  ###litter is array with dimensions:(nSites, nYears, nLayers, 3) !!!fourth dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  ####species is a matrix dims= nSites,nLayers
  #### initSoilC is array dim=nSites,5,3,nLayers;;!!! third dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  #### weatherYasso dims=nClimID, nYears, 3; dim 3 are weather inputs Tmean precip Tampl
  ###### climIDs vector of climIDs(nSites)
  ##### litterSize dimensions of litterfall matrix (nrow=3(stem,branch,fineLit), ncol=nSp) 
  nSites <- prebOut$nSites
  nYears <- dim(prebOut$multiOut)[2]
  nLayers <- dim(prebOut$multiOut)[3]
  litter <- array(NA,dim=c(nSites, nYears, nLayers, 3))
  litter[,,,1] <- prebOut$multiOut[,,26,,1] + prebOut$multiOut[,,27,,1]
  litter[,,,2] <- prebOut$multiOut[,,28,,1]
  litter[,,,3] <- prebOut$multiOut[,,29,,1]
  species <- prebOut$multiOut[,1,4,,1]
  climIDs <- prebOut$siteInfo[,2]
  litterSize <- prebOut$litterSize
  nSp <- ncol(parsAWEN)
  nClimID <- dim(prebOut$weatherYasso)[1]
  soilC <- array(0., dim=c(nSites,(nYears+1),5,3,nLayers))
  soilC[,1,,,] <- initSoilC
  
  xx <- .Fortran("runYasso",
                 litter = as.array(litter),
                 litterSize = as.array(litterSize),
                 nYears = as.integer(nYears),
                 nLayers = as.integer(nLayers), 
                 nSites = as.integer(nSites), 
                 nSp = as.integer(nSp),
                 species = as.matrix(species),
                 nClimID = as.integer(nClimID),
                 climIDs = as.integer(climIDs),
                 pAWEN = as.matrix(pAWEN),
                 pYASSO=as.double(pYASSO),
                 weatherYasso = as.array(prebOut$weatherYasso),
                 soilC = as.array(soilC)
  )
  return(xx)
}



stXX <- function(PrebOut){
  Tmean <- apply(PrebOut$weatherYasso[,,1],1,mean)
  Pre <- apply(PrebOut$weatherYasso[,,2],1,mean)
  Tamp <- apply(PrebOut$weatherYasso[,,3],1,mean)
  
  nLayers <- dim(PrebOut$multiOut)[4]
  nSites <- dim(PrebOut$multiOut)[1]
  species <- PrebOut$multiOut[,1,4,,1]
  if(nLayers==1){
    litFol <- apply(PrebOut$multiOut[,,26,1,1],1,mean) +
      apply(PrebOut$multiOut[,,27,1,1],1,mean)
    litBra <- apply(PrebOut$multiOut[,,28,1,1],1,mean)
    litSte <- apply(PrebOut$multiOut[,,29,1,1],1,mean)
    AWENf <- AWENb <- AWENs <- matrix(NA,nSites,5)
    soilC <- array(NA,dim=c(nSites,5,3))
    for(sitx in 1:nSites){
      AWENf[sitx,] <- compAWENH(litFol[sitx],parsAWEN = parsAWEN,spec = species[sitx],litType = 1)
      AWENb[sitx,] <- compAWENH(litBra[sitx],parsAWEN = parsAWEN,spec = species[sitx],litType = 2)
      AWENs[sitx,] <- compAWENH(litSte[sitx],parsAWEN = parsAWEN,spec = species[sitx],litType = 3)
      
      soilC[sitx,,1] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                          Tamp[sitx], Pre[sitx],rep(0,5),
                                          AWENf[sitx,], litterSizeDef[3,species[sitx]],
                                          Steadystate_pred=1)
      soilC[sitx,,2] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                          Tamp[sitx], Pre[sitx],rep(0,5),
                                          AWENb[sitx,], litterSizeDef[2,species[sitx]],
                                          Steadystate_pred=1)
      soilC[sitx,,3] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                          Tamp[sitx], Pre[sitx],rep(0,5),
                                          AWENs[sitx,],litterSizeDef[1,species[sitx]],
                                          Steadystate_pred=1)
    }
  }else{
    litFol <- apply(PrebOut$multiOut[,,26,,1],c(1,3),mean) +
      apply(PrebOut$multiOut[,,27,,1],c(1,3),mean)
    litBra <- apply(PrebOut$multiOut[,,28,,1],c(1,3),mean)
    litSte <- apply(PrebOut$multiOut[,,29,,1],c(1,3),mean)
    AWENf <- AWENb <- AWENs <- array(NA,dim=c(nSites,5,nLayers))
    soilC <- array(NA,dim=c(nSites,5,3,nLayers))
    for(layx in 1:nLayers){
      for(sitx in 1:nSites){
        AWENf[sitx,,layx] <- compAWENH(litFol[sitx,layx],parsAWEN = parsAWEN,spec = species[sitx,layx],litType = 1)
        AWENb[sitx,,layx] <- compAWENH(litBra[sitx,layx],parsAWEN = parsAWEN,spec = species[sitx,layx],litType = 2)
        AWENs[sitx,,layx] <- compAWENH(litSte[sitx,layx],parsAWEN = parsAWEN,spec = species[sitx,layx],litType = 3)
        
        soilC[sitx,,1,layx] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                                 Tamp[sitx], Pre[sitx],rep(0,5),
                                                 AWENf[sitx,,layx], litterSizeDef[3,species[sitx,layx]],
                                                 Steadystate_pred=1)
        soilC[sitx,,2,layx] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                                 Tamp[sitx], Pre[sitx],rep(0,5),
                                                 AWENb[sitx,,layx], litterSizeDef[2,species[sitx,layx]],
                                                 Steadystate_pred=1)
        soilC[sitx,,3,layx] <- Yasso15_R_version(pYAS, SimulationTime=1, Tmean[sitx],
                                                 Tamp[sitx], Pre[sitx],rep(0,5),
                                                 AWENs[sitx,,layx],litterSizeDef[1,species[sitx,layx]],
                                                 Steadystate_pred=1)
      }
    }
  }
  return(soilC)
}





#### calculate steady state C using prebas output.
####included option for ground vegetation 
stXX_GV <- function(prebOut, GVrun,pYASSO = pYAS, litterSize = NA, pAWEN=parsAWEN){
  Tmean <- apply(prebOut$weatherYasso[,,1],1,mean)
  Pre <- apply(prebOut$weatherYasso[,,2],1,mean)
  Tamp <- apply(prebOut$weatherYasso[,,3],1,mean)
  
  nLayers <- dim(prebOut$multiOut)[4]
  nSites <- dim(prebOut$multiOut)[1]
  species <- prebOut$multiOut[,1,4,,1]
  climIDs <- prebOut$siteInfo[,2]
  nYears <- dim(prebOut$multiOut)[2]
  
  if(GVrun==1){
    ###calculate steady state C for gv
    fAPAR <- prebOut$fAPAR
    fAPAR[which(is.na(prebOut$fAPAR),arr.ind = T)] <- 0.
    ###calculate soil C for gv
    fAPAR <- prebOut$fAPAR
    fAPARgv <- litGV <- matrix(0,nSites,nYears)
    fAPAR[which(is.na(prebOut$fAPAR),arr.ind = T)] <- 0.
    fAPAR[which(prebOut$fAPAR>1,arr.ind = T)] <- 1.
    fAPAR[which(prebOut$fAPAR<0,arr.ind = T)] <- 0.
    AWENgv <- array(0.,dim=c(dim(prebOut$fAPAR),4))
    soilCgv <- array(0.,dim=c(nSites,(nYears+1),5))
    p0 = prebOut$multiOut[,,6,1,1]
    ETSy = prebOut$multiOut[,,5,1,1]

    AWENgv <- .Fortran("multiGV",
                       fAPAR=as.matrix(fAPAR),
                       ETS=as.matrix(ETSy),
                       siteType = as.double(prebOut$siteInfo[,3]),
                       fAPARgv=as.matrix(fAPARgv),
                       litGV=as.matrix(litGV),
                       p0=as.matrix(p0),
                       AWENgv=as.array(AWENgv[,,1:4]),
                       nYears = as.integer(nYears),
                       nSites = as.integer(nSites))$AWENgv
    
    # for(ij in 1:nYears){
    #   AWENgv[,ij,] <- t(sapply(1:nrow(fAPAR), function(i) .Fortran("fAPARgv",fAPAR[i,ij],
    #                                               prebOut$multiOut[i,ij,5,1,1],
    #                                               prebOut$siteInfo[i,3],
    #                                               0,
    #                                               0,
    #                                               prebOut$multiOut[i,ij,6,1,1],
    #                                               rep(0,4))[[7]]))
    # }
    AWENgv2 <- apply(AWENgv,c(1,3),mean,na.rm=T)
    
    
    ###calculate steady state soil C per GV
    # ststGV <- matrix(NA,nSites,5)
    climIDs <- prebOut$siteInfo[,2]
    ststGV <- t(sapply(1:nSites, function(ij) .Fortran("mod5c",
                                                       pYAS,1.,colMeans(prebOut$weatherYasso[climIDs[ij],,]),
                                                       rep(0,5),c(AWENgv2[ij,],0),litterSize=0,leac=0.,rep(0,5),
                                                       stSt=1.)[[8]]))
  }
  if(nLayers>1){
    litter <- array(0.,dim=c(nSites, nLayers, 3))
    litter[,,1] <- apply(prebOut$multiOut[,,26,,1],c(1,3),mean,na.rm=T) + 
      apply(prebOut$multiOut[,,27,,1],c(1,3),mean,na.rm=T)
    litter[,,2] <- apply(prebOut$multiOut[,,28,,1],c(1,3),mean,na.rm=T)
    litter[,,3] <- apply(prebOut$multiOut[,,29,,1],c(1,3),mean,na.rm=T)
    species <- prebOut$multiOut[,1,4,,1]
    climIDs <- prebOut$siteInfo[,2]
    litterSize <- prebOut$litterSize
    nSp <- ncol(parsAWEN)
    weatherYasso <- apply(prebOut$weatherYasso,c(1,3),mean)
    nClimID <- dim(prebOut$weatherYasso)[1]
    soilC <- array(0.,dim=c(nSites,5,3,nLayers))
    
    soilC <- .Fortran("StstYasso",
                      litter = as.array(litter),
                      litterSize = as.matrix(litterSize),
                      nLayers = as.integer(nLayers), 
                      nSites = as.integer(nSites), 
                      nSp = as.integer(nSp),
                      species = as.matrix(species),
                      nClimID = as.integer(nClimID),
                      climIDs = as.integer(climIDs),
                      pAWEN = as.matrix(pAWEN),
                      pYASSO=as.double(pYASSO),
                      weatherYasso = as.matrix(weatherYasso),
                      soilC = as.array(soilC)
    )$soilC
    
    ####add gvsoilc to first layer foliage soilC
    # check in normal runs where ground vegetation soilC is calculated
    if(GVrun==1){
      soilC[,,1,1] <- soilC[,,1,1] + ststGV
    }  
  }else{
    # print("stXX_GV function needs to be coded for one layer prebas outut")
    litter <- array(NA,dim=c(nSites, 3))
    litter[,1] <- apply(prebOut$multiOut[,,26,1,1],1,mean,na.rm=T) + 
      apply(prebOut$multiOut[,,27,1,1],1,mean,na.rm=T)
    litter[,2] <- apply(prebOut$multiOut[,,28,1,1],1,mean,na.rm=T)
    litter[,3] <- apply(prebOut$multiOut[,,29,1,1],1,mean,na.rm=T)
    species <- prebOut$multiOut[,1,4,1,1]
    climIDs <- prebOut$siteInfo[,2]
    litterSize <- prebOut$litterSize
    nSp <- ncol(parsAWEN)
    weatherYasso <- apply(prebOut$weatherYasso,c(1,3),mean)
    nClimID <- prebOut$nClimID
    soilC <- array(0.,dim=c(nSites,5,3,nLayers))
    
    soilC <- .Fortran("StstYasso",
                      litter = as.array(litter),
                      litterSize = as.matrix(litterSize),
                      nLayers = as.integer(nLayers), 
                      nSites = as.integer(nSites), 
                      nSp = as.integer(nSp),
                      species = as.matrix(species),
                      nClimID = as.integer(nClimID),
                      climIDs = as.integer(climIDs),
                      pAWEN = as.matrix(pAWEN),
                      pYASSO=as.double(pYASSO),
                      weatherYasso = as.matrix(weatherYasso),
                      soilC = as.array(soilC)
    )$soilC
    
    ####add gvsoilc to first layer foliage soilC
    # check in normal runs where ground vegetation soilC is calculated
    if(GVrun==1){
      soilC[,,1,1] <- soilC[,,1,1] + ststGV
    }  
  }
  return(soilC)
}





####Wrapper function for YASSO runs (fortran version) with PREBAS inputs
yassoPREBASin <- function(prebOut,initSoilC,pYASSO = pYAS, litterSize = NA, pAWEN=parsAWEN){
  ###litter is array with dimensions:(nSites, nYears, nLayers, 3) !!!fourth dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  ####species is a matrix dims= nSites,nLayers
  #### initSoilC is array dim=nSites,5,3,nLayers;;!!! third dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  #### weatherYasso dims=nClimID, nYears, 3; dim 3 are weather inputs Tmean precip Tampl
  ###### climIDs vector of climIDs(nSites)
  ##### litterSize dimensions of litterfall matrix (nrow=3(stem,branch,fineLit), ncol=nSp) 
  nSites <- prebOut$nSites
  nYears <- dim(prebOut$multiOut)[2]
  nLayers <- dim(prebOut$multiOut)[4]
  litter <- array(NA,dim=c(nSites, nYears, nLayers, 3))
  litter[,,,1] <- prebOut$multiOut[,,26,,1] + prebOut$multiOut[,,27,,1]
  litter[,,,2] <- prebOut$multiOut[,,28,,1]
  litter[,,,3] <- prebOut$multiOut[,,29,,1]
  litter[which(is.na(litter))] <- 0.
  species <- prebOut$multiOut[,1,4,,1]
  climIDs <- prebOut$siteInfo[,2]
  if(all(is.na(litterSize))) litterSize <- prebOut$litterSize
  nSp <- ncol(parsAWEN)
  nClimID <- prebOut$nClimID
  soilC <- array(0., dim=c(nSites,(nYears+1),5,3,nLayers))
  soilC[,1,,,] <- initSoilC
  
  soilC <- .Fortran("runYasso",
                    litter = as.array(litter),
                    litterSize = as.array(litterSize),
                    nYears = as.integer(nYears),
                    nLayers = as.integer(nLayers), 
                    nSites = as.integer(nSites), 
                    nSp = as.integer(nSp),
                    species = as.matrix(species),
                    nClimID = as.integer(nClimID),
                    climIDs = as.integer(climIDs),
                    pAWEN = as.matrix(pAWEN),
                    pYASSO=as.double(pYASSO),
                    weatherYasso = as.array(prebOut$weatherYasso),
                    soilC = as.array(soilC)
  )$soilC
  
  
  ###calculate soil C for gv
  fAPAR <- prebOut$fAPAR
  fAPARgv <- litGV <- matrix(0,nSites,nYears)
  fAPAR[which(is.na(prebOut$fAPAR),arr.ind = T)] <- 0.
  fAPAR[which(prebOut$fAPAR>1,arr.ind = T)] <- 1.
  fAPAR[which(prebOut$fAPAR<0,arr.ind = T)] <- 0.
  AWENgv <- array(0.,dim=c(dim(prebOut$fAPAR),5))
  soilCgv <- array(0.,dim=c(nSites,(nYears+1),5))
  p0 = prebOut$multiOut[,,6,1,1]
  ETSy = prebOut$multiOut[,,5,1,1]
  
  AWENgv <- .Fortran("multiGV",
                     fAPAR=as.matrix(fAPAR),
                     ETS=as.matrix(ETSy),
                     siteType = as.double(prebOut$siteInfo[,3]),
                     fAPARgv=as.matrix(fAPARgv),
                     litGV=as.matrix(litGV),
                     p0=as.matrix(p0),
                     AWENgv=as.array(AWENgv[,,1:4]),
                     nYears = as.integer(nYears),
                     nSites = as.integer(nSites))$AWENgv
  
  soilCgv <- .Fortran("runYassoAWENin",
                      mAWEN=as.array(AWENgv),
                      nYears=as.integer(nYears),  
                      nSites=as.integer(nSites), 
                      litSize=as.double(0.),
                      nClimID=as.integer(nClimID),
                      climIDs=as.integer(climIDs),
                      pYasso=as.double(pYASSO),
                      climate=as.matrix(prebOut$weatherYasso),
                      soilCgv=as.array(soilCgv))$soilCgv
  
  
  soilC[,,,1,1] = soilC[,,,1,1] + soilCgv
  
  ###update model output fluxes
  soilByLay <- apply(soilC,c(1:2,5),sum)
  
  prebOut$soilC <- soilC[,2:(nYears+1),,,]
  prebOut$soilCtot <- apply(soilC[,2:(nYears+1),,,],1:2,sum)
  
  if(nLayers == 1){
    margins <- c(1:2)
  }else{
    margins <- c(1:2,5)
  }
  prebOut$multiOut[,,39,,1] <- apply(prebOut$soilC,margins,sum)
  
  if(nLayers == 1){
    margins <- c(1:2)
  }else{
    margins <- c(1:2,4)
  }
  #Rh calculations
  prebOut$multiOut[,,45,,1] <- (soilByLay[,1:nYears,] - soilByLay[,2:(nYears+1),] +
                                  apply(prebOut$multiOut[,1:nYears,26:29,,1],margins,sum))/10 ###/10 converts kgC/ha to gC/m2
  ###add GVlitter to first layer of multiout of Rh
  prebOut$multiOut[,,45,1,1] <- apply(AWENgv,1:2,sum)/10 + prebOut$multiOut[,,45,1,1] 
  ##NEP calculations
  prebOut$multiOut[,,46,,1] <-  prebOut$multiOut[,,44,,1] - prebOut$multiOut[,,9,,1] - 
    prebOut$multiOut[,,45,,1]
  ###add GVnpp to first layer of multiout of NEP  
  prebOut$multiOut[,,46,1,1] <- prebOut$multiOut[,,46,1,1] - prebOut$GVout[,,3]*0.5
  
  
  return(prebOut)
}
