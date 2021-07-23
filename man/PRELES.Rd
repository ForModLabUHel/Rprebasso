\name{PRELES}
\alias{PRELES}
\alias{preles}
\alias{Preles}
\title{A simple semi-empirical ecosystem carbon and water balance model.
}
\description{The model predicts gross primary production and evapotranpiration (and soil water balance) based on embedded empirical relationships and basic meteorological inputs.
}
\usage{ 
  PRELES(PAR, TAir, VPD, Precip, CO2, fAPAR, 
    GPPmeas = NA, ETmeas = NA, SWmeas = NA, 
      p = rep(NA, 30), DOY = NA, LOGFLAG = 0, control = 0,   
        returncols = c("GPP", "ET", "SW")) 
}
\arguments{
  \item{PAR}{A numeric vector of daily sums of photosynthetically active radiation, mmol/m2.  }
  \item{TAir}{A numeric vector of daily mean temperature, degrees C. }
  \item{VPD}{A numeric vector of daily mean vapour pressure deficits}
  \item{Precip}{A numeric vector of daily rainfall, mm}
  \item{CO2}{A numeric vector of air CO2}
  \item{fAPAR}{A numeric vector of fractions of absorbed PAR by the canopy, 0-1 unitless}
  
  
  ## OPTIONAL FOR BYPASSING PREDICTION. FOR DEVELOPMENT ONLY
  \item{GPPmeas}{NA}
  \item{ETmeas}{NA}
  \item{SWmeas}{NA}
  
  
  \item{p}{parameter vector of length 30. If parameter has value NA it is replaced with a default corresponding to Hyytiälä calibration. Vector p has following defaults:
    \itemize{
    \item SITE AND SOIL RELATED
    \item 1 soildepth = 413.0, soildepth
    \item 2 Effective field capacity = 0.450, ThetaFC
    \item 3 Permanent wilting point = 0.118, ThetaPWP
    \item 4 Drainage delay = 3 days, tauDrainage
    
    \item GPP_MODEL_PARAMETERS
    \item 5 LUE = 0.748018, betaGPP
    \item 6 tauGPP = 13.23383,  tauGPP
    \item 7 S0GPP = -3.9657867,   S0GPP
    \item 8 SmaxGPP = 18.76696,  SmaxGPP
    \item 9 kappaGPP = -0.130473,   kappaGPP
    \item 10 gammaGPP = 0.034459,   gammaGPP
    \item 11 soilthresGPP = 0.450828,  soilthresGPP
    \item 12 bCO2 = 0.5, CO2 modifier in GPP model
    \item 13 xCO2 = -0.364, CO2 modifier in ET model
    
    \item EVAPOTRANSPIRATION_PARAMETERS
    \item 14 transpiration efficiency = 0.324463, betaET
    \item 15 kappaET = 0.874151, kappaET
    \item 16 chiET = 0.075601, chiET
    \item 17 soilthresET = 0.541605, soilthresET
    \item 18 nuET = 0.273584, nu ET
    
    \item ## SNOW_RAIN_PARAMETERS
    \item 19 Metlcoeff = 1.2, Meltcoef
    \item 20 I_0 = 0.33, I_0
    \item 21 CWmax = 4.970496, CWmax, i.e. max canopy water
    \item 22 SnowThreshold = 0, SnowThreshold, 
    \item 23 T_0 = 0, T_0,
    
    \item ## START INITIALISATION PARAMETERS 
    \item 24 SWinit = 200, SWinit, Initial soil water
    \item 25 CWinit = 0, CWinit, Initial Canopy water
    \item 26 SOGinit = 0, SOGinit, Initial Snow on Ground 
    \item 27 Sinit = 20, Sinit, Initial temperature acclimation state
    \item 28 t0 = -999, t0 fPheno_start_date_Tsum_accumulation; conif -999, for birch 57
    \item 29 tcrit =  -999, tcrit, fPheno_start_date_Tsum_Tthreshold, 1.5 birch
    \item 30 tsumcrit = -999, tsumcrit, fPheno_budburst_Tsum, 134 birch}
    }
  
  \item{DOY}{Day of year integer vector. Needed for the prediction of deciduous species phenology. If DOY is not provided and deciduous species phenology parameters are not set to -999, it is assumed that the first values of all input vectors are from 1st Jan, and year has 365 day.}

  \item{LOGFLAG}{levels 0 (default), 1, 2. Generates increasing logging to preles.log file in the run directory.}
  
  \item{control}{
    Parameter that selects the transpiration model. Equals etmodel in c-code.
    
    if (etmodel == 0) 
    et = D * ET_par.beta*A/pow(D, ET_par.kappa) *
    pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
    fCO2mean +  // Mean effect of CO2 on transpiration
    ET_par.chi *  s / (s + psychom) * (1-fAPAR) *  fWsub * ppfd;
    
    if (etmodel == 1) 
    et = D * ET_par.beta*A/pow(D, ET_par.kappa) *
    pow(fWgpp, ET_par.nu) * 
    fCO2mean +  
    ET_par.chi * (1-fAPAR) *  fWsub * ppfd;
    
    if (etmodel == 2) 
    et = D * (1 + ET_par.beta/pow(D, ET_par.kappa)) * A / CO2 * 
    pow(fWgpp, ET_par.nu) *
    fCO2mean +  
    ET_par.chi * (1-fAPAR) *  fWsub * ppfd;
    
  }

  \item{returncols}{What columns are returned, defaults to
    c('GPP','ET','SW'). Other possible options are fD for vapour
    pressure deficit modifier, fW for soil water modifier, fE (minimum
    of fD, fW), fS (temperature/season modifier), SOG (snow on ground,
    mm weq.),  Throughfall, Interception, Snowmelt (mm per day),
    Drainage, i.e. water melted or precipitated that is above field capacity runs off with a small delay (see parameter tau), Canopywater for surfacial water storage (upper limit set by parameter Cmax and fAPAR), S season status (C) for fS calculation. }
}


\seealso{
  Package information: \code{\link{Rpreles}}
}

\examples{
## Run model with artificial inputs
CO2 <- 280:(2*380)
T=rep(18, length(CO2))
fAPAR=rep(1, length(CO2))
Precip=rep(3, length(CO2))
PAR=rep(20, length(CO2))

## Plot CO2 effect on GPP, ET, and SW. Feedbacks through soil
## eliminated with high precipitation
pdf('testCO2.pdf', hei=10, wid=10)
op <- par(mfrow=c(4,4), mar=c(1,1,1,1), oma=c(4,4,4,4))
for (D in c(0, 0.5, 1, 1.5)) {
    D <- rep(D, length(CO2))
    o1 <- PRELES(PAR, T, D, Precip, CO2, fAPAR,
  returncols=c("GPP", "ET", "SW", "fW", "fE"), LOGFLAG=0)
    plot(CO2, o1$GPP)
abline(v=380)
    plot(CO2, o1$ET)
abline(v=380)
    plot(CO2, o1$GPP/o1$ET)
abline(v=380)
    plot(CO2, o1$SW)
abline(v=380)
}
}
