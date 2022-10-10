\name{prebas}
\alias{prebas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PREBAS forest model
%%  ~~function to do ... ~~
}
\description{
This function runs PREBAS model.
}
\usage{
prebas(nYears, pCROBAS = pCROB, pPRELES = pPREL, etmodel = 0, pYASSO = pYAS, pAWEN = parsAWEN, siteInfo = NA, initVar = NA,
  thinning=NA,initClearcut = c(1.5,0.5,0.0431969,0.,0.),PAR,TAir,VPD,Precip,CO2, P0=NA,  soilC = NA,
  weatherYasso = NA, litterSize = NA, soilCtot = numeric(nYears), defaultThin = 1.,
  ClCut = 1., inDclct = NA, inAclct = NA, yassoRun = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nYears}{
  Number of years to run the model
%%     ~~Describe \code{nYears} here~~
}
  \item{pCROBAS}{
  Matrix of parameter sets, each column corresponds to a species. Default values pCROBAS = pCROB are the parameter sets for Scots pine, Norway spruce and Silver birch. print(pCROB) to see the parameter values and names.
%%     ~~Describe \code{pCROBAS} here~~
}
  \item{pPRELES}{
  Parameter vector for PRELES model. Default values (pPRELES = pPREL) is the calibration for European Boreal forests (Minunno et al. 2016)
%%     ~~Describe \code{pPRELES} here~~
}
  \item{etmodel}{
  Evapotranspiration model for PRELES. Default etmodel = 0. Possible values -1, 0, 1, 2
%%     ~~Describe \code{pPRELES} here~~
}
  \item{pYASSO}{
  Parameter vector for YASSO model. Default pYASSO = pYAS (Liski et al., 2017)
%%     ~~Describe \code{pYASSO} here~~
}
  \item{pAWEN}{
  Matrix of parameter sets for partitioning tree organs (foliage, branch and woody) in AWEN pools. Columns referes to different species. Note that the columns of pAWEN must correspond to the species of pPRELES. Default pAWEN = parsAWEN are for Scots pine, Norway spruce and Silver birch.
%%     ~~Describe \code{pYASSO} here~~
}
  \item{siteInfo}{
  Vector of site info SiteID, climID, siteType, SWinit (initial soil water), CWinit (initial crown water), SOGinit (initial snow on ground), Sinit (initial temperature acclimation state), soildepth, effective field capacity, permanent wilthing point. Default = c(1,1,3, 160, 0, 0, 20, 413.,0.45, 0.118), i.e. siteType = 3.
%%     ~~Describe \code{siteInfo} here~~
}
  \item{initVar}{
  Matrix with initial stand values for all the tree strata. Columns correspond to the layers in the stand. Initial information needed are: SpeciesID (a number corresponding to the species parameter values of pPRELES columns), Age (years), average height of the layer (H, m), average diameter at breast height of the layer (D, cm), basal area of the layer (BA, m2 ha-1), average height of the crown base of the layer (Hc, m). If initVar is not provided the model is initialized from plantation using default planting parameters (see initClearcut) and assuming that Pine, Spruce and Birch are equally present at plantation.
%%     ~~Describe \code{initVar} here~~
}
  \item{thinning}{
  A matrix with thinnig inputs. Rows correspond to a thinning event. Column 1 year from the start of the simulation; column 2 is siteID; column 3 layer where thinnings are carried out; column 4 to 7 stand variables (H, D, B, Hc); column 8 parameter that indicates if the stand variables (column 4:7) are provided as fraction of the actual model outputs (value=1 means that fraction is used); column 9 is the stand density after thinning if its value is not -999; colum 10 is Sapwood area of average tree at crown base (m2) if its value is not -999 (see examples).
%%     ~~Describe \code{thinning} here~~
}
  \item{initClearcut}{
  A numeric vector with initial stand variables after clearcut: H, D, BA, Hc, Ainit. Ainit is the year when the stand reaches the default values = c(1.5,0.5,0.0431969,0.,0.) and is automatically computed using air temperature.
%%     ~~Describe \code{initClearcut} here~~
}
  \item{PAR}{
  A numeric vector of daily sums of photosynthetically active radiation, mmol/m2.
%%     ~~Describe \code{PAR} here~~
}
  \item{TAir}{
  A numeric vector of daily mean temperature, degrees C.
%%     ~~Describe \code{TAir} here~~
}
  \item{VPD}{
  A numeric vector of daily mean vapour pressure deficits, kPa.
%%     ~~Describe \code{VPD} here~~
}
  \item{Precip}{
  A numeric vector of daily rainfall, mm
%%     ~~Describe \code{Precip} here~~
}
  \item{CO2}{
  A numeric vector of air CO2, ppm
%%     ~~Describe \code{CO2} here~~
}
  \item{P0}{
  A numeric vector with the annual potential photosynthesis (gC m-2 y-1). If P0 is not provided PRELES is used to compute P0 using fAPAR = 1.
%%     ~~Describe \code{P0} here~~
}
  \item{soilC}{
   Initial soil carbon compartments for each layer. Array with dimentions = c(nYears,5,3,nLayers). The second dimention (5) corresponds to the AWENH pools; the third dimention (3) corresponds to the tree organs (foliage, branch and stem).
%%     ~~Describe \code{soilC} here~~
}
  \item{weatherYasso}{
  Annual weather inputs for Yasso.
%%     ~~Describe \code{weatherYasso} here~~
}
  \item{litterSize}{
  Marix with litter inputs for YASSO. Rows are tree organs, columns correspond to the layers.
%%     ~~Describe \code{litterSize} here~~
}
  \item{soilCtot}{
  Vector with total initial soil carbon
%%     ~~Describe \code{soilCtot} here~~
}
  \item{defaultThin}{
  If defaultThin = 1 (default) Finnish standard managment practices are applied (ref).
%%     ~~Describe \code{defaultThin} here~~
}
  \item{ClCut}{
  If ClCut = 1 clearcuts are applied. If inDclct = NA and inAclct = NA Finnish standard clearcut practices are applied (ref).
%%     ~~Describe \code{ClCut} here~~
}
  \item{inDclct}{
  Vector of Diameter (cm) threshold for clearcut. Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
%%     ~~Describe \code{inDclct} here~~
}
  \item{inAclct}{
  Vector of Age (year) threshold for clearcut.  Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
%%     ~~Describe \code{inAclct} here~~
}
  \item{yassoRun}{
  If yassoRun=1 the YASSO model is run to compute the carbon balance of the soil.
%%     ~~Describe \code{yassoRun} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{output}{
  An array with annual model outputs. 1st dimension corresponds to the number of years of the simulation (nYears); 2nd dimension corresponds to the output variables (see list below); 3rd dimension corresponds to the number of layers in the stand (nLayers); 4th dimensions reports the state of the stand (1) and (2) the variables of the harvested trees (2).

  Output variables: \cr
  1."siteID" \cr
  2."gammaC" internal parameter \cr
  3."sitetype" site fertility class \cr
  4."species" \cr
  5."ETS" effective temperature sums \cr
  6."P0" Potential annual gross primary production (gC m-2 y-1) \cr
  7."age" Age of the layer (years) \cr
  8."DeadWoodVolume" Dead wood volume (m3 ha-1) \cr
  9."Respi_tot" Autotrophic respiration (gC m-2 y-1) \cr
  10."GPP/1000" Total tree GPP  (kgC m-2 y-1) \cr
  11."H" Layer average height (m) \cr
  12."D" Layer average diameter at breast height (cm) \cr
  13."BA" Layer basal area (m-2 ha-1) \cr
  14."Hc_base" Layer Base of crown height (m) \cr
  15."Cw" Crown width (m) \cr
  16."A" Sapwood area of average tree at crown base (m2) \cr
  17."N" Layer density \cr
  18."npp" net primary production (gC m-2 y-1) \cr
  19."leff" Effective leaf area \cr
  20."keff" Effective light extintion coefficient \cr
  21."lproj" Projected leaf area \cr
  22."ET_preles" Annual evapotranspiration (mm y-1) \cr
  23."weight" Layer weight on photosynthesis \cr
  24."Wbranch" Living Branch biomass (kgC ha-1) \cr
  25."WfineRoots" Fine roots biomass (kgC ha-1) \cr
  26."Litter_fol" Foliage litter (kgC ha-1) \cr
  27."Litter_fr" Fine root litter (kgC ha-1) \cr
  28."Litter_fWoody" fine woody litter (kgC ha-1) \cr
  29."Litter_cWood" coarse woody litter (kgC ha-1) \cr
  30."V" Layer volume (m3 ha-1) \cr
  31."Wstem" Stem Biomass (kgC ha-1) \cr
  32."W_croot" Course root Biomass (kgC ha-1) \cr
  33."wf_STKG" Foliage biomass (kgC ha-1) \cr
  34."wf_treeKG" Foliage biomass of the average tree (kgC ha-1) \cr
  35."B_tree" Basal area of average tree (m2) \cr
  36."Light" The proportion of light that has not been intercepted by canopy (meanlight)  \cr
  37."VroundWood" harvested round wood volume (m3 ha-1) \cr
  38."WroundWood" haversted round wood biomass (kgC ha-1) \cr
  39."soilC" totaal soil carbon (kgC ha-1) \cr
  40."aSW" average available soil water (mm) \cr
  41."dH" height growth (m) \cr
  42."Vmort" volume of dead trees (m3 ha-1) \cr
  43."grossGrowth" gross growth (m3 ha-1 y-1) \cr
  44."GPPtrees" Gross primary production per tree layer (gC m-2 y-1) \cr
  45."Rh" heterotrophic respiration (gC m-2 y-1) \cr
  46."NEP" Net ecosystem exchange (gC m-2 y-1), note 1st layer include fluxes from ground vegetation \cr
  47." W_wsap" sapwood biomass (kgC ha-1) \cr
  48."W_c" sapwood stem below Crown (kgC ha-1) \cr
  49."W_s" sapwood stem within crown (kgC ha-1) \cr
  50."Wsh" biomass of stem heartwood  (kgC ha-1) \cr
  51."Wdb" biomass of dead branches on living trees (kgC ha-1) \cr
  52."dHc" Height of the crown base change (m) \cr
  53."Wbh" biomass of branches heartwood (kgC ha-1) \cr
  54."Wcrh" biomass of coarse root heartwood (kgC ha-1)
  }
  \item{dailyPRELES}{
    Matrix with daily output from preles. 1st column is the daily GPP (gC m-2 d-1), 2nd column daily evapotranspiration (mm m-2), 3rd column daily soil water (mm -2).
  }%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francesco Minunno
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nYears, pCROBAS = pCROB, pPRELES = pPREL, pYASSO = pYAS,
    siteInfo = NA, thinning = matrix(0, 2, 10), nThinning = max(2,
        nrow(thinning)), initClearcut = c(1.5, 0.5, 0.0431969,
        0, 0), PAR, TAir, VPD, Precip, CO2, P0 = NA, initVar = NA,
    soilC = array(0, dim = c(nYears, 5, 3, 3)), weatherYasso = matrix(0,
        nYears, 3), litterSize = matrix(0, 3, 3), soilCtot = numeric(nYears),
    defaultThin = 1, ClCut = 1, inDclct = NA, inAclct = NA, yassoRun = 0)
{
    nLayers <- ifelse(!is.na(initVar), ncol(initVar), 3)
    nSp = ncol(pCROBAS)
    if (is.na(siteInfo))
        siteInfo = c(1, 1, 3)
    varNam <- getVarNam()
    nVar <- length(varNam)
    nLayers <- nSp <- ncol(pCROBAS)
    layerNam <- paste("layer", 1:nSp)
    output <- array(0, dim = c((nYears), nVar, nSp, 2), dimnames = list(NULL,
        varNam, layerNam, c("stand", "thinned")))
    fAPAR <- rep(0.7, nYears)
    Temp <- TAir[1:(365 * nYears)] - 5
    ETS <- pmax(0, Temp, na.rm = T)
    ETS <- matrix(ETS, 365, nYears)
    ETS <- colSums(ETS)
    if (is.na(P0)) {
        P0 <- PRELES(DOY = rep(1:365, nYears), PAR = PAR, TAir = TAir,
            VPD = VPD, Precip = Precip, CO2 = CO2, fAPAR = rep(1,
                length(PAR)), LOGFLAG = 0, p = pPRELES)$GPP
        P0 <- matrix(P0, 365, nYears)
        P0 <- colSums(P0)
    }
    ETSthres <- 1000
    ETSmean <- mean(ETS)
    if (any(!is.na(c(inDclct, inAclct)))) {
        if (is.na(inDclct))
            inDclct <- 9999999.99
        if (is.na(inAclct))
            inAclct <- 9999999.99
    }
    if (ClCut == 1 & is.na(initVar) & is.na(inDclct))
        inDclct <- c(ClCutD_Pine(ETSmean, ETSthres, siteInfo[3]),
            ClCutD_Spruce(ETSmean, ETSthres, siteInfo[3]), ClCutD_Birch(ETSmean,
                ETSthres, siteInfo[3]))
    if (ClCut == 1 & is.na(initVar) & is.na(inAclct))
        inAclct <- c(ClCutA_Pine(ETSmean, ETSthres, siteInfo[3]),
            ClCutA_Spruce(ETSmean, ETSthres, siteInfo[3]), ClCutA_Birch(ETSmean,
                ETSthres, siteInfo[3]))
    if (any(is.na(inDclct)))
        inDclct[is.na(inDclct)] <- 9999999.99
    if (length(inDclct) == 1)
        inDclct <- rep(inDclct, nSp)
    if (any(is.na(inAclct)))
        inAclct[is.na(inAclct)] <- 9999999.99
    if (length(inAclct) == 1)
        inAclct <- rep(inAclct, nSp)
    if (is.na(initVar)) {
        initVar <- matrix(NA, 6, nSp)
        initVar[1, ] <- 1:nSp
        initVar[3, ] <- initClearcut[1]
        initVar[4, ] <- initClearcut[2]
        initVar[5, ] <- initClearcut[3]/nSp
        initVar[6, ] <- initClearcut[4]
    }
    xx <- min(10, nYears)
    Ainit = 6 + 2 * 3.5 - 0.005 * (sum(ETS[1:xx])/xx) + 2.25
    initVar[2, which(is.na(initVar[2, ]))] <- initClearcut[5] <- round(Ainit)
    weatherPreles <- array(c(PAR, TAir, VPD, Precip, CO2), dim = c(365,
        nYears, 5))
    weatherPreles <- aperm(weatherPreles, c(2, 1, 3))
    prebas <- .Fortran("prebas", nYears = as.integer(nYears),
        nLayers = as.integer(nLayers), nSp = as.integer(nSp),
        siteInfo = as.numeric(siteInfo), pCROBAS = as.matrix(pCROBAS),
        initVar = as.matrix(initVar), thinning = as.matrix(thinning),
        output = as.array(output), nThinning = as.integer(nThinning),
        maxYearSite = as.integer(nYears), fAPAR = as.numeric(fAPAR),
        initClearcut = as.numeric(initClearcut), ETS = as.numeric(ETS),
        P0 = as.numeric(P0), weather = as.array(weatherPreles),
        DOY = as.integer(1:365), pPRELES = as.numeric(pPRELES),
        soilC = as.array(soilC), pYASSO = as.numeric(pYASSO),
        weatherYasso = as.matrix(weatherYasso), litterSize = as.matrix(litterSize),
        soilCtot = as.numeric(soilCtot), defaultThin = as.double(defaultThin),
        ClCut = as.double(ClCut), inDclct = as.double(inDclct),
        inAclct = as.double(inAclct), yassoRun = as.double(yassoRun))
    class(prebas) <- "prebas"
    return(prebas)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
