\name{InitMultiSite}

\alias{InitMultiSite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialization of multiPrebas forest model
}
\description{
This function runs InitMultiSite model, which initializes the needed variables for function multiPrebas().
}
\usage{
InitMultiSite(nYearsMS, pCROBAS = pCROB, pHcMod = pHcM, pPRELES = pPREL, etmodel = 0, pYASSO =pYAS, pAWEN = parsAWEN, siteInfo = NA, multiInitVar = NA, multiThin = NA, multiNthin = NA, multiInitClearCut = NA, fixBAinitClarcut = 1, initCLcutRatio = NA, areas = NA, PAR, TAir, VPD, Precip, CO2, multiP0=NA, soilC = NA, weatherYasso = NA, litterSize = litterSizeDef, soilCtot = NA, defaultThin = 1, ClCut = 1, energyCut = 0, inDclct = NA, inAclct = NA, yassoRun = 0, smoothP0 = 1, smoothETS = 1, smoothYear=5, HcModV= 2, tapioPars=pTapio, thdPer = NA, limPer = NA, ftTapioPar = ftTapio, tTapioPar = tTapio, GVrun = 1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{nYearsMS}{
  (1 x nSites vector) Number of years to run the model in each site.
%%     ~~Describe \code{nYears} here~~
}
  \item{pCROBAS}{
  (47 x nSpecies matrix) Matrix of parameter sets, each column corresponds to a species. Default values pCROBAS = pCROB are the parameter sets for Scots pine (Pinus sylvestris), Norway spruce (Picea abies), Silver birch (Betula pendula), European beech (Fagus sylvatica), Maritime pine (Pinus pinaster), Blue gum (Eucalyptus globulus), Black locust (Robinia pseudoacacia), Populus(in Romania), Eucalyptus grandis x Eucalyptus urophylla (in Paraguay), Norway spruce(in Germany), and holm oak (Quercus ilex). print(pCROB) to see the parameter values and names. 
%%     ~~Describe \code{pCROBAS} here~~
}

\item{pHcMod}{
(nSites x nSpecies matrix) Default values pHcMod = pHcM with columns 1: Scots pine (pisy), 2: Norway Spruce (piab), 3: Silver birch (beal), 4: European beech (fasy), 5: Maritime pine (pipi), 6: Blue gum (eugl), 7: Black locust (rops), 8: Poplar(popu), 9: Eucalyptus grandis x Eucalyptus urophylla (eugrur), 10:Norway spruce in Germany(piab(DE)). 
}
  
\item{pPRELES}{
  (1 x 30 vector) Parameter vector for PRELES model. Default values (pPRELES = pPREL) is the calibration for European Boreal forests (Minunno et al. 2016)
}

\item{etmodel}{
  (Value -1, 0, 1 or 2) Evapotranspiration model for PRELES. Default etmodel = 0. (Package Rpreles, Peltoniemi et al., 2015)
}

\item{pYASSO}{
  (1 x 35 vector) Parameter vector for YASSO model. Default pYASSO = pYAS (Liski et al., 2017)
}

\item{pAWEN}{
  (12 x nSpecies matrix) Matrix of parameter sets for partitioning tree organs (foliage, branch and woody) in AWEN pools. Columns referes to different species. Note that the columns of pAWEN must correspond to the species of pPRELES. Default pAWEN = parsAWEN are for 1: Scots pine (pisy), 2: Norway Spruce (piab), 3: Silver birch (beal), 4: European beech (fasy), 5: Maritime pine (pipi), 6: Blue gum (eugl), 7: Black locust (rops), 8: Poplar(popu), 9: Eucalyptus grandis x Eucalyptus urophylla (eugrur), 10:Norway spruce in Germany(piab(DE)),11:Quercus ilex (quil). 
}

\item{siteInfo}{%%     ~~Describe \code{siteInfo} here~~
  (nSites x 12 Matrix) Matrix of site info: SiteID, climID, siteType, SWinit (initial soil water), CWinit (initial crown water), SOGinit (initial snow on ground), Sinit (initial temperature acclimation state), nLayers, nSpecies, soildepth, effective field capacity, permanent wilthing point. Default = c(1,1,3, 160, 0, 0, 20, 413.,0.45, 0.118), i.e. siteType = 3.
  
Give for homogeneous sites as 
siteInfo <- matrix(c(NA,NA,3,160,0,0,20,nLayers,3,413,0.45,0.118),nSites,12,byrow = T)
# Set site index values to siteInfo columns 1-2
siteInfos[,2] <- siteInfos[,1] <- 1:nSites
}

\item{multiInitVar}{
  (nSite x 7 x nLayer matrix) Matrix with initial stand values for all the tree strata in each site. Third dimension corresponds to the layers in the stand. Initial information needed are: SpeciesID (a number corresponding to the species parameter values of pPRELES columns), Age (years), average height of the layer (H, m), average diameter at breast height of the layer (D, cm), basal area of the layer (BA, m2 ha-1), average height of the crown base of the layer (Hc, m). 8th column is updated automatically to Ac. 

If initVar is not provided the model is initialized from plantation using default planting parameters (see initClearcut) and assuming that Pine, Spruce and Birch are equally present at plantation.
}

\item{multiThin }{
  A array with thinnig inputs. Three dimensions: nSites x maxThin x 9. The first dimension is the number of sites. The second dimension, maxThin, is the number of thinnings. For the third demention, element 1 is year from the start of the simulation; element 2 is siteID; element 3 layer where thinnings are carried out; element 4 to 7 stand variables (H, D, B, Hc); element 8 parameter that indicates if the stand variables (column 4:7) are provided as fraction of the actual model outputs (value=1 means that fraction is used); element 9 is the stand density after thinning if its value is not -999;colum 10 is Sapwood area of average tree at crown base (m2) if its value is not -999.
}

\item{multiNthin}{?}

\item{multiInitClearCut}{?}

\item{fixBAinitClarcut}{?}

\item{initCLcutRatio}{?}

\item{areas}{(1 x nSites vector) Areas of the sites (Which unit)?}

\item{PAR}{
  (nSites x nYears*365 matrix) A numeric matrix of daily sums of photosynthetically active radiation, mmol/m2, for each site.}

\item{TAir}{
  (nSites x nYears*365 matrix) A numeric matrix of daily mean temperature, degrees C, for each site.}

\item{VPD}{
  (nSites x nYears*365 matrix) A numeric matrix of daily mean vapour pressure deficits, kPa, for each site.}
  
\item{Precip}{
  (nSites x nYears*365 matrix) A numeric matrix of daily rainfall, mm, for each site.}

\item{CO2}{
  (nSites x nYears*365 matrix) A numeric matrix of air CO2, ppm, for each site.}

\item{multiP0}{??
  (? x ? matrix/vector)  A numeric vector with the annual potential photosynthesis (gC m-2 y-1). If P0 is not provided PRELES is used to compute P0 using fAPAR = 1. (ADD REFERENCE!)
}

\item{soilC}{???
   (nSites x 5 x 3 x nLayers array) Initial soil carbon compartments for each layer. Array with dimensions = c(nYears,5,3,nLayers). The second dimension (5) corresponds to the AWENH pools; the third dimension (3) corresponds to the tree organs (foliage, branch and stem). (ADD REFERENCE?)
%%     ~~Describe \code{soilC} here~~
}

\item{weatherYasso}{
  (nSites x nYears matrix/vector) Annual weather inputs for Yasso. (Which variables, which units?)
  
  -> Can these be calculated from the other weather inputs??
%%     ~~Describe \code{weatherYasso} here~~
}

\item{litterSize}{
  (nSites x ? x ? matrix)
  Marix with litter inputs for YASSO. Rows are tree organs, columns correspond to the layers. (REFERENCE??)
%%     ~~Describe \code{litterSize} here~~
}

\item{soilCtot}{ ??
(? x ? matrix/vector/array)
  Vector with total initial soil carbon
%%     ~~Describe \code{soilCtot} here~~
}

\item{defaultThin}{
  If defaultThin = 1 (default) Finnish standard management practices are applied (REFERENCE?), otherwise if defaultThin = 0, ???. (REFERENCE)
%%     ~~Describe \code{defaultThin} here~~
}

\item{ClCut}{
  If ClCut = 1 clearcuts are applied. If inDclct = NA and inAclct = NA Finnish standard clearcut practices are applied (REFERENCE).
%%     ~~Describe \code{ClCut} here~~
}
  
\item{energyCut}{Energy cutting strategy is applied when set to 1. Default value is 0.}
  
\item{inDclct}{ (nSites x nLayer matrix??)
  Vector of Diameter (cm) threshold for clearcut. Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
%%     ~~Describe \code{inDclct} here~~
}

\item{inAclct}{(nSites x nLayer matrix??)
  Vector of Age (year) threshold for clearcut.  Each element correspond to a layer of the stand, if only one value is provided the same value is applied to all the layers. The different elements of the vector are for the different layers. The dominant species (highest basal area) is considered for clearcut.
%%     ~~Describe \code{inAclct} here~~
}

\item{yassoRun}{
  If yassoRun=1 the YASSO model is run to compute the carbon balance of the soil. Default value = 0.
%%     ~~Describe \code{yassoRun} here~~
}


\item{smoothP0}{?? Default = 1}

\item{smoothETS}{?? Default value = 1}

\item{smoothYear}{?? Default value smoothYear = 5}

\item{HcModV}{Default HcModV=2. Version of model to compute 
Hc. 1: use the version of based on ksi parameter, 2: use the empirical model. (REFERENCE??)}

\item{tapioPars}{(? x ? matrix) Thinning parameters. Default tapioPars= pTapio for conifers and decisuous trees (What is given in parameter values??) (REFERENCE?)}

\item{thdPer}{??}

\item{limPer}{??}

\item{ftTapioPar}{(? x ? matrix) parameters for the first thinning. Default ftTapioPar = ftTapio for pine, spruce and birch (betula pendula).}

\item{tTapioPar}{Default tTapioPar = tTapio}
                        
\item{GVrun}{1 if ground vegetation model is run, otherwise 0. Default GVrun = 1}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{output}{
  An array of initial values and matrices needed for function multiPrebas().

  These outputs need to be checked, they are now directly the outputs of prebas()-function!
  
  An array with annual model outputs. 1st dimension corresponds to the number of years of the simulation (nYears); 2nd dimension corresponds to the output variables (see list below); 3rd dimension corresponds to the number of layers in the stand (nLayers); 4th dimensions reports the state of the stand (1) and (2) the variables of the harvested trees (2).

  Output variables:
  1."siteID" \cr
  2."climID" \cr
  3."sitetype" \cr
  4."species" \cr
  5."ETS" effective temperature sums \cr
  6."P0" Potential annual gross primary production (gC m-2 y-1) \cr
  7."age" Age of the layer (years) \cr
  8."DeadWoodVolume" Dead wood volume (m3 ha-1) \cr
  9."Respi_tot" Autotrophic respiration (gC m-2 y-1) \cr
  10."GPP/1000" Total GPP  (kgC m-2 y-1) \cr
  11."H" Layer average height (m) \cr
  12."D" Layer average diameter at breast height (cm) \cr
  13."BA" Layer basal area (m-2 ha-1) \cr
  14."Hc_base" Base of crown height (m) \cr
  15."Cw" Crorn width (m) \cr
  16."A" Sapwood area of average tree at crown base (m2) \cr
  17."N" Layer density \cr
  18."npp" net primary production (gC m-2 y-1) \cr
  19."leff" Effective leaf area \cr
  20."keff" Effective light extintion coefficient \cr
  21."lproj" Projected leaf area \cr
  22."ET_preles" Annual evapotranspiration (mm y-1) \cr
  23."weight" Layer weight on photosynthesis \cr
  24."Wbranch" Branch biomass (kgC ha-1) \cr
  25."WfineRoots" Fine roots biomass (kgC ha-1) \cr
  26."Litter_fol" Foliage litter (kgC ha-1) \cr
  27."Litter_fr" Fine root litter (kgC ha-1) \cr
  28."Litter_branch" Branch litter (kgC ha-1) \cr
  29."Litter_wood" Woody litter (kgC ha-1) \cr
  30."V" Layer volume (m3 ha-1) \cr
  31."Wstem" Stem Biomass (kgC ha-1) \cr
  32."W_croot" Course root Biomass (kgC ha-1) \cr
  33."wf_STKG" Foliage biomass (kgC ha-1) \cr
  34."wf_treeKG" Foliage biomass of the average tree (kgC ha-1) \cr
  35."B_tree" Basal area of average tree (m2) \cr
  36."Light" light interseption \cr
  37."Vharvested" harvested volume (m3 ha-1) \cr
  38."Wharvested" haversted biomass (kgC ha-1) \cr
  39."soilC" totaal soil carbon (kgC ha-1) \cr
  40."aSW" average available soil water (mm) \cr
  41."summerSW" summer soil water (mm) \cr
  42."Vmort" volume of dead trees (m3 ha-1) \cr
  43."gross growth" (m3 ha-1 y-1) \cr
  44."GPPspecies" Gross primary production per layer (gC m-2 y-1) \cr
  45."Rh species" (gC m-2 y-1) \cr
  46."NEP sp" Net ecosystem exchange (gC m-2 y-1)}
  \item{dailyPRELES}{
    Matrix with daily output from preles. 1st column is the daily GPP (gC m-2 d-1), 2nd column daily evapotranspiration (mm m-2), 3rd column daily soil water (mm -2).
  }%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Minunno et al. 2016

Liski et al., 2017
%% ~put references to the literature/web site here ~
}
\author{
Virpi Junttila
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

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
