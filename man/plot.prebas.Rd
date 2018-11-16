\name{plot.prebas}
\alias{plot.prebas}

\title{
Plot PREBAS outputs
}
\description{
This function plots the outputs of PREBAS model
}
\usage{
plot.prebas(x, variableIDs=NA, siteIDs=NA, layerIDs=NA, leg=T, layerNam = NA, obsData=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a prebas or multiSitePrebas object.
}
  \item{variableIDs}{
 a numeric vector indicating the output variables IDs to plot. See list below.
}
  \item{siteIDs}{
a numeric vector indicating the site IDs (only for a multiSitePrebas object).
}
  \item{layerIDs}{
a numeric vector indicating the layer IDs to plot.
}
  \item{leg}{
logical. If TRUE the legend is included in the plots.
}
  \item{layerNam}{
a character vector with the layer names.
}
  \item{obsData}{
a matrix with observed data. Rows are observations; Columns are 7: site IDs, layer IDs,	the year of observation from the start of the simulations, variable IDs, the actual value observed,	sd of the observation (if available).
}
}
\details{
List of output variables:

  1."siteID" \cr
  2."climID" \cr
  3."sitetype" \cr
  4."species" \cr
  5."ETS" effective temperature sums \cr
  6."P0" Potential annual gross primary production (gC m-2 y-1) \cr
  7."age" Age of the layer (years) \cr
  8."Respi_m /10000." (kgC m-2 y-1) \cr
  9."Respi_tot" Autotrophic respiration (gC m-2 y-1) \cr
  10."GPP/1000" Total GPP  (kgC ha-1 y-1) \cr
  11."H" Layer average height (m) \cr
  12."D" Layer average diameter at breast height (cm) \cr
  13."BA" Layer basal area (m-2 ha-1) \cr
  14."Hc_base" Base of crown height (m) \cr
  15."Cw" Crorn width (m) \cr
  16."Lc" Length of the crown (m) \cr
  17."N" Layer density \cr
  18."npp" net primary production (gC m-2 y-1) \cr
  19."leff" Effective leaf area \cr
  20."keff" Effective light extintion coefficient \cr
  21."lproj" Projected leaf area \cr
  22."ET_preles" Annual evapotranspiration (mm m-2) \cr
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
  35."B_tree" Basal area of average tree (m2 ha-1) \cr
  36."Light" light interseption \cr
  37."Vharvested" harvested volume (m3 ha-1) \cr
  38."Vtot" total volume of the Layer considering also dead trees and harvested volume (m3 ha-1) \cr
  39."soilC" totaal soil carbon (kgC ha-1) \cr
  40."aSW" average available soil water (mm m-2) \cr
  41."summerSW" summer soil water (mm m-2) \cr
  42."Vmort" volume of dead trees (m3 ha-1) \cr
  43."gross growth" (m3 ha-1 y-1) \cr
  44."GPPspecies" Gross primary production per layer (gC m-2 y-1) \cr
  45."Rh species" (gC m-2 y-1) \cr
  46."NEP sp" Net ecosystem exchange (gC m-2 y-1)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
function (x, variables = NA, sites = NA, speciesIDs = NA, leg = T,
    speciesNam = NA, data = NA)
{
    varNam <- getVarNam()
    if (any(variables == "all") | anyNA(variables))
        variables <- c(5, 6, 8:18, 22, 24:34, 37:46)
    if (inherits(x, "prebas")) {
        if (anyNA(speciesIDs))
            speciesIDs <- 1:dim(x$output)[3]
        nSp <- length(speciesIDs)
        if (anyNA(speciesNam))
            speciesNam <- as.character(paste("sp", 1:nSp))
        count <- 0
        if (length(variables) > 1)
            par(mfrow = c(2, 3))
        else par(mfrow = c(1, 1))
        for (var in variables) {
            count <- count + 1
            plot(x$output[, var, speciesIDs[1], 1], type = "l",
                xaxt = "n", main = varNam[var], ylab = "units",
                xlab = "age (y)", col = speciesIDs[1], ylim = c(min(x$output[,
                  var, , 1]), max(x$output[, var, , 1])))
            if (nSp > 1)
                for (ij in speciesIDs[2:nSp]) lines(x$output[,
                  var, ij, 1], col = ij)
            axis(1, at = seq(1, (dim(x$output)[1]), length.out = 6),
                labels = x$output[seq(1, (dim(x$output)[1]),
                  length.out = 6), 7, 1, 1])
            if (leg == TRUE)
                legend("topleft", c(speciesNam[speciesIDs]),
                  lty = 1, col = 1:nSp)
            if (count\%\%6 == 0 & var != tail(variables, n = 1))
                pause()
        }
    }
    if (inherits(x, "multiPrebas")) {
        if (anyNA(speciesIDs))
            speciesIDs <- 1:dim(x$multiOut)[4]
        nSp <- length(speciesIDs)
        if (anyNA(speciesNam))
            speciesNam <- as.character(paste("sp", 1:nSp))
        if (anyNA(sites))
            sites <- 1:dim(x$multiOut)[1]
        for (iz in sites) {
            count <- 0
            if (length(variables) > 1)
                par(mfrow = c(2, 3))
            else par(mfrow = c(1, 1))
            for (var in variables) {
                plot(x$multiOut[iz, , var, speciesIDs[1], 1],
                  type = "l", xaxt = "n", main = varNam[var],
                  ylab = "units", xlab = "age (y)", col = speciesIDs[1],
                  ylim = c(min(x$multiOut[iz, , var, , 1]), max(x$multiOut[iz,
                    , var, , 1])))
                if (nSp > 1)
                  for (ij in speciesIDs[2:nSp]) lines(x$multiOut[iz,
                    , var, ij, 1], col = ij)
                axis(1, at = seq(1, (dim(x$multiOut)[2]), length.out = 6),
                  labels = x$multiOut[iz, seq(1, (dim(x$multiOut)[2]),
                    length.out = 6), 7, 1, 1])
                if (leg == TRUE)
                  legend("topleft", c(speciesNam[speciesIDs]),
                    lty = 1, col = speciesIDs)
                if (count\%\%6 == 0)
                  title(paste("Site:", x$multiOut[iz, 1, 1, 1,
                    1]), line = -22, outer = TRUE, cex.main = 2)
                count <- count + 1
                if (count\%\%6 == 0 & var != tail(variables, n = 1))
                  pause()
            }
            if (length(sites) > 1 & iz != tail(sites, n = 1))
                pause()
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
