#' R wrapper for Fortran subroutine `preles_fortran`
#'
#' @param PAR,TAir,VPD,Precip,CO2,fAPAR Numeric vectors of length N (double).
#' @param Site_par Numeric length 10.
#' @param GPP_par Numeric length 13.
#' @param ET_par Numeric length 5.
#' @param SnowRain_par Numeric length 5.
#' @param Genuchten_par Numeric length 8.
#' @param Ksat Numeric length 15.
#' @param WL Numeric vector length N (inout); initial groundwater level or state.
#' @param SW,ST,SR,SOG,S Numeric vectors length N (inout); initial soil/snow/canopy states.
#' @param etmodel Integer scalar (model choice).
#' @param dimTable Integer scalar; lookup table dimension.
#' @param LOGFLAG,multisiteNday Integer scalars.
#' @param day Integer vector length N (inout). Defaults to 1:N if NULL.
#' @param CO2model,soilmodel,REWmodel Integer scalars (model switches).
#' @param lib Path to compiled shared library (e.g., 'preles.so' or 'preles.dll').
#' @param load_library Logical: auto dyn.load the shared library if not already loaded.
#' @return A list with all outputs (vectors of length N) and updated inout states.
#' @examples
#' # See example section below for a working call template
preles_r <- function(
    PAR, TAir, VPD, Precip, CO2, fAPAR,
    Site_par, GPP_par, ET_par, SnowRain_par,
    Genuchten_par, Ksat,
    WL,
    dimTable,
    # inout state vectors (provide initial conditions; default to zeros)
    SW = NULL, ST = NULL, SR = NULL, SOG = NULL, S = NULL,
    # switches / flags
    etmodel = 1L, LOGFLAG = 0L, multisiteNday = 0L,
    day = NULL, CO2model = 1L, soilmodel = 1L, REWmodel = 1L
    # shared library handling
    # lib = "preles_fortran.so", load_library = TRUE
) {
  # ---------- Basic checks ----------
  N <- length(PAR)
  stopifnot(
    length(TAir) == N, length(VPD) == N, length(Precip) == N,
    length(CO2) == N, length(fAPAR) == N,
    length(WL) == N
  )
  # Required parameter lengths
  stopifnot(
    length(Site_par) == 10L,
    length(GPP_par) == 13L,
    length(ET_par) == 5L,
    length(SnowRain_par) == 5L,
    length(Genuchten_par) == 8L,
    length(Ksat) == 15L
  )
  
  # Initialize inout state vectors if not supplied
  if (is.null(SW))  SW  <- rep(0, N)
  if (is.null(ST))  ST  <- rep(0, N)
  if (is.null(SR))  SR  <- rep(0, N)
  if (is.null(SOG)) SOG <- rep(0, N)
  if (is.null(S))   S   <- rep(0, N)
  
  stopifnot(length(SW) == N, length(ST) == N, length(SR) == N,
            length(SOG) == N, length(S) == N)
  
  if (is.null(day)) day <- seq_len(N)
  
  # Fortran cannot handle NA/NaN gracefully; warn if present
  has_na <- function(x) any(!is.finite(x))
  if (any(vapply(list(PAR, TAir, VPD, Precip, CO2, fAPAR,
                      Site_par, GPP_par, ET_par, SnowRain_par,
                      Genuchten_par, Ksat, WL, SW, ST, SR, SOG, S),
                 has_na, logical(1)))) {
    warning("Some inputs contain NA/NaN/Inf; Fortran code may fail or propagate invalid values.")
  }
  
  # ---------- Load shared library if requested ----------
  # if (load_library) {
  #   if (!is.loaded("preles_fortran")) {
  #     dyn.load(lib)
  #     on.exit({
  #       try(dyn.unload(lib), silent = TRUE)
  #     }, add = TRUE)
  #   }
  # }
  
  # ---------- Coerce storage modes ----------
  NofDays       <- as.integer(N)
  PAR           <- as.double(PAR)
  TAir          <- as.double(TAir)
  VPD           <- as.double(VPD)
  Precip        <- as.double(Precip)
  CO2           <- as.double(CO2)
  fAPAR         <- as.double(fAPAR)
  
  Site_par      <- as.double(Site_par)
  GPP_par       <- as.double(GPP_par)
  ET_par        <- as.double(ET_par)
  SnowRain_par  <- as.double(SnowRain_par)
  Genuchten_par <- as.double(Genuchten_par)
  Ksat          <- as.double(Ksat)
  
  etmodel       <- as.integer(etmodel)
  dimTable      <- as.integer(dimTable)
  
  WL            <- as.double(WL)
  SW            <- as.double(SW)
  ST            <- as.double(ST)
  SR            <- as.double(SR)
  SOG           <- as.double(SOG)
  S             <- as.double(S)
  
  LOGFLAG       <- as.integer(LOGFLAG)
  multisiteNday <- as.integer(multisiteNday)
  day           <- as.integer(day)
  CO2model      <- as.integer(CO2model)
  soilmodel     <- as.integer(soilmodel)
  REWmodel      <- as.integer(REWmodel)
  
  # ---------- Allocate outputs (all length N) ----------
  GPP           <- double(N)
  ET            <- double(N)
  fS            <- double(N)
  fD            <- double(N)
  fW            <- double(N)
  fE            <- double(N)
  fL            <- double(N)
  Throughfall   <- double(N)
  Interception  <- double(N)
  Snowmelt      <- double(N)
  Drainage      <- double(N)
  Canopywater   <- double(N)
  GPPmeas       <- double(N)
  ETmeas        <- double(N)
  SWmeas        <- double(N)
  transp        <- double(N)
  evap          <- double(N)
  fWE           <- double(N)
  fOrg          <- double(N)
  
  # ---------- Call Fortran subroutine ----------
  # IMPORTANT: the argument order MUST exactly match the Fortran signature.
  out <- .Fortran(
    "preles_fortran",
    NofDays      = NofDays,
    PAR          = PAR,
    TAir         = TAir,
    VPD          = VPD,
    Precip       = Precip,
    CO2          = CO2,
    fAPAR        = fAPAR,
    Site_par     = Site_par,
    GPP_par      = GPP_par,
    ET_par       = ET_par,
    SnowRain_par = SnowRain_par,
    etmodel      = etmodel,
    Genuchten_par = Genuchten_par,
    Ksat         = Ksat,
    WL           = WL,
    dimTable     = dimTable,
    GPP          = GPP,
    ET           = ET,
    SW           = SW,
    ST           = ST,
    SR           = SR,
    SOG          = SOG,
    fL           = fL,
    fS           = fS,
    fD           = fD,
    fW           = fW,
    fE           = fE,
    Throughfall  = Throughfall,
    Interception = Interception,
    Snowmelt     = Snowmelt,
    Drainage     = Drainage,
    Canopywater  = Canopywater,
    GPPmeas      = GPPmeas,
    ETmeas       = ETmeas,
    SWmeas       = SWmeas,
    S            = S,
    LOGFLAG      = LOGFLAG,
    multisiteNday = multisiteNday,
    day          = day,
    transp       = transp,
    evap         = evap,
    fWE          = fWE,
    fOrg         = fOrg,
    CO2model     = CO2model,
    soilmodel    = soilmodel,
    REWmodel     = REWmodel#,
    # PACKAGE      = "Rprebasso"  # set to your package name if you embed in a package
  )
  
  # ---------- Collect results ----------
  # `.Fortran` returns a named list of all arguments after execution.
  res <- list(
    # Outputs
    GPP          = out$GPP,
    ET           = out$ET,
    SW           = out$SW,
    ST           = out$ST,
    SR           = out$SR,
    SOG          = out$SOG,
    fL           = out$fL,
    fS           = out$fS,
    fD           = out$fD,
    fW           = out$fW,
    fE           = out$fE,
    Throughfall  = out$Throughfall,
    Interception = out$Interception,
    Snowmelt     = out$Snowmelt,
    Drainage     = out$Drainage,
    Canopywater  = out$Canopywater,
    GPPmeas      = out$GPPmeas,
    ETmeas       = out$ETmeas,
    SWmeas       = out$SWmeas,
    transp       = out$transp,
    evap         = out$evap,
    fWE          = out$fWE,
    fOrg         = out$fOrg,
    day          = out$day,
    # Inout (updated states)
    WL           = out$WL,
    S            = out$S
  )
  
  res
}
