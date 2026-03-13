!=====================================================================
!  ETfun - compute evapotranspiration components
!=====================================================================
module water_module

  use gpp_module
  use prelesglobalsF_module
 use initruns_module
 
 
  implicit none
 

    contains

   !===================================================================
  !  ETfun – returns ET (et) = loss of water from soil through 
  !          evap and transp; and also provides transpiration,
  !          evaporation, and water‑stress factors via arguments.
  !===================================================================
  function ETfun(D, theta, ppfd, fAPAR, T, ET_par, Site_par,                   &
                 canw, theta_top, fE, A, fWgpp, fOrg, GPP_par, CO2, LOGFLAG,   &
                 etmodel, transp, evap, fWE, CO2model, soilmodel, REWmodel) result(et) 
    !             bind(C, name="ETfun")
     
    ! Arguments
    real(8), intent(in)    :: D
    real(8), intent(in)    :: theta
    real(8), intent(in)    :: ppfd
    real(8), intent(in)    :: fAPAR
    real(8), intent(in)    :: T
    real(8), intent(in)    :: ET_par(5)
    real(8), intent(in)    :: Site_par(10)
    real(8), intent(in)    :: canw
    real(8), intent(in)    :: theta_top
    real(8), intent(inout)   :: fE
    real(8), intent(in)    :: A
    real(8), intent(inout)    :: fWgpp
    real(8), intent(inout)   :: fOrg
    real(8), intent(in)    :: GPP_par(13)
    real(8), intent(in)    :: CO2
    integer,  intent(in)    :: LOGFLAG
    integer,  intent(in)    :: etmodel
    real(8), intent(out)   :: transp
    real(8), intent(out)   :: evap
    real(8), intent(inout)   :: fWE
    integer, intent(in)   :: CO2model
    integer, intent(in)   :: soilmodel
    integer, intent(in)   :: REWmodel
    

    ! Local variables
    real(8) :: et
    real(8) :: thetavol, REW, fWsub
    real(8) :: thetavol_top, REW_top
    real(8) :: fCO2mean, lambda, psychom, s
    real(8), parameter :: cp      = 1003.5_8      ! J/(kg·K)
    real(8), parameter :: MWratio = 0.622_8       ! water/air molecular weight ratio
    real(8), parameter :: pressure = 101300.0_8  ! Pa (default)

    real(8) :: D_loc   ! mutable copy of D (C code may modify D)
	real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
	real(8) :: GPP_par_beta,GPP_par_tau, GPP_par_S0,GPP_par_Smax,GPP_par_kappa, GPP_par_gamma,GPP_par_soilthres
	real(8) :: GPP_par_bCO2,GPP_par_xCO2,GPP_par_t0, GPP_par_tcrit,GPP_par_tsumcrit,GPP_par_soils
	real(8) :: ET_par_beta,ET_par_kappa,ET_par_chi,ET_par_soilthres,ET_par_nu

   	include 'init_gppPar_preles.h'
	include 'init_sitePar_preles.h'
	include 'init_etPar_preles.h'

    !-----------------------------------------------------------------
    !  contents of subroutine follow
    !-----------------------------------------------------------------
    D_loc = D
    thetavol = theta / Site_par_soildepth
    REW = (thetavol - Site_par_ThetaPWP) / (Site_par_ThetaFC - Site_par_ThetaPWP)
    
    thetavol_top = theta_top / Site_par_topdepth
    REW_top = thetavol_top

    fWsub = 1.0

    !--- CO2 mean effect ------------------------------------------------
    fCO2mean = fCO2_ET_model_mean(CO2model,CO2, GPP_par)

    !--- Physical constants ---------------------------------------------
    lambda = (-0.0000614342  * T**3 + 0.00158927  * T**2 -          &
               2.36418  * T + 2500.79 ) * 1000.0                ! J/kg
    psychom = cp * pressure / (lambda * MWratio)                     ! Pa/°C
    s = 1000.0  * 4098.0  * (0.6109  * exp( (17.27  * T) / (T + 237.3 ) )) &
        / (T + 237.3 )**2                                          ! Pa/°C

    !--- Soil water constraint -----------------------------------------
    if (ET_par_soilthres < -998.0 ) then
       fWsub = 1.0 
    else
        fWsub = fREW_fun(REW, GPP_par_soils, ET_par_soilthres, REWmodel)
    end if

    !--- Top soil water constraint -----------------------------------------
    if (Site_par_orgthres < -998.0 ) then
       fOrg = 1.0 
    else
        fOrg = fREW_fun(REW_top, GPP_par_soils, Site_par_orgthres, REWmodel)
    end if
    
    !--- Canopy water overrides soil constraint ------------------------
    if (canw > 1.0e-8 ) fWsub = 1.0 

    !--- Return water‑stress factors ------------------------------------
    fE  = fWsub
    fWE = fWsub
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! for drained peatlands use top water content as the basis of evaporation
    ! otherwise use total REW
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    select case (soilmodel)
    case (1)
        fOrg = fWsub
    case (2)
        fOrg = fOrg
    case default
        fOrg = fWsub
    end select

    !--- Minimum aerodynamic conductance --------------------------------
    if (D_loc < 0.01 ) D_loc = 0.01 

    !--- Model selection ------------------------------------------------
    select case (etmodel)
    case (-1)                         ! simple Penman–Monteith style
       transp = D_loc * ET_par_beta * A / D_loc**ET_par_kappa * &
                fWgpp**ET_par_nu * fCO2mean
       evap   = ET_par_chi * (1.0  - fAPAR) * fOrg * ppfd
       et = (transp + evap) * s / (s + psychom)

    case (0)                          ! same as -1 but evap term includes s/(s+psychom)
       transp = D_loc * ET_par_beta * A / D_loc**ET_par_kappa * &
                fWgpp**ET_par_nu * fCO2mean
       evap   = ET_par_chi * s / (s + psychom) * (1.0  - fAPAR) * &
                fOrg * ppfd
       et = transp + evap

    case (1)                          ! ET not limited by psychrometric term
       transp = D_loc * ET_par_beta * A / D_loc**ET_par_kappa * &
                fWgpp**ET_par_nu * fCO2mean
       evap   = ET_par_chi * (1.0  - fAPAR) * fOrg * ppfd
       et = transp + evap

    case (2)                          ! alternative formulation
       et = D_loc * (1.0  + ET_par_beta / D_loc**ET_par_kappa) * &
            A / CO2 * fWgpp**ET_par_nu * fCO2mean + &
            ET_par_chi * (1.0  - fAPAR) * fOrg * ppfd

    case default
       et = 0.0 
    end select

                 end function ETfun
                 



  !=================================================================
  !  Interception of rainfall by the canopy
  !=================================================================
  subroutine interception_fun(rain, intercepted, Temp, SnowRain_par, fAPAR)
    real(kind=8), intent(inout) :: rain        ! Incoming rainfall (mm)
    real(kind=8), intent(out)   :: intercepted ! Rain intercepted (mm)
    real(kind=8), intent(in)    :: Temp        ! Air temperature (°C)
    real(kind=8), intent(in)   :: SnowRain_par(5)
    real(kind=8), intent(in)    :: fAPAR       ! Fraction of absorbed PAR
	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0

!read parameters
	include 'init_SnowRainPar_preles.h'

    if (Temp > SnowRain_par_SnowThreshold) then
       intercepted = rain * (SnowRain_par_I0 * fAPAR / 0.75)
       rain        = rain - intercepted
    else
       intercepted = 0.0
    end if
  end subroutine interception_fun


  !=================================================================
  !  Soil water balance update mineral soils
  !=================================================================
  subroutine sw_balance(theta, throughfall, snowmelt, et, Site_par, &
                        drainage, snow, canw, SnowRain_par)
    real(kind=8), intent(inout) :: theta        ! Soil moisture (mm)
    real(kind=8), intent(in)    :: throughfall  ! Canopy throughfall (mm)
    real(kind=8), intent(in)    :: snowmelt     ! Snow melt (mm)
    real(kind=8), intent(inout) :: et           ! Evapotranspiration (mm)
    real(8), intent(in)    :: Site_par(10)
    real(kind=8), intent(out)   :: drainage     ! Drainage (mm)
    real(kind=8), intent(inout) :: snow         ! Snow water equivalent (mm)
    real(kind=8), intent(inout) :: canw         ! Canopy water (mm)
    real(kind=8), intent(in)   :: SnowRain_par(5)

    real(kind=8) :: st0
    real(kind=8) :: et_from_veg_and_soil
	real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0

!read parameters
	include 'init_sitePar_preles.h'
	include 'init_SnowRainPar_preles.h'

    et_from_veg_and_soil = 0.0 

    !--------------------------------------------------------------
    !  Evaporation first from wet canopy and then from snow
    !--------------------------------------------------------------
    if (SnowRain_par_CWmax > 1.0e-8 ) then
       if ( (canw + snow - et) > 0.0  ) then
          if ( (canw - et) > 0.0  ) then
             canw = canw - et
             et_from_veg_and_soil = 0.0 
          else if (canw - et < 0.0 ) then
             snow = snow + canw - et
             canw = 0.0 
             et_from_veg_and_soil = 0.0 
          end if
       else
          et_from_veg_and_soil = et - canw - snow
          canw = 0.0 
          snow = 0.0 
       end if
    else
       if ( (snow - et) > 0.0  ) then
          snow = snow - et
          et_from_veg_and_soil = 0.0 
       else if (snow - et < 0.0 ) then
          et_from_veg_and_soil = et - snow
          snow = 0.0 
       else
          snow = 0.0 
       end if
    end if

    et = et_from_veg_and_soil

    !--------------------------------------------------------------
    !  Water balance without drainage
    !--------------------------------------------------------------
    st0 = theta + throughfall + snowmelt - et
    if (st0 <= 0.0 ) st0 = 1.0e-4 

    !--------------------------------------------------------------
    !  Drainage (simple time‑delay model)
    !--------------------------------------------------------------
    if (Site_par_tauDrainage > 0.0 ) then
       if (st0 > Site_par_ThetaFC * Site_par_soildepth) then
          drainage = (st0 - Site_par_ThetaFC * Site_par_soildepth) / &
                     Site_par_tauDrainage
       else
          drainage = 0.0 
       end if
       theta = st0 - drainage
    else
       drainage = 0.0 
       theta    = st0
    end if

end subroutine sw_balance

   
  !=================================================================
  !  Soil water balance update drained peatlands
  !=================================================================
  subroutine swbal_peat(theta, theta_top, pond,  throughfall, snowmelt, transp, evap,  &
                        et, Site_par, Qdrain, snow, canw, SnowRain_par)
    real(kind=8), intent(inout) :: theta        ! Soil moisture (mm)
    real(kind=8), intent(inout) :: theta_top    ! Soil moisture top layer (mm)
    real(kind=8), intent(inout) :: pond         ! water accumulating on peat surface (mm)
    real(kind=8), intent(in)    :: throughfall  ! Canopy throughfall (mm)
    real(kind=8), intent(in)    :: snowmelt     ! Snow melt (mm)
    real(kind=8), intent(inout) :: evap         ! Evaporation (mm)
    real(kind=8), intent(inout) :: transp       ! Transpiration (mm)
    real(kind=8), intent(inout) :: et           ! Evaporanspiration (mm)
    real(kind=8), intent(in)        :: Site_par(10)
    real(kind=8), intent(in)    :: Qdrain       ! Drainage in ditches (mm)
    real(kind=8), intent(inout) :: snow         ! Snow water equivalent (mm)
    real(kind=8), intent(inout) :: canw         ! Canopy water (mm)
    real(kind=8), intent(in)        :: SnowRain_par(5)

    real(kind=8) :: st0, st1, st2               ! temp variables for water storage (mm)
    real(kind=8) :: et_from_veg_and_soil        ! evap after subtracting canopy and snow evap (mm/day)
    real(kind=8) :: Dr_pond                     ! drainage from pond (mm/day)
    real(kind=8) :: Interc_top                  ! water intercepted in top layer (mm/day)
    real(kind=8) :: Dr_top                      ! drainage from top layer (mm/day)
    real(kind=8) :: Exf_SW                      ! max exfiltration from soil pool to upper layers (mm/day)
    real(kind=8) :: Exf_top                     ! exfiltration to top layer (mm/day)
    real(kind=8) :: Exf_pond                    ! exfiltration to pond (mm/day)
    real(kind=8) :: SurfRunoff                  ! surface runoff: remainder of exfiltration (mm/day)
    real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0

!read parameters
	include 'init_sitePar_preles.h'
	include 'init_SnowRainPar_preles.h'

    et_from_veg_and_soil = 0.0 

    !--------------------------------------------------------------
    !  Evaporation first from wet canopy and then from snow
    !--------------------------------------------------------------
    if (SnowRain_par_CWmax > 1.0e-8 ) then
       if ( (canw + snow - evap) > 0.0  ) then
          if ( (canw - evap) > 0.0  ) then
             canw = canw - evap
             et_from_veg_and_soil = 0.0 
          else if (canw - et < 0.0 ) then
             snow = snow + canw - evap
             canw = 0.0 
             et_from_veg_and_soil = 0.0 
          end if
       else
          et_from_veg_and_soil = evap - canw - snow
          canw = 0.0 
          snow = 0.0 
       end if
    else
       if ( (snow - evap) > 0.0  ) then
          snow = snow - evap
          et_from_veg_and_soil = 0.0 
       else if (snow - et < 0.0 ) then
          et_from_veg_and_soil = evap - snow
          snow = 0.0 
       else
          snow = 0.0 
       end if
    end if

    et = et_from_veg_and_soil + transp

    !-----------------------------------------------------------------
    !  Calculate flux rates in soil including pond, top, and main pool
    !  uses Daily time step
    !-----------------------------------------------------------------
    
    Dr_pond = pond + throughfall + snowmelt
    Interc_top = (Site_par_topdepth - theta_top) * (1 - EXP(-Dr_pond/Site_par_topdepth))
    Dr_top = Dr_pond - Interc_top
    Exf_SW = MAX(theta + (Dr_top - transp - Qdrain) - Site_par_peatdepth, 0.0)
    Exf_top = MIN(Exf_SW, (Site_par_topdepth - theta_top))
    Exf_pond = MIN(Exf_SW - Exf_top , Site_par_MaxPond - pond)
    SurfRunoff = Exf_SW - Exf_top - Exf_pond
    
    
    !--------------------------------------------------------------
    !  Update water pools
    !--------------------------------------------------------------
   
     st0 = pond + throughfall + snowmelt + Exf_pond - Dr_pond
     st1 = theta_top + Dr_pond + Exf_top - evap - Dr_top
     st2 = theta + Dr_top - transp - Qdrain -  Exf_top - Exf_pond
     
     if (st0 <= 0.0 ) st0 = 1.0e-4 
     if (st1 <= 0.0 ) st1 = 1.0e-4 
     if (st2 <= 0.0 ) st2 = 1.0e-4 
     
     pond = st0
     theta_top = st1
     theta = st2

    end subroutine swbal_peat   
                        
                        
  !=================================================================
  !  Snow accumulation and melt
  !=================================================================
  subroutine snow_process(T, rain, snow, SnowRain_par, SnowMelt)
    real(kind=8), intent(in)    :: T            ! Air temperature (°C)
    real(kind=8), intent(inout) :: rain         ! Rainfall (mm)
    real(kind=8), intent(inout) :: snow         ! Snow water equivalent (mm)
    real(kind=8), intent(in)   :: SnowRain_par(5)
    real(kind=8), intent(out)   :: SnowMelt     ! Snow melt (mm)

    real(kind=8) :: new_snow
	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0

!read parameters
	include 'init_SnowRainPar_preles.h'


    if (T < SnowRain_par_SnowThreshold) then
       new_snow = rain
       rain     = 0.0 
    else
       new_snow = 0.0 
    end if

    if (T > SnowRain_par_T_0) then
       SnowMelt = SnowRain_par_MeltCoef * (T - SnowRain_par_T_0)
    else
       SnowMelt = 0.0 
    end if

    if (snow + new_snow - SnowMelt < 0.0 ) then
       SnowMelt = new_snow + snow
       snow     = 0.0 
    else
       snow = snow + new_snow - SnowMelt
    end if
  end subroutine snow_process
  
  !=================================================================
  ! Water retention curve
  !=================================================================
  function thetaFun(psi, Genuchten_par, layerswitch) result(theta)
    real(kind=8), intent(in)    :: psi           ! suction pressure (cm of water);
    real(kind=8), intent(in)    :: Genuchten_par(8) ! eight parameters (two layers)
    integer, intent(in)         :: layerswitch
 
    real(kind=8)                :: theta         ! water retention curve [m3m−3]; 
    real(kind=8)                :: thetaR        ! saturated water content [m3m−3];
    real(kind=8)                :: thetaS        ! residual water content [m3m−3];
    real(kind=8)                :: alpha         ! parameter related to the inverse of the air entry suction,  (cm−1); 
    real(kind=8)                :: n             ! measure of the pore-size distribution,  (dimensionless)
    real(kind=8)                :: m             ! derived parameter (dimensionless)
    real(8) :: Genuchten_par_thetaS,Genuchten_par_thetaR,Genuchten_par_alpha,Genuchten_par_n
	real(8) :: Genuchten_par_thetaST,Genuchten_par_thetaRT,Genuchten_par_alphaT,Genuchten_par_nT

!read parameters
	include 'init_GenuchtenPar_preles.h'

    select case (layerswitch)
    case (1)  ! top
    thetaR = Genuchten_par_thetaRT
    thetaS = Genuchten_par_thetaST
    alpha  = Genuchten_par_alphaT
    n      = Genuchten_par_nT
    
    case (2)  ! lower part
    thetaR = Genuchten_par_thetaR
    thetaS = Genuchten_par_thetaS
    alpha  = Genuchten_par_alpha
    n      = Genuchten_par_n
    
    case default
    thetaR = Genuchten_par_thetaR
    thetaS = Genuchten_par_thetaS
    alpha  = Genuchten_par_alpha
    n      = Genuchten_par_n
    
    end select 
    
    
    if(n>0.0) then
        m = 1 - 1/n
    else
        m = 1
    end if
        
    
    if((ABS(alpha*psi))>0.0) then
         theta = thetaR + (thetaS - thetaR) / (1 + (ABS(alpha*psi))**n)**m
    else
         theta = 0.00
    end if

    
  end function thetaFun
  
    
  !=================================================================
  ! Calculate tables for estimating water table depth
  !=================================================================
  subroutine waterTable(Genuchten_par, Site_par, Ksat, dimTable, SWTable) 
    integer, intent(in)         :: dimTable      ! water retention curve [m3m−3]; 
    real(kind=8), intent(in)        :: Genuchten_par(8) ! eight parameters 
    real(kind=8), intent(in)        :: Site_par(10)      ! water retention parameters for site
    real(kind=8), intent(in)    :: Ksat(12)     ! saturated conductivity for peat profile, Ksat_par%ksat(dimTable) (m/d)
    real(kind=8), intent(out)   :: SWTable(dimTable+1, dimTable + 3)            ! water retention table
    real(kind=8)                :: gwl           ! water table depth, running (m);
    real(kind=8)                :: thetaS        ! saturated water content [m3m−3];
    real(kind=8)                :: thetaST       ! saturated water content in topsoil [m3m−3];
    real(kind=8)                :: z             ! soil depth (cm)
    real(kind=8)                :: half1         ! half of peat slice (cm)
    real(kind=8)                :: half2         ! half of peat slice below ditch depth (cm)
    real(kind=8)                :: soilD         ! depth of peat soil considered (mm)
    real(kind=8)                :: ditchD        ! ditch depth (m)
    integer                     :: NToDitch      ! nr of slices to ditch bottom
    integer                     :: NToBot        ! nr of slices to soil bottom
    real(kind=8)                :: slice1        ! thickness of slice above ditch bottom (cm)
    real(kind=8)                :: slice2        ! thickness of slice below ditch bottom (cm)
    real(kind=8)                :: sum_var
    integer                     :: kk
    integer                     :: jj
    real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
	real(8) :: Genuchten_par_thetaS,Genuchten_par_thetaR,Genuchten_par_alpha,Genuchten_par_n
	real(8) :: Genuchten_par_thetaST,Genuchten_par_thetaRT,Genuchten_par_alphaT,Genuchten_par_nT

!read parameters
	include 'init_sitePar_preles.h'
	include 'init_GenuchtenPar_preles.h'
   
    ditchD = Site_par_ditchdepth * 100   ! convert from m to cm
    soilD = Site_par_peatdepth / 10.     ! convert from mm to cm
    
    ! dimTable must be divisible by 3
    
    
    NToDitch = int(2.0 * dimTable/3.0 + 0.1)     ! added 0.1 to make sure that value truncated to correct integer
    NToBot   = int(dimTable/3.0 + 0.1)
    
    slice1 = ditchD / NToDitch   
    slice2 = (soilD - ditchD ) / NToBot
    half1 = slice1/2.0
    half2 = slice2/2.0
    
    ! fill table with saturated values to begin with
    do kk = 1, dimTable + 1
        SWTable(kk,1) = Genuchten_par_thetaST
        do jj = 2, dimTable 
             SWTable(kk,jj) = Genuchten_par_thetaS
        end do
    end do
    
    ! fill table for non-saturated cases, first above ditch bottom
    
    do kk = 2 , NToDitch +1
        gwl = - (kk-1) * slice1
        do jj = 1 , kk-1
            z = - jj * slice1
            if (jj < 2) then
                SWTable(kk,jj) = thetaFun(gwl - z - half1, Genuchten_par,1 )
            else
                SWTable(kk,jj) = thetaFun(gwl - z - half1, Genuchten_par,2 )
            endif
        end do
    end do
    
    ! fill table for non-saturated cases, now below ditch bottom
    
     do kk = 1 , NToBot                                       
        gwl = - kk * slice2 - ditchD
        do jj = 1 , NToDitch 
            z = - jj * slice1 
            if(jj < 2) then
                SWTable(NToDitch+kk+1, jj) = thetaFun(gwl - z - half1, Genuchten_par,1 )
            else
                SWTable(NToDitch+kk+1, jj) = thetaFun(gwl - z - half1, Genuchten_par,2 )
            end if
        end do
        do jj = 1, kk 
            z = - ditchD - jj * slice2 
            SWTable(NToDitch+kk+1, NToDitch+jj) = thetaFun(gwl - z - half2, Genuchten_par,2 )
        end do
     end do
     
     ! Water content, insert in SWTable(kk,dimTable + 1)
     do kk = 1, dimTable + 1
         sum_var = 0
         do jj = 1,NtoDitch
             sum_var = sum_var + SWTable(kk,jj) * slice1
         end do
         do jj = NtoDitch + 1, dimTable
             sum_var = sum_var + SWTable(kk,jj) * slice2
         end do
         SWTable(kk, dimTable + 1) = sum_var * 10.  ! from cm to mm
     end do
     
     ! Ksat, insert in SWTable(kk,dimTable + 2)
     ! first insert zeros to all
     do kk = 1, dimTable + 1
         SWTable(kk,dimTable + 2) = 0.
     end do
     ! start insertinf Ksat_av
     do kk = 1, NtoDitch 
          sum_var = 0
          do jj =  kk , NtoDitch
             sum_var = sum_var + Ksat(jj) 
          end do
          if (NtoDitch - kk + 1 > 0) then
            SWTable(kk, dimTable + 2) = sum_var / (NtoDitch - kk + 1)
          else
             SWTable(kk, dimTable + 2) = 0.
          end if   
     end do
     
     do kk = 1, NtoDitch + 1
          SWTable(kk,dimTable + 3) = (kk-1) * slice1 * 0.01 ! conversion from cm to m
     end do
     do kk = 1 , NtoBot
         SWTable(kk+ NtoDItch + 1,dimTable + 3) = (DitchD + kk * slice2) * 0.01 ! conversion from cm to m
     end do
     
         
         
     
    
  end subroutine waterTable
  
   !--------------------------------------------------------------
      ! Peatland soils: Calculate water level and update drainage 
      !--------------------------------------------------------------
      subroutine water_level(theta, SWTable, gwl, Qdrain, theta_root, Site_par, dimTable)
      
    real(kind=8), intent(in)        :: Site_par(10)      ! For ditch depth and distance
    integer, intent(in)         :: dimTable      ! soil water table dimension variable 
    real(kind=8), intent(in)    :: SWTable(dimTable+1, dimTable + 3)            ! water retention table
    real(kind=8), intent(inout) :: gwl           ! water table depth from previous timestep to be updated (m)
    real(kind=8), intent(in)    :: theta         ! soil water
    real(kind=8), intent(inout) :: Qdrain        ! ditch drainage
    real(kind=8), intent(inout) :: theta_root    ! volumetric water in root zone (mm/mm)
    
    real(kind=8)                :: ditchD        ! ditch depth (m)
    real(kind=8)                :: ditchL        ! distance between ditches (m)
    real(kind=8)                :: slice         ! thickness of peat slice (m)
    real(kind=8)                :: sum           ! for calculating theta_root
    integer                     :: kk            ! index
    integer                     :: k0            ! index
    integer                     :: NtoRoot       ! nr of slices occupied by roots 
    real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth

!read parameters
	include 'init_sitePar_preles.h'
   
    
    ditchD = Site_par_ditchdepth   ! m
    ditchL = Site_par_ditchDist        ! m
    slice = SWTable(2,dimTable + 3)
    NtoRoot = ceiling(Site_par_soildepth * 0.001 / slice) ! nr of slices occupied by roots
        
    k0 = 0    
    
    if(theta >= SWTable(1,dimTable + 1)) then
        gwl = 0.
        k0 = 1
      else          
       do kk = 2, dimTable+1
                     
       if(theta > SWTable(kk,dimTable+1)) then
           if(k0 < 1) then
               k0 = kk
           else
               k0 = k0
           endif
       endif
       
      end do
    endif
      
       
       if(theta < SWTable(dimTable+1,dimTable+1) )then
        k0 = dimTable + 1
      endif
      
      ! calculate drainage
      if(k0 > 1) then
          gwl = (SWTable(k0,dimTable+3) + SWTable(k0 - 1 , dimTable+3))*0.5
      end if
           
      Qdrain = 345600000.0 * SWTable(k0,dimTable+2) * ((-gwl + ditchD)/ditchL)**2  ! Qdrain is in mm / d
      
      ! calculate theta in root zone
      sum = 0.
      do kk = 1,NtoRoot
          sum = sum + SWTable(k0,kk)
      end do
      theta_root = sum  * slice * 1000.   ! theta in mm
          
      
end subroutine water_level

end module water_module
