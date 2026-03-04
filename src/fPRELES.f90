!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PRELES MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Model function:
 !   1. Makes initialisations
 !   2. Estimates GPP
 !   3. Estimates snow and interception
 !   4. Estimates Evapotranspiration
 !   5. Updates soil water balance
    
    
  
    subroutine fpreles(NofDays, PAR, TAir, VPD, Precip, CO2, fAPAR, NpPRELES, pPRELES,etmodel, &
                  WL, dimTable, GPP, ET, SW, ST, SR, SOG, fL, fS, fD,        &                              
                   fW, fE, Throughfall, Interception, Snowmelt, Drainage, Canopywater, GPPmeas, ETmeas, SWmeas, &
                  S, LOGFLAG, multisiteNday, day, transp, evap, fWE, fOrg, CO2model, soilmodel, REWmodel)
	  
	  use gpp_module
	  use water_module
	  use prelesglobalsF_module
	  use initruns_module
    
    integer, intent(in) :: NofDays,NpPRELES
    real(8), intent(inout) :: PAR(NofDays), TAir(NofDays), VPD(NofDays), Precip(NofDays), CO2(NofDays), fAPAR(NofDays)
	real(8), intent(inout) :: pPRELES(NpPRELES)
    integer, intent(in)  :: etmodel
    real(8), intent(out) :: GPP(NofDays), ET(NofDays)
    real(8), intent(inout) :: ST(NofDays), SR(NofDays), WL(NofDays), SW(NofDays), SOG(NofDays), S(NofDays)
    real(8), intent(out) :: fS(NofDays), fD(NofDays), fW(NofDays), fE(NofDays), fL(NofDays)
    real(8), intent(out) :: Throughfall(NofDays), Interception(NofDays), Snowmelt(NofDays)
    real(8), intent(out) :: Drainage(NofDays), Canopywater(NofDays)
    real(8), intent(out) :: GPPmeas(NofDays), ETmeas(NofDays), SWmeas(NofDays)
    integer, intent(in) :: LOGFLAG
    integer, intent(in) :: multisiteNday
    integer, intent(inout) :: day(NofDays)
    real(8), intent(out) :: transp(NofDays), evap(NofDays), fWE(NofDays), fOrg(NofDays)
    integer, intent(in) ::  CO2model   ! 1 for Kolari, 2 for Launiainen, default Kolari
    integer, intent(in) ::  soilmodel   ! 1 for mineral soil, 2 for drained peatland, default mineral
    integer, intent(in) ::  REWmodel   ! 1 for smooth function, 2 for piecewise linear, default smooth
    integer, intent(in) ::  dimTable   ! dimension of look-up table
    integer :: i, kk,  jj
    real(8) :: SWTable(dimTable+1, dimTable+3)
    real(8) :: II, T, D, P, theta, theta_snow, theta_canopy, theta_top, theta_root, pond, S_state
    real(8) :: PhenoS, fPheno
    real(8) :: fEgpp, gpp380
    real(8) Qdrain,gwl
    type(p1) :: Site_par
    type(p2) :: GPP_par
    type(p3) :: ET_par
    type(p4) :: SnowRain_par
    type(p5) :: Init_par
    type(p7) :: Genuchten_par
    real(8) :: Ksat(15)

    
   ! Initialize variables
    Site_par%soildepth = pPRELES(1)
    Site_par%ThetaFC = pPRELES(2)
    Site_par%ThetaPWP = pPRELES(3)
    Site_par%tauDrainage = pPRELES(4)
    Site_par%topdepth = pPRELES(5)   !new  (mm)
    Site_par%orgthres = pPRELES(6)   !new  (unitless)
    Site_par%MaxPond = pPRELES(7)    !new  (mm)
    Site_par%ditchDepth = pPRELES(8) !new  (m)
    Site_par%ditchDist = pPRELES(9)  !new  (m)
    Site_par%peatdepth = pPRELES(10)
	
	GPP_par%beta = pPRELES(11)
    GPP_par%tau = pPRELES(12)
    GPP_par%S0 = pPRELES(13)
    GPP_par%Smax = pPRELES(14)
    GPP_par%kappa = pPRELES(15)
    GPP_par%gamma = pPRELES(16)
    GPP_par%soilthres = pPRELES(17)
    GPP_par%bCO2 = pPRELES(18)
    GPP_par%xCO2 = pPRELES(19)
    GPP_par%t0 = pPRELES(20)
    GPP_par%tcrit = pPRELES(21)
    GPP_par%tsumcrit = pPRELES(22)
    GPP_par%soils = pPRELES(23) !new  - parameter for soil threshold function

    ET_par%beta = pPRELES(24)
    ET_par%kappa = pPRELES(25)
    ET_par%chi = pPRELES(26)
    ET_par%soilthres = pPRELES(27)
    ET_par%nu = pPRELES(28)

    SnowRain_par%MeltCoef = pPRELES(29)
    SnowRain_par%I0 = pPRELES(30)
    SnowRain_par%CWmax = pPRELES(31)
    SnowRain_par%SnowThreshold = pPRELES(32)
    SnowRain_par%T_0 = pPRELES(33)

    Init_par%SW = pPRELES(34)
    Init_par%CW = pPRELES(35)
    Init_par%SOG = pPRELES(36)
    Init_par%S = pPRELES(37)
    Init_par%ST = pPRELES(38)  !new  (initial water level in top layer, mm)
    Init_par%WL = pPRELES(39)   !new (initial water table depth, mm, needs to be synchronised with water storage SW, to be consistent with SWTable)
  
    Genuchten_par%thetaS = pPRELES(40)
    Genuchten_par%thetaR = pPRELES(41)
    Genuchten_par%alpha = pPRELES(42)
    Genuchten_par%n = pPRELES(43)
    Genuchten_par%thetaST = pPRELES(44)
    Genuchten_par%thetaRT = pPRELES(45)
    Genuchten_par%alphaT = pPRELES(46)
    Genuchten_par%nT = pPRELES(47)
   
    Ksat = pPRELES(48:62)
	
    SW(1) = Init_par%SW
    Canopywater(1) = Init_par%CW
    SOG(1) = Init_par%SOG
	S(1) = Init_par%S
	ST(1) = Init_par%ST
	WL(1) = Init_par%WL / 1000.
	SR(1) = Site_par%ThetaFC * 200.


    PhenoS = 0.0d0
    fPheno = 0.0d0
    fEgpp = 0.0d0
    gpp380 = 0.0d0

    call init_conditions(PAR(1), TAir(1), VPD(1), Precip(1), CO2(1))

    theta = SW(1)
    theta_canopy = Canopywater(1)
    theta_snow = SOG(1)
    S_state = S(1)
    theta_top = ST(1)
    theta_root = SR(1)
    pond = Site_par%MaxPond
    Drainage(1) = 0.
    gwl = WL(1)
  
    
    ! calculate look-up table for estimating water table depth
    ! SWTable is output from this subroutine
    call waterTable(Genuchten_par, Site_par, Ksat, dimTable, SWTable)

    
! ---------------------------------------------------------------------
! START LOOPING DAYS
! ---------------------------------------------------------------------
do i = 1, NofDays

  ! Use previous day environment for prediction, if current values are missing,
  ! or suspicious
  if (i > 1) then
    if (PAR(i) < -900.0) PAR(i) = PAR(i-1)
    if (TAir(i) < -900.0) TAir(i) = TAir(i-1)
    if (VPD(i) < 0.0 .or. VPD(i) > 6.0) VPD(i) = VPD(i-1)
    if (Precip(i) < 0.0) Precip(i) = Precip(i-1) * 0.3
    ! On avg. P+1=0.315*P (in Sodis & Hyde)
    if (CO2(i) < 0.0) CO2(i) = CO2(i-1)
    if (GPPmeas(i) < -990.0) GPPmeas(i) = GPPmeas(i-1)
    if (ETmeas(i) < -990.0) ETmeas(i) = ETmeas(i-1)
    if (SWmeas(i) < 0.0) SWmeas(i) = SWmeas(i-1)
    if (SW(i) < -900.0) SW(i) = SW(i-1)
    if (SOG(i) < -900.0) SOG(i) = SOG(i-1)
  end if

  ! Assign current values for environment
  II = PAR(i)
  T = TAir(i)
  D = VPD(i)
  P = Precip(i)
  
        select case (soilmodel)
      case(1) 
          theta_top = theta
          theta_root = theta
          fOrg(i) = 0
      case(2)
          theta_top = theta_top
          theta_root = theta_root
          fOrg(i) = fOrg(i)
      case default
          theta_top = theta
          theta_root = theta
          fOrg(i) = 0
      end select


  ! Update temperature state that tells about seasonality -
  ! for GPP and through GPP to ET
  fS(i) = fS_model(S_state, T, GPP_par)

  ! Deciduous phenology - don't use if this information is inputted in fAPAR
  ! Note also that fapar is multiplied by 0 or 1 (i.e. leaf development is not accounted for)
  ! Model predicts budbreak based on critical threshold temperature sum
  ! Note that this implementation works only if data starts before t0-date of fPheno-model
  fPheno = fPheno_model(GPP_par, T, PhenoS, day(i), fS(i))

  fAPAR(i) = fAPAR(i) * fPheno

  call GPPfun(GPP(i), gpp380, II, D, CO2(i), theta_root, fAPAR(i), fS(i), &
              GPP_par, Site_par, fL(i), fD(i), fW(i), fEgpp, LOGFLAG, CO2model, REWmodel)

  ! Calculate amount of snow and snowmelt at the end of the day
  call snow_process(T, P, theta_snow, SnowRain_par, Snowmelt(i))

  ! NOTE: interception model could be better
  Throughfall(i) = P
  call interception_fun(Throughfall(i), Interception(i), T, SnowRain_par, fAPAR(i))

  ! Excess water from canopy will drip down to soil if not evaporated
  ! during the day, rest remains in canopy for the next day
  if (SnowRain_par%CWmax <= 1.0e-8) then
    Throughfall(i) = Throughfall(i) + Interception(i)
  else
    if (Interception(i) + theta_canopy > SnowRain_par%CWmax * fAPAR(i)) then
      Throughfall(i) = Throughfall(i) + Interception(i) + theta_canopy - SnowRain_par%CWmax * fAPAR(i)
      theta_canopy = SnowRain_par%CWmax * fAPAR(i)
    else
      theta_canopy = Interception(i) + theta_canopy
    end if
  end if


      !--------------------------------------------------------------
      ! Calculate daily evapotranspiration using the ET function.
      ! ET is throughfall to soil water after canopy interception
      ! and snow dynamics       
      !--------------------------------------------------------------
  
          
      ET(i) = ETfun( D, theta_root, II, fAPAR(i), T,                 &
                    ET_par, Site_par,                         &
                    theta_canopy,                             & ! input as scalar
                    theta_top,                             & ! input as scalar
                    fE(i),                                    & ! soil‑water constraint on evaporation
                    GPP(i),                                   & ! changed from gpp380
                    fW(i),                                    & ! soil‑water constraint of GPP at 380 ppm
                    fOrg(i),                                  & ! soil‑water constraint for evaporation in top layer
                    GPP_par,                                  & ! CO₂ influence on GPP
                    CO2(i), LOGFLAG, etmodel,                 &
                    transp(i), evap(i), fWE(i),               & ! outputs
                    CO2model, soilmodel, REWmodel)              ! alternative model switches
      

   
! Choose different subroutine for mineral soils and peatland
      
    if(soilmodel < 2) then 

      !--------------------------------------------------------------
      ! Mineral soils: Update soil‑water balance, drainage and 
      ! related variables.
      !--------------------------------------------------------------
      call sw_balance( theta, Throughfall(i), Snowmelt(i), ET(i), &
                      Site_par, Drainage(i),                     &
                      theta_snow, theta_canopy, SnowRain_par )

      !--------------------------------------------------------------
      ! Record the state variables that contain the water storages.
      !--------------------------------------------------------------
      SOG(i)         = theta_snow
      SW(i)          = theta
      CanopyWater(i) = theta_canopy
      S(i)           = S_state
      
    else
        
      !--------------------------------------------------------------
      ! Peatland soils: Update soil‑water balance 
      ! and related variables. Water balance for entire profile
      !--------------------------------------------------------------
      call swbal_peat(theta, theta_top, pond,  throughfall(i), snowmelt(i),  &
                        transp(i), evap(i), et(i), Site_par, Drainage(i),      &
                        theta_snow, theta_canopy, SnowRain_par)
      
      !--------------------------------------------------------------
      ! Peatland soils: Calculate water level and update drainage 
      ! and theta_root (depends on water level)
      !--------------------------------------------------------------
      if(i+1 < NofDays) then
      call water_level(theta, SWTable, gwl, Drainage(i+1), theta_root, Site_par, dimTable)    
      endif

      !--------------------------------------------------------------
      ! Record the state variables that contain the water storages.
      !--------------------------------------------------------------
      SOG(i)         = theta_snow
      SW(i)          = theta
      CanopyWater(i) = theta_canopy
      S(i)           = S_state
      ST(i)          = theta_top
      SR(i)          = theta_root
      WL(i)          = gwl
      
    end if
    
        

   end do   ! day loop


    end subroutine fpreles

    ! end module preles_module
    