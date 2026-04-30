!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PRELES MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Model function:
 !   1. Makes initialisations
 !   2. Estimates GPP
 !   3. Estimates snow and interception
 !   4. Estimates Evapotranspiration
 !   5. Updates soil water balance
    
! module preles_module
	  ! use gpp_module
	  ! use water_module
	  ! use prelesglobalsF_module
	  ! use initruns_module


      ! implicit none
      
    ! contains
    subroutine preles_f_crobas(NofDays,weather,DOY,fAPAR,prelesOut,pars,&
				GPP,ET,SW,etmodel,CO2model,soilmodel,REWmodel,SWTable, calc_SWtable)!,p0)

     implicit none


		integer, parameter :: dimTable=200
		integer, intent(inout) :: NofDays, calc_SWtable
		real (kind=8), intent(inout) :: weather(NofDays,5),fAPAR(NofDays)
		real (kind=8), intent(inout) :: prelesOut(19)!,p0
		real (kind=8), intent(inout) :: pars(60)
		integer, intent(inout):: DOY(NofDays), etmodel,CO2model
		real(8), intent(out) :: GPP(NofDays), ET(NofDays),SW(NofDays)
		integer, intent(in) :: soilmodel   ! 1 for mineral soil, 2 for drained peatland, default mineral
		integer, intent(in) :: REWmodel   ! 1 for smooth function, 2 for piecewise linear, default smooth
		real(8), intent(inout) :: SWTable(dimTable, dimTable+3)
    
		real(8) :: PAR(NofDays), TAir(NofDays), VPD(NofDays), Precip(NofDays), CO2(NofDays)
		real(8) :: Site_par(10)
		real(8) :: GPP_par(13)
		real(8) :: Init_par(7)
		real(8) :: ET_par(5)
		real(8) :: SnowRain_par(5)
		real(8) :: Genuchten_par(8)
		real(8) :: Ksat(12)
		real(8) :: ST(NofDays), SR(NofDays), WL(NofDays), SOG(NofDays), S(NofDays)
		real(8) :: fS(NofDays), fD(NofDays), fW(NofDays), fE(NofDays), fL(NofDays)
		real(8) :: Throughfall(NofDays), Interception(NofDays), Snowmelt(NofDays)
		real(8) :: Drainage(NofDays), Canopywater(NofDays)
		real(8) :: GPPmeas(NofDays), ETmeas(NofDays), SWmeas(NofDays)
		integer :: LOGFLAG = 1
		integer :: multisiteNday = 1
		real(8) :: transp(NofDays), evap(NofDays), fWE(NofDays), fOrg(NofDays)



		Ksat=0.
		PAR = weather(:,1)
		TAir = weather(:,2)
		VPD = weather(:,3)
		Precip = weather(:,4)
		CO2 = weather(:,5)
		Site_par=pars(1:10)
		GPP_par(1:9)=pars(11:19)
		GPP_par(10:12)=pars(38:40)
		GPP_par(13) = pars(20)
		ET_par=pars(21:25)
		SnowRain_par=pars(26:30)
		Init_par=pars(31:37)
		

		Genuchten_par=pars(41:48)
		Ksat=pars(49:60)

		SW  = 0.
		ST  = 0.
		SR  = 0.
		SOG = 0.
		S   = 0.
  
  
! LOGFLAG = 1
! multisiteNday = 1

!N of simulation days
! N = length(PAR)

! day = 1:N


		SW(1) = Init_par(1)
		Canopywater(1) = Init_par(2)
		SOG(1) = Init_par(3)
		S(1) = Init_par(4)
		ST(1) = Init_par(5)
		WL(1) = Init_par(6)!/1000.
		if(Init_par(7)<0.) then
			SR(1) = Site_par(2)*200.
		else
			SR(1) = Init_par(7)
		endif

		call preles_fortran(NofDays, PAR, TAir, VPD, Precip, CO2, fAPAR, Site_par, GPP_par, ET_par, SnowRain_par, etmodel, &
                  Genuchten_par, Ksat, WL, dimTable, GPP, ET, SW, ST, SR, SOG, fL, fS, fD,                                                             &
                   fW, fE, Throughfall, Interception, Snowmelt, Drainage, Canopywater, GPPmeas, ETmeas, SWmeas, &
                  S, LOGFLAG, multisiteNday, doy, transp, evap, fWE, fOrg, CO2model, soilmodel, REWmodel,SWTable, calc_SWtable)

		ET = transp + evap
		prelesOut(1) = sum(GPP(1:NofDays))
		prelesOut(2) = sum(ET(1:NofDays))
		prelesOut(3) = SW(NofDays)
		prelesOut(4) = SOG(NofDays)
		prelesOut(5) = fS(NofDays)
		prelesOut(6) = fD(NofDays)
		! prelesOut(7) = fW(NofDays)
		prelesOut(8) = fE(NofDays)
		prelesOut(9) = Throughfall(NofDays)
		prelesOut(10) = Interception(NofDays)
		prelesOut(11) = Snowmelt(NofDays)
		prelesOut(12) = Drainage(NofDays)
		prelesOut(13) = Canopywater(NofDays)
		prelesOut(14) = S(NofDays)
		prelesOut(15) = sum(SW(1:NofDays))/NofDays
		prelesOut(16) = sum(SW(152:243))/92
		prelesOut(17) = ST(NofDays)
		prelesOut(18) = SR(NofDays)
		prelesOut(19) = WL(NofDays)
		
		! ET = WL
	end subroutine preles_f_crobas



  
    subroutine preles_fortran(NofDays, PAR, TAir, VPD, Precip, CO2, fAPAR, Site_par, GPP_par, ET_par, SnowRain_par, etmodel, &
                  Genuchten_par, Ksat, WL, dimTable, GPP, ET, SW, ST, SR, SOG, fL, fS, fD,      &
                   fW, fE, Throughfall, Interception, Snowmelt, Drainage, Canopywater, GPPmeas, ETmeas, SWmeas, &
                  S, LOGFLAG, multisiteNday, day, transp, evap, fWE, fOrg, CO2model, soilmodel, REWmodel,SWTable,calc_SWtable)
	  
	  use gpp_module
	  use water_module
	  use prelesglobalsF_module
	  use initruns_module

      implicit none
      
    
    integer, intent(in) :: NofDays,dimTable, calc_SWtable
    real(8), intent(inout) :: PAR(NofDays), TAir(NofDays), VPD(NofDays), Precip(NofDays), CO2(NofDays), fAPAR(NofDays)
    real(8), intent(in) :: Site_par(10)
    real(8), intent(in) :: GPP_par(13)
    real(8), intent(in) :: ET_par(5)
    real(8), intent(in) :: SnowRain_par(5)
    real(8), intent(in) :: Genuchten_par(8)
    real(8), intent(in) :: Ksat(12)
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
    real(8), intent(inout) :: SWTable(dimTable, dimTable+3)
    ! integer, intent(in) ::  dimTable   ! dimension of look-up table
    integer :: i, kk,  jj
    real(8) :: II, T, D, P, theta, theta_snow, theta_canopy, theta_top, theta_root, pond, S_state
    real(8) :: PhenoS, fPheno
    real(8) :: fEgpp, gpp380
    real(8) :: Qdrain,gwl
	real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
	real(8) :: GPP_par_beta,GPP_par_tau, GPP_par_S0,GPP_par_Smax,GPP_par_kappa, GPP_par_gamma,GPP_par_soilthres
	real(8) :: GPP_par_bCO2,GPP_par_xCO2,GPP_par_t0, GPP_par_tcrit,GPP_par_tsumcrit,GPP_par_soils
	real(8) :: ET_par_beta,ET_par_kappa,ET_par_chi,ET_par_soilthres,ET_par_nu
	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0
	real(8) :: Genuchten_par_thetaS,Genuchten_par_thetaR,Genuchten_par_alpha,Genuchten_par_n
	real(8) :: Genuchten_par_thetaST,Genuchten_par_thetaRT,Genuchten_par_alphaT,Genuchten_par_nT

!read parameters
	include 'init_gppPar_preles.h'
	include 'init_sitePar_preles.h'
	include 'init_etPar_preles.h'
	include 'init_SnowRainPar_preles.h'
	include 'init_GenuchtenPar_preles.h'

    
   ! Initialize variables
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
    pond = Site_par_MaxPond
    Drainage(1) = 0.
    gwl = WL(1)
  
    
    ! calculate look-up table for estimating water table depth
    ! SWTable is output from this subroutine
	!the subroutine is used only if 
	! a peat site is simulated and if the model runs independently from crobas
    if(soilmodel==2 .and. calc_SWtable==1) call waterTable(Genuchten_par, Site_par, Ksat, dimTable, SWTable)

	! open(1,file="test1.txt")
    ! write(1,*) dimTable
	! write(1,*) Genuchten_par
	! write(1,*) Site_par
	! write(1,*) Ksat
	! write(1,*) SWTable
	! close(1)
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
  if (SnowRain_par_CWmax <= 1.0e-8) then
    Throughfall(i) = Throughfall(i) + Interception(i)
  else
    if (Interception(i) + theta_canopy > SnowRain_par_CWmax * fAPAR(i)) then
      Throughfall(i) = Throughfall(i) + Interception(i) + theta_canopy - SnowRain_par_CWmax * fAPAR(i)
      theta_canopy = SnowRain_par_CWmax * fAPAR(i)
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

ET = evap+transp
    end subroutine preles_fortran

    ! end module preles_module
    