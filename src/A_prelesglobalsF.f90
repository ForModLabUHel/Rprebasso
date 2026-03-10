module prelesglobalsF_module
  implicit none
  
  ! public:: p1,p2,p3,p4,p5,p7
  public:: TARGACCEPT, MAXK, MYINFINITY, PI, NUMBER_OF_MODEL_PARAMETERS, K, vectorlength

  integer, parameter :: TARGACCEPT = 44 ! 0.44 scaled by 100 for integer parameter
  integer, parameter :: MAXK = 1000
  real(8), parameter :: MYINFINITY = 999999999.9d0
  real(8), parameter :: PI = 3.1415926535d0
  integer, parameter :: NUMBER_OF_MODEL_PARAMETERS = 38

  integer :: K
  integer :: vectorlength
!  real(8) :: Site_par(10),GPP_par(13),ET_par(5),SnowRain_par(5),Genuchten_par(8)
!	real(8) :: Site_par_soildepth,Site_par_ThetaFC,Site_par_ThetaPWP,Site_par_tauDrainage,Site_par_topdepth
!	real(8) :: Site_par_orgthres,Site_par_MaxPond, Site_par_ditchDepth, Site_par_ditchDist, Site_par_peatdepth
!	real(8) :: GPP_par_beta,GPP_par_tau, GPP_par_S0,GPP_par_Smax,GPP_par_kappa, GPP_par_gamma,GPP_par_soilthres
!	real(8) :: GPP_par_bCO2,GPP_par_xCO2,GPP_par_t0, GPP_par_tcrit,GPP_par_tsumcrit,GPP_par_soils
!	real(8) :: ET_par_beta,ET_par_kappa,ET_par_chi,ET_par_soilthres,ET_par_nu
!	real(8) :: SnowRain_par_MeltCoef, SnowRain_par_I0, SnowRain_par_CWmax, SnowRain_par_SnowThreshold, SnowRain_par_T_0
!	real(8) :: Genuchten_par_thetaS,Genuchten_par_thetaR,Genuchten_par_alpha,Genuchten_par_n
!	real(8) :: Genuchten_par_thetaST,Genuchten_par_thetaRT,Genuchten_par_alphaT,Genuchten_par_nT


  ! type :: p1   !Site_par
    ! real(8) :: soildepth
    ! real(8) :: ThetaFC
    ! real(8) :: ThetaPWP
    ! real(8) :: tauDrainage
    ! real(8) :: topdepth   !new  (mm)
    ! real(8) :: orgthres   !new  (unitless)
    ! real(8) :: MaxPond    !new  (mm)
    ! real(8) :: ditchDepth !new  (m)
    ! real(8) :: ditchDist  !new  (m)
    ! real(8) :: peatdepth  !new  (mm)
  ! end type p1

  ! type :: p2  ! GPP_par
    ! real(8) :: beta
    ! real(8) :: tau
    ! real(8) :: S0
    ! real(8) :: Smax
    ! real(8) :: kappa
    ! real(8) :: gamma
    ! real(8) :: soilthres
    ! real(8) :: bCO2
    ! real(8) :: xCO2
    ! real(8) :: t0
    ! real(8) :: tcrit
    ! real(8) :: tsumcrit
    ! real(8) :: soils !new  - parameter for soil threshold function
  ! end type p2

  ! type :: p3  ! ET_par
    ! real(8) :: beta
    ! real(8) :: kappa
    ! real(8) :: chi
    ! real(8) :: soilthres
    ! real(8) :: nu

  ! end type p3

  ! type :: p4   ! SnowRain_par
    ! real(8) :: MeltCoef
    ! real(8) :: I0
    ! real(8) :: CWmax
    ! real(8) :: SnowThreshold
    ! real(8) :: T_0
  ! end type p4

  ! type :: p5   !Init_par
    ! real(8) :: SW
    ! real(8) :: CW
    ! real(8) :: SOG
    ! real(8) :: S
    ! real(8) :: ST  !new  (initial water level in top layer, mm)
    ! real(8) :: WL   !new (initial water table depth, mm, needs to be synchronised with water storage SW, to be consistent with SWTable)
    ! end type p5
  ! !
  ! !type :: p6
  ! !  real(8) :: cvGPP
  ! !  real(8) :: cvET
  ! !  real(8) :: cvSW
  ! !end type p6
  
   ! type :: p7   !Genuchten_par (all new, for one peat type)
    ! real(8) :: thetaS
    ! real(8) :: thetaR
    ! real(8) :: alpha
    ! real(8) :: n
    ! real(8) :: thetaST
    ! real(8) :: thetaRT
    ! real(8) :: alphaT
    ! real(8) :: nT
   ! end type p7
   
  

end module prelesglobalsF_module