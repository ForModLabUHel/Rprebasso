module prelesglobalsF_module
  implicit none
  
  public:: p1,p2,p3,p4,p5,p7
  public:: TARGACCEPT, MAXK, MYINFINITY, PI, NUMBER_OF_MODEL_PARAMETERS, K, vectorlength

  integer, parameter :: TARGACCEPT = 44 ! 0.44 scaled by 100 for integer parameter
  integer, parameter :: MAXK = 1000
  real(8), parameter :: MYINFINITY = 999999999.9d0
  real(8), parameter :: PI = 3.1415926535d0
  integer, parameter :: NUMBER_OF_MODEL_PARAMETERS = 38

  integer :: K
  integer :: vectorlength

  type :: p1   !Site_par
    real(8) :: soildepth
    real(8) :: ThetaFC
    real(8) :: ThetaPWP
    real(8) :: tauDrainage
    real(8) :: topdepth   !new  (mm)
    real(8) :: orgthres   !new  (unitless)
    real(8) :: MaxPond    !new  (mm)
    real(8) :: ditchDepth !new  (m)
    real(8) :: ditchDist  !new  (m)
    real(8) :: peatdepth  !new  (mm)
  end type p1

  type :: p2  ! GPP_par
    real(8) :: beta
    real(8) :: tau
    real(8) :: S0
    real(8) :: Smax
    real(8) :: kappa
    real(8) :: gamma
    real(8) :: soilthres
    real(8) :: bCO2
    real(8) :: xCO2
    real(8) :: t0
    real(8) :: tcrit
    real(8) :: tsumcrit
    real(8) :: soils !new  - parameter for soil threshold function
  end type p2

  type :: p3  ! ET_par
    real(8) :: beta
    real(8) :: kappa
    real(8) :: chi
    real(8) :: soilthres
    real(8) :: nu

  end type p3

  type :: p4   ! SnowRain_par
    real(8) :: MeltCoef
    real(8) :: I0
    real(8) :: CWmax
    real(8) :: SnowThreshold
    real(8) :: T_0
  end type p4

  type :: p5   !Init_par
    real(8) :: SW
    real(8) :: CW
    real(8) :: SOG
    real(8) :: S
    real(8) :: ST  !new  (initial water level in top layer, mm)
    real(8) :: WL   !new (initial water table depth, mm, needs to be synchronised with water storage SW, to be consistent with SWTable)
    end type p5
  !
  !type :: p6
  !  real(8) :: cvGPP
  !  real(8) :: cvET
  !  real(8) :: cvSW
  !end type p6
  
   type :: p7   !Genuchten_par (all new, for one peat type)
    real(8) :: thetaS
    real(8) :: thetaR
    real(8) :: alpha
    real(8) :: n
    real(8) :: thetaST
    real(8) :: thetaRT
    real(8) :: alphaT
    real(8) :: nT
   end type p7
   
  

end module prelesglobalsF_module