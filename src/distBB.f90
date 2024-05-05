subroutine TsumSBBfun(lat,temp,TsumSBB) 

  implicit none
  
  integer, parameter :: nDays=365
  real (kind=8), intent(in) :: lat,temp(365)
  real (kind=8), intent(inout) :: TsumSBB
  real (kind=8) :: pi = 3.1415927
  real (kind=8) :: two_pi,sin_dayl_deg_to_rad
  real (kind=8) :: deg_to_rad, pi_1, seconds_per_hour
  real (kind=8) :: doy(nDays),dec(nDays), mult(nDays),dayl_seconds(nDays),daylights(nDays)
  real (kind=8) :: sinld(nDays), cosld(nDays),aob(nDays),dayl_hours(nDays),effBTS(nDays)
  integer :: i
  real (kind=8) :: Tpar_alpha, Tpar_beta, Tpar_gamma, Tpar_Tmax, Tpar_T0, Tpar_DTL
  
  two_pi = 2. * pi
  sin_dayl_deg_to_rad = 0.3979486
  deg_to_rad = 0.01745329
  pi_1 = 0.3183099
  seconds_per_hour = 3600.
  ! TsumSBB model parameters from publication:
  Tpar_alpha=0.02876507
  Tpar_beta=3.5922336
  Tpar_gamma=1.24657367
  Tpar_Tmax=40.9958913
  Tpar_T0=30.4
  Tpar_DTL=8.3

  !inputs
  do i = 1,nDays
   doy(i) = real(i,8)
  enddo
  ! declination
  dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10 ) * 0.002739726 ) )
  !latitude in radians
  mult = lat * deg_to_rad
  !day length is estimated as the ratio of sin and cos of the product of declination an latitude in radians
  sinld = sin( mult ) * sin( dec )
  cosld = cos( mult ) * cos( dec )
  !cosld = cos( mult ) * cos( dec )
  aob = sinld / cosld
  do i = 1,ndays
   if(aob(i) > 1.) aob(i) = 1.
   if(aob(i) < -1.) aob(i) = -1.
  enddo
  !estimate day length in hours and seconds and upload to module variables
  dayl_hours = 12. * ( 1. + 2. * asin( aob ) * pi_1 )
  dayl_seconds = dayl_hours * seconds_per_hour
  
  do i=1,nDays
   if(temp(i)> Tpar_DTL .and. temp(i) < 38.9) then
    effBTS(i) = (Tpar_T0-Tpar_DTL)*(exp(Tpar_alpha*temp(i))- exp(Tpar_alpha*Tpar_Tmax-(Tpar_Tmax-temp(i))/Tpar_beta)-Tpar_gamma)
   else
    effBTS(i) = 0.
   endif
   if(dayl_hours(i)>= 14.5) then
	daylights(i)=1.
   else
    daylights(i)=0.
   endif
  enddo
  effBTS = effBTS*daylights
  TsumSBB = sum(effBTS)/557. ! this is used as input in SBB model
end subroutine


!#####
! From outputs, calculate:
! BA_spruce = Basal area of spruce
! BAspruceFract = Spruce fractions of basal area
! age_spruce = Spruce age (maximum among spruce layers)
! SMI_Tprev = SMI from the previous year
! TsumSBB_Tprev2, TsumSBB_Tprev, TsumSBB_T: TsumSBB from the years t-2, t-1 and t (the current year)
subroutine riskBB(pBB,TsumSBBs,BA_spruce,BAtot,age_spruce,SMI) 

  implicit none
  
  real (kind=8), intent(in) :: TsumSBBs(3),BA_spruce,BAtot,age_spruce,SMI !TsumSBBs vector of three elements: two years ago, previous year and current year
  real (kind=8), intent(inout) :: pBB(5)
  real (kind=8) :: BAspruceFract,PI_agespruce,PI_BAspruce
  real (kind=8) :: x0, k, PI_spruceFract,PI_SMITprev
  real (kind=8) :: TsumsMean(3), C(3,3),TsumSBBcurr,TsumSBBprev,TsumSBBprev2
  real (kind=8) :: TsumsMat, Tsums,x1,x2,gen
  real (kind=8) :: aspruceshare, aage, aBA, adrought,PI

!!parameters
  TsumsMean(:) = (/1.418036, 1.632891, 1.408584/)
  C(1,:) = (/10.16699, -22.40190, -13.48438/)
  C(2,:) = (/0.00000,  20.23612, -17.89729/)
  C(3,:) = (/0.00000,   0.00000,  31.11091/)
! Predisposition PI is a weighted sum of four PI_components:
  aspruceshare = 0.3
  aage = 0.25
  aBA = 0.15
  adrought = max(1.-aspruceshare-aage-aBA,0.)
  
  BAspruceFract = BA_spruce/BAtot

! PI for BA spruceFract
  x0 = 0.4
  k = -10.
  PI_spruceFract = 1./(1.+exp(k* (BAspruceFract - x0)))

! PI for age_spruce
  x0 = 60.
  k = -0.1
  PI_agespruce = 0.2 + 0.8/(1.+exp(k* (age_spruce - x0)))

! BA_spruce
  x0 = 1.
  k = -1.
  PI_BAspruce = 1./(1.+exp(k* (BA_spruce - x0)))

! SMI_Tprev
  x0 = 0.89
  k = 100.
  PI_SMITprev = 1./(1.+exp(k* (SMI - x0)))

! Zero probability for SBB if the long term temperature is not high enough 
  TsumsMat = sum(MATMUL((TsumSBBs-TsumsMean),C))
  Tsums = exp(-TsumsMat**2.)/2.
  if(Tsums < 0.0001) then
   Tsums =  0.
  else
   Tsums =  1.
  endif

  PI = (aspruceshare*PI_spruceFract + aage*PI_agespruce + &
         aBA*PI_BAspruce + adrought*PI_SMITprev)

! GEN is the bark beetle generation index, which depends on temperature
   ! The previous year TsumSBB is used for bark beetle generations
  gen = 0.  
  if(TsumSBBs(2)<1.) gen = 0.
  if(TsumSBBs(2) >= 1. .and. TsumSBBs(2) < 1.5) gen = 0.1
  if(TsumSBBs(2) >= 1.5 .and. TsumSBBs(2) < 2.) gen = 0.2
  if(TsumSBBs(2) >= 2. .and. TsumSBBs(2) < 2.5) gen = 0.6
  if(TsumSBBs(2) >= 2.5) gen = 1

! probability function coefficients from Seild et al. 2007
  x1 = -1.51
  x2 = 1.65

! SBB probability
  pBB(1) =(1-exp(Tsums*x1*PI**x2)**gen)
  pBB(2) = PI_spruceFract
  pBB(3) = PI_agespruce 
  pBB(4) = PI_BAspruce 
  pBB(5) = PI_SMITprev

end subroutine




!!!Calculate spruce variables   
subroutine spruceVars(standInfo,nLayers,spruceIDs,nSpIDs,spruceStandVars)
  implicit none
  logical :: spruceLayer(nLayers)
  integer, intent(in) :: nLayers,nSpIDs
  real(8), intent(in) :: standInfo(3,nLayers),spruceIDs(nSpIDs) !standInfo first argument is species, age and BA by layer
  real(8), intent(inout) :: spruceStandVars(3)
  integer :: i
  real(8) :: ba(nLayers), age(nLayers),species(nLayers)
  real(8) :: BAspruce,ageMaxSpruce,BAspruceShare,BAtot
  real(8) :: par_BAshare,par_PIba,par_PIage,par_PIdrought,PI_age,PI_ba,PI_drought

ba = standInfo(3,:)
age = standInfo(2,:)
species = standInfo(1,:)
BAtot = sum(ba)
BAspruce = 0.d0
ageMaxSpruce = 0.d0
BAspruceShare = 0.d0

do i =  1,nLayers
  spruceLayer(i) = any(spruceIDs .eq. int(species(i)))
  if(spruceLayer(i)) then
    BAspruce = BAspruce + BA(i)
    ageMaxSpruce = max(0.d0,age(i))
    BAspruceShare = BAspruceShare + BA(i)/BAtot
  endif
end do

 spruceStandVars(1) = BAspruce
 spruceStandVars(2) = ageMaxSpruce
 spruceStandVars(3) = BAtot
 
end subroutine

