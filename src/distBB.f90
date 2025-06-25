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
subroutine riskBB(pBB,TsumSBBs,BA_spruce,BAtot,age_spruce,SMI,sitetype) 

  implicit none
  
  real (kind=8), intent(in) :: TsumSBBs(4),BA_spruce,BAtot,age_spruce,SMI,sitetype !TsumSBBs vector of four elements: three years ago,two years ago, previous year and current year
  real (kind=8), intent(inout) :: pBB(5)
  real (kind=8) :: BAspruceFract,PI_agespruce,PI_BAspruce
  real (kind=8) :: x0, k, PI_spruceFract,PI_SMITprev,x,f0, PI_sitetype, PI_SMIprev
  real (kind=8) :: x1,x2,gen,n
  real (kind=8) :: aspruceshare, aage, aBA, PI

!!parameters
! Predisposition PI is a weighted sum of four PI_components:
  aspruceshare = 0.45
  aage = 0.35
  aBA = 0.2
  ! adrought = 0.3!max(1.-aspruceshare-aage-aBA,0.)
  BAspruceFract=0.
  pBB = 0.
  if(BAtot>0.) BAspruceFract = BA_spruce/BAtot

! PI for BA spruceFract
  n = 2.
  k = 1.
  PI_spruceFract = min(1., k*BAspruceFract**n)

! PI for age_spruce
  x0 = 50.
  k = -0.09
  PI_agespruce = 1.0/(1.+exp(k* (age_spruce - x0)))

! BA_spruce
  n = 1.4
  k = 0.0068
  PI_BAspruce = min(1., k* BA_spruce**n)
 
! site type
  x0 = 3.5
  k = -0.7
  PI_sitetype = 1./(1.+exp(k* (sitetype - x0)))

 ! SMI_Tprev
  x0 = 0.88
  k = 50.
  PI_SMIprev = 1./(1.+exp(k* (SMI - x0)))

! Zero probability for SBB if the long term temperature is not high enough 
  x = sum(TsumSBBs)/4.
  x0 = 1.58
  k = -11. 
  f0 = 1./(1.+exp(k* (x - x0)))
  ! if(f0 < 0.0001) then
   ! f0 =  0.
  ! else
   ! f0 =  1.
  ! endif

  PI = (aspruceshare*PI_spruceFract + aage*PI_agespruce + &
         aBA*PI_BAspruce) * PI_sitetype * PI_SMIprev

! GEN is the bark beetle generation index, which depends on temperature
   ! The previous year TsumSBB is used for bark beetle generations
  gen = 0.0d0  
  if(TsumSBBs(3)<1.) gen = 0.0d0
  if(TsumSBBs(3) >= 1.0d0 .and. TsumSBBs(3) < 1.5d0) gen = 0.1d0
  if(TsumSBBs(3) >= 1.5d0 .and. TsumSBBs(3) < 2.0d0) gen = 0.2d0
  if(TsumSBBs(3) >= 2.0d0 .and. TsumSBBs(3) < 2.5d0) gen = 0.6d0
  if(TsumSBBs(3) >= 2.5d0) gen = 1.0d0

! probability function coefficients from Seild et al. 2007
  x1 = -1.51
  x2 = 1.65

! SBB probability
  pBB(1) =(1.0d0-exp(x1*PI**x2)**gen) * f0
  pBB(2) = PI_spruceFract
  pBB(3) = PI_agespruce 
  pBB(4) = PI_BAspruce 
  pBB(5) = PI_SMITprev
end subroutine




!!!Calculate spruce variables   
subroutine spruceVars(standInfo,nLayers,spruceIDs,nSpIDs,spruceStandVars,rBAspruce)
  implicit none
  logical :: spruceLayer(nLayers)
  integer, intent(in) :: nLayers,nSpIDs,spruceIDs(nSpIDs)
  real(8), intent(in) :: standInfo(3,nLayers) !standInfo first argument is species, age and BA by layer
  real(8), intent(inout) :: spruceStandVars(3), rBAspruce(nLayers)
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
rBAspruce(:) = 0.d0
spruceLayer(:) = .false.
if(BAtot>0.) then
 do i =  1,nLayers
  spruceLayer(i) = any(spruceIDs .eq. int(species(i)))
  if(spruceLayer(i)) then
    rBAspruce(i) = BA(i)
    BAspruce = BAspruce + BA(i)
    ageMaxSpruce = max(ageMaxSpruce,age(i))
    BAspruceShare = BAspruceShare + BA(i)/BAtot
  endif
 end do
 if(BAspruce>0.d0) rBAspruce = rBAspruce/BAspruce
endif

 spruceStandVars(1) = BAspruce
 spruceStandVars(2) = ageMaxSpruce
 spruceStandVars(3) = BAtot

end subroutine

!!! bark beetle impact model   
! returns the share of damaged spruce
subroutine bb_imp_mod(SMI,BA_spruceshare,dam_int)
  implicit none
  real(8), intent(in) :: SMI,BA_spruceshare
  real(8), intent(out) :: dam_int 
  real(8) :: a, SHI,x0,k,n

!!initialize parameters
  a = 1.
  n = 1.3
  k = 11.6566

  SHI = min(1.,(1.-SMI)/a)* BA_spruceshare
  dam_int = min(1., k*SHI**n)

end subroutine



!!!Calculate spruce ba mortality by layer based on BB disturbance
!calculate the BA mortality buy each layer starting from the bigger trees based on H   
subroutine baBBdist_bylay_fun(standInfo,nLayers,spruceIDs,nSpIDs,BAdist_bylay,BA_dist_tot)
  implicit none
  logical :: spruceLayer(nLayers)
  integer, intent(in) :: nLayers,nSpIDs,spruceIDs(nSpIDs)
  real(8), intent(inout) :: standInfo(3,nLayers), BA_dist_tot !standInfo first argument is species, age and BA by layer
  real(8), intent(inout) :: BAdist_bylay(nLayers)
  integer :: i, orderStand(nLayers),layerx
  real(8) :: ba(nLayers), h(nLayers),species(nLayers),batot
  
  
ba = standInfo(3,:)
h = standInfo(2,:)
species = standInfo(1,:)
BAdist_bylay = 0.

call order_desc(nLayers,h,orderStand)

BAtot = sum(ba)
spruceLayer(:) = .false.

if(BAtot>0.) then
 do i =  1,nLayers
  layerx = orderStand(i)
  spruceLayer(layerx) = any(spruceIDs .eq. int(species(layerx)))
  if(spruceLayer(layerx) .and. BA_dist_tot>0.) then
    if(ba(layerx) <= BA_dist_tot) then
     BAdist_bylay(layerx) = ba(layerx)
     BA_dist_tot = BA_dist_tot - ba(layerx) 
	else
     BAdist_bylay(layerx) = BA_dist_tot
     BA_dist_tot = 0. 	 
    endif	 
  endif
 enddo
endif

end subroutine
