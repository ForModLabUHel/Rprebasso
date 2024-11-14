

subroutine fireDist(Cpool_litter_woodIn,Cpool_litter_greenIn,livegrass,soil_moisture, & 
			TAir,NI,Precip,FDI) 
! livegrass = gv biomass
! Cpool_litter_wood = soilC woody component AWEN
! Cpool_litter_green = soilC non-woody component AWEN

  implicit none
  
  integer, parameter :: nDays=365
  real (kind=8), intent(in) :: Cpool_litter_woodIn,Cpool_litter_greenIn,livegrass
  real (kind=8), intent(in) :: soil_moisture, TAir(nDays), NI(nDays), Precip(nDays)
  real (kind=8), intent(inout) :: FDI(nDays)
  real (kind=8) :: alpha_livegrass(nDays),alpha_fuel(nDays)
  real (kind=8) :: char_alpha_fuel(nDays),rel_fuel_moisture(nDays)
  real (kind=8) :: leaf_moisture,Cpool_litter_wood(nDays),Cpool_litter_green(nDays)
  real (kind=8) :: SurfArea2Vol(3), moisture_scaling, frac_green_active, frac_1hr_wood, frac_10hr_wood
  real (kind=8) :: frac_100hr_wood, frac_1000hr_wood, moistfactor, moistfactor_livegrass
  real (kind=8) :: fuel_1hr(nDays),fuel_10hr(nDays),fuel_100hr(nDays),fuel_1000hr(nDays),fuel_1to100hr_sum(nDays)
  real (kind=8) :: ratio_dead_fuel(nDays),ratio_live_fuel(nDays),char_moistfactor(nDays)
  real (kind=8) :: alpha_fuel_1hr,alpha_fuel_10hr,alpha_fuel_100hr
  integer :: i

!initialize
 Cpool_litter_wood(:) = Cpool_litter_woodIn
 Cpool_litter_green(:) = Cpool_litter_greenIn
 ! NI(:) = 0.0
 alpha_livegrass(:) = 0.0
 rel_fuel_moisture(:) = 0.0
 FDI(:) = 0.0
 
 ! Parameters for the fire model
 ! Surface area to volume of fuels
 SurfArea2Vol = (/66.,3.58,0.98/)
 moisture_scaling = 3.*10.**4.
! Allocate C pools to fuel classes
 frac_green_active = 0.3
 frac_1hr_wood = 0.045
 frac_10hr_wood = 0.12
 frac_100hr_wood = 0.21
 frac_1000hr_wood = 0.67
! Moisture of extinction
 moistfactor = 0.3
 moistfactor_livegrass = 0.2

 !Fire risk modelling
! Fuel load in g/m2
 fuel_1hr = Cpool_litter_green*frac_green_active+Cpool_litter_wood*frac_1hr_wood
 fuel_10hr = Cpool_litter_wood*frac_10hr_wood
 fuel_100hr = frac_100hr_wood*frac_100hr_wood
 fuel_1000hr = frac_1000hr_wood*frac_1000hr_wood

 fuel_1to100hr_sum = fuel_1hr+fuel_10hr+fuel_100hr
 ratio_dead_fuel = fuel_1to100hr_sum / (fuel_1to100hr_sum + livegrass)
 ratio_live_fuel = livegrass / (fuel_1to100hr_sum + livegrass)

! ! Nesterov Index
! ! A cumulative function of daily Tmax and dew-point temperature Tdew, eq. 5 in TH2010
 ! do i = 2,nDays
	 ! if(Precip(i)<3. .and. (Tmin(i)-4)>=0.) NI(i) = (Tmax(i)*(Tmax(i)-Tmin(i)-4.))+(Tmax(i-1.)*(Tmax(i-1.)-Tmin(i-1.)-4.))
 ! end do

 ! Calculate moisture scaling of fuels
 alpha_fuel_1hr = SurfArea2Vol(1)/moisture_scaling
 alpha_fuel_10hr = SurfArea2Vol(2)/moisture_scaling
 alpha_fuel_100hr = SurfArea2Vol(3)/moisture_scaling

 ! This should be daily
 alpha_fuel = (alpha_fuel_1hr*fuel_1hr+alpha_fuel_10hr*fuel_10hr+alpha_fuel_100hr*fuel_100hr)/fuel_1to100hr_sum

 leaf_moisture = max(0., (10./9.*soil_moisture-1./9.)) ! Eq. B2 in TH2010
 alpha_livegrass = -log(leaf_moisture)/NI
 char_alpha_fuel = alpha_fuel*ratio_dead_fuel+alpha_livegrass*ratio_live_fuel

! Moisture of extinction
 char_moistfactor = moistfactor*ratio_dead_fuel+moistfactor_livegrass*ratio_live_fuel

! A weighted-average of the relative moisture content of fuels, eq. 6 in TH2010
! Pspread decreases linearly as litter moisture increases towards its moisture of extinction
do i = 1, nDays
! for (d in 1:length(TAir)){
  rel_fuel_moisture(i) = exp(-char_alpha_fuel(i)*NI(i)) !Fuel moisture
  FDI(i) = max(0., (1.-(1./char_moistfactor(i)*rel_fuel_moisture(i)))) !Fire Danger Index
enddo
  
end subroutine
