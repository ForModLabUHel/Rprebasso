!--------------------------------------------------------------------
! Subroutine: init_conditions
! Purpose   : Replace missing first‑day values with reasonable defaults.
!--------------------------------------------------------------------
    
module initruns_module
implicit none

    contains
    
subroutine init_conditions(par, tair, vpd, precip, co2)

   ! Arguments – all are scalars that will be modified (intent inout)
   real(kind=8), intent(inout) :: par      ! Photosynthetically active radiation
   real(kind=8), intent(inout) :: tair    ! Air temperature (°C)
   real(kind=8), intent(inout) :: vpd     ! Vapor pressure deficit (kPa)
   real(kind=8), intent(inout) :: precip  ! Precipitation (mm)
   real(kind=8), intent(inout) :: co2     ! Atmospheric CO₂ concentration (ppm)

   !-----------------------------------------------------------------
   ! Apply default values when input values are implausible or missing.
   ! The thresholds follow the original C logic.
   !-----------------------------------------------------------------
   if (par < -900.0d0) then
      par = 5.0d0               ! Dark winter day assumption
   end if

   if (tair < -900.0d0) then
      tair = 0.0d0              ! Reasonable baseline temperature
   end if

   if (vpd < 0.0d0 .or. vpd > 6.0d0) then
      vpd = 0.5d0               ! Typical moderate VPD; >6 kPa is implausible
   end if

   if (precip < 0.0d0) then
      precip = 0.0d0            ! No negative precipitation
   end if

   if (co2 < 0.1d0) then
      co2 = 380.0d0             ! Approximate pre‑industrial CO₂ level
   end if

end subroutine init_conditions

end module initruns_module