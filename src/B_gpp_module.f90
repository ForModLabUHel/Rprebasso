!=====================================================================
!  Translated from the original C source (Kalliokoski version)
!  Implements the seasonality, phenology, CO₂ modifiers and GPP model
!=====================================================================
module gpp_module
    
  use prelesglobalsF_module
  
  implicit none

contains

  !=================================================================
  !  fS_model – seasonality model (Mäkelä et al. 2004)
  !=================================================================
  function fS_model(S, T, GPP_par) result(fS)
    real(8), intent(inout) :: S          
    real(8), intent(in)    :: T
    type(p2), intent(in)    :: GPP_par
    real(8)                :: fS
    real(8)                :: temp

    ! *S = *S + (T-*S)/GPP_par.tau;
    S = S + (T - S) / GPP_par%tau

    ! if (0 > *S-GPP_par.S0) fS=0; else fS= *S-GPP_par.S0;
    temp = S - GPP_par%S0
    if (temp < 0.0) then
       fS = 0.0
    else
       fS = temp
    end if

    ! if (1 < fS/GPP_par.Smax) fS=1; else fS=fS/GPP_par.Smax;
    temp = fS / GPP_par%Smax
    if (temp > 1.0) then
       fS = 1.0
    else
       fS = temp
    end if
  end function fS_model


  !=================================================================
  !  fPheno_model – phenology model
  !=================================================================
  function fPheno_model(GPP_par, T, PhenoS, DOY, fS) result(fPheno)
    type(p2),    intent(in)    :: GPP_par
    real(8),    intent(in)    :: T
    real(8),    intent(inout) :: PhenoS   ! updated in place
    integer,     intent(in)    :: DOY
    real(8),    intent(in)    :: fS
    real(8)                   :: fPheno
    real(8)                   :: m

    if (GPP_par%t0 > -998.0) then
       ! Budbreak must occur between specified min. date and end of July
       if ( (DOY > (GPP_par%t0 - 0.5)) .and. (DOY < 213) ) then
          m = T - GPP_par%tcrit
          if (m < 0.0) m = 0.0
          PhenoS = PhenoS + m
       else
          PhenoS = 0.0
       end if

       if (PhenoS > GPP_par%tsumcrit - 0.005) then
          fPheno = 1.0
       else
          fPheno = 0.0
       end if

       ! Quick solution to leaf‑out after end of July
       if (DOY > 212) then
          fPheno = fS * fS
          if (fPheno < 0.5) fPheno = 0.0
       end if
    else
       ! No t0 → evergreen
       fPheno = 1.0
    end if
  end function fPheno_model


  !=================================================================
  !  CO₂ modifier – mean GPP response - select
  !=================================================================
 
  function fCO2_model_mean(CO2model, CO2, GPP_par) result(fCO2)
    real(8), intent(in) :: CO2
    type(p2), intent(in) :: GPP_par
    integer,  intent(in) :: CO2model
    real(8)             :: fCO2
 
    select case (CO2model)
    case (1)      ! Kolari
    fCO2 = 1.0 + (CO2 - 380.0) / (CO2 - 380.0 + GPP_par%bCO2)
 
    case (2)      ! Launiainen
    fCO2 = 1.0 + (CO2-380)/(CO2-380+GPP_par%bCO2)
     
   case default   ! Kolari
   fCO2 = 1.0 + (CO2 - 380.0) / (CO2 - 380.0 + GPP_par%bCO2)
  
    end select
    
end function fCO2_model_mean

  !=================================================================
  !  CO₂ modifier – ET (evapotranspiration) response - select
  !=================================================================
  
  
      
  function fCO2_ET_model_mean(CO2model, CO2, GPP_par) result(fCO2_ET)
    real(8), intent(in) :: CO2
    type(p2), intent(in) :: GPP_par
    integer,  intent(in) :: CO2model
    real(8)             :: fCO2_ET
    REAL(8)             :: f_C_P
    REAL(8)             :: f_C_ET
    
    select case (CO2model)
    case (1)   ! Kolari
    fCO2_ET = 1.0 - (CO2 - 380.0) / (CO2 - 380.0 + GPP_par%bCO2)
  
   case (2)   ! Launiainen
   f_C_P = 1 + GPP_par%bCO2 * log(CO2/380);
   f_C_ET = 1 + GPP_par%xCO2 * log(CO2/380);
   fCO2_ET = (1/f_C_P)*f_C_ET;
 
   case default   ! Kolari
   fCO2_ET = 1.0 - (CO2 - 380.0) / (CO2 - 380.0 + GPP_par%bCO2)
   
   endselect
   
  end function fCO2_ET_model_mean
                 
 !=================================================================
  !  Soil water sensitivity function - new in gpp module
  !=================================================================
  function fREW_fun(REW, s, rho, REWmodel) result(fREW)
    real(kind=8), intent(in)    :: REW        ! relative extractable water
    real(kind=8), intent(in)    :: s          ! exponent parameter
    real(kind=8), intent(in)    :: rho        ! threshold parameter
    integer, intent(in)         :: REWmodel   ! REWmodel = 1: continuous REWmodel = 2: piecewise linear
   
    real(kind=8)                :: fREW       ! result function (soil water modifier)

   select case (REWmodel)
   case(1)
       fREW = (1 + EXP(s * REW)) / (1 + exp(s * (REW - rho/2.0)))
       
   case(2)
       if (REW < rho) then
          if (REW > 0.01 ) then
             fREW = REW / rho
          else
             fREW = 0.0 
          end if
       else
          fREW = 1.0 
       end if
       
   case default
    fREW = (1 + EXP(s * REW)) / (1 + exp(s * (REW - rho/2.0)))
       
   end select
   
  end function fREW_fun


     
  !=================================================================
  !  GPPfun – Gross Primary Production model (modified Mäkelä et al. 2008)
  !=================================================================
  subroutine GPPfun(gpp, gpp380, ppfd, D, CO2, theta, fAPAR, fSsub, &
                    GPP_par, Site_par, fL, fD, fW, fE, LOGFLAG, CO2model, REWmodel)
    real(8), intent(out)   :: gpp          ! actual GPP
    real(8), intent(out)   :: gpp380       ! reference‑CO₂ GPP
    real(8), intent(in)    :: ppfd
    real(8), intent(in)    :: D
    real(8), intent(in)    :: CO2
    real(8), intent(in)    :: theta
    real(8), intent(in)    :: fAPAR
    real(8), intent(in)    :: fSsub
    type(p2), intent(in)    :: GPP_par
    type(p1), intent(in)    :: Site_par
    real(8), intent(out)   :: fL           ! light stress factor
    real(8), intent(out)   :: fD           ! drought stress factor
    real(8), intent(out)   :: fW           ! water stress factor
    real(8), intent(out)   :: fE           ! combined stress factor (unused here)
    integer,  intent(in)    :: LOGFLAG      ! kept for compatibility
    integer, intent(in)    :: CO2model
    integer, intent(in)    :: REWmodel

    real(8) :: thetavol, REW
    real(8) :: fCO2, fDsub, fWsub, fLsub, fEsub
    real(8) :: pow_term

    ! --------------------------------------------------------------
    ! Soil water potential (relative to field capacity)
    ! --------------------------------------------------------------
    thetavol = theta / Site_par%soildepth
    REW = (thetavol - Site_par%ThetaPWP) / &
          (Site_par%ThetaFC - Site_par%ThetaPWP)

    ! --------------------------------------------------------------
    ! Drought stress – exponential response : select
    ! --------------------------------------------------------------
    
    select case (CO2model)
    case (1)   ! Kolari
    
    pow_term = (380.0 / CO2) ** GPP_par%xCO2
    fDsub = exp( GPP_par%kappa * pow_term * D )
    if (fDsub > 1.0) fDsub = 1.0
    
    case (2)   ! Launiainen
    
    pow_term = (380.0 / CO2) ** GPP_par%xCO2
    fDsub = exp( GPP_par%kappa * D )
    if (fDsub > 1.0) fDsub = 1.0
    
    case default   ! Kolari
    
    pow_term = (380.0 / CO2) ** GPP_par%xCO2
    fDsub = exp( GPP_par%kappa * pow_term * D )
    if (fDsub > 1.0) fDsub = 1.0
    
    endselect

    ! --------------------------------------------------------------
    ! Water stress – based on relative extractable water (REW)
    ! --------------------------------------------------------------
    if (GPP_par%soilthres < -998.0) then
       fWsub = 1.0                     ! no water control
    else
        fWsub = fREW_fun(REW, GPP_par%soils, GPP_par%soilthres, REWmodel)
    end if

    ! --------------------------------------------------------------
    ! Light limitation factor (L)
    ! --------------------------------------------------------------
    fLsub = 1.0 / ( GPP_par%gamma * ppfd + 1.0 )

    ! --------------------------------------------------------------
    ! Combined stress factor (minimum of drought and water stress)
    ! --------------------------------------------------------------
    if (fDsub > fWsub) then
       fEsub = fWsub
    else
       fEsub = fDsub
    end if

    ! Return stress factors to caller
    fW = fWsub
    fD = fEsub
    fE = fEsub   ! kept for compatibility with original interface
    fL = fLsub   ! added for consistency

    ! --------------------------------------------------------------
    ! Reference GPP (CO₂ = 380 ppm)
    ! --------------------------------------------------------------
    gpp380 = GPP_par%beta * ppfd * fAPAR * fSsub * fLsub * fEsub

    ! --------------------------------------------------------------
    ! Apply CO₂ effect
    ! --------------------------------------------------------------
    fCO2 = fCO2_model_mean(CO2model, CO2, GPP_par)
    gpp = gpp380 * fCO2
  end subroutine GPPfun

end module gpp_module
