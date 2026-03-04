 ! Seasonality model of Mäkelä et al 2004 
subroutine fS_model(T,tau,S0,Smax,S,fS)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: S, fS
   !inputs
   real(8), intent(in) :: T
   !parameters
   real(8), intent(in) :: tau, S0, Smax

  S = S + (T-S)/tau
  fS = max(0.d0, S-S0)
  fS = min(1.0d0, fS/Smax)
end subroutine


!fenology model
subroutine fPheno_model(T, DOY, fS, t0,tcrit,tsumcrit,PhenoS,fPheno) 

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: PhenoS, fPheno 
   !inputs
   real(8), intent(in) :: T, fS
   integer, intent(in) :: DOY
   !parameters
   real(8), intent(in) ::  t0, tcrit, tsumcrit
   !others
   real(8) :: m 


  if(t0 > -998) then!{ // ie not -999 
    !/* Budbreak must occur between specified min. date and end of July */
    if((DOY > (t0 - 0.5)) .and. (DOY < 213) ) then
      m = (T - tcrit)
      if(m < 0.d0) m = 0.d0
      PhenoS = PhenoS + m
    else
      PhenoS = 0.d0
    endif
    
    if(PhenoS > tsumcrit - 0.005d0) then
    fPheno = 1.d0
  else 
    fPheno = 0.d0
  endif
    ! /* Quick solution to leaf out: 
     ! * After end of July we just apply season prediction based on conifer fS 
     ! *  for gradual leaf out. Assume leaves drop much faster that fS.
     ! *  ...essentially this should be light driven process...i think. */
    if(DOY > 212) then
    fPheno = fS * fS
     if(fPheno < 0.5) fPheno = 0.d0
    endif
    
    ! /* If there is no t0 parameter, it is an evergreen */
  else
    fPheno = 1.d0
  endif
  
endsubroutine




! /* *****************************************************************/
! /* f-modifiers for increasing CO2 prepared by P. Kolari, pers. comm.*/
! /*double fCO2_model_mean(double CO2, p2 GPP_par ) {
  ! return(1 + (CO2-380)/(CO2-380+GPP_par.bCO2));
! }
! double fCO2_VPD_exponent(double CO2, p2 GPP_par ) {
  ! return(pow(CO2/380, GPP_par.xCO2));
! }
! */
! /*
! double fCO2_model_mean(double CO2, double b ) {
  ! return(1 + (CO2-380)/(CO2-380+b));
! }
! double fCO2_VPD_exponent(double CO2, double xCO2 ) {
  ! return(pow(380/CO2, xCO2));
! }
! */

! /* Note: 'ET_par.bC02' is the same as GPP_par.bCO2 */
! /*
! double fCO2_ET_model_mean(double CO2, p2 GPP_par ) {
  ! return(1 - 1.95*(CO2-380)/(CO2-380+(GPP_par.bCO2)));
! }
! */
! /* *****************************************************************/
! /* New CO2 modifiers based on APES simulator (Launiainen et al.)
   ! which account for the energy balance of the forest. Fitted responses
   ! to model predicted
   ! bCO2 = 0.5; xCO2 = -0.364
! */
subroutine fCO2_model_mean(CO2,bCO2,fCO2)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: fCO2
   !inputs
   real(8), intent(in) :: CO2
   !parameters
   real(8), intent(in) :: bCO2

  fCO2 = 1.d0 + bCO2 * log(CO2/380.d0) 
end subroutine


subroutine fCO2_ET_model_mean(CO2,xCO2,fCO2_ET)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: fCO2_ET
   !inputs
   real(8), intent(in) :: CO2
   !parameters
   real(8), intent(in) :: xCO2

  fCO2_ET = 1.d0 + xCO2 * log(CO2/380.d0) 
end subroutine


! /* GPP model, modified from Mäkelä et al 2008 */
subroutine GPPfun(gpp, gpp380, ppfd, D, CO2, theta, &
         fAPAR, fSsub, fD, fW, & 
         fE,  pars) 
 !      // double *fE, FILE *flog, int LOGFLAG) {
 !  extern double fCO2_model_mean(double CO2, p2 b ) ;
 !   //    extern double fCO2_VPD_exponent(double CO2, double xCO2 ) ;


 implicit none

   !state variables and outputs
   real(8), intent(inout) :: gpp, gpp380, fSsub, fD , fW, fE
   !inputs
   real(8), intent(in) :: CO2,ppfd,theta,fAPAR,D
   !parameters
   real(8), intent(in) :: pars(30)!preles parameters
   !other variables

real(8) :: thetavol,fCO2,REW!,GPPsub,GPP380sub
real(8) :: fEsub,fWsub,fLsub,fDsub
real(8) :: soildepth,ThetaPWP,ThetaFC
real(8) :: kappa,soilthres, gamma,beta,bCO2
!initialise
  soildepth = pars(1)
  ThetaFC = pars(2)
  ThetaPWP = pars(3)
  kappa = pars(9)
  soilthres = pars(11)
  gamma = pars(10)
  beta = pars(5)
    bCO2 = pars(12)
  
  thetavol = theta/soildepth
    REW=(thetavol-ThetaPWP)/(ThetaFC-ThetaPWP)
  

!    /* Calculate first the reference condition (ca=380 ppm) effect */
    fDsub = min(exp(kappa * D),1.0d0)
    ! fDsub = fDsub > 1 ? 1 : fDsub;

    if (soilthres < -998.0d0) then 
      fWsub = 1.0   ! /* e.g. -999 means no water control of GPP*/
    else 
      if (REW < soilthres) then
        if (REW > 0.01) then
      fWsub = REW/soilthres
    else
      fWsub = 0.0d0
    endif
      else
        fWsub = 1.0d0
      endif
    endif

    fLsub = 1.0d0 / (gamma * ppfd + 1.0d0)

    fEsub = min(fWsub,fDsub)
  fW = fWsub
    fD = fEsub  !!!!here it should be reported fDsub 
  fE = fEsub

    gpp380 = beta *  ppfd *  fAPAR * fSsub * fLsub * fEsub
    call fCO2_model_mean(CO2,bCO2,fCO2)
    gpp = gpp380 * fCO2

    
    ! // if (LOGFLAG > 1.5) 
      ! // fprintf(flog, 
        ! // "   gpp(): Modifiers: fAPAR %lf\tfSsub %lf\t fLsub %lf\t fDsub %lf\t fWsub %lf\tfEsub %lf\t fCO2 %lf\n                    gpp380 %lf\t gpp %lf\n",
        ! // fAPAR, fSsub, fLsub, fDsub, fWsub, fEsub, fCO2, *gpp380, *gpp);



    ! /* This has been removed, and simpler multiplicative CO2 modifier to gpp380 is used.
    ! * CO2 effect not only influences fD but also fW, due to stomatal action
    
    ! fDCO2sub = fDsub * pow(exp(-0.4 * D),
         ! fCO2_VPD_exponent(CO2, GPP_par.xCO2)) / exp(-0.4 * D) ;
    ! fWCO2sub = fWsub * pow(fWsub, fCO2_VPD_exponent(CO2, GPP_par.xCO2));

    ! if (LOGFLAG > 1.5) 
      ! fprintf(flog, 
        ! "   gpp(): Modifier values for GPP at %lf\n      fD=%lf\tfW=%lf\tfCO2mean=%lf\n", 
        ! CO2, fDCO2sub, fWCO2sub, fCO2_model_mean(CO2, GPP_par.bCO2));
    
    ! if (fDCO2sub > fWCO2sub) fECO2sub=fWCO2sub; else fECO2sub = fDCO2sub;
    
    ! *fECO2 = fECO2sub;
    
    ! *gpp = GPP_par.beta *  ppfd *  fAPAR * fSsub * fLsub  * fECO2sub * 
      ! fCO2_model_mean(CO2, GPP_par.bCO2);
! */

end subroutine


!/* Rain is snow below T > 0 C, and snow melts above O C. */
subroutine SnowFun(T, rain, snow, SnowMelt, &
  MeltCoeff, T_0 , CWmax,SnowThreshold)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: rain, snow, SnowMelt
   !inputs
   real(8), intent(in) :: T
   !parameters
   real(8), intent(in) :: MeltCoeff, T_0, CWmax,SnowThreshold
   !other variables
   real(8) :: NewSnow

   if (T < SnowThreshold) then
    NewSnow = rain 
    rain = 0
   else 
    NewSnow=0
   endif
  
  if (T > T_0) then  
    SnowMelt = MeltCoeff*(T-T_0)
  else 
  SnowMelt=0
  endif
  
  if(snow + NewSnow - SnowMelt < 0) then
    SnowMelt=NewSnow + snow
    snow =0;
  else
    snow = snow + NewSnow - SnowMelt
  endif
  
end subroutine

!/*Interception is a fraction of daily rainfall, fraction depending on fAPAR*/
subroutine interceptionfun(rain, intercepted, Temp, fAPAR, &
  I_0 , SnowThreshold)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: rain, intercepted
   !inputs
   real(8), intent(in) :: Temp, fAPAR
   !parameters
   real(8), intent(in) :: I_0, SnowThreshold
   
   if(Temp > SnowThreshold)  then
    intercepted = rain * (I_0 * fAPAR / 0.75) 
    rain = rain - intercepted
   else
    intercepted = 0
  endif
endsubroutine



! /* Soil water balance is updated with snowmelt and canopy throughfall
 ! * and evapotranspiration. No drainage occurs below field capacity */
subroutine swbalance(theta, throughfall, snowmelt, et, &
               tauDrainage,ThetaFC,soildepth, &
         drainage, snow, canw, &
         MeltCoeff, T_0 , CWmax,SnowThreshold)

 implicit none

   !state variables and outputs
   real(8), intent(inout) :: theta, throughfall, SnowMelt,et
   !inputs
   real(8), intent(inout) :: drainage,snow,canw
   !parameters
   real(8), intent(in) :: tauDrainage, ThetaFC, soildepth
   real(8), intent(in) :: MeltCoeff, T_0, CWmax,SnowThreshold
   !other variables
   real(8) :: st0, etfromvegandsoil=0

  ! /* Evaporate first from wet canopy and snow on ground */

  if (CWmax > 0.00000001) then 
    if ( (canw + snow - et) > 0 ) then             
      if ( (canw - et) > 0 ) then 
    canw = canw -et
    etfromvegandsoil = 0
      elseif (canw - et < 0) then !// in this case, there's enough snow left
    snow = snow + canw - et
    canw = 0
    etfromvegandsoil = 0
      endif    
    else
      etfromvegandsoil = et - canw - snow
      canw=0.0
      snow = 0.0
    endif

  else  
    if ( (snow - et) > 0 ) then             
      snow = snow - et
      etfromvegandsoil = 0
    elseif (snow - et < 0) then !// in this case, there's enough snow left
      etfromvegandsoil = et - snow
      snow = 0
    else
      snow = 0.0
    endif
  endif

  et = etfromvegandsoil

  ! /* Water balance without drainage */
  st0 = theta + throughfall + snowmelt  - et
  if (st0 <= 0) st0 = 0.0001

  ! /* Calculate what is left to drainage after partial balance update above: */
  if (tauDrainage > 0) then  


    ! // Simple time delay drainage above FC:
    if (st0 > ThetaFC * soildepth) then 
      drainage = (st0 - ThetaFC * soildepth) /tauDrainage      
    else
      drainage = 0
    endif
    theta = st0 - drainage

    ! /* Include marginal drainage below FC.
     ! * This was needed for model calibration only, below FC drainage
     ! * was practically zero, but important for convergence */
    ! /*
! if (st0 > sitepar.ThetaFC * sitepar.soildepth) {
      ! *drainage = (st0 - sitepar.ThetaFC * sitepar.soildepth) / 
  ! sitepar.tauDrainage;      
    ! }
    ! if (*drainage < (sitepar.ThetaFC * sitepar.soildepth - 
         ! sitepar.ThetaPWP * sitepar.soildepth) / 
  ! 10000) //pow(sitepar.tauDrainage, 5)) 
  ! *drainage = (sitepar.ThetaFC * sitepar.soildepth - 
         ! sitepar.ThetaPWP * sitepar.soildepth) / 
    ! 10000; //pow(sitepar.tauDrainage, 5);
    
    ! if (st0 <= sitepar.ThetaFC * sitepar.soildepth && 
  ! st0 > sitepar.ThetaPWP * sitepar.soildepth) { 
      ! *drainage = (st0 - sitepar.ThetaPWP * sitepar.soildepth) / 
  ! 10000; //pow(sitepar.tauDrainage, 5);      
      ! *theta = st0 - *drainage;
    ! }
  
    ! if (st0 <= sitepar.ThetaPWP * sitepar.soildepth) {
      ! *drainage = 0;
      ! *theta = st0;
    ! }
    ! *theta = st0 - *drainage;
    ! */
    ! //****************************************************** */

  endif
  
endsubroutine





! /* Estimate Evapotranspiration according to a simple empirical model
 ! * that uses GPP prediction to calculate transpiration, as driven
 ! * by VPD. Evaporation is estimated with PPFD, which is a surrogate
 ! * for Rnet */
subroutine ETfun(D, theta, ppfd, fAPAR, T, &
             beta, kappa,chi,soilthres,nu, &
       soildepth, ThetaFC,ThetaPWP, &
             canw, fE, A, fWgpp, fCO2mean, & 
       xCO2, & !//double fCO2mean, 
       CO2, etmodel, transp, evap, fWE, et) 

  
 implicit none

   !state variables,inputs and outputs
   real(8), intent(inout) :: D, theta, ppfd, fAPAR, T
   real(8), intent(inout) :: canw, fE, A, fWgpp, fCO2mean
   real(8), intent(inout) :: CO2, transp, evap, fWE, et
   !inputs
   integer, intent(in) :: etmodel
   !parameters
   real(8), intent(in) :: beta, kappa,chi,soilthres,nu
   real(8), intent(in) :: soildepth, ThetaFC,ThetaPWP
   real(8), intent(in) :: xCO2
   !other variables
   real(8) :: st0, etfromvegandsoil=0.d0,thetavol,REW,fWsub
   real(8) :: lambda, psychom, s,cp,MWratio,pressure
  
  thetavol = theta/soildepth
  REW=(thetavol-ThetaPWP)/(ThetaFC-ThetaPWP)
  ! //  double fEsub = -999; /* Minimum of fW and fD returned if ET-model 
  ! //      * flag indicates similar modifier as for GPP */
  fWsub=1.d0
  ! //  double fDsub=1;
  cp = 1003.5 ! // J/(kg K) (nearly constant, this is dry air on sea level)
  MWratio = 0.622 ! // Ratio of molecular weigths of water vapor and dry air;
  ! // double R = 287.058; // J/(kg K) Specific gas constant for dry air, wiki
  ! // double zh, zm, d, zom, zoh;
  ! /*If pressure is not inputted use default */
  pressure = 101300. ! // Pa  

  call fCO2_ET_model_mean(CO2,xCO2,fCO2mean)

  ! // rho=pressure/(R * (T+273.15) ); // Dry air density, kg/m3
  lambda = (-0.0000614342 * T**3 + 0.00158927 * T**2 - & 
      2.36418 * T +  2500.79) * 1000. ! // J/kg
  psychom= cp * pressure / (lambda * MWratio) ! // Pa/C, wiki
  s = 1000. * 4098.0 * (0.6109 * exp((17.27 * T)/(T+237.3))) / & 
    ((T+237.3)**2)  !  // Pa/C! (Ice has nearly the same slope)


  ! /* Calculate soil constraint, simple linear following Granier 1987*/
  if (soilthres < -998) then ! /*-999 omits water control*/
    fWsub = 1
  else
    if (REW < soilthres) then
      if (REW > 0.01) then
     fWsub = REW/soilthres
    else
      fWsub = 0.0
    endif
    else
      fWsub = 1.0;
    endif
  endif

  
  ! /* If there is any water in canopy, evaporation is not reduced by
   ! * low soil water */
  if (canw > 0.00000001) fWsub = 1

  !//  if (fDsub > fWsub) fEsub=fWsub; else fEsub = fDsub;     
  
  fE = fWsub   
  fWE = fWsub

  if (D < 0.01) D=0.01

  ! // if (LOGFLAG > 1.5)
    ! // fprintf(flog, "   ETfun(): CO2mean=%lf\tat CO2=%lf\n", 
      ! // fCO2mean, CO2);
  
  if (etmodel == -1) then
    transp = D * beta*A/(D**kappa) * &
      (fWgpp**nu) * fCO2mean!// ET differently sensitive to soil water than GPP
    evap = chi *  (1-fAPAR) *  fWsub * ppfd
    et = (transp + evap) * s / (s + psychom) 
  endif

  if (etmodel == 0) then
    transp = D * beta*A/(D**kappa) * &
      (fWgpp**nu) * &!// ET differently sensitive to soil water than GPP
      fCO2mean
    evap = chi * s / (s + psychom) * (1-fAPAR) *  fWsub * ppfd
    ! //    et = D * ET_par.beta*A/pow(D, ET_par.kappa) *
    ! //  pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
    ! //  fCO2mean +  // Mean effect of CO2 on transpiration
    ! //  ET_par.chi *  s / (s + psychom) * (1-fAPAR) *  fWsub * ppfd;
    et = transp + evap
  endif
  if (etmodel == 1) then
    transp = D * beta*A/(D**kappa) * &
      (fWgpp**nu) * &!// ET differently sensitive to soil water than GPP
      fCO2mean
    evap = chi * (1-fAPAR) *  fWsub * ppfd
    ! //et = D * ET_par.beta*A/pow(D, ET_par.kappa) *
    ! //  pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
    ! //  fCO2mean +  // Mean effect of CO2 on transpiration
    ! //  ET_par.chi * (1-fAPAR) *  fWsub * ppfd;
    et = transp + evap
  endif
  if (etmodel == 2) then
      et = D * (1 + beta/(D**kappa)) * A / CO2 * & 
  (fWgpp**nu) * & !// ET differently sensitive to soil water than GPP
  fCO2mean + & !// Mean effect of CO2 on transpiration      
  chi * (1-fAPAR) *  fWsub * ppfd
  endif

  ! // if (LOGFLAG > 2.5)
    ! // fprintf(flog, "      ETfun(): Model=%d\nD\t%lf\nET_par.beta\t%lf\nA\t%lf\npow(D, ET_par.kappa)\t%lf\npow(fWgpp, ET_par.nu)\t%lf\nfWgpp\t%lf\nET_par.nu\t%lf\nfCO2mean\t%lf\nCO2\t%lf\nET_par.chi\t%lf\ns/(s+psychom)\t%lf\n1-fAPAR\t%lf\nfWsum\t%lf\nppfd\t%lf\n-->et\t%lf\n",      
      ! // etmodel, D, ET_par.beta, A, pow(D, ET_par.kappa), 
      ! // pow(fWgpp, ET_par.nu), fWgpp, ET_par.nu,
      ! // fCO2mean, 
      ! // CO2,
      ! // ET_par.chi ,   s / (s + psychom), 1-fAPAR, fWsub,  ppfd, et);

endsubroutine