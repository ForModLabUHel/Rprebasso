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
 !		  // double *fE, FILE *flog, int LOGFLAG) {
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
