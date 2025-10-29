

subroutine initBiomasses(pCrobas,initVar,siteType,biomasses,nVar,nPar) 
  implicit none
  
  integer, intent(in) :: nVar,nPar
  real (kind=8), parameter :: pi = 3.1415927
  real (kind=8), intent(in) :: pCrobas(npar),initVar(7), siteType
  real (kind=8), intent(inout) :: biomasses(nvar)
  real (kind=8) :: par_betab, par_x, par_beta0, par_betas, par_mf, par_mr
  real (kind=8) :: par_mw, par_alfar, par_c, par_rhof, par_rhor , par_rhow
  real (kind=8) :: par_S_branchMod, gammaC, Tbd
  real (kind=8) :: A, ba, d, N, h, hc, B, Lc, betab, beta0, beta1, beta2, betaC, V
  real (kind=8) :: wf_STKG, W_froot, W_wsap, W_c, W_s, W_branch, W_croot, Wdb, W_stem, Wsh
  real (kind=8) :: W_crs, W_crh
    real (kind=8) :: age_factor, par_fAa, par_fAb, par_fAc
   
  !initBiomasses = function(pCrobas,initVarX){
  ! initVarX<-as.matrix(initVarX) change vector to matrix when maxlayer=1
  ! siteType = siteInfo(3)
  !##set parameters
  par_betab = pCrobas(13)
  par_x = pCrobas(19)
  par_beta0 = pCrobas(12)
  par_betas = pCrobas(14)
  par_mf = pCrobas(8)
  par_mr = pCrobas(9)
  par_mw = pCrobas(10)
  par_c = pCrobas(7)
  par_rhof = pCrobas(15)
  par_rhow = pCrobas(2)
  par_S_branchMod = pCrobas(27)
  par_fAa = pCrobas(45)
  par_fAb = pCrobas(46)
  par_fAc = pCrobas(47)
  gammaC = 0. !#initVarX(8,)
  Tbd = 10 !#####to include in the parameters
  
  !###set variables
  A = initVar(7)
  ba = initVar(5); d = initVar(4)
  N = ba/(pi*((d/2/100)**2))
  h = initVar(3); hc = initVar(6)
  B = ba/N
  Lc = h - hc
  
  age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - h)/par_fAc)))/par_fAa
  par_alfar = pCrobas(int(20+min(siteType,5.))) * age_factor
  par_rhor = par_alfar * par_rhof
  beta0 = par_beta0 * age_factor
  
  betab =  par_betab * Lc**(par_x-1)
  beta1 = (beta0 + betab + par_betas) 
  beta2 = 1. - betab - par_betas     
  betaC = (beta1 + gammaC * beta2) / par_betas
  wf_STKG = par_rhof * A * N
  W_froot = par_rhor * A * N  
  W_wsap = par_rhow * A * N * (beta1 * h + beta2 * hc) 
  W_c = par_rhow * A * N * hc 
  W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
  W_branch =  par_rhow * A * N * betab * Lc !branches biomass
  Wsh = max((A+B+sqrt(A*B)) * hc * par_rhow * N/2.9 - W_c,0.) !#initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
  !#initialize Wdb dead branches biomass
  if(par_S_branchMod == 1.) then
  Wdb = Tbd * W_branch * ((0.0337+0.000009749*N)*exp(-0.00456*d**2)+0.00723)
  else
    Wdb = Tbd * W_branch *((-0.00513+0.000012*N)*exp((0.00000732-0.000000764*N)*d**2)+0.00467)
  endif
  W_stem = W_c + W_s + Wsh
  W_crs = par_rhow * beta0 * A * H * N
  W_crh = Wsh * beta0
  W_croot = W_crh + W_crs! max(0.,(par_rhow *Lc * beta0 * A * N + (W_c + Wsh) * beta0)) !#coarse root biomass

  V = W_stem / par_rhow
  
  biomasses(33) = wf_STKG
  biomasses(25) = W_froot
  biomasses(47) = W_wsap
  biomasses(48) = W_c
  biomasses(49) = W_s
  biomasses(24) = W_branch
  biomasses(32) = W_croot
  biomasses(50) = Wsh
  biomasses(51) = Wdb
  biomasses(31) = W_stem
  biomasses(30) = V
  biomasses(54) = W_crh
  
end subroutine

SUBROUTINE Ffotos2(STAND_all,nClass,nSp,pCrobas,nVar,nPar,MeanLight,coeff,qcTOT)
implicit none

 integer, intent(in) :: nclass,nSp,nVar,nPar

!*****************************************************************************************
 real (kind=8), intent(in) :: pCrobas(npar,nSp)
 real (kind=8), intent(inout) :: STAND_all(nVar,nclass)
 real (kind=8), intent(inout) :: coeff(nclass) , qcTOT
!****************************************************************************************
 integer  :: ki
 real (kind=8) :: param(nPar), ln2 = 0.693147181
 real (kind=8) :: ht(nclass),hc(nclass),h(nclass)
 real (kind=8) :: LAIe(nclass),qc(nclass),btc(nclass),LAI(nclass),N(nclass)
 real (kind=8) :: l(2*nclass),vrel(2*nclass,nclass)
 real (kind=8) :: lpt(2*nclass,nclass),lt(2*nclass)
 real (kind=8) :: bt(2*nclass), k(nclass), par_betab(nclass), rc(nclass)
 real (kind=8) :: kLAIetot, kLAItot, Atot
 real (kind=8), intent(inout) :: MeanLight(nclass)
 real (kind=8) :: x1,x2,apuJ,apuI,par_sla,par_sla0,par_tsla,age
     integer :: iclass,i2,i1,species,nv        !!**!! nv defined as integer
       integer :: i, j, ii(2*nclass), iapu
 real (kind=8) :: apu, b1,  qctot0, qctot1, wwx, dc, e1
!****************************************************************************************
 real (kind=8) :: pi = acos(-1.)
 
 MeanLight = 0.
 coeff = 0.
 qcTOT = 0.
 
 do i = 1,nclass
   species = int(stand_all(4,i))
   if(species==0) then
     ht(i) = 0.   ! H
     hc(i) = 0.   ! Hc
     h(i) = ht(i) - hc(i)        ! Lc
     LAIe(i) = 0. ! leff
     k(i) = 0.               ! k 
     LAI(i) = 0. * par_sla / 10000.   ! WF_stand * sla
     ! par_betab(i) = PARAM(17)   ! betab
     rc(i) = 0./2.         ! rc crown radius
     N(i) = 0. / 10000.   ! N per m2
   else
     param = pCrobas(:,species)
     qc(i) = 0.

     par_sla = param(3)
     par_sla0 = param(39)
     par_tsla = param(40)
     age = STAND_all(7,i)
     if(par_tsla .gt. 0.) then
      par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)
     else
      par_sla = par_sla
     endif

   
     ht(i) = STAND_all(11,i)   ! H
     hc(i) = STAND_all(14,i)   ! Hc
     h(i) = ht(i) - hc(i)        ! Lc
     LAIe(i) = STAND_all(19,i) ! leff
     k(i) = PARAM(4)               ! k 
     LAI(i) = STAND_all(33,i) * par_sla / 10000.   ! WF_stand * sla
     ! par_betab(i) = PARAM(17)   ! betab
     rc(i) = STAND_all(15,i)/2.         ! rc crown radius
     N(i) = STAND_all(17,i) / 10000.   ! N per m2
   endif
 end do
       
       nv= 2*nclass

do  i = 1, nv
    ii(i) = i
end do
    

! **  sort tree tops and crown bases in descending order into vector l

do  i=1,nclass
  l(i) = ht(i)
  l(i+nclass) = hc(i)
end do 
        

do  i=1,nv-1
  do j=i+1,nv
          if(l(i).lt.l(j)) then 
             apu = l(i)
             l(i) = l(j)
             l(j) = apu

!  ii-table sorts the l-table indeces so that later the corresponding "locations" for hc and ht values can be located
          iapu = ii(i)
          ii(i) = ii(j)
          ii(j) = iapu
        endif
  end do
end do
        
! ** end sort
! ** calculate effective leaf area for each species in canopy layers determined by heights and hc:s
! ** use function wwx to calculate foliage distribution, defined by species

lt(1) = 0.
bt(1) = 0.
do i=1,nv-1
        
  lt(i+1) = 0.
   do j=1,nclass
    species = j
    apuJ = wwx(0.0d+0,1.0d+0,ht(j)-hc(j),species)
    if(l(i).gt.hc(j).and.l(i+1).lt.ht(j)) then
    if((ht(j)-hc(j)).gt.0.) then
                   x1 = (ht(j)-l(i))/(ht(j)-hc(j))
                   x2 = (ht(j)-l(i+1))/(ht(j)-hc(j))
                  else 
                   x1=0.
                   x2=0.
                  endif
                  apuI = wwx(x1,x2,ht(j)-hc(j),species)
             vrel(i+1,j) = apuI / apuJ
                else
             vrel(i+1,j) = 0.
                endif     
                lpt(i+1,j) = k(j) * LAIe(j) * vrel(i+1,j)
                lt(i+1) = lt(i+1) + lpt(i+1,j)
  end do
   bt(i+1) = bt(i) + lt(i+1)
end do
    

do j=1,nclass
       dc = 0.
        i1 = 0
        i2 = 0
   do i=1,nv
        if(ht(j).eq.l(i)) i1=i
      if(hc(j).eq.l(i)) i2=i

      if(ii(i)==j) i1 = i
      if(ii(i) == j+nclass) i2 = i

!        if(ht(j).gt.l(i).and.hc(j).le.l(i)) dc=dc+lt(i)
  end do
    e1 = exp(-bt(i1)) - exp(-bt(i2))
              b1 = bt(i2) - bt(i1)

     if (b1 .ne. 0)  qc(j) = k(j) * laie(j) * e1  / b1 

       btc(j) = bt(i2)

!           MeanLight(j) = 0.5 * (exp(-bt(i1)) + exp(-bt(i2)))
            MeanLight(j) = exp(-bt(i2))
end do
        



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  **  add a new more simple way of dividing between size classes
!       
!      Here the fAPAR for the whole stand is calculated 
!      (stored in qc(in), same for all classes), and trees 
!      utilise this in proportion to their foliage mass
! 
!      AM 2.7.2008
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
kLAIetot = 0.
kLAItot = 0.
Atot = 0.
qcTOT1 = 0
do  i =1,nclass
  kLAIetot = kLAIetot + k(i) * LAIe(i)
    kLAItot = kLAItot + k(i) * LAI(i)
    Atot = Atot + N(i) * pi*(rc(i)**2 )
  qcTOT1 = qcTOT1 + qc(i)
end do
    
! calculate LPJ style fAPAR and use the smaller of the two
  
     if(Atot > 0. ) then
         qcTOT0 = (1. - exp(-kLAItot/ Atot)) * Atot
     else
         qctot = 0.
     endif
     
     qcTOT = min(qcTOT0,qctot1)
   qcTOT = min(qcTOT,1.)
   qcTOT = max(qcTOT,0.)
!    qctot = qctot1
     
!     if(stand_P(7) > 150) then
!        continue
!     endif
     
     
! calculate weights - on the basis of qcTOT1 but all downscaled if qcTOT0 < qcTOT1
!         
!do i = 1,nclass
     
   if(qcTOT1.gt.0.) then
    coeff  = qc / qcTOT1 * qcTOT / qcTOT1 ! weight

!        coeff_SP  = qc(2) / qcTOT1 * qcTOT / qcTOT1 ! 

!        coeff_B  = qc(3) / qcTOT1 * qcTOT / qcTOT1  ! 



!!!!FMadded
!        qcTOT = qcTOT * (coeff_P+coeff_SP+coeff_B)

!        coeff_P = coeff_P/(coeff_P+coeff_SP+coeff_B)
!        coeff_SP = coeff_SP/(coeff_P+coeff_SP+coeff_B)
!        coeff_B = coeff_B/(coeff_P+coeff_SP+coeff_B)
!!!!


!    if(kLAIetot.gt.0.) then
!          coeff_P  = k(1)*LAIe(1) / kLAIetot ! weight

!          coeff_SP  = k(2)*LAIe(2) / kLAIetot ! 
           
!             coeff_B  = k(3)*LAIe(3) / kLAIetot ! 


!   else
!       coeff_P = 1./3.
!       coeff_SP = 1./3.
!       coeff_B = 1./3.
   
! end do  
    endif      



 !     write(60,*)qcTOT, qcTOT1, qc(1), qc(2), qc(3)

!81      continue


  end subroutine Ffotos2



!***************************************************************
!  WWX
!
!  Alkuper?inen Annikki M?kel?
!
!  MUUTOKSIA (Sanna H?rk?nen, 14.8.2002):
!  -Lis?tty muuttujien m??rittelyt
!
!  MUUTOKSIA (Annikki M?kel? 19.11.2004):
!  -Muutettu funktiokutsun argumentteja (hc, species)
!  -M?nnyn malli ennallaan kuuselle uusi formulointi
!  -Kuusella latvan maksimi aina annetulla kohdalla, <= 5m
!
!***************************************************************
    real(kind=8)  function wwx(x1,x2,Lc,species)

  implicit none
    
    real (kind=8), parameter :: hmax0 = 5., p = 1., q = 1.

  real (kind=8) :: x1,x2,Lc
  integer :: species

  real (kind=8) :: pp,qq
  integer :: N,i
  real (kind=8) :: x,dx,apu,a,b,c,d,w



!  if(species==1 .OR. species==3)then
    pp = p
    qq = q
!  endif


!!!check with Annikki
!  if(species==4)then
!!     hmax = amin1(0.9*Lc,hmax0+0.3*Lc)
!    pp = p
!! *** If-lause lis?tty 2011/10/14 by TL
!    if (Lc .gt. hmax0) then
!           qq = 0.18 * Lc - 0.6
!          else
!          qq = 0.18 * hmax0 - 0.6
!          endif
!  endif

  N = max(1,int((x2-x1)*10 + 0.5))
      dx = (x2-x1)/float(N)

      w = 0.
      x = x1

      do  i = 1,N
      if(x+dx>1.)x=1.-dx
    A = x**pp * (1.-x)**qq
    B = (x+dx/2.)**pp * (1.-x-dx/2.)**qq
    C = B
    apu = max(0.,1.-x-dx)
    D = (x+dx)**pp * apu**qq
    w = w + (A+2.*B+2.*c+D)/6.
    x = x+dx
      end do
      
      wwx = w*dx
      end function wwx
!*************************************************************



subroutine preles(weather,DOY,fAPAR,prelesOut,pars,GPP,ET,SW,etmodel,CO2model)!,p0)

implicit none

 INTERFACE
   SUBROUTINE call_preles( &
        PAR, TAir, VPD, Precip, CO2, fAPARc, &  !inputs
         GPPmeas, ETmeas, SWmeas, &
                GPP, ET, SW, SOG, fS, fD, fW, fE, & !outputs
    Throughfall, Interception, Snowmelt, Drainage, &
    Canopywater, S, &
      soildepth,ThetaFC, ThetaPWP, tauDrainage, beta, & !parameters
    tau,S0,Smax,kappa,gamma, soilthres, bCO2, xCO2, &
    ETbeta, ETkappa, ETchi,ETsoilthres,ETnu, MeltCoef, &
    I0,CWmax,SnowThreshold,T_0,SWinit,CWinit,SOGinit, &
    Sinit,t0,tcrit,tsumcrit, &  
    etmodel, LOGFLAG,NofDays, &
    day, &!!!!this is DOY
    transp, evap, fWE,CO2model) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_CHAR, C_PTR, C_DOUBLE  
     real ( C_DOUBLE ) :: PAR(365), TAir(365), VPD(365), Precip(365), CO2(365), fAPARc(365)
     real ( c_double ) :: GPPmeas(365), ETmeas(365), SWmeas(365)
     real ( c_double ) :: GPP(365), ET(365), SW(365), SOG(365), fS(365), fD(365), fW(365), fE(365)
     real ( c_double ) :: Throughfall(365), Interception(365), Snowmelt(365), Drainage(365)
     real ( c_double ) :: Canopywater(365), S(365)
     real ( c_double ) :: soildepth,ThetaFC, ThetaPWP, tauDrainage, beta !parameters
     real ( c_double ) :: tau,S0,Smax,kappa,gamma, soilthres, bCO2, xCO2
     real ( c_double ) :: ETbeta, ETkappa, ETchi,ETsoilthres,ETnu, MeltCoef
     real ( c_double ) :: I0,CWmax,SnowThreshold,T_0,SWinit,CWinit,SOGinit
     real ( c_double ) :: Sinit,t0,tcrit,tsumcrit
     integer(c_int) :: etmodel, LOGFLAG,NofDays,CO2model
     integer(c_int) :: day(365)
     real ( c_double ) :: transp(365), evap(365), fWE(365)
   END SUBROUTINE call_preles
 END INTERFACE

 real (kind=8), intent(inout) :: weather(365,5),fAPAR(365)
 real (kind=8), intent(out) :: prelesOut(16)!,p0
 real (kind=8), intent(inout) :: pars(30)
 integer, intent(in):: DOY(365), etmodel,CO2model

     real (kind=8) PAR(365), TAir(365), VPD(365), Precip(365), CO2(365), fAPARc(365)
     real (kind=8) GPPmeas(365), ETmeas(365), SWmeas(365)
     real (kind=8), intent(inout) :: GPP(365), ET(365),SW(365)
     real (kind=8) SOG(365), fS(365), fD(365), fW(365), fE(365)
     real (kind=8) Throughfall(365), Interception(365), Snowmelt(365), Drainage(365)
     real (kind=8) Canopywater(365), S(365)
     real (kind=8) :: soildepth,ThetaFC, ThetaPWP, tauDrainage, beta !parameters
     real (kind=8) :: tau,S0,Smax,kappa,gamma, soilthres, bCO2, xCO2
     real (kind=8) :: ETbeta, ETkappa, ETchi,ETsoilthres,ETnu, MeltCoef
     real (kind=8) :: I0,CWmax,SnowThreshold,T_0,SWinit,CWinit,SOGinit
     real (kind=8) :: Sinit,t0,tcrit,tsumcrit
     integer LOGFLAG,NofDays
     integer day(365),nDays
     real (kind=8) transp(365), evap(365), fWE(365)!,fAPAR0

!init inputs
PAR = weather(:,1)
TAir = weather(:,2)
VPD = weather(:,3)
Precip = weather(:,4)
CO2 = weather(:,5)

day = DOY
NofDays = 365
!fAPAR0 = 1

GPPmeas(:) = 0.
ETmeas(:) = 0.
SWmeas(:) = 0.
!etmodel=0.
LOGFLAG=0

!if(PAR(366)==-999) then
 nDays = 365 
!else
! nDays = 366
!endif

!init preles parameters
soildepth = pars(1)
ThetaFC = pars(2)
ThetaPWP = pars(3)
tauDrainage = pars(4)
beta = pars(5)
tau = pars(6)
S0 = pars(7)
Smax = pars(8)
kappa = pars(9)
gamma = pars(10)
soilthres = pars(11)
bCO2 = pars(12)
xCO2 = pars(13)
ETbeta = pars(14)
ETkappa = pars(15)
ETchi = pars(16)
ETsoilthres = pars(17)
ETnu = pars(18)
MeltCoef = pars(19)
I0 = pars(20)
CWmax = pars(21)
SnowThreshold = pars(22)
T_0 = pars(23)
SWinit = pars(24)
CWinit = pars(25)
SOGinit = pars(26)
Sinit = pars(27)
t0 = pars(28)
tcrit = pars(29)
tsumcrit = pars(30)


fAPARc = fAPAR
 call call_preles( &
        PAR, TAir, VPD, Precip, CO2, fAPARc, &  !inputs
         GPPmeas, ETmeas, SWmeas, &!end inputs
                GPP, ET, SW, SOG, fS, fD, fW, fE, & !outputs
    Throughfall, Interception, Snowmelt, Drainage, &
    Canopywater, S, & !end outputs
      soildepth,ThetaFC, ThetaPWP, tauDrainage, beta, & !parameters
    tau,S0,Smax,kappa,gamma, soilthres, bCO2, xCO2, &
    ETbeta, ETkappa, ETchi,ETsoilthres,ETnu, MeltCoef, &
    I0,CWmax,SnowThreshold,T_0,SWinit,CWinit,SOGinit, &
    Sinit,t0,tcrit,tsumcrit, & !end parameters  
    etmodel, LOGFLAG, NofDays, &
    day, &!!!!this is DOY
    transp, evap, fWE,CO2model)

 call SMIfromPRELES(GPP,fW,prelesOut(7),sum(fAPARc)/365.)

prelesOut(1) = sum(GPP(1:nDays))
prelesOut(2) = sum(ET(1:nDays))
prelesOut(3) = SW(nDays)
prelesOut(4) = SOG(nDays)
prelesOut(5) = fS(nDays)
prelesOut(6) = fD(nDays)
! prelesOut(7) = fW(nDays)
prelesOut(8) = fE(nDays)
prelesOut(9) = Throughfall(nDays)
prelesOut(10) = Interception(nDays)
prelesOut(11) = Snowmelt(nDays)
prelesOut(12) = Drainage(nDays)
prelesOut(13) = Canopywater(nDays)
prelesOut(14) = S(nDays)
prelesOut(15) = sum(SW(1:nDays))/nDays
prelesOut(16) = sum(SW(152:243))/92

end subroutine


SUBROUTINE mod5c(theta,time,climate,init,b,d,leac,xt,steadystate_pred)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
    !********************************************* &
    ! returns the model prediction xt for the given parameters
    ! 1-16 matrix A entries: 4*alpha, 12*p

    ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

    ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

    ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

    ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

    ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

    ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

    ! 33-35 Woody parameters: theta_1, theta_2, r

    REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
    REAL (kind=8),INTENT(IN) :: time,d,leac ! time,size,leaching
    REAL (kind=8),DIMENSION(3),INTENT(IN) :: climate ! climatic conditions
    REAL (kind=8),DIMENSION(5),INTENT(IN) :: init ! initial state
    REAL (kind=8),DIMENSION(5),INTENT(IN) :: b ! infall
    REAL (kind=8),DIMENSION(5),INTENT(OUT) :: xt ! the result i.e. x(t)
    REAL (kind=8),INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution 
    ! LOGICAL,OPTIONAL,INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution 
    ! in steady-state conditions (which sould give equal solution as if time is set large enough)
    REAL (kind=8),DIMENSION(5,5) :: A,At,mexpAt
    INTEGER :: i
    REAL (kind=8),PARAMETER :: pi = 3.141592653589793
    REAL (kind=8) :: tem,temN,temH,size_dep
    REAL (kind=8),DIMENSION(5) :: te
    REAL (kind=8),DIMENSION(5) :: z1,z2
    REAL (kind=8),PARAMETER :: tol = 1E-12
    LOGICAL :: ss_pred = .FALSE.

    ! IF(PRESENT(steadystate_pred)) THEN
        ! ss_pred = steadystate_pred
    ! ENDIF
    IF(steadystate_pred == 1.) THEN
        ss_pred = .true. 
  else
    ss_pred = .false.
    ENDIF

    !#########################################################################
    ! Compute the coefficient matrix A for the differential equation

    ! temperature annual cycle approximation
    te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
    te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
    te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
    te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

    tem = 0.0
    temN = 0.0
    temH = 0.0
    DO i = 1,4 ! Average temperature dependence
        tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
    END DO

    ! Precipitation dependence
    tem = tem*(1.0-EXP(theta(28)*climate(2)/1000.0))
    temN = temN*(1.0-EXP(theta(29)*climate(2)/1000.0))
    temH = temH*(1.0-EXP(theta(30)*climate(2)/1000.0))

    ! Size class dependence -- no effect if d == 0.0
    size_dep = MIN(1.0,(1.0+theta(33)*d+theta(34)*d**2.0)**(-ABS(theta(35))))

    ! check rare case where no decomposition happens for some compartments 
    ! (basically, if no rain)
    IF (tem <= tol) THEN
        xt = init + b*time
        return
    END IF

    ! Calculating matrix A (will work ok despite the sign of alphas)
    DO i = 1,3
        A(i,i) = -ABS(theta(i))*tem*size_dep
    END DO
    A(4,4) = -ABS(theta(4))*temN*size_dep
    
    A(1,2) = theta(5)*ABS(A(2,2))
    A(1,3) = theta(6)*ABS(A(3,3))
    A(1,4) = theta(7)*ABS(A(4,4))
    A(1,5) = 0.0 ! no mass flows from H -> AWEN
    A(2,1) = theta(8)*ABS(A(1,1))
    A(2,3) = theta(9)*ABS(A(3,3))
    A(2,4) = theta(10)*ABS(A(4,4))
    A(2,5) = 0.0
    A(3,1) = theta(11)*ABS(A(1,1))
    A(3,2) = theta(12)*ABS(A(2,2))
    A(3,4) = theta(13)*ABS(A(4,4))
    A(3,5) = 0.0
    A(4,1) = theta(14)*ABS(A(1,1))
    A(4,2) = theta(15)*ABS(A(2,2))
    A(4,3) = theta(16)*ABS(A(3,3))
    A(4,5) = 0.0
    A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
    DO i = 1,4
        A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
    END DO

    ! Leaching (no leaching for humus)
    DO i = 1,4
        A(i,i) = A(i,i)+leac*climate(2)/1000.0
    END DO

    !#########################################################################
    ! Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init
  IF(ss_pred) THEN
    ! Solve DE directly in steady state conditions (time = infinity)
    ! using the formula 0 = x'(t) = A*x + b => x = -A**-1*b
    CALL solve(-A, b, xt)
  ELSE
    ! Solve DE in given time
    z1 = MATMUL(A,init) + b
    At = A*time !At = A*t
    CALL matrixexp(At,mexpAt)
    z2 = MATMUL(mexpAt,z1) - b
    CALL solve(A,z2,xt) ! now it can be assumed A is non-singular
    ENDIF

    END SUBROUTINE mod5c

    !#########################################################################
    ! Functions for solving the diff. equation, adapted for the Yasso case
    SUBROUTINE matrixexp(A,B)
        IMPLICIT NONE
        ! Approximated matrix exponential using Taylor series with scaling & squaring
        ! Accurate enough for the Yasso case
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: B
        REAL (kind=8),DIMENSION(n,n) :: C,D
        REAL (kind=8) :: p,normiter
        INTEGER :: i,q,j
        q = 10 ! #terms in Taylor
        B = 0.0
        DO i = 1,n
            B(i,i) = 1.0
        END DO
        normiter = 2.0 ! Amount of scaling & squaring
        j = 1
        CALL matrixnorm(A, p)
        DO
            IF (p<normiter) THEN
                EXIT
            END IF
            normiter = normiter*2.0
            j = j+1
        END DO
        !write(*,*) normiter
        C = A/normiter ! scale
        B = B+C
        D = C
        DO i = 2,q ! compute Taylor expansion
            D = MATMUL(C,D)/REAL(i)
            B = B+D
        END DO
        DO i = 1,j ! square
            B = MATMUL(B,B)
        END DO
    END SUBROUTINE matrixexp

    SUBROUTINE matrixnorm(A,B)
        !returns elementwise (i.e. Frobenius) norm of a square matrix
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),INTENT(OUT) :: b
        INTEGER :: i
        b = 0.0
        DO i = 1,n
            b = b+SUM(A(:,i)**2.0)
        END DO
        b = SQRT(b)
    END SUBROUTINE matrixnorm


    SUBROUTINE solve(A, b, x)
        ! Solve linear system A*x = b
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: x
        REAL (kind=8),DIMENSION(n,n) :: U
        REAL (kind=8),DIMENSION(n) :: c
        INTEGER :: i

        ! transform the problem to upper diagonal form
        CALL pgauss(A, b, U, c)

        ! solve U*x = c via back substitution
        x(n) = c(n)/U(n,n)
        DO i = n-1,1,-1
            x(i) = (c(i) - DOT_PRODUCT(U(i,i+1:n),x(i+1:n)))/U(i,i)
        END DO
    END SUBROUTINE solve

    SUBROUTINE pgauss(A, b, U, c)
        ! Transform the lin. system to upper diagonal form using gaussian elimination
        ! with pivoting
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: U
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: c
        INTEGER :: k, j
        REAL,PARAMETER :: tol = 1E-12

        U = A
        c = b
        DO k = 1,n-1
            CALL pivot(U,c,k) ! do pivoting (though may not be necessary in our case)
            IF (ABS(U(k,k)) <= tol) THEN
                write(*,*) 'Warning!!! Matrix is singular to working precision!'
            END IF
            U(k+1:n,k) = U(k+1:n,k)/U(k,k)
            DO j = k+1,n
                U(j,k+1:n) = U(j,k+1:n) - U(j,k)*U(k,k+1:n)
            END DO
            c(k+1:n) = c(k+1:n) - c(k)*U(k+1:n,k)
        END DO
    END SUBROUTINE pgauss

    SUBROUTINE pivot(A, b, k)
        ! perform pivoting to matrix A and vector b at row k
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(INOUT) :: A
        REAL (kind=8),DIMENSION(n),INTENT(INOUT) :: b
        INTEGER,INTENT(IN) :: k
        INTEGER :: q, pk

        !write(*,*) 'Pivot elements are: ', A(k:n,k)
        q = MAXLOC(ABS(A(k:n,k)),1)
        !write(*,*) q
        IF (q > 1) THEN
            pk = k-1+q
            A(k:pk:pk-k,:) = A(pk:k:k-pk,:)
            b(k:pk:pk-k) = b(pk:k:k-pk)
        END IF
        !write(*,*) 'Pivot elements are: ', A(k:n,k)
    END SUBROUTINE pivot


  
    ! SUBROUTINE deadWoodV(y,nY,deadVol,dbh, pars)
        ! ! calculating deadwood volume decay
        ! IMPLICIT NONE
        ! INTEGER,intent(in) :: nY
    ! REAL (kind=8),intent(in) :: y(nY),dbh,pars(4)
    ! REAL (kind=8),intent(inout) :: deadVol(nY)
        ! !parameters
! !    REAL (kind=8) :: p1 = -2.653,p2 = -2.948,p3 = -3.324,p4 = .055,p5 = .059,p6 = .135,p7 = -0.03

    ! !###Gomprtz models
    ! deadVol = exp(-exp(pars(1) + pars(2)*y + pars(3)*dbh + pars(4)))
  ! END SUBROUTINE deadWoodV




! Note for Birch Betula pubenscens and brown leaves is used
    SUBROUTINE compAWENH(Lit,AWENH,parsAWEN)
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: AWENH
        REAL (kind=8),INTENT(IN) :: Lit,parsAWEN(4)
  AWENH(1) = parsAWEN(1)*Lit
  AWENH(2) = parsAWEN(2)*Lit
  AWENH(3) = parsAWEN(3)*Lit
  AWENH(4) = parsAWEN(4)*Lit
  AWENH(5) = 0.
    END SUBROUTINE compAWENH

  
SUBROUTINE runYasso(litter,litterSize,nYears, nLayers, nSites, nSp,&
        species, nClimID,climIDs,pAWEN,pYasso, weatherYasso,soilC)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION 
    !********************************************* &
    ! run yasso for some years with litterfal inputs from prebas.

  integer, intent(in) :: nYears, nLayers, nSites, nSp,nClimID
  REAL (kind=8),INTENT(IN) :: litter(nSites, nYears, nLayers, 3) !!!fourth dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  REAL (kind=8),INTENT(IN) :: weatherYasso(nClimID, nYears, 3)
  REAL (kind=8),INTENT(IN) :: species(nSites, nLayers),litterSize(3,nSp)
  REAL (kind=8),INTENT(IN) :: pAWEN(12, nSp), pYasso(35)
  real (kind=8),INTENT(inout) :: soilC(nSites,(nYears+1),5,3,nLayers)
  integer,INTENT(IN) :: climIDs(nSites)
  INTEGER :: year, site, layer, spec
  real (kind=8) :: t=1.,Lst,Lb,Lf,leac=0.,stSt=0. !leaching parameter for Yasso
  real (kind=8),DIMENSION(5) :: fbAWENH,folAWENH,stAWENH
  

fbAWENH = 0.
folAWENH = 0.
stAWENH = 0.

!!!!run Yasso
do site = 1, nSites
 do layer = 1,nLayers
  do year = 1,nYears

   Lst = litter(site,year,layer,3)
   Lb = litter(site,year,layer,2)
   Lf = litter(site,year,layer,1)

   spec = int(species(site,layer))
   call compAWENH(Lf,folAWENH,pAWEN(1:4,spec))   !!!awen partitioning foliage
   call compAWENH(Lb,fbAWENH,pAWEN(5:8,spec))   !!!awen partitioning branches
   call compAWENH(Lst,stAWENH,pAWEN(9:12,spec))         !!!awen partitioning stems

   call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,1,layer),stAWENH,litterSize(1,spec), &
  leac,soilC(site,(year+1),:,1,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,2,layer),fbAWENH,litterSize(2,spec), &
  leac,soilC(site,(year+1),:,2,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,3,layer),folAWENH,litterSize(3,spec), &
  leac,soilC(site,(year+1),:,3,layer),stSt)
  
  enddo
 enddo
enddo

END SUBROUTINE runYasso  

subroutine calWf(pars,Wf,inputs,nData,As)
 IMPLICIT NONE
 integer, intent(in) :: nData
 REAL (kind=8),INTENT(inOUT) :: Wf(nData,2),As(ndata,2) !!!Wf(:,1) Wf as function of As; 
                  !Wf(:,2) as function of Lc
 REAL (kind=8),INTENT(INout) :: pars(3),inputs(nData,3) !inputs col#1 = basal area; 
                !col#2=height; col#3 = height of crown base 
 REAL (kind=8) ba(ndata), h(ndata), hc(ndata), Lc(ndata) !!variables
 REAL (kind=8) par_rhof, par_ksi, par_z !!parameters
 
  ba = inputs(:,1)
  h = inputs(:,2)
  hc = inputs(:,3)
  Lc = h-hc
  As(:,1) = ba * Lc/(H-1.3)
  par_z = pars(1)
  par_rhof = pars(2)
  par_ksi = pars(3)
  As(:,2) = par_ksi/par_rhof * Lc ** par_z 
  Wf(:,1) = par_rhof * As(:,1)
  Wf(:,2) = par_ksi * Lc ** par_z 
END SUBROUTINE calWf


SUBROUTINE StstYasso(litter,litterSize, nLayers, nSites, nSp,species,nClimID,climIDs,pAWEN,pYasso,weatherYasso,soilC)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION 
    !********************************************* &
    ! run yasso for some years with litterfal inputs from prebas.

  integer, intent(in) ::  nLayers, nSites, nSp,nClimID
  REAL (kind=8),INTENT(IN) :: litter(nSites, nLayers, 3) !!!fourth dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  REAL (kind=8),INTENT(IN) :: weatherYasso(nClimID, 3)
  REAL (kind=8),INTENT(IN) :: species(nSites, nLayers),litterSize(3,nSp)
  REAL (kind=8),INTENT(IN) :: pAWEN(12, nSp), pYasso(35)
  real (kind=8),INTENT(inout) :: soilC(nSites,5,3,nLayers)
  integer,INTENT(IN) :: climIDs(nSites)
  INTEGER :: year, site, layer, spec
  real (kind=8) :: t=1.,Lst=0.,Lb=0.,Lf=0.,leac=0.,stSt=1. !leaching parameter for Yasso
  real (kind=8),DIMENSION(5) :: fbAWENH,folAWENH,stAWENH


fbAWENH = 0.
folAWENH = 0.
stAWENH = 0.
soilC = 0.

!!!!run Yasso
do site = 1, nSites
 do layer = 1,nLayers
  
   Lst = litter(site,layer,3)
   Lb = litter(site,layer,2)
   Lf = litter(site,layer,1)
   spec = int(species(site,layer))
  if(Lst>0) then
    call compAWENH(Lst,stAWENH,pAWEN(9:12,spec))         !!!awen partitioning stems
    call mod5c(pYasso,t,weatherYasso(climIDs(site),:),soilC(site,:,1,layer),stAWENH,litterSize(1,spec), &
      leac,soilC(site,:,1,layer),stSt)
  else
    soilC(site,:,1,layer) = 0.
   endif
   if(Lf>0) then
    call compAWENH(Lf,folAWENH,pAWEN(1:4,spec))   !!!awen partitioning foliage
    call mod5c(pYasso,t,weatherYasso(climIDs(site),:),soilC(site,:,3,layer),folAWENH,litterSize(3,spec), &
      leac,soilC(site,:,3,layer),stSt)
  else
    soilC(site,:,3,layer) = 0.
   endif
   if(Lb>0) then
    call compAWENH(Lb,fbAWENH,pAWEN(5:8,spec))   !!!awen partitioning branches
    call mod5c(pYasso,t,weatherYasso(climIDs(site),:),soilC(site,:,2,layer),fbAWENH,litterSize(2,spec), &
      leac,soilC(site,:,2,layer),stSt)
   else
    soilC(site,:,2,layer) = 0.
   endif
 enddo
enddo

END SUBROUTINE StstYasso

subroutine calW(pars,Wf,Wbr,Wstem,inputs,nData)
 IMPLICIT NONE
 integer, intent(in) :: nData
 REAL (kind=8),INTENT(OUT) :: Wf(nData,2) !!!Wf(:,1) Wf as function of As; Wf(:,2) as function of Lc
 REAL (kind=8),INTENT(OUT) :: Wbr(nData),Wstem(nData)
 REAL (kind=8) W_c(nData),W_s(nData),Wsh(nData)
 REAL (kind=8),INTENT(IN) :: pars(6),inputs(nData,3) !inputs col#1 = basal area;
!col#2=height; col#3 = height of crown base
 REAL (kind=8) ba(ndata), h(ndata), hc(ndata), Lc(ndata), As(ndata) !!variables
 REAL (kind=8) par_rhow, par_z, par_betab, par_betas, par_rhof, par_ksi!!parameters

ba = inputs(:,1)
h = inputs(:,2)
hc = inputs(:,3)
Lc = h-hc
As = ba * Lc/(H-1.3)

par_rhow = pars(1)
par_z = pars(2)
par_betab = pars(3)
par_betas = pars(4)
par_rhof = pars(5)
par_ksi = pars(6)

Wf(:,1) = par_rhof * As
Wf(:,2) = par_ksi * Lc ** par_z
Wbr =  par_rhow * As *  par_betab * Lc !branches biomass

W_c = par_rhow * As * hc !sapwood stem below Crown
  W_s = par_rhow * As * par_betas * Lc !sapwood stem within crown
Wsh = max((As+ba+sqrt(As*ba)) * hc * par_rhow /2.9 - W_c,0.0) !initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
Wstem = W_c + W_s + Wsh

END SUBROUTINE calW

!***************************************************************
!  tapioThin
!
!  subroutine to calculate the BA limits (ba_lim) and the to apply thinnings
!  and the BA after thinnings are applied (ba_thd)
!***************************************************************
subroutine tapioThin(forType,siteType,ETSmean,Hdom,tapioPars,baThin,BAthdPer, BAlimPer)

  implicit none
    real (kind=8),dimension(2) :: baThin
    real (kind=8) :: forType !1 for conifers; 2 for deciduous
  real (kind=8) :: siteType,ETSmean, Hdom !siteType; average ETS of the site, average height of the stand before thinning 
    real (kind=8) :: BA_lim, BA_thd, BA_limLow, BA_limUp, BA_thdLow, BA_thdUp
  real (kind=8) :: HthinStart,HthinLim, ETSlim, tapioPars(5,2,3,20) !!dimensions are: 1st=SiteType; 2nd = ForType; 3rd= ETS; 4th=nTapioPars
  real (kind=8) :: pX(3,20) !pX(1) = ETS threshold; pX(2)= Hlim;  pX(3:20) equation parameters
    real (kind=8) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16
    real (kind=8) :: BAthdPer, BAlimPer ! 1 for the upper limit, 0 for the lower limit

 pX = tapioPars(int(siteType), int(ForType),:,:)
 if(ETSmean > pX(1,1)) then !if we are in South Finland
  HthinStart =  pX(1,3)
  HthinLim =  pX(1,4)
    p1 = pX(1,5)
   p2 = pX(1,6)
   p3 = pX(1,7)
    p4 = pX(1,8)
   p5 = pX(1,9)
   p6 = pX(1,10)
   p7 = pX(1,11)
   p8 = pX(1,12)
   p9 = pX(1,13)
   p10 = pX(1,14)
   p11 = pX(1,15)
    p12 = pX(1,16)
   p13 = pX(1,17)
   p14 = pX(1,18)
   p15 = pX(1,19)
   p16 = pX(1,20)
 elseif(ETSmean <= pX(1,1) .and. ETSmean >= pX(1,2)) then !if we are in Central Finland
  HthinStart =  pX(2,3)
  HthinLim =  pX(2,4)
    p1 = pX(2,5)
   p2 = pX(2,6)
   p3 = pX(2,7)
    p4 = pX(2,8)
   p5 = pX(2,9)
   p6 = pX(2,10)
   p7 = pX(2,11)
   p8 = pX(2,12)
   p9 = pX(2,13)
   p10 = pX(2,14)
   p11 = pX(2,15)
    p12 = pX(2,16)
   p13 = pX(2,17)
   p14 = pX(2,18)
   p15 = pX(2,19)
   p16 = pX(2,20)
 else !if we are in Northern Finland
   HthinStart =  pX(3,3)
  HthinLim =  pX(3,4)
    p1 = pX(3,5)
   p2 = pX(3,6)
   p3 = pX(3,7)
    p4 = pX(3,8)
   p5 = pX(3,9)
   p6 = pX(3,10)
   p7 = pX(3,11)
   p8 = pX(3,12)
   p9 = pX(3,13)
   p10 = pX(3,14)
   p11 = pX(3,15)
    p12 = pX(3,16)
   p13 = pX(3,17)
   p14 = pX(3,18)
   p15 = pX(3,19)
   p16 = pX(3,20)
 endif


 if(Hdom>HthinStart .and. Hdom<HthinLim) then !!!first check if height is above 12 meters
    ! BA_lim = p1*H**3. + p2*H**2. + p3*H + p4
    ! BA_thd = p5*H**3. + p6*H**2. + p7*H + p8
    BA_limLow = p1*Hdom**3. + p2*Hdom**2. + p3*Hdom + p4
    BA_limUp = p5*Hdom**3. + p6*Hdom**2. + p7*Hdom + p8
    BA_lim = BA_limLow + (BA_limUp - BA_limLow) * BAlimPer !(0:1)
    
    BA_thdLow = p9*Hdom**3. + p10*Hdom**2. + p11*Hdom + p12
    BA_thdUp = p13*Hdom**3. + p14*Hdom**2. + p15*Hdom + p16
    BA_thd = BA_thdLow + (BA_thdUp - BA_thdLow) * BAthdPer !(0:1)
  
  baThin(1) = BA_lim
  baThin(2) = BA_thd
 else
  baThin(1) = 9999999999.9
  baThin(2) = 0.
 endif
end subroutine tapioThin
!*************************************************************

!tapioClearcut

!subroutine to check out if it's time for clearcut
!*************************************************************
subroutine tapioClearcut(species, siteType, ETSmean, dbh, age, ccTapio, ccPer, ccMature)
  implicit none
    LOGICAL :: ccMature 
    real (kind=8) :: species !1 for pine; 2 for spruce; 3 for betula pendula
  real (kind=8) :: siteType, ETSmean, dbh, age !siteType; average ETS of the site, average dbh of the stand 
  real (kind=8) :: ccTapio(5,3,3,5) !!dimensions are: 1st=SiteType; 2nd = species; 3rd= ETS; 4th=nTapioCC
  real (kind=8) :: pX(3,5) !pX(1) = ETS threshold; pX(2)= Hlim;  pX(3:5) equation parameters
    real (kind=8) :: dbhLim, dbhLimL, dbhLimU, ageLim
    real (kind=8) :: ccPer ! 0 = clearcut is done as soon as the first dbh limit is reached, 1 = clearcut is done at the upper dbh limit
  
pX = ccTapio(int(siteType), int(species),:,:)
 if(ETSmean > pX(1,1)) then !if we are in South Finland
  dbhLimL = pX(1,3)
  dbhLimU = pX(1,4)
  ageLim =  pX(1,5)
 elseif(ETSmean <= pX(1,1) .and. ETSmean >= pX(1,2)) then !if we are in Central Finland
  dbhLimL = pX(2,3)
  dbhLimU = pX(2,4)
  ageLim =  pX(2,5)
 else !if we are in Northern Finland
  dbhLimL = pX(3,3)
  dbhLimU = pX(3,4)
  ageLim =  pX(3,5)
 endif
 
 dbhLim = dbhLimL + (dbhLimU - dbhLimL) * ccPer 
 if(dbh>dbhLim .or. age>ageLim) then ! if the limits are met, the stand is ready for clearcut
  ccMature = .TRUE.
 else
  ccMature = .FALSE.
 endif
end subroutine tapioClearcut
!*************************************************************


!*************************************************************
!tapioFirstThin

!subroutine to do the first commercial thinning according to the Tapio rules
!selection thinning for pine or one thinning model for spruce has to be done manually!
!right now function does not check if stand is too old for first thinning!
!*******************************************

subroutine tapioFirstThin(species, siteType, ETSmean, ftTapio, hPer, densPer, early, output,nSp)
  implicit none
  integer,intent(in) :: nSp
    LOGICAL :: early ! if first thinning is done early; matters only for northern finland 
  real (kind=8),dimension(3) :: output ! outputs: 1 for height limit, 2 for density limit and 3 for result
    real (kind=8) :: species !1 for pine; 2 for spruce; 3 for betula pendula
  real (kind=8) :: siteType, ETSmean !siteType; average ETS of the site
  real (kind=8) :: ftTapio(5,nSp,3,7) !!dimensions are: 1st=SiteType; 2nd = species; 3rd = ETS; 4th = nftTapio
  real (kind=8) :: pX(3,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning
    real (kind=8) :: hLimL, hLimU, hLim, densityL, densityU, densityNew, densThd
    real (kind=8) :: hPer, densPer, densityMin ! hPer=0 first thinning is done as soon as the first height limit is reached, hPer=1 clearcut is done at the upper height limit
! densPer for adjusting the thinning result the same way

!open(1, file="tapioFTlog.txt")
  
pX = ftTapio(int(siteType), int(species),:,:)
 if(ETSmean > pX(1,1)) then !if we are in South or Central Finland
  densThd = pX(1,2)
  hLimL = pX(1,3)
  hLimU = pX(1,4)
  densityL = pX(1,5)
  densityU = pX(1,6)
 else if(early) then !parameters for early first thinning North Finland
    densThd = pX(2,2)
    hLimL = pX(2,3)
    hLimU = pX(2,4)
    densityL = pX(2,5)
    densityU = pX(2,6)
 else !parameters for late first thinning North Finland 
    densThd = pX(3,2)
    hLimL = pX(3,3)
    hLimU = pX(3,4)
    densityL = pX(3,5)
    densityU = pX(3,6)
 endif
 
! write(1, *) "densThd:", densThd, "densityU:", densityU
 hLim = hLimL + (hLimU - hLimL) * hPer
 densityMin = densityU*(1+densThd)
 densityNew = densityL + (densityU-densityL)*densPer
!  write(1, *) "hLim:", hLim, "densMin:", densityMin, "densNew:", densityNew

 ! if height is over the limit and density high enough, thinning is made to the new density
 !if(Hdom>hLim .and. density>(densityU*(1+densThd))) then 
  output(1) = hLim
  output(2) = densityMin
  output(3) = densityNew

!  close(1)
 !endif
end subroutine tapioFirstThin


!*************************************************************
!tapioTend

!subroutine to do the non-commercial thinning (tending of seedling stand) according to the Tapio rules
!tending for dense sown pine stands has to be done manually!
!right now function does not check if stand is too old for tending!
!*******************************************

subroutine tapioTend(species, siteType, ETSmean, tTapio, hPer, densPer, output,nSp)
  implicit none
  integer, intent(in) :: nSp
  real (kind=8),dimension(3) :: output ! outputs: 1 for height limit, 2 for density limit and 3 for result
    real (kind=8) :: species !1 for pine; 2 for spruce; 3 for betula pendula
  real (kind=8) :: siteType, ETSmean !siteType; average ETS of the site
  real (kind=8) :: tTapio(5,nSp,2,7) !!dimensions are: 1st=SiteType; 2nd = species; 3rd = ETS; 4th = nftTapio
  real (kind=8) :: pX(2,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning
    real (kind=8) :: hLimL, hLimU, hLim, densityL, densityU, densityNew, densThd
    real (kind=8) :: hPer, densPer, densityMin ! hPer=0 first thinning is done as soon as the first height limit is reached, hPer=1 clearcut is done at the upper height limit
! densPer for adjusting the thinning result the same way
  
pX = tTapio(int(siteType), int(species),:,:)
 if(ETSmean > pX(1,1)) then !if we are in South or Central Finland
  densThd = pX(1,2)
  hLimL = pX(1,3)
  hLimU = pX(1,4)
  densityL = pX(1,5)
  densityU = pX(1,6)
 else !if we are in North Finland
  densThd = pX(2,2)
  hLimL = pX(2,3)
  hLimU = pX(2,4)
  densityL = pX(2,5)
  densityU = pX(2,6)
 endif
 
! open(1, file="tapioTendLog.txt")
 
 hLim = hLimL + (hLimU - hLimL) * hPer
 densityMin = densityU*(1+densThd)
 densityNew = densityL + (densityU-densityL)*densPer
 
! write(1, *) "hLim:", hLim, "densMin:", densityMin, "densNew:", densityNew

 ! if height is over the limit and density high enough, thinning is made to the new density
 !if(Hdom>hLim .and. density>(densityU*(1+densThd))) then 
!  density = densityNew
  
output(1) = hLim
output(2) = densityMin
output(3) = densityNew

!write(1, *) output
!close(1)
 !endif
end subroutine tapioTend

!*****************************************
! chooseThin
! subroutine to pick the right thinning function: returns thinning = 1 for tapioTend, 2 for tapioFirstThin and 3 for tapioThin
!*****************************************

subroutine chooseThin(species, siteType, ETSmean, density, Hdom, tTapio, ftTapio, thinning,nSp)
  implicit none
    integer,intent(in) :: species,nSp !1 for pine; 2 for spruce; 3 for betula pendula
  real (kind=8),intent(out) :: thinning ! 1 for tapioTend, 2 for tapioFirstThin, 3 for tapioThin
  real (kind=8),intent(in) :: siteType, ETSmean, density, Hdom !siteType; average ETS of the site
  real (kind=8),intent(in) :: tTapio(5,nSp,2,7), ftTapio(5,nSp,3,7) 
  real (kind=8) :: pX1(3,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning; pX(7) height limit to move on to next thinning function
  real (kind=8) :: pX2(2,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning; pX(7)  height limit to move on to next thinning function
    real (kind=8) :: densityU1, hNext1, densityU2, hNext2

! log for testing 
!open(1, file = "chooseThinLog.txt")

! parameters for first thinning
pX1 = ftTapio(int(siteType), species,:,:)
 if(ETSmean > pX1(1,1)) then !if we are in South or Central Finland
  densityU1 = pX1(1,6)
  hNext1 = pX1(1,7)
 else 
  densityU1 = pX1(2,6)
  hNext1 = pX1(2,7)
 endif
 
! parameters for tending
pX2 = tTapio(int(siteType), species,:,:)
 if(ETSmean > pX2(1,1)) then !if we are in South or Central Finland
  densityU2 = pX2(1,6)
  hNext2 = pX2(1,7)
 else 
  densityU2 = pX2(2,6)
  hNext2 = pX2(2,7)
 endif

! writing test log 
!write(1, *) "first thinning: density under", densityU1, "stems/ha or dheight over", hNext1, &
!  "m; tending: density under", densityU2, "stems/ha or dheight over", hNext2, "m"

! if the stand is already thinner than the thinning result or the dominant height is over the limit, we move on to tapioThin subroutine
 if(density < densityU1 .or. Hdom > hNext1) then
  thinning = 3. 
 else if(density < densityU2 .or. Hdom > hNext2) then
  thinning = 2.
 else 
  thinning = 1.
 endif

 ! closing test log
!close(1) 
  
end subroutine chooseThin


!  function to calculate fAPAR of ground vegetation
!***************************************************************
! subroutine fAPARgv(fAPARstand,ets,siteType,agW,bgW,fAPAR_gv,litAG,litBG)
subroutine fAPARgv(fAPARstand,ets,siteType,totfAPAR_gv,totlitGV,p0,AWENs,totWGV) !reduced input output  
  implicit none
    real (kind=8) :: fAPARstand,ets,siteType, litAG(3), litBG(2),p0
  real (kind=8) :: totfAPAR_gv,totlitGV,totWGV
  real (kind=8) :: bgW(2),agW(3), xx(3),lai_gv(3),fAPAR_gv(3)!x_g, x_s, x_m !%cover grass&herbs, shrubs and mosses&lichens
    real (kind=8) :: b_g,a_g,a_s,a_m,b_m,alpha_ag(3),beta_ag(3),alpha_bg(2),beta_bg(2),laB(3)
  real (kind=8) :: turnAG(3),turnBG(2),p0ref
  real (kind=8) :: AWENs(4), AWENsh(4), AWENghAG(4), AWENml(4), AWENghBG(4)
 
 
 p0ref=1400.0
 !!!set parameters %cover
 if(siteType == 1.) then
  a_g = 0.4; a_s = 0.1; a_m = 0.1; b_m = 0.3; b_g = 0.4
 elseif(siteType == 2.) then
  a_g = 0.3; a_s = 0.3; a_m = 0.1; b_m = 0.4; b_g = 0.3
 elseif(siteType == 3.) then
  a_g = 0.3; a_s = 0.55; a_m = 0.3; b_m = 0.6; b_g = 0.1
 elseif(siteType == 4.) then
  a_g = 0.1; a_s = 0.65; a_m = 0.6; b_m = 0.8; b_g = 0.05
 elseif(siteType >4.5) then
  a_g = 0.05; a_s = 0.7; a_m = 0.6; b_m = 0.8; b_g = 0.
 endif
 
 alpha_ag = (/3152.4,1390.1,1637.1/)
 beta_ag = (/1.107,0.9496,0.8289/)
 alpha_bg = (/5016.0,1278.3/)
 beta_bg = (/0.831,0.6628/)
 if(ets<700.) then
  alpha_ag(1) = 3627.1
  beta_ag(1) = 0.948
  alpha_bg(1) = 5587.7
  beta_bg(1) = 0.45
 endif
 
 laB = (/6.,7.,1./)
 turnAG = (/0.37,0.54,0.2/) !!!turnovers above ground srbs,g&h, m&l
 turnBG = (/0.08,0.59/) !!!turnovers below ground srbs,g&h, m&l
 
 AWENsh = (/0.557,0.225264,0.086736,0.131/) !!!!Awen parameters for shrubs
 AWENghAG = (/0.273,0.427518,0.274482,0.025/) !!!!Awen parameters for grass&herbs aboveground
 AWENghBG = (/0.273,0.506844,0.195156,0.025/) !!!!Awen parameters for grass&herbs belowground
 AWENml = (/0.786,0.088445,0.034055,0.0915/) !!!!Awen parameters for lichens and mosses

 !% cover calculations
 ! if(fAPARstand==0.) then
  ! fAPARstand = 0.04
 ! endif
 xx(1) = a_s * (fAPARstand+0.2) * (log((1.04)/(fAPARstand+0.04))) **0.5 !!1.04 0.04 is to avoid division by 0 and to make the argument of log 1 at fAPAR=1
 xx(2) = a_g * (1-fAPARstand) + b_g
 xx(3) = a_m * (1-fAPARstand) + b_m * fAPARstand
 
 !! calculate biomasses
 ! agW = P0/p0ref * alpha_ag * xx ** (beta_ag)*0.5
 ! bgW = P0/p0ref * alpha_bg * xx(1:2) ** (beta_bg)*0.5  !!!!0.5 converts DW to carbon
 agW = alpha_ag * xx ** (beta_ag)*0.5
 bgW = alpha_bg * xx(1:2) ** (beta_bg)*0.5  !!!!0.5 converts DW to carbon

 !! calculate litterfal
 litAG = agW * turnAG
 litBG = bgW * turnBG
 !! calculate AWEN
 AWENs =  litAG(1) * AWENsh + litAG(2)*AWENghAG + litAG(3) * AWENml + &
  litBG(1)*AWENsh + litBG(2) *AWENghBG
 
 ! !calculate LAI
 lai_gv = agW * laB / (10000 * 0.5)   !!!!0.5 converts SLA (m2/kg DW) to (m2/kg C) 
  
 fAPAR_gv(1) = (1-fAPARstand) * (1-exp(-0.5*(lai_gv(1))))
 fAPAR_gv(2) = (1-fAPARstand-fAPAR_gv(1)) * (1-exp(-0.5*(lai_gv(2))))
 fAPAR_gv(3) = (1-fAPARstand-fAPAR_gv(1)-fAPAR_gv(2)) * (1-exp(-0.5*(lai_gv(3))))
 totfAPAR_gv = sum(fAPAR_gv)   !!!alternatively 0.6354*fAPARstand + 0.3716
 totlitGV = sum(litAG) + sum(litBG)
 totWGV = sum(agW) + sum(bgW)
 
end subroutine fAPARgv


!  function to calculate fAPAR of ground vegetation and report GVoutputs splitted by vegetation type (grass&herbs, srhubs, mosses&lichens)
!***************************************************************
! subroutine fAPARgv(fAPARstand,ets,siteType,agW,bgW,fAPAR_gv,litAG,litBG)
subroutine fAPARgvByVtypes(inputs, output) !reduced input output  
  implicit none
  real (kind=8) :: output(13) !fAPAR(3), Lit_ab(3),Lit_bg(2),W_ab(3),W_bg(2)
  real (kind=8) :: inputs(3) !fAPARstand,ets,siteType
    real (kind=8) :: fAPARstand,ets,siteType, litAG(3), litBG(2)
  real (kind=8) :: totfAPAR_gv,totlitGV,totWGV
  real (kind=8) :: bgW(2),agW(3), xx(3),lai_gv(3),fAPAR_gv(3)!x_g, x_s, x_m !%cover grass&herbs, shrubs and mosses&lichens
    real (kind=8) :: b_g,a_g,a_s,a_m,b_m,alpha_ag(3),beta_ag(3),alpha_bg(2),beta_bg(2),laB(3)
  real (kind=8) :: turnAG(3),turnBG(2),p0ref
  real (kind=8) :: AWENs(4), AWENsh(4), AWENghAG(4), AWENml(4), AWENghBG(4)
 
 
 fAPARstand = inputs(1)
 ets = inputs(2)
 siteType = inputs(3)
 p0ref=1400.0
 !!!set parameters %cover
 if(siteType == 1.) then
  a_g = 0.4; a_s = 0.1; a_m = 0.1; b_m = 0.3; b_g = 0.4
 elseif(siteType == 2.) then
  a_g = 0.3; a_s = 0.3; a_m = 0.1; b_m = 0.4; b_g = 0.3
 elseif(siteType == 3.) then
  a_g = 0.3; a_s = 0.55; a_m = 0.3; b_m = 0.6; b_g = 0.1
 elseif(siteType == 4.) then
  a_g = 0.1; a_s = 0.65; a_m = 0.6; b_m = 0.8; b_g = 0.05
 elseif(siteType >4.5) then
  a_g = 0.05; a_s = 0.7; a_m = 0.6; b_m = 0.8; b_g = 0.
 endif
 
 alpha_ag = (/3152.4,1390.1,1637.1/)
 beta_ag = (/1.107,0.9496,0.8289/)
 alpha_bg = (/5016.0,1278.3/)
 beta_bg = (/0.831,0.6628/)
 if(ets<700.) then
  alpha_ag(1) = 3627.1
  beta_ag(1) = 0.948
  alpha_bg(1) = 5587.7
  beta_bg(1) = 0.45
 endif
 
 laB = (/6.,7.,1./)
 turnAG = (/0.37,0.54,0.2/) !!!turnovers above ground srbs,g&h, m&l
 turnBG = (/0.08,0.59/) !!!turnovers below ground srbs,g&h, m&l
 
 AWENsh = (/0.557,0.225264,0.086736,0.131/) !!!!Awen parameters for shrubs
 AWENghAG = (/0.273,0.427518,0.274482,0.025/) !!!!Awen parameters for grass&herbs aboveground
 AWENghBG = (/0.273,0.506844,0.195156,0.025/) !!!!Awen parameters for grass&herbs belowground
 AWENml = (/0.786,0.088445,0.034055,0.0915/) !!!!Awen parameters for lichens and mosses

 !% cover calculations
 ! if(fAPARstand==0.) then
  ! fAPARstand = 0.04
 ! endif
 xx(1) = a_s * (fAPARstand+0.2) * (log((1.04)/(fAPARstand+0.04))) **0.5 !!1.04 0.04 is to avoid division by 0 and to make the argument of log 1 at fAPAR=1
 xx(2) = a_g * (1-fAPARstand) + b_g
 xx(3) = a_m * (1-fAPARstand) + b_m * fAPARstand
 
 !! calculate biomasses
 ! agW = P0/p0ref * alpha_ag * xx ** (beta_ag)*0.5
 ! bgW = P0/p0ref * alpha_bg * xx(1:2) ** (beta_bg)*0.5  !!!!0.5 converts DW to carbon
 agW = alpha_ag * xx ** (beta_ag)*0.5
 bgW = alpha_bg * xx(1:2) ** (beta_bg)*0.5  !!!!0.5 converts DW to carbon

 !! calculate litterfal
 litAG = agW * turnAG
 litBG = bgW * turnBG
 !! calculate AWEN
 AWENs =  litAG(1) * AWENsh + litAG(2)*AWENghAG + litAG(3) * AWENml + &
  litBG(1)*AWENsh + litBG(2) *AWENghBG
 
 ! !calculate LAI
 lai_gv = agW * laB / (10000 * 0.5)   !!!!0.5 converts SLA (m2/kg DW) to (m2/kg C) 
  
 fAPAR_gv(1) = (1-fAPARstand) * (1-exp(-0.5*(lai_gv(1))))
 fAPAR_gv(2) = (1-fAPARstand-fAPAR_gv(1)) * (1-exp(-0.5*(lai_gv(2))))
 fAPAR_gv(3) = (1-fAPARstand-fAPAR_gv(1)-fAPAR_gv(2)) * (1-exp(-0.5*(lai_gv(3))))
 output(1:3) = fAPAR_gv
 output(4:6) = litAG
 output(7:8) = litBG
 output(9:11) = agW
 output(12:13) = bgW

 
end subroutine fAPARgvByVtypes




SUBROUTINE runYassoAWENin(AWENin,nYears, nSites, litSize,nClimID,climIDs,pYasso,weatherYasso,soilC)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION 
    !********************************************* &
    ! run yasso for some years with litterfal inputs from prebas.

  integer, intent(in) :: nYears, nSites, nClimID
  REAL (kind=8),INTENT(IN) :: AWENin(nSites, nYears, 5) 
  REAL (kind=8),INTENT(IN) :: weatherYasso(nClimID, nYears, 3)
  ! REAL (kind=8),INTENT(IN) :: species(nSites, nLayers)
  REAL (kind=8),INTENT(IN) :: pYasso(35),litSize
  real (kind=8),INTENT(inout) :: soilC(nSites,(nYears+1),5)
  integer,INTENT(IN) :: climIDs(nSites)
  INTEGER :: year, site, layer, spec
  real (kind=8) :: t=1.,Lst,Lb,Lf,leac=0.,stSt=0. !leaching parameter for Yasso
  real (kind=8),DIMENSION(5) :: AWENH
  

! fbAWENH = 0.
! folAWENH = 0.
! stAWENH = 0.

!!!!run Yasso
do site = 1, nSites
  do year = 1,nYears

   ! Lst = litter(site,year,layer,3)
   ! Lb = litter(site,year,layer,2)
   AWENH = AWENin(site,year,:)

   ! spec = int(species(site,layer))
   ! call compAWENH(Lf,folAWENH,pAWEN(1:4,spec))   !!!awen partitioning foliage
   ! call compAWENH(Lb,fbAWENH,pAWEN(5:8,spec))   !!!awen partitioning branches
   ! call compAWENH(Lst,stAWENH,pAWEN(9:12,spec))         !!!awen partitioning stems

   ! call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,1,layer),stAWENH,litterSize(1,spec), &
  ! leac,soilC(site,(year+1),:,1,layer),stSt)
   ! call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,2,layer),fbAWENH,litterSize(2,spec), &
  ! leac,soilC(site,(year+1),:,2,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:),AWENH,litSize, &
  leac,soilC(site,(year+1),:),stSt)
  
 enddo
enddo

END SUBROUTINE runYassoAWENin  




! subroutine multiGV(fAPARstand,ets,siteType,agW,bgW,fAPAR_gv,litAG,litBG)
subroutine multiGV(fAPARstand,ets,siteType,totfAPAR_gv,totlitGV,p0,AWENs,nYears,nSites,totWGV) !reduced input output  
  implicit none
  integer,INTENT(IN) :: nYears,nSites
    real (kind=8) :: fAPARstand(nSites,nYears),ets(nSites,nYears),siteType(nSites),p0(nSites,nYears)
  real (kind=8) :: totfAPAR_gv(nSites,nYears),totlitGV(nSites,nYears),totWGV(nSites,nYears)
  real (kind=8) :: AWENs(nSites,nYears,4)
  integer :: site,year
  
 do site = 1, nSites
  do year = 1,nYears
  call fAPARgv(fAPARstand(site,year),ets(site,year),siteType(site),totfAPAR_gv(site,year), &
      totlitGV(site,year),p0(site,year),AWENs(site,year,:),totWGV(site,year))
  enddo !year
 enddo !site
 
 
end subroutine multiGV



SUBROUTINE runYassoAWENinMonthly(AWENin,nYears, nSites, litSize,nClimID,climIDs,pYasso,weatherYasso,soilC)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION 
    !********************************************* &
    ! run yasso for some years with litterfal inputs from prebas.

  integer, intent(in) :: nYears, nSites, nClimID
  REAL (kind=8),INTENT(IN) :: AWENin(nSites, nYears, 5) 
  REAL (kind=8),INTENT(IN) :: weatherYasso(nClimID, nYears, 3)
  ! REAL (kind=8),INTENT(IN) :: species(nSites, nLayers)
  REAL (kind=8),INTENT(IN) :: pYasso(35),litSize
  real (kind=8),INTENT(inout) :: soilC(nSites,(nYears+1),5)
  integer,INTENT(IN) :: climIDs(nSites)
  INTEGER :: year, site, layer, spec
  real (kind=8) :: t=1./12,Lst,Lb,Lf,leac=0.,stSt=0. !leaching parameter for Yasso
  real (kind=8),DIMENSION(5) :: AWENH
  

! fbAWENH = 0.
! folAWENH = 0.
! stAWENH = 0.

!!!!run Yasso
do site = 1, nSites
  do year = 1,nYears

   ! Lst = litter(site,year,layer,3)
   ! Lb = litter(site,year,layer,2)
   AWENH = AWENin(site,year,:)

   ! spec = int(species(site,layer))
   ! call compAWENH(Lf,folAWENH,pAWEN(1:4,spec))   !!!awen partitioning foliage
   ! call compAWENH(Lb,fbAWENH,pAWEN(5:8,spec))   !!!awen partitioning branches
   ! call compAWENH(Lst,stAWENH,pAWEN(9:12,spec))         !!!awen partitioning stems

   ! call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,1,layer),stAWENH,litterSize(1,spec), &
  ! leac,soilC(site,(year+1),:,1,layer),stSt)
   ! call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:,2,layer),fbAWENH,litterSize(2,spec), &
  ! leac,soilC(site,(year+1),:,2,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),year,:),soilC(site,year,:),AWENH,litSize, &
  leac,soilC(site,(year+1),:),stSt)

  
 enddo
enddo

END SUBROUTINE runYassoAWENinMonthly  



SUBROUTINE runYassoMonthly(litter,litterSize,nMonths, nLayers, nSites, nSp,species,nClimID,climIDs,pAWEN,pYasso, &
      weatherYasso,soilC)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION 
    !********************************************* &
    ! run yasso for some years with litterfal inputs from prebas.

  integer, intent(inout) :: nMonths, nLayers, nSites, nSp,nClimID
  REAL (kind=8),INTENT(inout) :: litter(nSites, nMonths, nLayers, 3) !!!fourth dimension (3) 1 is fine litter, 2 = branch litter, 3=stemLitter
  REAL (kind=8),INTENT(inout) :: weatherYasso(nClimID, nMonths, 3)
  REAL (kind=8),INTENT(inout) :: species(nSites, nLayers),litterSize(3,nSp)
  REAL (kind=8),INTENT(inout) :: pAWEN(12, nSp), pYasso(35)
  real (kind=8),INTENT(inout) :: soilC(nSites,(nMonths+1),5,3,nLayers)
  integer,INTENT(inout) :: climIDs(nSites)
  INTEGER :: month, site, layer, spec
  real (kind=8) :: t=1./12,Lst,Lb,Lf,leac=0.,stSt=0. !leaching parameter for Yasso
  real (kind=8),DIMENSION(5) :: fbAWENH,folAWENH,stAWENH
  

fbAWENH = 0.
folAWENH = 0.
stAWENH = 0.

!!!!run Yasso
do site = 1, nSites
 do layer = 1,nLayers
  do month = 1,nMonths

   Lst = litter(site,month,layer,3)
   Lb = litter(site,month,layer,2)
   Lf = litter(site,month,layer,1)

   spec = int(species(site,layer))
   call compAWENH(Lf,folAWENH,pAWEN(1:4,spec))   !!!awen partitioning foliage
   call compAWENH(Lb,fbAWENH,pAWEN(5:8,spec))   !!!awen partitioning branches
   call compAWENH(Lst,stAWENH,pAWEN(9:12,spec))         !!!awen partitioning stems

   call mod5c(pYasso,t,weatherYasso(climIDs(site),month,:),soilC(site,month,:,1,layer),stAWENH,litterSize(1,spec), &
  leac,soilC(site,(month+1),:,1,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),month,:),soilC(site,month,:,2,layer),fbAWENH,litterSize(2,spec), &
  leac,soilC(site,(month+1),:,2,layer),stSt)
   call mod5c(pYasso,t,weatherYasso(climIDs(site),month,:),soilC(site,month,:,3,layer),folAWENH,litterSize(3,spec), &
  leac,soilC(site,(month+1),:,3,layer),stSt)
  
  ! soilC(site,(month+1),:,:,layer) = soilC(site,month,:,:,layer) + (soilC(site,(month+1),:,:,layer) -soilC(site,month,:,:,layer))/12

  enddo
 enddo
enddo

END SUBROUTINE runYassoMonthly  



!*****************************************************************************************
!    SUBROUTINE FMortality
!
!    This subroutine updates the number of trees in size !lasses
!    according to alternative mortality assumptions
!     
!
!    Mortality has been removed from the FGrowth2 subroutine (December 2010 AM)
!    W: biomass matrix, components Wf, Wfr, Wb, Wcr,Wsapw, Whw,
!    rf: relative growth rate of foliage ,
!    rf: in Prebas we can look at the length of the crown
!****************************************************************************************
      SUBROUTINE FMortality(D13, BA, N, h, dN, dBA,rf, kokoluokka, ind)

  implicit none
    
  ! integer Maxhakkuu, maxkokoluokka, Maxvuodet, MaxOksat, VolWithbark
  ! Parameter (Maxhakkuu = 5, maxkokoluokka = 451, Maxvuodet = 451)
  ! Parameter (MaxOksat = 50, VolWithBark = 1)

  integer,INTENT(inout) :: kokoluokka, ind
  REAL (kind=8),INTENT(in) :: D13(kokoluokka), N(kokoluokka)
  REAL (kind=8),INTENT(in) :: BA(kokoluokka), rf(kokoluokka) 
  REAL (kind=8),INTENT(out) :: dN(kokoluokka)
  REAL (kind=8),INTENT(in) :: h(5,kokoluokka),dBA(kokoluokka)
    
!**********************************************************
    REAL (kind=8) :: delta1, delta2, Diam, dDiam, rNs, rN0
  REAL (kind=8) :: phi, xij, dbhT, alpha, beta1, beta2
  REAL (kind=8) :: a1, a2, CI, a0, mprob, sigmau
    REAL (kind=8) :: alpha_d, D13_d,k_mort1,k_mort2,greff,coeff
    REAL (kind=8) :: m_intr, m_greff, BAtot, Ntot, Reineke, cr, deltaN, dev
  integer :: j
! ********************************************************************************* 
! these data commands relate to mortality in the model by Peltoniemi and Mäkipää
! choose first line for model 8 (no dbh effect) and second line for model 9
! third line is common for both. these are individual tree model, different layers
! *********************************************************************************
!  data a0,a1,a2, sigmau /-5.9,0.00606,0, 0.207/
   data a0,a1,a2, sigmau /-6.9,0.00854,3.55, 0.204/
  data alpha,beta1,beta2/ 800., 0.035, 0.04/
      data alpha_d,D13_d,k_mort1,k_mort2/2.3,210.,0.005,0.3/
      data delta1,delta2/30000.,2000./


      Diam = sqrt(1.2732 * BA(ind)) 

      if(Diam>0.) then
      dDiam = 0.6366 * dBA(ind) / Diam
      else
      dDiam = sqrt(1.2732 * dBA(ind))
      endif
      
! !**********************************************************************
! !
! ! now apply Mikko's mortality model. result in mortality in each layer, we need to take account all the layers
! !
! !**********************************************************************

  
  CI = 0.
      Ntot = 0
  dbhT = beta1 + beta2 * D13(ind) / 100.
      
      do j = 1, kokoluokka
          Ntot = Ntot + N(j)
      enddo
      

      do j = 1, kokoluokka

          if(N(j) > 0.) then

    if(j.ne.ind) then 

! *** Mortality due to other size classes
        
      xij = (D13(j) - D13(ind)) / 100.
          if((alpha*(xij-dbhT)) < 10.) then
       phi = exp(alpha*(xij-dbhT))
       phi = phi / (1. + phi)
          else
       phi = 1.
      endif


    CI = CI + N(j)*phi * sqrt(D13(j)/100.)  

    else

! ***  account for the same size class assuming 5% deviation
!      (has hardly any impact)

          dev = 0.05 
!          dev = N(ind) / Ntot

          xij = dev * D13(ind)/100.

      if((alpha*(xij-dbhT)) < 10.) then
       phi = exp(alpha*(xij-dbhT))
       phi = phi / (1. + phi)
          else
       phi = 1.
      endif

    CI = CI + dev * N(j)*phi * sqrt(D13(j)/100.)
      
    if((alpha*(-xij-dbhT)) < 10.) then
       phi = exp(alpha*(-xij-dbhT))
       phi = phi / (1. + phi)
          else
       phi = 1.
      endif

    CI = CI + dev * N(j)*phi * sqrt(D13(j)/100.)

     
    endif     
          endif
          
      if(h(5,ind)>2.)then
    if ((1.+ delta1*dDiam + delta2*Diam*Diam)>0 &
       .and.(1.-1./(1.+ delta1*dDiam + delta2*Diam*Diam)>0)) then
        rNs = -log(1.-1./(1.+ delta1*dDiam + delta2*Diam*Diam) )
    endif
      else
     rNs = 0.
      endif
  end do

  mprob = a1 * CI + a2 * D13(ind) / 100. + a0
!  write(6,*) ind, CI, D13(ind), mprob

  mprob = exp(mprob + sigmau/2.)
      
      if(D13(ind)>0.) then

  dN(ind) = - mprob * N(ind) - rNs*N(ind)
      if(rf(ind)< 0.0) then
      dN(ind) = dN(ind) + rf(ind) * N(ind)
      endif
      else
          
          dN(ind) = 0.
      endif
      
end subroutine Fmortality




subroutine calWf_fA(par_rhof,Wf,nData,As)
 IMPLICIT NONE
 integer, intent(in) :: nData
 REAL (kind=8),INTENT(inOUT) :: Wf(nData),As(ndata)         
 REAL (kind=8),INTENT(INout) :: par_rhof !!parameters
 
   Wf = par_rhof * As

END SUBROUTINE calWf_fA

subroutine calWf_fLc(pars,Wf,nData,Lc)
 IMPLICIT NONE
 integer, intent(in) :: nData
 REAL (kind=8),INTENT(inOUT) :: Wf(nData),Lc(ndata)         
 REAL (kind=8),INTENT(INout) :: pars(2) !!parameters
 REAL (kind=8) ksi, z
    
  ksi = pars(1)
  z = pars(2)
   Wf = ksi * Lc ** z 
  
END SUBROUTINE calWf_fLc

subroutine calAs_fLc(pars,As,nData,Lc)
 IMPLICIT NONE
 integer, intent(in) :: nData
 REAL (kind=8),INTENT(inOUT) :: As(nData),Lc(ndata)         
 REAL (kind=8),INTENT(INout) :: pars(3) !!parameters
 REAL (kind=8) ksi, z, rhof
 
  ksi = pars(1)
  z = pars(2)
  rhof = pars(3)
   As = ksi/rhof * Lc ** z 
  
END SUBROUTINE calAs_fLc



function MatVecMult(A, v) result (w)
   implicit none

   real, dimension(:,:), intent(in) :: A
   real, dimension(:), intent(in)   :: v
   real, dimension( SIZE(A,1) )     :: w

   integer :: i, j
   integer :: N

   N = size(v)

   w = 0.0       !! clear whole vector                      
   DO i = 1, N
      w = w + v(i) * A( :, i )
   END DO
  end function




subroutine calRein(outputs,nLayers,pRein,nVar,nSp,reinX)
 IMPLICIT NONE
 real (kind=8), parameter :: pi = 3.1415927
 integer, intent(in) :: nLayers, nVar,nSp
 real (kind=8), intent(in) :: outputs(nVar,nLayers) !!takes multiOut(maxYears,nVar,maxNlayers,2)
                            !!dim1 nVar, dim2 maxNlayers 1 year and 1standing trees
 real (kind=8), intent(in) :: pRein(nSp)
 real (kind=8), intent(inout) ::reinX
 real (kind=8) :: valX(nLayers), domSp(1),  Ntot, B, Reineke(nLayers)
 integer :: i, ijx,layerX(nLayers)
 
 ! nLayers = integer(size(outputs,3))
  
 valX = outputs(11,:)
 do ijx = 1, nLayers
  domSp = maxloc(valX)
  layerX(ijx) = int(domSp(1))
  valX(layerX(ijx)) = -999.
  
  Ntot = sum(outputs(17,layerX(1:ijx)))
  B = sum(outputs(35,layerX(1:ijx))*outputs(17,layerX(1:ijx)))/Ntot   !!!!!!!!!#####changed

   if(Ntot>0. .and. (outputs(17,layerX(ijx))>0.)) then
     Reineke(layerX(ijx)) = Ntot*(sqrt(B*4/pi)*100./25.)**(1.66) &
    /pRein(int(outputs(4,ijx)))*(outputs(13,ijx)/sum(outputs(13,:)))
   else
     Reineke(layerX(ijx)) = 0.
   endif
 enddo
  
 reinX = sum(Reineke)
END subroutine calRein


  subroutine changeOrder(oldOrd,age,newOrd,nSites,ageX)
  implicit none
  ! integer, allocatable :: indices1(:)
  ! integer, allocatable :: indices2(:)
  integer,intent(in) :: nSites
  real(8),intent(in) :: ageX
  real(8),intent(inout) :: oldOrd(nSites),newOrd(nSites)
  real(8), intent(inout) :: age(nSites)
  integer i, nX, index1(nSites), index2(nSites),indexX,nY
  
  nX=0
  nY=0
  
  do i = 1, nSites
  indexX = int(oldOrd(i))
   if(age(indexX) <= ageX) then
    nX=nX+1
    index1(nX) = oldOrd(i)
   else
    nY=nY+1
    index2(nY) = oldOrd(i)
   endif
  enddo
  newOrd(1:nX) = index1(1:nX)
  newOrd((nX+1):nSites) = index2(1:nY)
  
  !!!!old version 
    ! indices1 = PACK([(i, i=1,nSites)], age<=ageX)
    ! ! if(ij==1) write(1,*) indices1
    ! nX = size(indices1)
    ! newOrd(1:nX) = oldOrd(indices1)
    ! indices2 = PACK([(i, i=1,nSites)], age>ageX)
    ! newOrd((nX+1):nSites) = oldOrd(indices2)
  !!!!end old version 
  end subroutine changeOrder 
  
  
subroutine testOption(a,b,c,valX1,valX2,valX3)
  implicit none
  real(8),intent(inout) :: a, b, c
  real(8), optional, intent(in) :: valX1,valX2,valX3
  real(8) :: v1,v2,v3
  v1=1.
  v2=2.
  v3=3.
  if(present(valX1)) then 
   v1 = valX1
  endif
  c=a+b+v1+v2+v3
end subroutine testOption 


subroutine calcAlfar(siteTypeOriginal,species,pCrobas,nLayers,siteTalfar,nSp,nYearsFert,npar,deltaSiteTypeFert)
  implicit none
  integer, intent(in) :: nLayers,nSp,nYearsFert,npar
  real(8),intent(inout) :: siteTypeOriginal,pCrobas(npar,nSp)
  real(8),intent(inout) :: siteTalfar(nYearsFert,nLayers,2),species(nLayers),deltaSiteTypeFert
  ! integer,intent(inout) :: year!,ind(3,2)
  integer :: i
  real(8):: alfarUnfert(nLayers), alfarFert(nLayers),slope,interc
  
    ! ind(:,1) = int(max(20+min(modOut(1,year,3,:,1),5.)-1.,21.))
  ! ind(:,2) = int(modOut(1,year,4,:,1))
  do i = 1,nLayers
  
    alfarUnfert(i) = pCrobas(int(max(20+min(siteTypeOriginal,5.),21.)),int(max(species(i),1.)))
    if(deltaSiteTypeFert<1.) then !!!this is not used in fertilization at thinning
      alfarFert(i) = alfarUnfert(i) * deltaSiteTypeFert
    else
      alfarFert(i) = pCrobas(int(max(20+min(siteTypeOriginal,5.)-deltaSiteTypeFert,21.)), &
        int(max(species(i),1.)))
    endif
    ! write(*,*) i,siteTAlpha(i,1),species(i),alfarFert(i), alfarUnfert(i)
    siteTalfar(:,i,1) = max(0.,siteTypeOriginal-1)
    siteTalfar(1:(nYearsFert/2),i,2) = alfarFert(i)
    slope = (alfarUnfert(i) - alfarFert(i))/(nYearsFert/2+1.-0.)
    interc = alfarFert(i) - slope*1.
    siteTalfar((nYearsFert/2+1):nYearsFert,i,2) = slope* real((/(i, i = 2, nYearsFert/2+1)/), kind=8) + interc
  enddo
  
end subroutine calcAlfar  
  

!!!random mortality calculations based on Siilipehto et 2020  
subroutine randMort(age,d,ba,N,rBApine,rBAbrd,slope,pSize,pMort,perBAmort,step,baDead)
  implicit none
  real(8),intent(inout) :: age,d,ba,N,rBApine,rBAbrd,slope,pSize,pMort,perBAmort,step,baDead
  !parameters
  ! real(8) :: m=0.0004d0, Interc=-0.01d0  !parameters of probability of mortality obtained fitting a linear model as function of tree density (Fig 3A siipilehto et al.2020)
  ! real(8) :: a=0.2643402d0, b= 0.9987199d0  !parameters of proportion of dead BA obtained fitting an exponential model as function of tree density (Fig 6C siipilehto et al.2020)
  real(8) :: randX, nX, baX, XB1,rp1,rp=314.16d0, XB2,rp2
  
  !parameters mod1&2 siilipehto et al. 2020
  real(8) :: int_m1, Age_m1, DQ_m1, sqrtBA_m1, BApineRelxAge_m1
  real(8) :: BAbrdRelxAge_m1,Slope_100_m1, west_m1, TH5xBA_m1, peat_m1
  real(8) :: int_m2, lnN_m2, ba_m2, sqrtN_m2, lnDQ_m2, BApineRelxlnAge_m2
  real(8) :: BAbrdRelxlnAge_m2,Slope_100west_m2, Norway_m2

!initParameters 
    int_m1 = 3.1751d0
  Age_m1 = -0.00417d0
  DQ_m1 = 0.01382d0
  sqrtBA_m1 = -0.6693d0
  BApineRelxAge_m1 = 0.004356d0
  BAbrdRelxAge_m1 = -0.01112d0
  Slope_100_m1 = -1.0542d0
  west_m1 = 0.1354d0
  TH5xBA_m1 = 0.005803d0
  peat_m1 = -0.2201d0
  int_m2 = -7.3745d0
  lnN_m2=  1.4010d0
  sqrtN_m2 = -0.0296d0
  ba_m2 = -0.0117d0
  lnDQ_m2 = 0.5633d0
  BApineRelxlnAge_m2 = 0.0790d0
  BAbrdRelxlnAge_m2 = -0.0500d0
  Slope_100west_m2 = 0.2438d0
  Norway_m2 = 0.1657d0

  rp1=(step/5.d0) * (pSize/rp)
  XB1 = int_m1 + Age_m1 * age + DQ_m1*d**(1.5d0) + sqrtBA_m1 * sqrt(ba) + BApineRelxAge_m1* rBApine * age + &
    BAbrdRelxAge_m1* rBAbrd * age + Slope_100_m1* slope + west_m1 * 0.25d0+ TH5xBA_m1 * ba*0.005d0 + peat_m1 * 0.001d0
  rp2=(rp/pSize)
  XB2 = int_m2 + lnN_m2 * log(N) + sqrtN_m2 *sqrt(N) + ba_m2*ba + lnDQ_m2 *log(d) + &
    BApineRelxlnAge_m2 * rBApine * log(age) + BAbrdRelxlnAge_m2 * rBAbrd * log(age) + &
    Slope_100west_m2 *slope*0.25d0 + Norway_m2*0.d0
  
  pMort = 1.d0 - (1.d0 + exp(-(XB1)))**(-(rp1))!probability of survival
  perBAmort = 1.d0 - (1.d0+exp(-XB2))**(-rp2)
  ! if(age<30) pMort = max(min(pMort*2,0.5),0.4)
  ! if(age<30) perBAmort = min(perBAmort*2,0.5)
  
  call random_number(randX)
   if (randX < pMort) then
    baDead = perBAmort * ba
  endif
   
end subroutine
  
  
  !!!!calculate C:N based on ets and site type 
!Note: parameters should be included as arguments inputs instead of being internally assigned
subroutine CNratio(CN, ETS, st,pars)
implicit none

  real(8),intent(in) :: ETS, st, pars(3)
  real(8),intent(out) :: CN
  !parameters
  real(8) :: int_CN, p_ETS, p_st
  
  !init parameters
  int_CN = pars(1)
  p_ETS = pars(2)
  p_st = pars(3)

 !calculate CN ratio
 CN = int_CN + p_ETS * ETS + p_st * st
endsubroutine  

!!!!calculate ration of ECM biomass to fine root biomass (from Ostonen et al. 2011)
!Note: parameters should be included as arguments inputs instead of being internally assigned
subroutine rhoMcalc(rho_M, CN,pars)
implicit none

  real(8),intent(in) :: CN, pars(4)
  real(8),intent(out) :: rho_M
  !parameters
  real(8) :: p1, p2, p3, p4
  
  !init parameters
  p1 = pars(1)
  p2 = pars(2)  
  p3 = pars(3)  
  p4 = pars(4)  

 !calculate rho_M
 rho_M = p1 + p2 / (1.d0 + exp(p3*(CN+p4)))

endsubroutine



subroutine CUEcalc(ETS, st,r_r,W_RT,r_RT,rm_aut_roots,litt_RT,exud,P_RT,normFactETS,pars)

implicit none
  real(8),intent(in) :: ETS, st,r_r,W_RT, pars(12), normFactETS
  real(8),intent(out) :: rm_aut_roots,litt_RT,exud,r_RT,P_RT

  real(8) :: r_F, s_F, rho_M,CN
  !parameters
  real(8) :: h_M, s_H, phi_M, ksi_M, gamma_M !, r_M  

!initialise parameters
h_M = pars(1) !extramatrical hyphal biomass parameter
s_H = pars(2) * normFactETS !hyphae specific turnover rate
phi_M = pars(3) !priming parameter 
ksi_M = pars(4) * normFactETS !rate of exudation
gamma_M = pars(5) * normFactETS !apparent hyphal respiration rate 

!!!!!!!!!! I do not find the value of r_M
! r_M  = 0.d0  !Maintenance respiration rate of fungal tips
!!!!!!!!!! I do not find the value of r_M

call CNratio(CN, ETS, st,pars(6:8))
call rhoMcalc(rho_M, CN,pars(9:12))


!  r_F=(∆_M+(1+c_H) h_M s_H+r_M+h_M r_H+ξ)   !original
r_F = h_M * s_H + (1-phi_M) * ksi_M+ gamma_M + r_r !r_M    !assuming d_M = 0.

! s_F=(ds_M+h_M * s_H+ksi_M) !original
s_F=(h_M * s_H + ksi_M) !assuming ds_M = 0


! Here, r_RT is the apparent maintenance respiration rate of fine roots when C input to the fungi has been taken into account.
!Fine root maintenance respiration is replaced with the following:
r_RT = (r_r+r_F * rho_M)/(1+rho_M)

! When calculating the related loss of C to the atmosphere, we need to subtract the losses from this, i.e., use 
!rm_aut_roots  losses of C from roots & ECM to 
! atmosphere, i.e., actual autotrophic maintenance respiration of roots + ECM (here 
! respiration of ECM is regarded autotrophic because it comes from photosynthates
rm_aut_roots = (r_r+(r_F-s_F) * rho_M)/(1 + rho_M)

!normalise parameters according to P and ETS factors
! s_H = s_H * normFactETS
! phi_M = phi_M * normFactETS
rm_aut_roots = rm_aut_roots

! To Yasso as root litter (or later maybe special composition as hyphal litter)
litt_RT = (rho_M * h_M * s_H)/(1+rho_M) * W_RT
! To Yasso as sugar (W) exudates
exud = (rho_M *ksi_M)/(1+rho_M) * W_RT

P_RT =  (rho_M *ksi_M*phi_M)/(1+rho_M) * W_RT

end subroutine

! subroutine test(ciao,species)
  
  ! implicit none
  ! integer, intent(inout) :: ciao(4)
  ! real(8) :: species(1)
  ! species(1) = 0
  ! ciao(2) = 1/species(1)
  ! species = maxloc(ciao)
  
  ! if(species(1) /= species(1)) species(1)=999
  
! end subroutine test

! calculate climate dependence on decomposition based on YASSO equations
SUBROUTINE fTyasso(theta,climate,fTaweNH)
IMPLICIT NONE
    !********************************************* &
    ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
    !********************************************* &
    ! returns the model prediction xt for the given parameters
    ! 1-16 matrix A entries: 4*alpha, 12*p

    ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

    ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

    ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

    ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

    ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

    ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

    ! 33-35 Woody parameters: theta_1, theta_2, r

    REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
    REAL (kind=8),DIMENSION(3),INTENT(IN) :: climate ! climatic conditions
    REAL (kind=8),DIMENSION(3),INTENT(out) :: fTaweNH ! climate dependence on decomposition 

    INTEGER :: i
    REAL (kind=8),PARAMETER :: pi = 3.141592653589793
    REAL (kind=8) :: tem,temN,temH
    REAL (kind=8),DIMENSION(5) :: te
    REAL (kind=8),DIMENSION(5) :: z1,z2
    REAL (kind=8),PARAMETER :: tol = 1E-12
    LOGICAL :: ss_pred = .FALSE.

    !#########################################################################
    ! Compute the coefficient matrix A for the differential equation

    ! temperature annual cycle approximation
    te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
    te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
    te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
    te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

    tem = 0.0
    temN = 0.0
    temH = 0.0
    DO i = 1,4 ! Average temperature dependence
        tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
    END DO

    ! Precipitation dependence
    tem = tem*(1.0-EXP(theta(28)*climate(2)/1000.0))
    temN = temN*(1.0-EXP(theta(29)*climate(2)/1000.0))
    temH = temH*(1.0-EXP(theta(30)*climate(2)/1000.0))

fTaweNH(1) = tem
fTaweNH(2) = temN
fTaweNH(3) = temH

    END SUBROUTINE

!!!toy model disturbance calculations   
subroutine pDistTest(ETS,parFire,pDistFireX)
  implicit none
  real(8), intent(inout) :: ETS,pDistFireX
  !parameters
  real(8), intent(in) :: parFire
  
  if (ETS > parFire) then
    call random_number(pDistFireX)
  endif
end subroutine

!!!percentage of BA killed by disturbance   
subroutine intTest(pDistFire,pBAmort)
  implicit none
  real(8), intent(inout) :: pDistFire,pBAmort
    
    call random_number(pBAmort) !just some random number for now

end subroutine


!!!matrix of interactive number of rows
subroutine testAlloc(pDistFire,mat1,mat2)
  implicit none
  real(8), allocatable, intent(inout) :: mat1(:,:),mat2(:,:)
  ! real(8), intent(inout) :: mat1x(10,3),mat2x(10,3)
  real(8), intent(inout) :: pDistFire
  real(8) :: pX
  integer :: i, count1,count2
  
open(1,file="test1.txt")
  count1 = 0 
  count2 = 0
  do i = 1,10 
  write(1,*) i
  call random_number(pX)
  if(pX > pDistFire) then 
  write(1,*) "1", i,pX
    if(count1 ==0 .and. count2==0) then
      count1 = 1
      allocate(mat1(count1,3))
      mat1(1,1) = count1
      mat1(1,2) = count2
      mat1(1,3) = pX
      write(1,*) "2", i,count1,count2,pX
    elseif(count1>count2) then
    count2 = count2+2
      write(1,*) "2a", i,count1,count2,pX
    !deallocate(mat2)
    allocate(mat2(count2,3))
    mat2(1:count1,:) = mat1
    mat2(count2,1) = count1
    mat2(count2,2) = count2
    mat2(count2,3) = pX
    write(3,*) "3", i,count1,count2,pX
    deallocate(mat1)
    elseif(count2>count1 .and. count1>0) then
    count1 = count1+2
    write(3,*) "3a", i,count1,count2,pX
    !deallocate(mat1)
    allocate(mat1(count1,3))
    mat1(1:count2,:) = mat2
    mat1(count1,1) = count1
    mat1(count1,2) = count2
    mat1(count1,3) = pX
    deallocate(mat2)
    write(4,*) "4", i,count1,count2,pX
    endif
  endif
  enddo
  if(allocated(mat1)) then
    mat2 = mat1
  else
    mat1 = mat2
  endif
  
  ! mat1x(1:5,:) = mat1(1:5,:)
  ! mat2x(1:5,:) = mat2(1:5,:)
  close(1)
end subroutine



! subroutine foo(px,mat1,mat2)
! implicit none
  ! real(8), intent(inout) :: px
  ! real(8), intent(inout) :: mat1(10,3),mat2(10,3)
  ! real(8),allocatable :: mm1x(:,:),mm2x(:,:)
  
  ! call testAlloc(px,mm1x,mm2x)
  ! ! mat1 = mat1x
  ! ! mat2 = mat2x

! end subroutine 


!update parameter value as linear function of sitetype
subroutine linearUpdateParam(pars,siteType,par_New) 
  implicit none
  real (kind=8), intent(in) :: pars(2), siteType
  real (kind=8), intent(inout) :: par_New
   !!!calculate reineke parameter as a function of siteType
   par_New = pars(1) + pars(2) * siteType 
endsubroutine


!!calculate the soil moisture index to be used in the bark beatle disturbance calculations
subroutine SMIfromPRELES(GPP,fW,SMI,fAPAR)
 real (kind=8), intent(in) :: GPP(365),fW(365),fAPAR
 real (kind=8), intent(out) :: SMI
 real(kind=8) :: gpp_threshold
 integer :: startSeason, endSeason

!gpp threshold that should be considered for start and end of the season
 gpp_threshold=5.0d0 *fAPAR !!! it is divided by fAPAR to aproximate to p0

 startSeason = 1
 do while (GPP(startSeason) <= gpp_threshold .and. startSeason < 366)
  startSeason = startSeason + 1
 enddo


 endSeason = 365
 do while (GPP(endSeason) <= gpp_threshold .and. endSeason > 1)
  endSeason = endSeason - 1
 enddo

 SMI = sum(fW(startSeason:endSeason))/(endSeason-startSeason+1)

endsubroutine

!!calculate minimum fAPAR of last 15 years
subroutine minFaparCalc(fAPARtrees,nYears,minFapar,fAparFactor)
   integer, intent(in) :: nYears
   real (kind=8), intent(in) :: fAPARtrees(nYears), fAparFactor
   real (kind=8), intent(out) :: minFapar
   integer :: lastYears = 5 !lastYears are the number of years before the current years that should not be considered in the minimum fAPAR calculations
   integer :: maxYears =15 !maxYears are the total number of years before the current years that should be considered in the minimum fAPAR calculations
   integer :: firstYear
   
   if(nYears <= lastYears) minFapar = fAPARtrees(1)*fAparFactor
   if(nYears > lastYears .and. nYears <= int((maxYears - lastYears)/2)) minFapar = minval(fAPARtrees(1:(nYears-5)))*fAparFactor
   if(nYears > lastYears .and. nYears > int((maxYears - lastYears)/2)) then
   if((nYears-maxYears+1) < 0) then
    firstYear=1
    minFapar = minval(fAPARtrees(firstYear:(nYears-5)))*fAparFactor
   else
    firstYear = nYears-maxYears+1
    minFapar = minval(fAPARtrees(firstYear:(nYears-5)))  
   endif
   endif
     
endsubroutine


subroutine calcAlfar_MultiSite(siteTAlpha,species,pCrobas,nLayers,nSp,&
    nYearsFert,npar,siteTypeOrig,deltaSiteTypeFert,nSites,nYears,yearFert)

  implicit none
  integer, intent(in) :: nSp,nYearsFert,npar,nSites,nLayers,nYears,yearFert(nSites)
  real(8),intent(inout) :: siteTAlpha(nSites,nYears,nLayers,2),pCrobas(npar,nSp)
  real(8),intent(inout) :: siteTypeOrig(nSites),deltaSiteTypeFert(nSites)
  real(8),intent(inout) :: species(nSites,nLayers)
  ! integer,intent(inout) :: year!,ind(3,2)
  integer :: i,it,maxYearSim
  real(8) :: alfarUnfert(nLayers), alfarFert(nLayers),slope,interc,SiteTypeFert,siteType(nYearsFert,nLayers),CN
  real(8) :: siteTAlphaX(nYearsFert,nLayers,2)
  
  
  do i = 1,nSites
    siteTAlphaX(:,:,:) = 1.
    maxYearSim = min((nYears-yearFert(i)+1),nYearsFert)
    siteTAlphaX(1:maxYearSim,:,:) = siteTAlpha(i,yearFert(i):(maxYearSim+yearFert(i)-1),:,:)
    call calcAlfar(siteTypeOrig(i),species(i,:),pCrobas,nLayers, & 
        siteTAlphaX,nSp,nYearsFert,npar,deltaSiteTypeFert(i))

    siteTAlpha(i,yearFert(i):(maxYearSim+yearFert(i)-1),:,:) = siteTAlphaX(1:maxYearSim,:,:)

  end do
endsubroutine



!order a vector in descendete order
subroutine order_desc(v_size,vector_x,v_descendente)
  integer, intent(in) :: v_size
  real(8),intent(in) :: vector_x(v_size)
  integer,intent(out) :: v_descendente(v_size)
  LOGICAL, DIMENSION(v_size) :: mk
  real(8) :: indx(1)
  integer ix
! open(1,file="test1.txt")
mk = .TRUE.
DO ix = 1, v_size
   ! write(1,*) MAXVAL(vector_x,mk)
   indx = MAXLOC(vector_x,mask=mk,dim=1)
   v_descendente(ix) = int(indx(1))
   ! write(1,*)  MAXLOC(vector_x,mask=mk,dim=1)
   mk(MAXLOC(vector_x,mk)) = .FALSE.
END DO
! close(1)
endsubroutine

subroutine alternative_chooseThin(species,H,age, BA, density, tTapio, ftTapio, thinning,BA_thd,dens_thd,doThin,nSp)
  implicit none

  integer, intent(in) :: nSp,species
  real (kind=8),intent(out) :: thinning,BA_thd,dens_thd ! 1 for tapioTend, 2 for tapioFirstThin, 3 for tapioThin
  real (kind=8),intent(in) :: H, density,BA,age !siteType; average ETS of the site
  real (kind=8),intent(in) :: tTapio(5,nSp,2,7), ftTapio(5,nSp,3,7) 
  logical,intent(out) :: doThin
  ! real (kind=8) :: pX1(3,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning; pX(7) height limit to move on to next thinning function
  ! real (kind=8) :: pX2(2,7) !pX(1) = ETS threshold; pX(2)= densThd;  pX(3:4) height limits; pX(5:6) density after thinning; pX(7)  height limit to move on to next thinning function
    real (kind=8) :: densityCT2, ageCT1,ageCT2,hCT2, hCT1, hPCT, densityCT1, densityPCT
	real (kind=8) :: dens_thd_CT1,BA_thd_CT2,dens_thd_PCT

! log for testing 
!open(1, file = "chooseThinLog.txt")

! parameters for first thinning
  densityPCT = ftTapio(1,species,1,1) 
  hPCT = ftTapio(1,species,1,2) 
  dens_thd_PCT = ftTapio(1,species,1,3) 

! parameters for commercial thinning basend on tree density
  densityCT1 = tTapio(1,species,1,1) 
  hCT1 = tTapio(1,species,1,2)
  ageCT1 = tTapio(1,species,2,2)
  dens_thd_CT1 = tTapio(1,species,1,3) 

! parameters for commercial thinning basend on basal area
  densityCT2 = tTapio(2,species,1,1) 
  hCT2 = tTapio(2,species,1,2) 
  ageCT2 = tTapio(2,species,2,2)
  BA_thd_CT2 = tTapio(2,species,1,3) 

  doThin = .false.
  thinning = 0.
! writing test log 
!write(1, *) "first thinning: density under", densityU1, "stems/ha or dheight over", hNext1, &
!  "m; tending: density under", densityU2, "stems/ha or dheight over", hNext2, "m"

! if the stand is already thinner than the thinning result or the dominant height is over the limit, we move on to tapioThin subroutine
 if(density > densityCT2 .and. (H > hCT2 .or. (age > ageCT2 .and. age<ageCT2+5))) then
  thinning = 3. 
  if(BA_thd_CT2 < 1.) BA_thd_CT2 = BA_thd_CT2 * BA 
  BA_thd = BA_thd_CT2
  doThin = .true.
 endif
  
 if(density > densityCT1 .and. (H > hCT1 .or. (age > ageCT1 .and. age<ageCT1+5))) then
  thinning = 2. 
  if(dens_thd_CT1 < 1.) dens_thd_CT1 = dens_thd_CT1 * density
  dens_thd = dens_thd_CT1
  doThin = .true.
 endif
  
 if(density > densityPCT .and. H > hPCT) then
  thinning = 1. 
  if(dens_thd_PCT < 1.) dens_thd_PCT = dens_thd_PCT * density
  dens_thd = dens_thd_PCT 
  doThin = .true.
 endif  

end subroutine alternative_chooseThin
