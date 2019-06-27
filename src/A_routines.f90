SUBROUTINE Ffotos2(STAND_all,nClass,nSp,pCrobas,nVar,nPar,MeanLight,coeff,qcTOT)
implicit none

 integer, intent(in) :: nclass,nSp,nVar,nPar

!*****************************************************************************************
 real (kind=8), intent(in) :: pCrobas(npar,nSp)
 real (kind=8), intent(inout) :: STAND_all(nVar,nclass)
 real (kind=8), intent(inout) :: coeff(nclass) , qcTOT
!****************************************************************************************
 integer  :: ki
 real (kind=8) :: param(nPar)
 real (kind=8) :: ht(nclass),hc(nclass),h(nclass)
 real (kind=8) :: LAIe(nclass),qc(nclass),btc(nclass),LAI(nclass),N(nclass)
 real (kind=8) :: l(2*nclass),vrel(2*nclass,nclass)
 real (kind=8) :: lpt(2*nclass,nclass),lt(2*nclass)
 real (kind=8) :: bt(2*nclass), k(nclass), par_betab(nclass), rc(nclass)
 real (kind=8) :: kLAIetot, kLAItot, Atot
 real (kind=8), intent(inout) :: MeanLight(nclass)
 real (kind=8) :: x1,x2,apuJ,apuI
	   integer :: iclass,i2,i1,species,nv				!!**!! nv defined as integer
       integer :: i, j, ii(2*nclass), iapu
 real (kind=8) :: apu, b1,  qctot0, qctot1, wwx, dc, e1
!****************************************************************************************
 real (kind=8) :: pi = acos(-1.)
 
 MeanLight = 0.
 coeff = 0.
 qcTOT = 0.
 
 do i = 1,nclass
	 species = int(stand_all(4,i))
     param = pCrobas(:,species)
     qc(i) = 0.
     
     ht(i) = STAND_all(11,i)   ! H
     hc(i) = STAND_all(14,i)   ! Hc
     h(i) = ht(i) - hc(i)        ! Lc
     LAIe(i) = STAND_all(19,i) ! leff
     k(i) = PARAM(4)               ! k 
     LAI(i) = STAND_all(33,i) * PARAM(3) / 10000.   ! WF_stand * sla
     ! par_betab(i) = PARAM(17)   ! betab
     rc(i) = STAND_all(15,i)/2.         ! rc
     N(i) = STAND_all(17,i) / 10000.   ! N per m2
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

!	ii-table sorts the l-table indeces so that later the corresponding "locations" for hc and ht values can be located
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

!				if(ht(j).gt.l(i).and.hc(j).le.l(i)) dc=dc+lt(i)
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
!    qctot = qctot1
     
!     if(stand_P(7) > 150) then
!        continue
!     endif
     
     
! calculate weights - on the basis of qcTOT1 but all downscaled if qcTOT0 < qcTOT1
!	       
!do i = 1,nclass
     
   if(qcTOT1.gt.0.) then
  	coeff  = qc / qcTOT1 * qcTOT / qcTOT1 ! weight

!      	coeff_SP  = qc(2) / qcTOT1 * qcTOT / qcTOT1 ! 

!        coeff_B  = qc(3) / qcTOT1 * qcTOT / qcTOT1  ! 



!!!!FMadded
!        qcTOT = qcTOT * (coeff_P+coeff_SP+coeff_B)

!        coeff_P = coeff_P/(coeff_P+coeff_SP+coeff_B)
!        coeff_SP = coeff_SP/(coeff_P+coeff_SP+coeff_B)
!        coeff_B = coeff_B/(coeff_P+coeff_SP+coeff_B)
!!!!


!    if(kLAIetot.gt.0.) then
!       	 coeff_P  = k(1)*LAIe(1) / kLAIetot ! weight

!        	coeff_SP  = k(2)*LAIe(2) / kLAIetot ! 
           
!             coeff_B  = k(3)*LAIe(3) / kLAIetot ! 


!   else
!       coeff_P = 1./3.
!       coeff_SP = 1./3.
!       coeff_B = 1./3.
   
! end do	
    endif      



 !     write(60,*)qcTOT, qcTOT1, qc(1), qc(2), qc(3)

!81    	continue


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



!	if(species==1 .OR. species==3)then
	  pp = p
	  qq = q
!	endif


!!!check with Annikki
!	if(species==4)then
!! 	  hmax = amin1(0.9*Lc,hmax0+0.3*Lc)
!	  pp = p
!! *** If-lause lis?tty 2011/10/14 by TL
!	  if (Lc .gt. hmax0) then
!       	  qq = 0.18 * Lc - 0.6
!          else
!          qq = 0.18 * hmax0 - 0.6
!          endif
!	endif

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



subroutine preles(weather,DOY,fAPAR,prelesOut,pars,GPP,ET,SW,etmodel)!,p0)

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
		transp, evap, fWE) BIND(C)
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
     integer(c_int) :: etmodel, LOGFLAG,NofDays
     integer(c_int) :: day(365)
     real ( c_double ) :: transp(365), evap(365), fWE(365)
   END SUBROUTINE call_preles
 END INTERFACE

 real (kind=8), intent(in) :: weather(365,5),fAPAR(365)
 real (kind=8), intent(out) :: prelesOut(16)!,p0
 real (kind=8), intent(inout) :: pars(30)
 integer, intent(in):: DOY(365), etmodel

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
		transp, evap, fWE)


prelesOut(1) = sum(GPP(1:nDays))
prelesOut(2) = sum(ET(1:nDays))
prelesOut(3) = SW(nDays)
prelesOut(4) = SOG(nDays)
prelesOut(5) = fS(nDays)
prelesOut(6) = fD(nDays)
prelesOut(7) = fW(nDays)
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
		! using the formula 0 = x'(t) = A*x + b => x = -A^-1*b
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
! !		REAL (kind=8) :: p1 = -2.653,p2 = -2.948,p3 = -3.324,p4 = .055,p5 = .059,p6 = .135,p7 = -0.03

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
	! AWENH(5) = 0.
    END SUBROUTINE compAWENH

!! Note for Birch Betula pubenscens and brown leaves is used
!    SUBROUTINE foliageAWENH(Lf,folAWENH)
!        IMPLICIT NONE
!        INTEGER,PARAMETER :: n = 5, nSp=3
!        REAL (kind=8),DIMENSION(nSp,n),INTENT(OUT) :: folAWENH
!        REAL (kind=8),DIMENSION(nSp),INTENT(IN) :: Lf
!!	folAWENH = 0.
!	folAWENH(1,1) = 0.518*Lf(1)
!	folAWENH(2,1) = 0.4826*Lf(2)
!	folAWENH(3,1) = 0.4079*Lf(3)
!	folAWENH(1,2) = 0.1773*Lf(1)
!	folAWENH(2,2) = 0.1317*Lf(2)
!	folAWENH(3,2) = 0.198*Lf(3)
!	folAWENH(1,3) = 0.0887*Lf(1)
!	folAWENH(2,3) = 0.0658*Lf(2)
!	folAWENH(3,3) = 0.099*Lf(3)
!	folAWENH(1,4) = 0.216*Lf(1)
!	folAWENH(2,4) = 0.3199*Lf(2)
!	folAWENH(3,4) = 0.2951*Lf(3)
!   END SUBROUTINE foliageAWENH

!! Branches are here
!! It seems that there is only valiues for pine (these are applied for others as well)
!    SUBROUTINE branchesAWENH(Lb, fbAWENH)
!        IMPLICIT NONE
!        INTEGER,PARAMETER :: n = 5, nSp=3
!        REAL (kind=8),DIMENSION(nSp,n),INTENT(OUT) :: fbAWENH
!        REAL (kind=8),DIMENSION(nSp),INTENT(IN) :: Lb
!!	fbAWENH = 0.
!	fbAWENH(:,1) = 0.47466*Lb
!	fbAWENH(:,2) = 0.019012*Lb
!	fbAWENH(:,3) = 0.078308*Lb
!	fbAWENH(:,4) = 0.430248*Lb
!   END SUBROUTINE branchesAWENH


!    SUBROUTINE stemAWENH(Lst, stAWENH)
!        IMPLICIT NONE
!        INTEGER,PARAMETER :: n = 5, nSp=3
!        REAL (kind=8),DIMENSION(nSp,n),INTENT(OUT) :: stAWENH
!        REAL (kind=8),DIMENSION(nSp),INTENT(IN) :: Lst
!!	stAWENH = 0.
!	stAWENH(1,1) = 0.5*(0.66+0.68)*Lst(1)
!  	stAWENH(2,1) = 0.5*(0.63+0.7)*Lst(2)
!  	stAWENH(3,1) = 0.5*(0.65+0.78)*Lst(3)
!  	stAWENH(1,2) = 0.5*(0.03+0.015)*Lst(1)
!  	stAWENH(2,2) = 0.5*(0.03+0.005)*Lst(2)
!  	stAWENH(3,2) = 0.5*(0.03+0)*Lst(3)
!  	stAWENH(1,3) = 0.5*(0+0.015)*Lst(1)
!  	stAWENH(2,3) = 0.5*(0+0.005)*Lst(2)
!  	stAWENH(3,3) = 0
!  	stAWENH(1,4) = 0.5*(0.28+0.29)*Lst(1)
!  	stAWENH(2,4) = 0.5*(0.33+0.28)*Lst(2)
!  	stAWENH(3,4) = 0.5*(0.22+0.33)*Lst(3)
!   END SUBROUTINE stemAWENH


