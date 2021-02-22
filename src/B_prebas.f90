
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine bridging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prebas(nYears,nLayers,nSp,siteInfo,pCrobas,initVar,thinning,output,nThinning,maxYearSite,fAPAR,initClearcut,&
		fixBAinitClarcut,initCLcutRatio,ETSy,P0y,weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,pAWEN,weatherYasso,&
		litterSize,soilCtotInOut,defaultThin,ClCut,energyCut,inDclct,&
		inAclct,dailyPRELES,yassoRun,energyWood,tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,GVrun) !energyCut


implicit none

 integer, parameter :: nVar=54,npar=43, inttimes = 1!, nSp=3
 real (kind=8), parameter :: pi = 3.1415927, t=1. , ln2 = 0.693147181
 real (kind=8), parameter :: energyRatio = 0.7, harvRatio = 0.9 !energyCut
 ! logical steadystate_pred= .false.
!define arguments
 integer, intent(in) :: nYears,nLayers,nSp
 real (kind=8), intent(in) :: weatherPRELES(nYears,365,5)
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(inout) :: pPRELES(30),tapioPars(5,2,3,20),thdPer,limPer
 real (kind=8), intent(inout) :: tTapio(5,3,2,7), ftTapio(5,3,3,7)
 real (kind=8), intent(inout) :: thinning(nThinning,9)
 real (kind=8), intent(inout) :: initClearcut(5)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
 real (kind=8), intent(in) :: pCrobas(npar,nSp),pAWEN(12,nSp)
 integer, intent(in) :: maxYearSite
 real (kind=8), intent(in) :: defaultThin,ClCut,energyCut,yassoRun,fixBAinitClarcut		!energyCut

 real (kind=8), intent(in) :: inDclct(nSp),inAclct(nSp)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning
 !!!ground vegetation variables
 real (kind=8) :: AWENgv(4)  !!!ground vegetation
 integer, intent(in) :: gvRun			!!!ground vegetation
 real (kind=8), intent(inout) :: fAPAR(nYears),GVout(nYears,3) !fAPAR_gv,litGV,photoGV			!!!ground vegetation
 real (kind=8), intent(inout) :: dailyPRELES((nYears*365),3)
 real (kind=8), intent(inout) :: initVar(7,nLayers),P0y(nYears,2),ETSy(nYears),initCLcutRatio(nLayers)!
 real (kind=8), intent(inout) :: siteInfo(10)
 real (kind=8), intent(inout) :: output(nYears,nVar,nLayers,2), energyWood(nYears,nLayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nYears,5,3,nLayers),soilCtotInOut(nYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(inout) :: pYasso(35), weatherYasso(nYears,3),litterSize(3,nSp) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: prelesOut(16),fAPARsite
 real (kind=8) :: leac=0 !leaching parameter for Yasso
 real (kind=8),DIMENSION(nLayers,5) :: fbAWENH,folAWENH,stAWENH
 real (kind=8),DIMENSION(nLayers) :: Lb,Lf,Lst
! real (kind=8),DIMENSION(nLayers) :: speciesIDs
 real (kind=8),DIMENSION(nLayers) :: valX
 integer,DIMENSION(nLayers) :: layerX
 real (kind=8) :: STAND(nVar),STAND_tot(nVar),param(npar)!, output(nYear,nSites,nVar)
 integer :: i, ij, ijj,species,layer,nSpec,ll! tree species 1,2,3 = scots pine, norway spruce, birch

 real (kind=8) :: p0_ref, ETS_ref,P0yX(nYears,2)
 integer :: time, ki, year,yearX,Ainit, countThinning,domSp(1)
 real (kind=8) :: step, totBA

 real (kind=8) :: stand_all(nVar,nLayers)
 real (kind=8) :: outt(nVar,nLayers,2)
 real (kind=8) :: modOut((nYears+1),nVar,nLayers,2)
 real (kind=8) :: soilC((nYears+1),5,3,nLayers),soilCtot((nYears+1))
 real (kind=8) :: par_phib,par_phic,par_alfat,par_alfar1,par_alfar2,par_alfar3,par_alfar4
 real (kind=8) :: par_alfar5,par_etab,par_k,par_vf,par_vr,par_sla,par_mf,par_mr,par_mw,par_vf0, mrFact,par_vr0
 real (kind=8) :: par_z,par_rhos,par_cR, par_x, Light,MeanLight(nLayers),par_mf0,par_mr0,par_mw0,par_zb!par_zb= ratio between dead branches l and mean pipe leanght
 real (kind=8) :: par_sarShp, par_S_branchMod
 real (kind=8) :: par_rhof, par_rhor, par_rhow, par_c, par_beta0, par_betab, par_betas
 real (kind=8) :: par_s1, par_p0, par_ksi, par_cr2,par_kRein,Rein, c_mort
 real (kind=8) :: BA, dA, dB, reineke(nLayers), dN, wf_test,par_thetaMax, par_H0max,par_kH, par_gamma,par_H0
 real (kind=8) :: par_rhof0, par_rhof1, par_rhof2, par_aETS,dHcCum,dHCum,pars(30), thinningType = 0.

!management routines
 real (kind=8) :: A_clearcut, D_clearcut, BAr(nLayers), BA_tot,BA_lim, BA_thd, ETSthres = 1000
 real (kind=8) :: dens_thd, dens_lim, Hdom_lim

!define varibles
 real (kind=8) :: LAT, LONG, sitetype, P0, age, meantemp, mintemp, maxtemp, rainfall, ETS
 real (kind=8) :: H, D, B, Hc, Cw, Lc, N, Ntree, Ntot,dNtot
 real (kind=8) :: wf_treeKG, wf_STKG, sar_con, sar_ell, rc, ppow, sar,W_stem
 real (kind=8) :: lproj, leff,laPer_sar, keff, slc
 real (kind=8) :: hb, A, B2,beta0, beta1,beta2, betab!,betas
 real (kind=8) :: c,dHc,dH,dLc,g0,g1,g2,g3,g4,g5
 real (kind=8) :: npp, p_eff_all, gammaC, betaC, W_c, W_s,Wdb,Wsh
 real (kind=8) :: p_eff, par_alfar,p,gpp_sp
 real (kind=8) :: s0,par_s0scale,par_sla0,par_tsla
 real (kind=8) :: weight, dNp,dNb,dNs
 real (kind=8) :: W_wsap, respi_m, respi_tot, V_scrown, V_bole, V,Vold
 real (kind=8) :: coeff(nLayers), denom,W_froot,W_croot, lit_wf,lit_froot
 real (kind=8) :: S_wood,Nold, Nthd, S_branch,S_fol,S_fr,W_branch
 real (kind=8) :: W_stem_old,wf_STKG_old,W_bh, W_crh,W_bs, W_crs,dW_bh,dW_crh,dWdb,dWsh

!fix parameters
 real (kind=8) :: qcTOT0,Atot,fAPARprel(365)
!v1 version definitions
 real (kind=8) :: theta,Tdb=10.,f1,f2, Gf, Gr,mort
 real (kind=8) :: ETSmean, BAtapio(2), tapioOut(3)
 logical :: doThin, early = .false.
 real (kind=8) :: Hdom
  ! open(1,file="test1.txt")
  ! open(2,file="test2.txt")

!###initialize model###!
energyWood = 0.
fbAWENH = 0.
folAWENH = 0.
stAWENH = 0.
yearX = 0
modOut = 0.
modOut(1,:,:,:) = output(1,:,:,:)
soilC = 0.
countThinning = 1
pars = pPRELES
pars(1:3) = siteInfo(8:10)
soilC(1,:,:,:) = soilCinout(1,:,:,:)
pars(24) = siteInfo(4)!SWinit
pars(25) = siteInfo(5)!CWinit
pars(26) = siteInfo(6) !SOGinit
pars(27) = siteInfo(7) !Sinit
Ainit = initClearcut(5)
P0yX = P0y
Reineke(:) = 0.
ETSmean = sum(ETSy)/nYears

 modOut(:,1,:,1) = siteInfo(1)  !! assign siteID 
 ! modOut(1,2,:,1) = initVar(8,:)	!! assign initial gammaC values !!newX
 modOut(1,11,:,1) = initVar(3,:)
 modOut(1,16,:,1) = initVar(7,:)
 ! modOut(1,12,:,1) = initVar(4,:)
 ! modOut(1,13,:,1) = initVar(5,:)
 modOut(1,14,:,1) = initVar(6,:)
 ! modOut(1,17,:,1) = modOut(1,13,:,1)/(pi*((modOut(1,12,:,1)/2/100)**2))
 ! modOut(1,35,:,1) =  modOut(1,13,:,1)/modOut(1,17,:,1)
 modOut(:,3,:,1) = siteInfo(3);sitetype = siteInfo(3)! assign site type
 soilCtot(1) = sum(soilC(1,:,:,:)) !assign initial soilC
 do i = 1,nLayers
  modOut(:,4,i,1) = initVar(1,i)  ! assign species
  modOut(:,7,i,1) = initVar(2,i) ! assign initAge !age can be made species specific assigning different ages to different species
  modOut(1,39,i,1) = sum(soilC(1,:,:,i)) !assign initial soilC
  modOut(:,5,i,1) = ETSy! assign ETS
  modOut(:,6,i,1) = P0yX(:,2)	! assign P0
  modOut(1,12,i,1) = initVar(4,i)
  modOut(1,13,i,1) = initVar(5,i)
  if(modOut(1,12,i,1) > 0.) then
	modOut(1,17,i,1) = modOut(1,13,i,1)/(pi*((modOut(1,12,i,1)/2/100)**2))
	modOut(1,35,i,1) = modOut(1,13,i,1)/modOut(1,17,i,1)
  else
	modOut(1,17,i,1) = 0. 
	modOut(1,35,i,1) = 0.
  endif
			! wf_STKG = N * wf_treeKG  !!newX
			! W_froot = par_rhor * A * N  !!to check  !!newX
			! W_wsap = par_rhow * A * N * (beta1 * H + beta2 * Hc) !!newX
			! W_c = par_rhow * A * N * Hc !sapwood stem below Crown
			! W_s = par_rhow * A * N * betas * Lc !sapwood stem within crown
			! W_branch =  par_rhow * A * N * betab * Lc !branches biomass
 	        ! W_croot = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)	!coarse root biomass
			! Wsh = (A+B+sqrt(A*B)) * Hc * par_rhow * N/2.9 - Wc !initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
			! !initialize Wdb dead branches biomass
		! if(par_S_branchMod .eq. 1.) then
			! Wdb = Tdb * W_branch * ((0.0337+0.000009749*N)*exp(-0.00456*D**2)+0.00723)
		! else
			! Wdb = Tdb * W_branch *((-0.00513+0.000012*N)*exp((0.00000732-0.000000764*N)*D**2)+0.00467)
		! endif
 enddo

!######!

do year = 1, (nYears)

! write(1,*) year,yearX,Ainit,initClearcut(5)
  if(year==int(min(yearX,nYears))) then
   Ainit = int(min(Ainit, Ainit + nYears - yearX))
      totBA = sum(modOut((year-Ainit-1),13,:,1))
   do ijj = 1,nLayers
     species = int(modOut(year,4,ijj,1))  ! read species
	 if(fixBAinitClarcut==1) then
	  modOut(year,13,ijj,1) = initClearcut(3) * initCLcutRatio(ijj)
	 else
      modOut(year,13,ijj,1) = initClearcut(3) * modOut((year-Ainit-1),13,ijj,1)/ totBA
     endif
	 modOut(year,11,ijj,1) = initClearcut(1)
	 modOut(year,2,ijj,1) = 0.
     modOut(year,12,ijj,1) = initClearcut(2)
     modOut(year,14,ijj,1) = initClearcut(4)
	 modOut(year,16,ijj,1) = pCrobas(38,species)/pCrobas(15,species) * (initClearcut(1) - &
		initClearcut(4))**pCrobas(11,species)!A = p_ksi/p_rhof * Lc**p_z
	if(modOut(year,12,ijj,1) > 0.) then
	  modOut(year,17,ijj,1) = modOut(year,13,ijj,1)/(pi*((modOut(year,12,ijj,1)/2/100)**2))
	  modOut(year,35,ijj,1) =  modOut(year,13,ijj,1)/modOut(year,17,ijj,1)
    else
	  modOut(year,17,ijj,1) = 0. 
	  modOut(year,35,ijj,1) = 0. 
    endif
     ! modOut(year,35,ijj,1) = modOut(year,13,ijj,1) / modOut(year,17,ijj,1)
	 
	   siteType = siteInfo(3)
  !!set parameters
  par_betab = pCrobas(13,int(initVar(1,ijj)))
  par_x = pCrobas(19,int(initVar(1,ijj)))
  par_beta0 = pCrobas(12,int(initVar(1,ijj)))
  par_betas = pCrobas(14,int(initVar(1,ijj)))
  par_mf = pCrobas(8,int(initVar(1,ijj)))
  par_mr = pCrobas(9,int(initVar(1,ijj)))
  par_mw = pCrobas(10,int(initVar(1,ijj)))
  par_alfar = pCrobas(int(20+min(siteType,5.)),int(initVar(1,ijj)))
  par_c = pCrobas(7,int(initVar(1,ijj)))
  par_rhof = pCrobas(15,int(initVar(1,ijj)))
  par_rhor = par_alfar * par_rhof
  par_rhow = pCrobas(2,int(initVar(1,ijj)))
  par_S_branchMod = pCrobas(27,int(initVar(1,ijj)))
  gammaC = 0. !initVarX(8,)
  ! Tbd = 10 !!!!to include in the parameters
  
  !!set variables
  A = modOut(year,16,ijj,1)
  ba = modOut(year,13,ijj,1)
  d = modOut(year,12,ijj,1)
  N = ba/(pi*((d/2/100)**2))
  h = modOut(year,11,ijj,1)
  hc = modOut(year,14,ijj,1)
  B = ba/N
  Lc = h - hc
  betab =  par_betab * Lc**(par_x-1)
  beta0 = par_beta0
  beta1 = (beta0 + betab + par_betas) 
  beta2 = 1. - betab - par_betas 		
  betaC = (beta1 + gammaC * beta2) / par_betas
  wf_STKG = par_rhof * A * N
  W_froot = par_rhor * A * N  !!to check  ##newX
  W_wsap = par_rhow * A * N * (beta1 * h + beta2 * hc) 
  W_c = par_rhow * A * N * hc !sapwood stem below Crown
  W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
  W_branch =  par_rhow * A * N * betab * Lc !branches biomass
  Wsh = 0.!max((A+B+sqrt(A*B)) * hc * par_rhow * N/2.9 - W_c,0.) !initialize heart wood, only stem considered. W_bole (total biomass below crown)  - Wc
!  W_croot = par_rhow * beta0 * A * h * N !coarse root biomass
  W_croot = Lc * beta0 * A / par_betas * N + (W_c + Wsh) * beta0 !coarse root biomass
  !initialize Wdb dead branches biomass
  Wdb = 0.
  W_stem = W_c + W_s + Wsh
  V = W_stem / par_rhow
  
  modOut(year,33,ijj,1) = wf_STKG
  modOut(year,25,ijj,1) = W_froot
  modOut(year,47,ijj,1) = W_wsap
  modOut(year,48,ijj,1) = W_c
  modOut(year,49,ijj,1) = W_s
  modOut(year,24,ijj,1) = W_branch
  modOut(year,32,ijj,1) = W_croot
  modOut(year,50,ijj,1) = Wsh
  modOut(year,51,ijj,1) = Wdb
  modOut(year,31,ijj,1) = W_stem
  modOut(year,30,ijj,1) = V
  
   enddo
   do ki = 1,int(Ainit)
    do ijj = 1,nLayers
     modOut((year-Ainit+ki),7,ijj,1) = ki !#!#
     modOut((year-Ainit+ki),4,ijj,1) = initVar(1,ijj) !#!#
    enddo
   enddo	
	yearX = 0
	
	
  endif

  stand_all = modOut(year,:,:,1)

  step = 1. / float(inttimes)
  outt(:,:,2) = 0.

 !----------------------------------
 !PHOTOSYNTHESIS MODEL PART 1
 !----------------------------------

  do time = 1, inttimes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! do ki = 1, nSites
  ! calculate self-thinning using all tree classes
 if(time==inttimes)then
    valX = STAND_all(11,:)
	do ij = 1, nLayers
	  domSp = maxloc(valX)
 	  layerX(ij) = int(domSp(1))
	  valX(layerX(ij)) = -999.
	
     Ntot = sum(STAND_all(17,layerX(1:ij)))
     B = sum(STAND_all(35,layerX(1:ij))*STAND_all(17,layerX(1:ij)))/Ntot   !!!!!!!!!#####changed
	 ! B = STAND_all(35,layerX(1:ij))/STAND_all(17,layerX(1:ij))   !!!!!!!!!#####changed
     if(Ntot>0.) then
         Reineke(layerX(ij)) = Ntot*(sqrt(B*4/pi)*100./25.)**(1.66)
     else
         Reineke(layerX(ij)) = 0.
     endif
	enddo
 endif

do ij = 1 , nLayers 		!loop Species

 STAND=STAND_all(:,ij)
 species = int(stand(4))
 param = pCrobas(:,species)

 par_cR= param(1)
 par_rhow= param(2)
 par_sla = param(3)
 par_k =param(4)
 par_vf0 =param(5)
 par_vr0 =param(6)
 par_c=param(7)
 par_mf0=param(8)
 par_mr0=param(9)
 par_mw0=param(10)
 par_z=param(11)
 par_beta0=param(12)
 par_betab=param(13)
 par_betas = param(14)
 par_rhof2 = param(15)
 par_s1 = param(16)
 par_kRein = param(17)
 par_s0scale = param(18)
 par_x = param(19)
 par_aETS = param(20)
 par_alfar1 =param(21)
 par_alfar2 =param(22)
 par_alfar3 =param(23)
 par_alfar4 = param(24)
 par_alfar5 = param(25)
 par_sarShp = param(26) !Shape surface area of the crown: 1.= cone; 2.=ellipsoide
 par_S_branchMod = param(27) !model for branch litter model
 p0_ref = param(29)
 ETS_ref = param(30)
 par_thetaMax = param(31)
 par_H0max = param(32)
 par_gamma = param(33)
 par_rhof1 = 0.!param(20)
 par_Cr2 = 0.!param(24)
 par_kH = param(34)
 par_sla0 = param(39)
 par_tsla = param(40)

! do siteNo = 1, nSites  !loop sites

if (year > maxYearSite) then
  STAND(2) = 0. !!newX
  STAND(8:21) = 0. !#!#
  STAND(23:37) = 0. !#!#
  STAND(42:44) = 0. !#!#
  STAND(47:nVar) = 0. !#!#
else
! initialize site variables
!  sitetype = STAND(3)

  age = STAND(7)
  H = STAND(11)
  D = STAND(12)
  BA = STAND(13)
  Hc = STAND(14)
  N = BA/(pi*((D/2/100)**2))
  B = BA/N! * par_ops2
  A = stand(16)
!  Cw = STAND(15)
  Lc = H - Hc
  hb = par_betab * Lc ** par_x
  Cw = 2. * hb
  STAND(15) = Cw
  ! STAND(16) = LC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ETS = STAND(5) !!##!!2
  Light = STAND(36)
  V = stand(30)
  mort = stand(40)
  par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)
  
if (N>0.) then

  par_rhof0 = par_rhof1 * ETS_ref + par_rhof2
  par_rhof = par_rhof1 * ETS + par_rhof2
  par_vf = par_vf0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref)
  par_vr = par_vr0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref) !!!new version

 !calculate derived variables
  rc = Lc / (H-1.3) !crown ratio
  ! A = rc * B
  wf_treeKG = par_rhof * A
  ! par_ksi = wf_treeKG / (Lc ** par_z)
  wf_STKG = wf_treeKG * N !needle mass per STAND in units C
  ppow=1.6075

! Mortality - use Reineke from above
!      if((Reineke(siteNo) > par_kRein .OR. Light < par_cR) .and. siteThinning(siteNo) == 0) then !
     if(time==inttimes) then
      Rein = Reineke(ij) / par_kRein

      if(Rein > 1.) then
           dN = - 0.02 * N * Rein
      else
           dN = 0.
      endif
	  if(mort == 888.) then
		dN = min(dN,-(0.03*N)) !!!!reduce try density of 3% if there is no growth
		mort=0.
		stand(40) = 0.
	  endif
      Vold = V
      Nold = N
      ! if(N < 5.) N = 0.0

      N = max(0.0, N + step*dN)

	  if (dN<0. .and. Nold>0.) then
			W_wsap = stand(47)
			W_froot = stand(25)
			W_c = stand(48) !sapwood stem below Crown
			W_s = stand(49) !sapwood stem within crown
			W_branch = stand(24) !branches biomass
			W_croot = stand(32) !coarse root biomass
			Wsh = stand(50)
			Wdb = stand(51)
			W_bh = stand(53)
			W_crh = stand(54)
			W_stem = stand(31)
	S_branch = (W_branch + W_croot*0.83 + Wdb) * min(1.,-dN*step/Nold) 
	S_wood = (W_croot*0.17 + W_stem) * min(1.,-dN*step/Nold)
    S_fol = wf_STKG * min(1.,-dN*step/Nold) !foliage litterfall
    S_fr  = W_froot * min(1.,-dN*step/Nold)  !fine root litter
			W_wsap = W_wsap * N/Nold
			W_froot = W_froot * N/Nold
			W_c = W_c * N/Nold
			W_s = W_s * N/Nold
			W_branch = W_branch * N/Nold
			W_croot = W_croot * N/Nold
			Wsh = Wsh * N/Nold
			W_bh = W_bh * N/Nold
			W_crh = W_crh * N/Nold
			Wdb = Wdb * N/Nold
			W_stem = W_stem * N/Nold
			V = V * N/Nold
			BA = BA * N/Nold
			wf_STKG = wf_STKG * N/Nold
  STAND(24) = W_branch
  STAND(25) = W_froot
  STAND(26) = S_fol
  STAND(27) = S_fr
  STAND(28) = S_branch
  STAND(29) = S_wood
  STAND(31) = W_stem
  STAND(32) = W_croot
  STAND(42) = Vold - V!* min(1.,-dN*step/Nold)
  STAND(47) = W_wsap
  STAND(48) = W_c
  STAND(49) = W_s
  STAND(50) = Wsh
  STAND(53) = W_bh
  STAND(54) = W_crh
  STAND(51) = Wdb

      else
  STAND(26) = 0.
  STAND(27) = 0.
  STAND(28) = 0.
  STAND(29) = 0.
  STAND(42) = 0.
      endif

	  !!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
	  if(dN<0.) then
	  modOut((year+1),8,ij,1) = modOut((year+1),8,ij,1) + V* min(1.,-dN*step/N)
	    do ijj = 1,(nyears-year)
			modOut((year+ijj+1),8,ij,1) = modOut((year+ijj+1),8,ij,1) + (V/N) * (-dN*step) * &
				exp(-exp(pCrobas(35,species) + pCrobas(36,species)*ijj + pCrobas(37,species)*D + 0.))
		enddo
	  end if
	 endif
!!keff calculations
  if(par_sarShp==0.) then
	keff = par_k
  else
	  !Surface area of the crown
	  sar_ell= 4. * pi *  (((((Lc/2)**ppow)*((Cw/2)**ppow)+((Lc/2)**ppow)*((Cw/2)**ppow)+((Cw/2)**ppow)*((Cw/2)**ppow))/3)**(1/ppow))!surface area per tree
	  sar_con = pi * ((0.8*hb)**2) * (1 + sqrt(1 + 1 / (((0.8*hb) / Lc)**2))) !surface area per tree
	  !Ellipsoid for pine and birch, cone for spruce
	  if(par_sarShp==1.) then
	   sar = sar_ell
	  else
	   sar = sar_con
	   slc = 0.005
	  end if

	  !specific leaf area ------------------------------------------------
	  laPer_sar = wf_treeKG * par_sla / sar !leaf area per tree  /  crown surface area
	  keff = 0.4 * (1. - exp( - par_k / 0.4 * laPer_sar)) / laPer_sar !effective extinction coefficient    }
  endif
  
  !projected leaf area on the STAND -----------------------------------
  if (wf_STKG>0.) then
   lproj = par_sla * wf_STKG / 10000.
  else
   lproj = 0.
  end if
  !weight per tree STAND$species stratum ------------------------------------
  leff= (keff/par_k)*(wf_STKG*par_sla) / 10000. !effective lai

  STAND(7) = age
  STAND(19) = leff
  STAND(20) = keff
  STAND(21) = lproj
  STAND(23) = weight
  ! STAND(24) = W_branch
  ! STAND(25) = W_froot
  STAND(11) = H
  STAND(12) = D
  STAND(13) = BA ! * par_ops2
  STAND(14) = Hc
  STAND(15) = Cw
  STAND(17) = N
  STAND(33) = wf_STKG
  STAND(34) = wf_treeKG
  STAND(35) = B
  STAND(30) = V
else
  STAND(2) = 0. !#!#
  STAND(8:21) = 0. !#!#
  STAND(23:37) = 0. !#!#
  STAND(42:44) = 0. !#!#
  STAND(47:nVar) = 0. !#!#
endif
endif
! end do !!!!!!!end loop sites

 STAND_all(:,ij)=STAND
end do !!!!!!!end loop layers

!!!calculate species weight for photosynthesis
!do siteNo = 1, nSites

if (year <= maxYearSite) then
	nSpec = nSp
	ll = nLayers
    call Ffotos2(STAND_all,nLayers,nSpec,pCrobas,&
		nVar,nPar,MeanLight,coeff,fAPARsite)
   STAND_all(36,:) = MeanLight
   STAND_all(23,:) = coeff
! fAPARsite=0.7
   if(sum(modOut(year,11,:,1)) == 0. .and. yearX == 0) then
	if((nYears-year)<10) then
		! if(initClearcut(5)>0.) then
			! Ainit = initClearcut(5)
		! else
			Ainit = nint(6. + 2*sitetype - 0.005*modOut(year,5,1,1) + 2.25)
		! endif
	else
		! if(initClearcut(5)>0.) then
			! Ainit = 11!initClearcut(5)
		! else
			Ainit = nint(6. + 2*sitetype - 0.005*(sum(modOut(year:(year+9),5,1,1))/10) + 2.25)
		! endif
	endif
	yearX = Ainit + year 
!	initClearcut(5) = Ainit
   endif

   fAPARprel(:) = fAPARsite
   fAPAR(year) = fAPARsite
      
   
   call preles(weatherPRELES(year,:,:),DOY,fAPARprel,prelesOut, pars, &
		dailyPRELES((1+((year-1)*365)):(365*year),1), &  !daily GPP
		dailyPRELES((1+((year-1)*365)):(365*year),2), &  !daily ET
		dailyPRELES((1+((year-1)*365)):(365*year),3), &  !daily SW
		etmodel)		!type of ET model

   STAND_all(22,:) = prelesOut(2)  	!ET
   ! STAND_all(40,:) = prelesOut(15)  !aSW
   ! STAND_all(41,:) = prelesOut(16)  !summerSW 

   pars(24) = prelesOut(3);siteInfo(4) = prelesOut(3)!SWinit
   pars(25) = prelesOut(13); siteInfo(5) = prelesOut(13) !CWinit
   pars(26) = prelesOut(4); siteInfo(6) = prelesOut(4) !SOGinit
   pars(27) = prelesOut(14); siteInfo(7) = prelesOut(14) !Sinit

   STAND_all(10,:) = prelesOut(1)/1000. ! Photosynthesis in g C m-2 (converted to kg C m-2)


   !!!ground vegetation
   !!!fapar_gv compute fapar, biomasses and litter of gv with routine
   if(gvRun==1) then
	if(fAPARsite>0.) then
     call fAPARgv(fAPARsite,ETSmean,siteType,GVout(year,1),GVout(year,2),sum(P0yX(:,1))/nYears,AWENgv) !reduced input output
     GVout(year,3) = prelesOut(1) * GVout(year,1)/fAPARsite! Photosynthesis in g C m-2 
     ! GVout(year,4) = GVout(year,3)*0.5 !where to put those two variables
   	 ! STAND_all(26,1) = STAND_all(26,1) + GVout(year,2)	!add !!!ground vegetation to the 1st layer
    elseif(fAPARsite==0.) then
	 call fAPARgv(fAPARsite,ETSmean,siteType,GVout(year,1),GVout(year,2),sum(P0yX(:,1))/nYears,AWENgv) !reduced input output
     fAPARprel(:) = GVout(year,1)
    !!!fapar_gv run preles for ground vegetation
     call preles(weatherPRELES(year,:,:),DOY,fAPARprel,prelesOut, pars, &
		dailyPRELES((1+((year-1)*365)):(365*year),1), &  !daily GPP
		dailyPRELES((1+((year-1)*365)):(365*year),2), &  !daily ET
		dailyPRELES((1+((year-1)*365)):(365*year),3), &  !daily SW
		etmodel)		!type of ET model
      GVout(year,3) = prelesOut(1) ! Photosynthesis in g C m-2 
      ! GVout(year,4) =  GVout(year,3)*0.5 
	endif
   endif

endif
!enddo !! end site loop

do ij = 1 , nLayers
 STAND=STAND_all(:,ij)
 species = int(stand(4))
 param = pCrobas(:,species)

 par_cR=param(1)
 par_rhow=param(2)
 par_sla =param(3)
 par_k =param(4)
 par_vf0 =param(5)
 par_vr0 =param(6)
 par_c=param(7)
 par_mf0=param(8)
 par_mr0=param(9)
 par_mw0=param(10)
 par_z=param(11)
 par_beta0=param(12)
 par_betab=param(13)
 par_betas = param(14)
 par_rhof2 = param(15)
 par_s1 = param(16)
 par_kRein = param(17)
 par_s0scale = param(18)
 par_x = param(19)
 par_aETS = param(20)
 par_alfar1 =param(21)
 par_alfar2 =param(22)
 par_alfar3 =param(23)
 par_alfar4 =param(24)
 par_alfar5 =param(25)
 par_sarShp = param(26) !Shape surface area of the crown: 1.= cone; 2.=ellipsoide
 par_S_branchMod = param(27) !model for branch litter model
 p0_ref = param(29)
 ETS_ref = param(30)
 par_thetaMax = param(31)
 par_H0max = param(32)
 par_gamma = param(33)
 par_kH = param(34)
 par_rhof1 = 0.!param(20)
 par_Cr2 = 0.!param(24)
 par_sla0 = param(39)
 par_tsla = param(40)
 par_zb = param(41)

! do siteNo = 1, nSites  !start site loop

if (year > maxYearSite) then
  STAND(2) = 0. !!newX
  STAND(8:21) = 0. !#!#
  STAND(23:37) = 0. !#!#
  STAND(42:44) = 0. !#!#
  STAND(47:nVar) = 0. !#!#
else

! initialize site variables
!  sitetype = STAND(3)
  ! gammaC = STAND(2)
  age = STAND(7)
  p0 = STAND(6)/1000.	! convert    g C m-2    to    kg C m-2  !#!#
  ETS = STAND(5)
  H = STAND(11)
  D = STAND(12)
  BA = STAND(13)! * par_ops2
  Hc = STAND(14)
  Cw = STAND(15)
  ! Lc = STAND(16)
  N = STAND(17)
  Lc = H - Hc
  leff = STAND(19)
  keff = STAND(20)
  lproj = STAND(21)
  p_eff_all = STAND(10)*P0yX(year,2)/P0yX(year,1) !!##!!2    smoothing PHOTOSYNTHESIS
  weight = STAND(23)
  A = STAND(16)
  dHc=stand(52)
  S_fol = STAND(26)
  S_fr = STAND(27)
  S_branch = STAND(28)
  par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)

  rc = Lc / (H-1.3) !crown ratio
  B = BA / N
  ! A = rc * B

  wf_STKG = STAND(33)
  wf_treeKG = STAND(34)
  ! B = STAND(35)
  Light = STAND(36)
  hb = par_betab * Lc ** par_x
  Cw = 2. * hb

if (N>0.) then

!!!!###here starts stand2 subroutine!!!!!!!!!!!#########
  if (sitetype <= 1.) then
   par_alfar = par_alfar1
  else if (sitetype==2.) then
   par_alfar = par_alfar2
  else if (sitetype==3.) then
   par_alfar = par_alfar3
  else if (sitetype==4.) then
   par_alfar = par_alfar4
  else
   par_alfar = par_alfar5
  end if

!!relate metabolic and structural parameters to site conditions

!  par_mf = par_mf0 * p0 / p0_ref
!  par_mr = par_mr0 * p0 / p0_ref
!  par_mw = par_mw0 * p0 / p0_ref

  par_H0 = par_H0max * (1 - exp(-par_kH * ETS/par_alfar)) !!!new version
  theta = par_thetaMax / (1. + exp(-(H - par_H0)/(par_H0*par_gamma)))   !!!!new version

  mrFact = max(0., par_aETS * (ETS_ref-ETS)/ETS_ref) !!!new version
  par_mr = par_mr0* p0 / p0_ref + (1+par_c) * mrFact / par_vr0    !!!new version !!newX
  par_mf = par_mf0* p0 / p0_ref !!newX
  par_mw = par_mw0* p0 / p0_ref !!newX

  par_rhof0 = par_rhof1 * ETS_ref + par_rhof2
  par_rhof = par_rhof1 * ETS + par_rhof2
  par_vf = par_vf0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref)
  par_vr = par_vr0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref) !!new version
  par_rhor = par_alfar * par_rhof

    ! -------------------------------------
    !GPP all STAND$species   UNITS: g C  /  m2
    ! -------------------------------------
        p_eff = weight * p_eff_all
		! gpp_sp = weight * STAND(10)

    if(wf_STKG > 0.) then
!if(H < 10.) then
 !       s0 = min(0.35 * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.)
!else
        s0 = min(par_s0scale * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.)
! endif
		! if(ij==1 .and. output(year,1,ij,1) == 3.) write(1,*) year, s0, &
	! par_s0scale * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.
  ! if(ij==2) write(2,*) s0,par_s0scale * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.
  ! if(ij==3) write(3,*) s0,par_s0scale * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.
    else
        s0 = 0.
    endif
	gpp_sp = max(0.,(s0 - par_s1 * H) * wf_STKG / 10000.)
	
        !---------------------------------------
        ! DYNAMIC GROWTH MODEL STARTS
        !Updating the tree H, D, Hc and Cw for the next year, according to the method by Valentine & Makela (2005)
        !Valentine & Makela 2005. Bridging process - based and empirical approaches to modeling tree growth.
        ! HERE the units are kg / ha
			gammaC = par_cR/Light
            betab = hb/Lc
            beta0 = par_beta0
            beta1 = (beta0 + betab + par_betas) !!newX
            beta2 = 1. - betab - par_betas 		!!newX
			betaC = (beta1 + gammaC * beta2) / par_betas
			W_wsap = stand(47)
            ! wf_treeKG = par_rhof * A  !!newX
			! wf_STKG = N * wf_treeKG  !!newX
			W_froot = stand(25)
			W_c = stand(48) !sapwood stem below Crown
			W_s = stand(49) !sapwood stem within crown
			W_branch = stand(24) !branches biomass
			W_croot = stand(32) !W_stem * (beta0 - 1.)	!coarse root biomass
			Wsh = stand(50)
			Wdb = stand(51)
			W_bh = stand(53)
			W_crh = stand(54)
			W_bs =  par_rhow * A * N * betab * Lc
			W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)	!coarse root biomass
			Respi_m = (par_mf + par_alfar*par_mr)* wf_STKG + par_mw * W_wsap  !!newX
			npp = (gpp_sp - Respi_m / 10000.) / (1.+par_c)  !!newX
			Respi_tot = gpp_sp - npp
				 ! ! litter fall in the absence of thinning
      S_fol = S_fol + wf_STKG / par_vf	!foliage litterfall
      S_fr  = S_fr + W_froot / par_vr	!fine root litter
	  S_branch = S_branch + Wdb/Tdb
		
	  ! S_branch = S_branch + N * par_rhow * betab * A * (dHc + theta*Lc)
			
        !Height growth-----------------------
			! if(ij==1 .and. stand(1)==13429) write(1,*) dH,H,Hc,npp,wf_STKG,par_vf,W_froot, &
				! par_vr,theta,W_wsap, par_z, W_wsap,gammaC, W_c,W_bs, betaC, W_s
		f1 = npp*10000 - (wf_STKG/par_vf) - (W_froot/par_vr) - (theta * W_wsap)
		f2 = (par_z* (wf_STKG + W_froot + W_wsap)* (1-gammaC) + par_z * gammaC * (W_c + &
				par_zb *W_bs + beta0 * W_c) + betaC * W_s)
		dH = max(0.,((H-Hc) * f1/f2))
		Gf = par_z * wf_STKG/(H-Hc) * (1-gammac)*dH
		Gr = par_z * W_froot/(H-Hc) * (1-gammac)*dH
		
		if(f1 < 0.) then
			dH = 0.
			mort = 888.
		elseif(f2<=0. .or. Gf<0. .or. Gr < 0.) then
			gammaC = 1.
			f2 = (par_z* (wf_STKG + W_froot + W_wsap)* (1-gammaC) + par_z * gammaC * (W_c + &
				par_zb *W_bs + beta0 * W_c) + betaC * W_s)
			dH = max(0.,((H-Hc) * f1/f2))
			mort = 888.
		endif

 ! if(stand_all(1,1)==4746. .and. ij==1) then
	 ! write(1,*) f1,f2,dH,gammaC, ij,&
	 ! H,Hc,npp,wf_STKG,par_vf,W_froot,par_vr,theta,W_wsap, & 
	 ! par_z,W_c,W_bs,betaC,W_s
 ! endif
 ! if(stand_all(1,1)==4746. .and. ij==2) then
	 ! write(2,*) f1,f2,dH,gammaC, ij,&
	 ! H,Hc,npp,wf_STKG,par_vf,W_froot,par_vr,theta,W_wsap, & 
	 ! par_z,W_c,W_bs,betaC,W_s
 ! endif
! if(stand_all(1,1)==4746. .and. ij==3) then
	 ! write(3,*) f1,f2,dH,gammaC, ij,&
	 ! H,Hc,npp,wf_STKG,par_vf,W_froot,par_vr,theta,W_wsap, & 
	 ! par_z,W_c,W_bs,betaC,W_s
 ! endif
  
 !-----------------------------------
        !crown rise
!         if(H - Hc > par_Cr2*100./sqrt(N)) then
!        if(2.*hb > 100./sqrt(N) ) then
        dHc = min(gammaC * dH,(H - Hc))
			! if(ij==1 .and. stand(1)==13429) write(2,*) dH,H,Hc,npp,wf_STKG,par_vf,W_froot, &
				! par_vr,theta,W_wsap, par_z, W_wsap,gammaC, W_c,W_bs, betaC, W_s
		
if(time==1)then
      dHcCum = 0.
      dHCum = 0.
endif
        dHcCum = dHcCum + dHc
        dHCum = dHCum + dH
!	    else
!		dHc = 0  !CAN BE DIFFERENT FROM THE PAPER HARKONEN ET AL. 2013 CANADIAN JOURNAL, SEE THE EQUATION THERE
!        endif
        if(dHc < 0. )dHc = 0.

            !----------------------------------
            !New values for H, Hc and Lc

       ! diameter growth
            if(Lc > 0.) then
                dA = max(0.,par_z*A*(dH-dHc)/Lc)
!                dB = par_z * (A / Lc) * dH
                dB = max(0.,par_z * A / Lc * dH + theta * A) !!!! v1
            else
                dA = 0.
                dB = 0.
            endif
            
! Calculate change rates for non-living wood   - DECLARE dWsh etc...          
            dWsh = max(0.,par_rhow * par_z * A/Lc * dHc * Hc * N	+ theta*(W_c + W_s))
            dW_bh = max(0.,W_bs*theta - W_bh * gammaC * dH / Lc) 
            dW_crh = max(0.,W_crs*theta + par_z * W_c * beta0 / Lc * gammaC * dH)
            dWdb = max(0.,W_branch/Lc * par_zb * gammaC * dH - Wdb/Tdb)

!!  Update state variables
          H = H + step * dH
          A = A + step * dA
          B = B + step * dB
		  Hc = Hc + step * dHc
          
          Wsh = Wsh + dWsh * step
          W_bh = W_bh + dW_bh * step
		  W_crh = W_crh + dW_crh * step
		  Wdb = Wdb + dWdb * step
          
! Update dependent variables
      wf_treeKG = par_rhof * A
      wf_STKG = N * wf_treeKG
      BA = N * B
      D = sqrt(B*4./pi)*100. ! * 100 converts meters in cm
      Lc = H - Hc
      rc = Lc / (H-1.3)
      if(rc > 0.) B2 = A / rc
      hb = par_betab * Lc**par_x
      ! hb = betab * Lc!!!to check one cause of not matching npp because betab was computed before growth

! more dependent variables (not used in calculation)
        betab = hb/Lc
        beta0 = par_beta0
        beta1 = (beta0 + betab + par_betas) !!newX
        beta2 = 1. - betab - par_betas 		!!newX
		betaC = (beta1 + gammaC * beta2) / par_betas
		W_froot = par_rhor * A * N  !!to check  !!newX
		W_wsap = par_rhow * A * N * (beta1 * H + beta2 * Hc) !!newX
		W_c = par_rhow * A * N * Hc !sapwood stem below Crown
		W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
		W_bs =  par_rhow * A * N * betab * Lc !branches biomass
        W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)	!coarse root sapwood biomass
		! Wsh = Wsh + par_z * gammaC * W_c/Lc * dH + theta*(W_c + W_s)
		W_stem = W_c + W_s + Wsh
		V = W_stem / par_rhow
		W_branch = W_bs + W_bh
		W_croot = W_crs + W_crh
		
  age = age + step
  ! if (dH==0.) mort=888.
  ! if(ij==1) write(11,*)gammaC , dH, dHc,H,Hc
  ! if(ij==2) write(12,*)gammaC , dH, dHc,H,Hc
  ! if(ij==3) write(13,*)gammaC , dH, dHc,H,Hc
! write(3,*) gammaC, dH,dHc
  STAND(2) = gammaC
  STAND(40) = mort
  STAND(41) = dH
  STAND(7) = age !#!#
  STAND(18) = npp
  !STAND(8) = Respi_m /10000.
  STAND(9) = Respi_tot
  STAND(11) = H
  STAND(12) = D
  STAND(13) = BA
  STAND(14) = Hc
  STAND(15) = Cw
  STAND(16) = A
  STAND(17) = N
  STAND(24) = W_branch
  STAND(25) = W_froot
  STAND(26) = S_fol
  STAND(27) = S_fr
  STAND(28) = S_branch
  ! STAND(29) = S_wood
  STAND(30) = V
  STAND(31) = W_stem
  STAND(32) = W_croot
  STAND(33) = wf_STKG
  STAND(34) = wf_treeKG
  STAND(35) = B
  STAND(36) = Light
  ! STAND(42) = Vold* min(1.,-dN*step/Nold)
  STAND(44) = gpp_sp
  STAND(47) = W_wsap
  STAND(48) = W_c
  STAND(49) = W_s
  STAND(50) = Wsh
  STAND(53) = W_bh
  STAND(54) = W_crh
  STAND(51) = Wdb
  STAND(52) = dHc
  ! stand(22) = theta
else
  STAND(2) = 0. !!newX
  STAND(8:21) = 0. !#!#
  STAND(23:37) = 0. !#!#
  STAND(42:44) = 0. !#!#
  STAND(47:nVar) = 0. !#!#
  STAND(7) = STAND(7) + step
endif
endif
 
  !Perform user defined thinning or defoliation events for this time period
  If (countThinning <= nThinning .and. time==inttimes) Then
   If (year == int(thinning(countThinning,1)) .and. ij == int(thinning(countThinning,3))) Then! .and. siteNo == thinning(countThinning,2)) Then
	STAND_tot = STAND
	if(thinning(countThinning,9) .NE. -999) then
	 thinning(countThinning,6) = thinning(countThinning,9) * (pi*((D/2./100.)**2.))
	endif
    if(thinning(countThinning,4)==0.) then
     STAND(2) = 0. !!newX
	 STAND(8:21) = 0. !#!#
     STAND(23:37) = 0. !#!#
     STAND(43:44) = 0. !#!#
	 STAND(47:nVar) = 0. !#!#
 !! calculate litter including residuals from thinned trees
  !energyCut
	 S_fol = wf_STKG + S_fol
     S_fr = W_froot + S_fr
	 if(energyCut==1.) then
	  energyWood(year,ij,2) = (W_branch + W_croot*0.3 + W_stem* (1-harvRatio)) * energyRatio
	  energyWood(year,ij,1) = energyWood(year,ij,2) / par_rhow
      S_branch = (W_branch + Wdb) * (1-energyRatio) + W_croot*0.83 + S_branch + &
				W_stem* (1-harvRatio)* (1-energyRatio) 
      S_wood = S_wood + W_croot*0.17 * (1-energyRatio)!(1-harvRatio) takes into account of the stem residuals after thinnings
	 else
      S_branch = W_branch + Wdb + W_croot*0.83 + S_branch + W_stem* (1-harvRatio) 
      S_wood = S_wood + W_croot*0.17!(1-harvRatio) takes into account of the stem residuals after thinnings
	 endif
  !energyCut
     STAND(26) = S_fol
     STAND(27) = S_fr
     STAND(28) = S_branch
     STAND(29) = S_wood
    else
     if(thinning(countThinning,8)==1.) then
	if(thinning(countThinning,4) < 2. .and. thinning(countThinning,4) > 0.) then
	 thinning(countThinning,4) = H * thinning(countThinning,4)
	endif
	if(thinning(countThinning,5) < 2. .and. thinning(countThinning,5) > 0.) then
	 thinning(countThinning,5) = D * thinning(countThinning,5)
	endif
	if(thinning(countThinning,6) < 1. .and. thinning(countThinning,6) > 0.) then
	 thinning(countThinning,6) = BA * thinning(countThinning,6)
	endif
	if(thinning(countThinning,7) < 2. .and. thinning(countThinning,7) > 0.) then
	 thinning(countThinning,7) = Hc * thinning(countThinning,7)
	endif
     endif
     if (thinning(countThinning,4) /= -999.) H = thinning(countThinning,4)
     if (thinning(countThinning,7) /= -999.) stand(14) = thinning(countThinning,7)
     if (thinning(countThinning,5) /= -999.) D = thinning(countThinning,5)
     BA = thinning(countThinning,6)
	 Hc=stand(14)
     Lc = H - Hc !Lc
     rc = Lc / (H-1.3) !crown ratio
     Nold = N
     wf_STKG_old = wf_STKG
     W_stem_old = W_stem
     N = BA/(pi*((D/2./100.)**2.)) ! N
     Nthd = max(0.,(Nold-N)) ! number of cutted trees
     B = BA/N!(pi*((D/2/100)**2))
     A = stand(16) * B/stand(35)

     ! Update dependent variables
	 hb = par_betab * Lc**par_x
	 gammaC = par_cR/stand(36)
	 par_rhof = par_rhof1 * ETS + par_rhof2
 	 par_rhor = par_alfar * par_rhof
     betab = hb/Lc
     beta0 = par_beta0
     beta1 = (beta0 + betab + par_betas) !!newX
     beta2 = 1. - betab - par_betas 		!!newX
	 betaC = (beta1 + gammaC * beta2) / par_betas


	    W_bh = stand(53)
		W_crh = stand(54)
		wf_treeKG = par_rhof * A
		wf_STKG = N * wf_treeKG
		W_froot = par_rhor * A * N  !!to check  !!newX
		W_wsap = par_rhow * A * N * (beta1 * H + beta2 * Hc) !!newX
		W_c = par_rhow * A * N * Hc !sapwood stem below Crown
		W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
		W_bs =  par_rhow * A * N * betab * Lc !branches biomass
        W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)	!coarse root biomass
		! Wsh = Wsh + par_z * gammaC * W_c/Lc * dH + theta*(W_c + W_s)
		Wsh = stand(50) * N/Nold !!this is not proportional to the size of the standing trees
		W_stem = W_c + W_s + Wsh
		W_crh = W_crh * N/Nold
		W_bh = W_bh * N/Nold
		V = W_stem / par_rhow
		W_croot = W_crs + W_crh
		W_branch = W_bs + W_bh
		Wdb = Wdb * N/Nold
!! calculate litter including residuals from thinned trees
  !energyCut
	S_fol = stand(26) + stand(33) - wf_STKG
    S_fr = stand(27) + stand(25) - W_froot
	if(energyCut==1.) then
	 energyWood(year,ij,2) = (stand(24) - W_branch + (stand(32) - W_croot)*0.3 + &
					(stand(31) - W_stem) * (1-harvRatio)) * energyRatio
	 energyWood(year,ij,1) = energyWood(year,ij,2) / par_rhow
     S_branch = stand(28) + (stand(24) - W_branch + stand(51) - Wdb) * (1-energyRatio) + &
				(stand(32) - W_croot)*0.83 + &
				(stand(31) - W_stem) * (1-harvRatio) * (1-energyRatio)
     S_wood = stand(29) +(stand(32) - W_croot)*0.17* (1-energyRatio)
	else
    S_branch = stand(28)+stand(24)-W_branch+stand(51)-Wdb+(stand(32)-W_croot)*0.83 + &
		(stand(31) - W_stem) * (1-harvRatio)
    S_wood = stand(29)  + (stand(32) - W_croot)*0.17
	endif
  !energyCut	
! !! calculate litter including residuals from thinned trees
    ! S_fol = stand_all(26,ij) + stand_all(33,ij) - wf_STKG
    ! S_fr = stand_all(27,ij) + stand_all(25,ij) - W_froot
    ! S_branch = stand_all(28,ij) + stand_all(24,ij) - W_branch
    ! S_wood = stand_all(29,ij) + (stand_all(31,ij) - W_stem) * (1-harvRatio) + stand_all(32,ij) - W_croot


     outt(11,ij,2) = STAND_tot(11)
     outt(12,ij,2) = STAND_tot(12)
     outt(13,ij,2) = STAND_tot(13) - BA
     outt(14,ij,2) = STAND_tot(14)
     outt(15,ij,2) = STAND_tot(15)
     outt(16,ij,2) = STAND_tot(16)
     outt(17,ij,2) = Nthd
     outt(18:23,ij,2) = -999.
     outt(24,ij,2) = STAND_tot(24) - W_branch
     outt(25,ij,2) = STAND_tot(25) - W_froot
     outt(26:29,ij,2) = -999.
     outt(30,ij,2) = STAND_tot(30) - V
     outt(31,ij,2) = STAND_tot(31) - W_stem
     outt(32,ij,2) = STAND_tot(32) - W_croot
     outt(33,ij,2) = STAND_tot(33) - wf_STKG
     outt(34,ij,2) = (STAND_tot(34)*Nold - wf_treeKG*N)/Nthd
     outt(35,ij,2) = -999.; outt(36,ij,2)= -999.

    stand(11) = H
    stand(12) = D
    stand(13) = BA
    stand(16) = A
    stand(17) = N
    stand(26) = S_fol
    stand(27) = S_fr
    stand(28) = S_branch
    stand(29) = S_wood
    stand(24) = W_branch
    stand(25) = W_froot
    stand(30) = V  !
    stand(31) = W_stem
    stand(32) = W_croot
    stand(33) = wf_STKG
    stand(34) = wf_treeKG
    stand(35) = B
	stand(47) = W_wsap
	stand(48) = W_c
	stand(49) = W_s
	stand(53) = W_bh
	stand(54) = W_crh
	stand(50) = Wsh
	stand(51) = Wdb
    endif

	countThinning = countThinning + 1

   End If
  End If

	STAND_all(:,ij)=STAND
end do !!!!end loop species
 end do !!!!end loop inttimes

!Perform thinning or defoliation events for this time period using standard management routines!!!!!!!!!!!!!!!!
!do siteNo = 1, nSites
 ! write(2,*) "before clcut"
 ! if(stand_all(1,1)==6944.) then
	! write(2,*) stand_all(30,1), stand_all(30,2), stand_all(30,3)
 ! endif
!!!!test for clearcut!!!!
 domSp = maxloc(STAND_all(13,:))
 layer = int(domSp(1))
if (ClCut == 1.) then
	species = int(stand_all(4,layer))
	D_clearcut = inDclct(species)
	A_clearcut = inAclct(species)
	D = stand_all(12,layer)
	age = stand_all(7,layer)
 if ((D > D_clearcut) .or. (age > A_clearcut)) then
  do ij = 1, nLayers
  ! if(stand_all(1,1)==6944. .and. ij==1) then
	! write(1,*) stand_all(30,1), stand_all(30,2), stand_all(30,3)
  ! endif
   outt(6:nVar,ij,2) = stand_all(6:nVar,ij)

  !energyCut
   S_fol = stand_all(33,ij) + stand_all(26,ij)
   S_fr = stand_all(25,ij) + stand_all(27,ij)
	if(energyCut==1.) then
	 energyWood(year,ij,2) = energyWood(year,ij,2) + (stand_all(24,ij) + &
					stand_all(32,ij)*0.3 + stand_all(31,ij) * (1-harvRatio)) * energyRatio
	 energyWood(year,ij,1) = energyWood(year,ij,2) / par_rhow
	 S_branch = stand_all(28,ij) + (stand_all(24,ij) + stand_all(51,ij)) * &
		(1-energyRatio) + stand_all(32,ij) * 0.83 +&
		stand_all(31,ij)* (1-harvRatio) * (1-energyRatio)
	 S_wood = stand_all(32,ij) *0.17 * (1-energyRatio) + stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts
	else
	 S_branch = stand_all(51,ij)+stand_all(24,ij)+stand_all(28,ij)+stand_all(32,ij)* 0.83 +&
			stand_all(31,ij)* (1-harvRatio)
	 S_wood = stand_all(32,ij) *0.17 + stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts
	endif
  !energyCut
   stand_all(2,ij) = 0. !!newX
   stand_all(8:21,ij) = 0.
   stand_all(23:38,ij) = 0.
   stand_all(41,ij) = 0.
   stand_all(43:44,ij) = 0.
   stand_all(47:nVar,ij) = 0.
   stand_all(26,ij) = S_fol
   stand_all(27,ij) = S_fr
   stand_all(28,ij) = S_branch
   stand_all(29,ij) = S_wood
!!update age
  ! do ki = 1, min(20,(nyears-year))
   ! modOut((year+ki),7,ij,1) = ki !#!#
   ! modOut((year+ki),4,ij,1) = initVar(1,ij) !#!#
  ! enddo
  enddo
 endif
endif

!!!!test for thinnings!!!!
 !!!!!!!for coniferous dominated stands!!!!!!
if(defaultThin == 1.) then
 sitetype = siteInfo(3)
 BA_tot = sum(stand_all(13,:))!+stand_all(13,2)+stand_all(13,3)
 BAr = stand_all(13,:)/BA_tot
 domSp = maxloc(STAND_all(13,:))
 layer = int(domSp(1))
 H = stand_all(11,layer)
 species = int(stand_all(4,layer))

 ! counting the dominant height of the dominant species
 Hdom = pCrobas(42,species)+pCrobas(43,species)*H 
 Ntot = sum(STAND_all(17,:))
	!! here we decide what thinning function to use; 3 = tapioThin, 2 = tapioFirstThin, 1 = tapioTend
 call chooseThin(species, siteType, ETSmean, Ntot, Hdom, tTapio, ftTapio, thinningType) 
 if(thinningType == 3.) then    
	call tapioThin(pCrobas(28,species),siteType,ETSmean,Hdom,tapioPars,BAtapio,thdPer,limPer)  
	BA_lim = BAtapio(1) ! BA limit to start thinning
	BA_thd = BAtapio(2) ! BA after thinning
	if(BA_tot > BA_lim) then 
		doThin = .true.
	else
		doThin = .false.
	endif
 else if(thinningType == 2.) then
	call tapioFirstThin(pCrobas(28,species),siteType,ETSmean,ftTapio,limPer,thdPer,early,tapioOut)
	Hdom_lim = tapioOut(1) ! Hdom limit to start thinning
	dens_lim = tapioOut(2) ! density limit to start thinning; both need to be reached
	dens_thd = tapioOut(3) ! density after thinning
	if(Hdom > Hdom_lim .and. Ntot > dens_lim) then 
		doThin = .true.
	else
		doThin = .false.
	endif
 else if(thinningType == 1.) then
	call tapioTend(pCrobas(28,species),siteType,ETSmean,tTapio,limPer,thdPer,tapioOut)
	Hdom_lim = tapioOut(1)! Hdom limit to start thinning
	dens_lim = tapioOut(2) ! density limit to start thinning; both need to be reached
	dens_thd = tapioOut(3) ! density after thinning
	if(Hdom > Hdom_lim .and. Ntot > dens_lim) then 
		doThin = .true.
	else
		doThin = .false.
	endif
 endif
 


 if(doThin) then

  do ij = 1, nLayers

   if(stand_all(17,ij)>0.) then
    STAND_tot = stand_all(:,ij)
	species = int(stand_all(4,ij))
    param = pCrobas(:,species)
    par_cR=param(1)
    par_rhow=param(2)
    par_sla =param(3)
    par_k =param(4)
    par_vf0 =param(5)
    par_vr0 =param(6)
    par_c=param(7)
    par_mf0=param(8)
    par_mr0=param(9)
    par_mw0=param(10)
    par_z=param(11)
    par_beta0=param(12)
    par_betab=param(13)
    par_betas = param(14)
    par_rhof2 = param(15)
    par_s1 = param(16)
    par_kRein = param(17)
    par_s0scale = param(18)
    par_x = param(19)
    par_aETS = param(20)
    par_alfar1 =param(21)
    par_alfar2 =param(22)
    par_alfar3 =param(23)
    par_alfar4 =param(24)
    par_alfar5 =param(25)
    par_sarShp = param(26) !Shape surface area of the crown: 1.= cone; 2.=ellipsoide
    par_S_branchMod = param(27) !model for branch litter model
    p0_ref = param(29)
    ETS_ref = param(30)
    par_thetaMax = param(31)
    par_H0max = param(32)
    par_gamma = param(33)
	par_kH = param(34)
    par_rhof1 = 0.!param(20)
    par_Cr2 = 0.!param(24)
    par_rhof = par_rhof1 * stand_all(5,ij) + par_rhof2
	Nold = stand_all(17,ij)
	
	if(thinningType == 1 .or. thinningType == 2) then
		! N = number of trees in the current layer after thinning
		N = (stand_all(17,ij)/Ntot) * dens_thd
		H = stand_all(11,ij)
		D = stand_all(12,ij)
		BA = N*pi*(D/2./100.)**2.
	else if(thinningType == 3) then 
		BA_tot = BA_thd
		BA = BAr(ij) * BA_thd

		if(par_sarShp==1.) then
			H = stand_all(11,ij) *  (1.2147-0.2086 * (BA/ stand_all(13,ij)))
			D = stand_all(12,ij) * (1.2192 -0.2173 * (BA/ stand_all(13,ij)))
		else
			H = stand_all(11,ij) *  (1.07386 -0.06553 * (BA/ stand_all(13,ij)))
			D = stand_all(12,ij) * (1.1779 -0.1379 * (BA/ stand_all(13,ij)))
		endif
		N = BA/(pi*((D/2./100.)**2.))
	endif

	stand_all(13,ij) = BA	
    Nthd = max(0.,(Nold - N))
    Hc = stand_all(14,ij)
	Wdb = stand_all(51,ij)
    Lc = H - Hc !Lc
    rc = Lc / (H-1.3) !crown ratio
    wf_STKG_old = stand_all(33,ij)
    W_stem_old = stand_all(31,ij)
    B = BA/N
    A = stand_all(16,ij) * B/stand_all(35,ij) !!! to check 
    hb = par_betab * Lc ** par_x
    Cw = 2. * hb


!!update biomasses
    if (sitetype <= 1.) then
     par_alfar = par_alfar1
    else if (sitetype==2.) then
     par_alfar = par_alfar2
    else if (sitetype==3.) then
     par_alfar = par_alfar3
    else if (sitetype==4.) then
     par_alfar = par_alfar4
    else
     par_alfar = par_alfar5
    end if

     ! Update dependent variables
	 gammaC = par_cR/stand_all(36,ij)
	 par_rhof = par_rhof1 * ETS + par_rhof2
 	 par_rhor = par_alfar * par_rhof
     betab = hb/Lc
     beta0 = par_beta0
     beta1 = (beta0 + betab + par_betas) !!newX
     beta2 = 1. - betab - par_betas 		!!newX
	 betaC = (beta1 + gammaC * beta2) / par_betas

		W_bh = stand_all(53,ij)
		W_crh = stand_all(54,ij)
		wf_treeKG = par_rhof * A
		wf_STKG = N * wf_treeKG
		W_froot = par_rhor * A * N  !!to check  !!newX
		W_wsap = par_rhow * A * N * (beta1 * H + beta2 * Hc) !!newX
		W_c = par_rhow * A * N * Hc !sapwood stem below Crown
		W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
		W_bs =  par_rhow * A * N * betab * Lc !branches biomass
        W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)	!coarse root biomass
		! Wsh = Wsh + par_z * gammaC * W_c/Lc * dH + theta*(W_c + W_s)
		Wsh = stand_all(50,ij) * N/Nold !!this is not proportional to the size of the standing trees
		W_stem = W_c + W_s + Wsh
		W_crh = W_crh * N/Nold
		W_bh = W_bh * N/Nold
		V = W_stem / par_rhow
		W_croot = W_crs + W_crh
		W_branch = W_bs + W_bh
		Wdb = stand_all(51,ij) * N/Nold

!! calculate litter including residuals from thinned trees
  !energyCut
    S_fol = stand_all(26,ij) + stand_all(33,ij) - wf_STKG
    S_fr = stand_all(27,ij) + stand_all(25,ij) - W_froot
    if(energyCut==1.) then
	 energyWood(year,ij,2) = energyWood(year,ij,2) + (stand_all(24,ij) - W_branch + &
		(stand_all(32,ij) - W_croot) * 0.3 + &
	    (stand_all(31,ij) - W_stem) * (1-harvRatio)) * energyRatio
	 energyWood(year,ij,1) = energyWood(year,ij,2) / par_rhow

     S_branch = stand_all(28,ij) + (stand_all(24,ij) - W_branch + stand_all(51,ij) - Wdb) * &
		(1-energyRatio) + (stand_all(32,ij) - W_croot) * 0.83 +&
		(stand_all(31,ij)-W_stem)*(1-harvRatio)*(1-energyRatio)
	 S_wood = stand_all(29,ij)+(stand_all(32,ij) - W_croot) * 0.17 * (1-energyRatio)
	else
     S_branch = stand_all(28,ij) + stand_all(24,ij) - W_branch + stand_all(51,ij) - Wdb + &
		(stand_all(32,ij) - W_croot) * 0.83+ (stand_all(31,ij) - W_stem) * (1-harvRatio)
     S_wood = stand_all(29,ij)  + (stand_all(32,ij) - W_croot) * 0.17
	endif
  !energyCut
	
    outt(11,ij,2) = STAND_tot(11)
    outt(12,ij,2) = STAND_tot(12)
    outt(13,ij,2) = STAND_tot(13) - BA
    outt(14,ij,2) = STAND_tot(14)
    outt(15,ij,2) = STAND_tot(15)
    outt(16,ij,2) = STAND_tot(16)
    outt(17,ij,2) = Nthd
    outt(18:23,ij,2) = -999.
    outt(24,ij,2) = STAND_tot(24) - W_branch
    outt(25,ij,2) = STAND_tot(25) - W_froot
    outt(26:29,ij,2) = -999.
    outt(30,ij,2) = STAND_tot(30) - V
    outt(31,ij,2) = STAND_tot(31) - W_stem
    outt(32,ij,2) = Nthd * W_croot/N
    outt(33,ij,2) = STAND_tot(33) - wf_STKG
    outt(34,ij,2) = outt(33,ij,2)/Nthd
    outt(35,ij,2) = -999.; outt(36,ij,2)= -999.

    stand_all(11,ij) = H
    stand_all(12,ij) = D
    stand_all(13,ij) = BA
    stand_all(16,ij) = A
    stand_all(17,ij) = N
    stand_all(26,ij) = S_fol
    stand_all(27,ij) = S_fr
    stand_all(28,ij) = S_branch
    stand_all(29,ij) = S_wood
    stand_all(24,ij) = W_branch
    stand_all(25,ij) = W_froot
    stand_all(30,ij) = V  !
    stand_all(31,ij) = W_stem
    stand_all(32,ij) = W_croot
    stand_all(33,ij) = wf_STKG
    stand_all(34,ij) = wf_treeKG
    stand_all(35,ij) = B
	stand_all(47,ij) = W_wsap
	stand_all(48,ij) = W_c
	stand_all(49,ij) = W_s
	stand_all(53,ij) = W_bh
	stand_all(54,ij) = W_crh
	stand_all(50,ij) = Wsh
	stand_all(51,ij) = Wdb

   endif
  enddo
 endif !default thin
endif
 ! write(2,*) "after thinnings"
outt(:,:,1) = STAND_all

modOut((year+1),2,:,:) = outt(2,:,:)
modOut((year+1),7,:,:) = outt(7,:,:)
modOut((year+1),9:nVar,:,:) = outt(9:nVar,:,:)

!!!!run Yasso
 if(yassoRun==1.) then
  do ijj = 1, nLayers
   Lst(ijj) = outt(29,ijj,1)
   Lb(ijj) =  outt(28,ijj,1)
   Lf(ijj) = outt(26,ijj,1)+outt(27,ijj,1)

   species = int(initVar(1,ijj))
   call compAWENH(Lf(ijj),folAWENH(ijj,:),pAWEN(1:4,species))   !!!awen partitioning foliage
   if(GVrun==1 .and. ijj==1) then 
    folAWENH(ijj,1:4) = folAWENH(ijj,1:4) + AWENgv			 !!!add AWEN gv to 1st layer
   endif
   call compAWENH(Lb(ijj),fbAWENH(ijj,:),pAWEN(5:8,species))   !!!awen partitioning branches
   call compAWENH(Lst(ijj),stAWENH(ijj,:),pAWEN(9:12,species))         !!!awen partitioning stems

   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,1,ijj),stAWENH(ijj,:),litterSize(1,species), &
	leac,soilC((year+1),:,1,ijj),0.)
   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,2,ijj),fbAWENH(ijj,:),litterSize(2,species), &
	leac,soilC((year+1),:,2,ijj),0.)
   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,3,ijj),folAWENH(ijj,:),litterSize(3,species), &
	leac,soilC((year+1),:,3,ijj),0.)
  enddo
 ! write(2,*) "after yasso"

  soilCtot(year+1) = sum(soilC(year+1,:,:,:))
  ! write(*,*) soilCtot(year+1)
 endif !end yassoRun if
enddo !end year loop

! write(2,*) "after loop years"

!soil and harvested volume outputs
modOut(:,37,:,1) = modOut(:,30,:,2) * harvRatio!! harvRatio takes into account the residuals left in the soil 
modOut(:,38,:,1) = modOut(:,31,:,2) * harvRatio!! harvRatio takes into account the residuals left in the soil 

do year = 1,(nYears+1)
  do ijj = 1, nLayers
	! modOut(year,38,ijj,1) = sum(modOut(1:year,30,ijj,2)) + &
		! sum(modOut(1:year,42,ijj,1)) + modOut(year,30,ijj,1)
	modOut(year,39,ijj,1) = sum(soilC(year,:,:,ijj))
	! modOut(year,38,ijj,1) = pCrobas(2,int(modOut(year,4,ijj,1))) * modOut(year,37,ijj,1)
	if(year > 1.5) then
	!compute gross growth
	  modOut(year,43,ijj,1) = modOut(year,30,ijj,1) - modOut((year-1),30,ijj,1) + &
				modOut(year,42,ijj,1) + modOut(year,37,ijj,1)/harvRatio
	endif
	! write(*,*) modOut(year,39,ijj,1)
  enddo
enddo



 ! write(2,*) "here2"

!compute fluxes in g C m−2 y−1
 modOut(:,44,:,1) = modOut(:,44,:,1)*1000. !*1000 coverts units to g C m−2 y−1
 modOut(:,9,:,1) = modOut(:,9,:,1)*1000.    !*1000 coverts units to g C m−2 y−1
 modOut(:,18,:,1) = modOut(:,18,:,1)*1000.    !*1000 coverts units to g C m−2 y−1

	modOut(2:(nYears+1),45,:,1) = modOut(1:(nYears),39,:,1)/10. - modOut(2:(nYears+1),39,:,1)/10. + &	!/10 coverts units to g C m−2 y−1
	modOut(2:(nYears+1),26,:,1)/10. + modOut(2:(nYears+1),27,:,1)/10. + &
	modOut(2:(nYears+1),28,:,1)/10. + modOut(2:(nYears+1),29,:,1)/10.
    if(GVrun==1) modOut(2:(nYears+1),45,1,1) = modOut(2:(nYears+1),45,1,1) + GVout(:,2)/10.  !/10 coverts units to g C m−2 y−1
	
modOut(:,46,:,1) = modOut(:,44,:,1) - modOut(:,9,:,1) - modOut(:,45,:,1) !!Gpp is not smoothed
!modOut(:,46,:,1) = modOut(:,18,:,1) - modOut(:,45,:,1) !!!everything smoothed

!!!!ground vegetation Add Npp ground vegetation to the NEE first layer
if(GVrun==1) modOut(2:(nYears+1),46,1,1) = modOut(2:(nYears+1),46,1,1) + GVout(:,3)*0.5 

 output = modOut(2:(nYears+1),:,:,:)
 output(:,5:6,:,:) = modOut(1:(nYears),5:6,:,:)
 soilCinOut = soilC(2:(nYears+1),:,:,:)
 soilCtotInOut = soilCtot(2:(nYears+1))

 ! write(2,*) "end"
 ! close(1)
 ! close(2)
 ! close(3)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
