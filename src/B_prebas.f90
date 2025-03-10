
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine bridging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prebas(nYears,nLayers,nSp,siteInfo,pCrobas,initVar,thinning,output, &
     nThinning,maxYearSite,fAPAR,initClearcut,&
     fixBAinitClarcut,initCLcutRatio,ETSy,P0y,weatherPRELES,DOY,pPRELES,&
     soilCinOut,pYasso,pAWEN,weatherYasso,&
     litterSize,soilCtotInOut,defaultThin,ClCut,energyCut,clct_pars,&
     dailyPRELES,yassoRun,energyWood,tapioPars,thdPer,limPer,&
     ftTapio,tTapio,GVout,thinInt, &
   flagFert,nYearsFert,mortMod,pECMmod, &
   layerPRELES,LUEtrees,LUEgv, siteInfoDist, outDist, prebasFlags, &
   latitude, TsumSBBs)




implicit none

!! Constants
 integer, parameter :: nVar=54, npar=53, inttimes = 1 ! no. of variables, parameters, simulation time-step (always 1)
 real (kind=8), parameter :: pi = 3.1415927, t=1. , ln2 = 0.693147181, fAparFactor=0.9
 real (kind=8), parameter :: energyRatio = 0.7, harvRatio = 0.9 !energyCut
 integer, intent(in) :: nYears, nLayers, nSp ! no of year, layers, species (only to select param.)
 real (kind=8), intent(in) :: weatherPRELES(nYears,365,5) ! R, T, VPD, P, CO2
 integer, intent(in) :: DOY(365)!, ECMmod, etmodel fvec
 real (kind=8), intent(inout) :: pPRELES(30), tapioPars(5,2,3,20), thdPer, limPer ! tapioPars(sitetype, conif/decid, south/center/north, thinning parameters), and parameters for modifying thinnig limits and thresholds
 real (kind=8), intent(inout) :: LUEtrees(nSp),LUEgv,latitude, TsumSBBs(4)!TsumSBB = temp sums bark beetle (1)= previous two years,(2)= previous year, (1)= current year
 real (kind=8), intent(inout) :: tTapio(5,3,2,7), ftTapio(5,3,3,7) ! Tending and first thinning parameter.
 real (kind=8), intent(inout) :: thinning(nThinning, 11) ! User defined thinnings, BA, height of remaining trees, year, etc. Both Tapio rules and user defined can act at the same time. Documented in R interface
 real (kind=8), intent(inout) :: initClearcut(5) !initial stand conditions after clear cut: (H, D, totBA, Hc, Ainit). If not given, defaults are applied. Ainit is the year new stand appears.
 real (kind=8), intent(inout) :: pCrobas(npar, nSp), pAWEN(12, nSp),mortMod,pECMmod(12)
 integer, intent(in) :: maxYearSite ! absolute maximum duration of simulation.
!disturbances

 logical :: disturbance_wind, disturbance_bb, disturbance_fire !fvec
 real (kind=8), intent(inout) :: siteInfoDist(10), outDist(nYears,10) !inputs(siteInfoDist) & outputs(outDist) of disturbance modules
 real (kind=8) :: rndm !random number for disturbance sampling
 integer :: distvloc, sevclasslength !integer for element of NFI sevclass list of directly damaged relative volumes; n of elements in that list
real (kind=8) :: sc1vols(87), sc2vols(15), sc3vols(6)
real (kind=8) :: wriskLayers(nLayers, 6)
!!! wind dist / salvlog development
real (kind=8) :: wdistproc(7) !to replace siteinfodist
REAL (kind=8)::  wrisk5, wrisk0, wrisk ! 5-year wind risk (suvanto output), pre-logit value, annual risk
REAL (kind=8)::  hthresh, htresh_ba !
REAL (kind=8):: wrisk5dd1, wrisk5dd2, wrisk5dd3 !5-year wind risk of each damage density class
REAL (kind=8)::  V_tot, vdam ! vol of all layers, site-level damaged vol
REAL (kind=8):: BAdist(nLayers) !disturbed BA per layer


 real (kind=8), intent(in) :: defaultThin, ClCut, energyCut, yassoRun, fixBAinitClarcut  ! flags. Energy cuts takes harvest residues out from the forest.
 !!oldLayer scenario
 integer, intent(in) :: layerPRELES !oldLayer, fvec
!!! fertilization parameters
 !integer, intent(inout) :: fertThin !!! flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 fvec
 integer, intent(inout) :: flagFert !!! flag that indicates if fertilization has already been applied along the rotation
 integer :: yearsFert !!actual number of years for fertilization (it depends if the thinning occur close to the end of the simulations)
 integer, intent(inout) :: nYearsFert !!number of years for which the fertilization is effective
 real(8) :: alfarFert(nYearsFert,nLayers,2)
!! define arguments, inputs and outputs
 real (kind=8), intent(in) :: clct_pars(nSp,3)! parameters for clearcut (dbh, age). For mixed species is identified according to BA fraction.
 real (kind=8), intent(in) :: thinInt !parameter that determines the thinning intensity; from below (thinInt>1) or above (thinInt<1);
                    !thinInt=999. uses the default value from tapio rules
 ! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning ! user defined n of thinnings.

!!!ground vegetation variables
 real (kind=8) :: AWENgv(4)  !!! ground vegetation, Yasso params.
 !integer :: gvRun !!! flag for including ground vegetation !fvec
 real (kind=8), intent(inout) :: fAPAR(nYears), GVout(nYears, 5) ! GVout contains: fAPAR_gv, litGV, photoGV, Wgv,GVnpp !!! ground vegetation
 real (kind=8), intent(inout) :: dailyPRELES((nYears*365), 3) ! GPP, ET, SW
 real (kind=8), intent(inout) :: initVar(7, nLayers), P0y(nYears,2), ETSy(nYears), initCLcutRatio(nLayers) ! initCLcutRatio sets the initial layer compositions after clearcut.
 real (kind=8), intent(inout) :: siteInfo(10)
 real (kind=8), intent(inout) :: output(nYears, nVar, nLayers, 2), energyWood(nYears, nLayers, 2) ! last dimension: 1 is for stand and 2 is for harvested sum of wood.
 real (kind=8), intent(inout) :: soilCinOut(nYears, 5, 3, nLayers), soilCtotInOut(nYears) ! dimensions: nyears, AWENH, woody/fineWoody/foliage, layers
 real (kind=8), intent(inout) :: pYasso(35), weatherYasso(nYears,3), litterSize(3, nSp) ! litterSize dimensions: treeOrgans, species

!! Parameters internal to the model
 real (kind=8) :: prelesOut(16), fAPARsite, fAPARgvX, fAPARtrees, lastGVout(5), minFapar  !!!state of the GV at the last year
 real (kind=8) :: fAPARlayers(1+nLayers), LUElayers(1+nLayers),LUEsite
 real (kind=8) :: leac=0 ! leaching parameter for Yasso, not used
 real (kind=8), DIMENSION(nLayers, 5) :: fbAWENH, folAWENH, stAWENH
 real (kind=8), DIMENSION(nLayers) :: Lb, Lf, Lst
! real (kind=8),DIMENSION(nLayers) :: speciesIDs
 real (kind=8), DIMENSION(nLayers) :: valX
 integer, DIMENSION(nLayers) :: layerX
 real (kind=8) :: STAND(nVar), STAND_tot(nVar), param(npar) ! temp variables fillled for each layer and fills  output(nYear, nSites, nVar).
 integer :: i, ij, ijj,dayx, species, layer, nSpec, ll ! tree species 1,2,3 = scots pine, norway spruce, birch
 integer, allocatable :: indices(:)
 real(kind=8) :: rPine, rBirch, perBAmort, pMort!,dailyPRELESnoStored((nYears*365), 3) ! GPP, ET, SW

 real (kind=8) :: p0_ref, ETS_ref, P0yX(nYears, 2)
 integer :: time, ki, year, yearX, Ainit, countThinning, domSp(1)
 real (kind=8) :: step, totBA,GVnpp(nYears)

 real (kind=8) :: stand_all(nVar, nLayers)
 real (kind=8) :: outt(nVar, nLayers, 2)
 real (kind=8) :: modOut((nYears+1), nVar, nLayers, 2)
 real (kind=8) :: soilC((nYears+1), 5, 3, nLayers), soilCtot((nYears+1))
 real (kind=8) :: par_phib, par_phic, par_alfat, par_alfar1, par_alfar2, par_alfar3, par_alfar4
 real (kind=8) :: par_alfar5, par_etab, par_k, par_vf, par_vr, par_sla, &
      par_mf, par_mr, par_mw, par_vf0, mrFact, par_vr0 ! parameters filled from Crobas-vector
 real (kind=8) :: par_z, par_rhos, par_cR, par_x, Light, MeanLight(nLayers), &
      par_mf0,par_mr0,par_mw0,par_zb !par_zb: ratio between dead branches l and mean pipe lenght
 ! Light is light at canopy bottom, MeanLight is layer-specific mean. Calculated by &
 ! Llight extinction subroutine.
 real (kind=8) :: par_sarShp, par_S_branchMod ! crobas param
 real (kind=8) :: par_rhof, par_rhor, par_rhow, par_c, par_beta0, par_betab, par_betas ! crobas param
 real (kind=8) :: par_s1, par_p0, par_ksi, par_cr2,par_kRein,Rein, c_mort ! crobas param
 real (kind=8) :: BA, dA, dB, reineke(nLayers), dN, wf_test, par_thetaMax, par_H0max, par_kH, par_gamma,par_H0 ! more variables and crobas params.
 real (kind=8) :: par_rhof0, par_rhof1, par_rhof2, par_aETS, dHcCum, dHCum,pars(30), thinningType = 0. ! thinningType initialization zero, see thinning subroutines. Determines type of thinning that will occur next.

!management routines
 real (kind=8) :: A_clearcut, D_clearcut,H_clearcut, BAr(nLayers), BA_tot, BA_lim, BA_thd, ETSthres = 1000
 real (kind=8) :: dens_thd, dens_lim, Hdom_lim ! thinning parameters

!define varibles
 real (kind=8) :: sitetype, P0, age, meantemp, mintemp, maxtemp, rainfall, ETS
 real (kind=8) :: H, D, B, Hc, Cw, Lc, N, Ntree, Ntot, dNtot
 real (kind=8) :: wf_treeKG, wf_STKG, sar_con, sar_ell, rc, ppow, sar, W_stem
 real (kind=8) :: lproj, leff, laPer_sar, keff, slc
 real (kind=8) :: hb, A, B2, beta0, beta1,beta2, betab !,betas
 real (kind=8) :: c,dHc,dH,dLc,g0,g1,g2,g3,g4,g5
 real (kind=8) :: npp, p_eff_all, gammaC, betaC, W_c, W_s, Wdb, Wsh,nppCost !nppCost is the npp used in growth that considers the cost of ECM to the tree
 real (kind=8) :: p_eff, par_alfar, p, gpp_sp
 real (kind=8) :: s0,par_s0scale,par_sla0,par_tsla
 real (kind=8) :: age_factor,par_fAa,par_fAb,par_fAc ! Changes in fine root allocation for young seedlings, affecting beta0 and alfar
 real (kind=8) :: weight, dNp,dNb,dNs,perVmort
 real (kind=8) :: W_wsap, respi_m, respi_tot, V_scrown, V_bole, V, Vold
 real (kind=8) :: coeff(nLayers), denom,W_froot,W_croot, lit_wf,lit_froot
 real (kind=8) :: S_wood,Nold, Nthd, S_branch,S_fol,S_fr,W_branch,Vmort
 real (kind=8) :: W_stem_old,wf_STKG_old,W_bh, W_crh,W_bs, W_crs,dW_bh,dW_crh,dWdb,dWsh
!variables for random mortality calculations & disturbances
real (kind=8) :: Nmort, BAmort, VmortDist(nLayers),deltaSiteTypeFert=1.
!!ECMmodelling
 real (kind=8) :: r_RT, rm_aut_roots, litt_RT, exud(nLayers), P_RT
 real (kind=8) :: Cost_m,normFactETS !normFactP,normFactETS,!!Cost_m is the "apparent maintenance respiration" rate of fine roots when C input to the fungi has been taken into account.

!fix parameters
 real (kind=8) :: qcTOT0,Atot,fAPARprel(365)
!v1 version definitions
 real (kind=8) :: theta,Tdb=10.,f1,f2, Gf, Gr,mort
 real (kind=8) :: ETSmean, BAtapio(2), tapioOut(3)
 logical :: doThin, early = .false., flagInitWithThin = .false.
 real (kind=8) :: Hdom,thinClx(nYears,2),pDomRem, randX
 !!user thinnings
real (kind=8) :: pHarvTrees, hW_branch, hW_croot, hW_stem, hWdb
real (kind=8) :: remhW_branch, remhW_croot,remhW_stem,remhWdb

integer :: CO2model, AinitFix,etmodel, gvRun, fertThin, ECMmod, oldLayer !not direct inputs anymore, but in prebasFlags fvec
integer, intent(inout) :: prebasFlags(9)
real (kind=8) :: dailySW(365)

!fire disturbances
real (kind=8) :: Cpool_litter_wood,Cpool_litter_green,livegrass,soil_moisture(365)
real (kind=8) :: Tmin(365),Tmax(365),FDI(365), NI((nYears*365)),n_fire_year!
!BB disturbances
real (kind=8) :: rBAspruce(nLAyers), spruceStandVars(3),pBB(5), SMI, SMIt0, intenSpruce, SHI !SMIt0 = SMI previous year

!!! 'un-vectorise' flags, fvec
etmodel = int(prebasFlags(1))
gvRun = int(prebasFlags(2))
fertThin = int(prebasFlags(3))
oldLayer = int(prebasFlags(4))
ECMmod = int(prebasFlags(5))
CO2model = int(prebasFlags(7))
AinitFix = int(prebasFlags(8))

!!set disturbance flags
! set all dist to 0 and then choose based on flag
!if(prebasFlags(6)==0) then
 disturbance_wind = .FALSE.
 disturbance_bb = .FALSE.
 disturbance_fire = .FALSE.
!endif
if(prebasFlags(6)==1) disturbance_wind = .TRUE.
if(prebasFlags(6)==2) disturbance_bb = .TRUE.
if(prebasFlags(6)==3) disturbance_fire = .TRUE.
if(prebasFlags(6)==12) then
 disturbance_wind = .TRUE.
 disturbance_bb = .TRUE.
endif
if(prebasFlags(6)==12) then
 disturbance_wind = .TRUE.
 disturbance_bb = .TRUE.
endif
if(prebasFlags(6)==13) then
 disturbance_wind = .TRUE.
 disturbance_fire = .TRUE.
endif
if(prebasFlags(6)==23) then
 disturbance_fire = .TRUE.
 disturbance_bb = .TRUE.
endif
if(prebasFlags(6)==123) then
 disturbance_wind = .TRUE.
 disturbance_bb = .TRUE.
 disturbance_fire = .TRUE.
endif


  ! open(1,file="test1.txt")
  ! open(2,file="test2.txt")

!###initialize model###!
NI(:) = dailyPRELES(:,3) !read nestorov index and reset to -999 the dailyPreles output
dailyPRELES(:,3) = -999.0
SMIt0 = output(1,46,1,2) !initialize SMI previous year
output(1,46,1,2) = 0.d0
lastGVout = 0.
thinClx = 0.
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

! Fill in model output structure with initial values and ids.
 modOut(:,1,:,1) = siteInfo(1)  !! assign siteID
 ! modOut(1,2,:,1) = initVar(8,:)  !! assign initial gammaC values !!newX
 modOut(1,11,:,1) = initVar(3,:) ! height of trees
 modOut(1,16,:,1) = initVar(7,:)
 ! modOut(1,12,:,1) = initVar(4,:)
 ! modOut(1,13,:,1) = initVar(5,:)
 modOut(1,14,:,1) = initVar(6,:)
 ! modOut(1,17,:,1) = modOut(1,13,:,1)/(pi*((modOut(1,12,:,1)/2/100)**2))
 ! modOut(1,35,:,1) =  modOut(1,13,:,1)/modOut(1,17,:,1)
 !init siteType
 modOut(1,3,:,:) = output(1,3,:,:) ! assign site type and alfar
 modOut(2:nYears,3,:,:) = output(:,3,:,:) ! assign site type and alfar
 soilCtot(1) = sum(soilC(1,:,:,:)) !assign initial soilC
 modOut(:,45,:,1) = 0. !set heterotrophic respiration to 0
 do i = 1,nLayers
  modOut(:,4,i,1) = initVar(1,i)  ! assign species
  modOut(:,7,i,1) = initVar(2,i) ! assign initAge !age can be made species specific assigning different ages to different species
  modOut(1,39,i,1) = sum(soilC(1,:,:,i)) !assign initial soilC
  modOut(:,5,i,1) = ETSy! assign ETS
  modOut(:,6,i,1) = P0yX(:,2)  ! assign P0
  modOut(1,12,i,1) = initVar(4,i)
  modOut(1,13,i,1) = initVar(5,i)
  if(modOut(1,12,i,1) > 0.) then
  modOut(1,17,i,1) = modOut(1,13,i,1)/(pi*((modOut(1,12,i,1)/2/100)**2))
  modOut(1,35,i,1) = modOut(1,13,i,1)/modOut(1,17,i,1)
  else
  modOut(1,17,i,1) = 0.
  modOut(1,35,i,1) = 0.
  endif
 enddo


!######! SIMULATION START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do year = 1, (nYears)
  VmortDist=0.

siteInfoDist(2) = siteInfoDist(2)+1 !counter for time since thinning (wind disturbance model predictor)

!!!! check if clearcut occured. If yes initialize forest (start)
  if (year == int(yearX)) then
  !if (year == int(min(yearX, nYears))) then ! yearX is the running simulation year when stand is initialized after clearcut
   !Ainit = int(min(Ainit, Ainit + nYears - yearX)) ! Ainit is the age stand is measureable. This is to avoid some special conditions that might occur in some simulation cases.
      ! Ainit = int(Ainit, Ainit + nYears - yearX))
    Ainit = int(Ainit)
    totBA = sum(modOut((year-Ainit-1), 13, :, 1)) ! BA structure before clearcut, used for estimating spec. proportions at initialization if fixBAratio is = 0 (1 is for user defined)
   do ijj = 1, nLayers
     species = int(max(1.,modOut(year, 4, ijj, 1)))  ! read species
   if (fixBAinitClarcut==1) then
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

  siteType = modOut(year,3,ijj,1) !siteInfo(3)
  !!set parameters
  par_betab = pCrobas(13,int(initVar(1,ijj)))
  par_x = pCrobas(19,int(initVar(1,ijj)))
  par_beta0 = pCrobas(12,int(initVar(1,ijj)))
  par_betas = pCrobas(14,int(initVar(1,ijj)))
  par_mf = pCrobas(8,int(initVar(1,ijj)))
  par_mr = pCrobas(9,int(initVar(1,ijj)))
  par_mw = pCrobas(10,int(initVar(1,ijj)))
  par_fAa = pCrobas(45,int(initVar(1,ijj)))
  par_fAb = pCrobas(46,int(initVar(1,ijj)))
  par_fAc = pCrobas(47,int(initVar(1,ijj)))
  par_c = pCrobas(7,int(initVar(1,ijj)))
  par_rhof = pCrobas(15,int(initVar(1,ijj)))
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

  age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - h)/par_fAc)))/par_fAa
  par_alfar = modOut(year,3,ijj,2) * age_factor
  ! par_alfar = pCrobas(int(20+min(siteType,5.)),int(initVar(1,ijj))) * age_factor
  par_rhor = par_alfar * par_rhof

  hc = modOut(year,14,ijj,1)
  B = ba/N
  Lc = h - hc
  betab =  par_betab * Lc**(par_x-1)
  beta0 = par_beta0 * age_factor
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
   modOut((year-Ainit):year,48,1,2) = 0.
   do ki = 1,int(Ainit)
    do ijj = 1,nLayers
     modOut((year-Ainit+ki),7,ijj,1) = ki !#!#
     modOut((year-Ainit+ki),4,ijj,1) = initVar(1,ijj) !#!#
    enddo
   enddo
  yearX = 0
  endif
!!!! check if clearcut occured. If yes initialize forest (end)

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

! !calculate reneike and random mortality
! include 'mortalityCalc.h'

do ij = 1 , nLayers     !loop Species
 STAND=STAND_all(:,ij)
 species = int(max(1.,stand(4)))
 param = pCrobas(:,species)
 sitetype=STAND(3)

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
 par_fAa = param(45)
 par_fAb = param(46)
 par_fAc = param(47)
! do siteNo = 1, nSites  !loop sites

    !!!!update kRein and cR
   !!!!update par_kRein as a function of sitetype if parameters (param(50>-999.))) are are provided
   if(param(50)>-999.d0) call linearUpdateParam(param(50:51),stand(3),par_kRein)
   !!!!update par_cR as a function of sitetype if parameters (param(52>-999.))) are are provided
   if(param(52)>-999.d0) call linearUpdateParam(param(52:53),stand(3),par_cR)

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
  ! mort = stand(40)
  if(par_tsla .gt. 0.) then
   par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)
  else
   par_sla = par_sla
  endif

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
 else
  STAND(2) = 0. !#!#
  STAND(8:21) = 0. !#!#
  STAND(23:37) = 0. !#!#
  STAND(42:44) = 0. !#!#
  STAND(47:nVar) = 0. !#!#
 endif

 STAND_all(:,ij)=STAND
end do !!!!!!!end loop layers

!!!calculate species weight for photosynthesis
!do siteNo = 1, nSites

if (year <= maxYearSite) then
  nSpec = nSp
  ll = nLayers

  domSp = maxloc(STAND_all(13,:))
  layer = int(domSp(1))
  siteType = modOut(year,3,layer,1) !siteInfo(3)

    call Ffotos2(STAND_all,nLayers,nSpec,pCrobas,&
    nVar,nPar,MeanLight,coeff,fAPARtrees)
   STAND_all(36,:) = MeanLight
   STAND_all(23,:) = coeff

!!calculate year of replanting after a clearcut
!if scenario = "oldLayer" do not consider the old layer
   if(oldLayer==1) then
    ij=max((nLayers-1),1)
   else
    ij=nLayers
   endif
   if(sum(modOut(year,11,1:ij,1)) > 0.) yearX = 0
   if(sum(modOut(year,11,1:ij,1)) == 0. .and. yearX == 0) then
  if(AinitFix < 1) then
	  if((nYears-year)<10) then
		  Ainit = max(nint(6. + 2*sitetype - 0.005*modOut(year,5,1,1) + 2.25 + 2.0),2)!! + 2.0 to account for the delay between planting and clearcut
	  else
		  Ainit = max(nint(6. + 2*sitetype - 0.005*(sum(modOut(year:(year+9),5,1,1))/10.) + 2.25 + 2.0),2)!! + 2.0 to account for the delay between planting and clearcut
	  endif
  endif 
  yearX = Ainit + year
   endif
   !!!ground vegetation
   !!!fapar_gv compute fapar, biomasses and litter of gv with routine
   if(gvRun==1) then

! if(isnan(siteType)) siteType = siteInfo(3)
! if(siteType==0.) siteType = siteInfo(3)

    call fAPARgv(fAPARtrees, ETSmean, siteInfo(3), fAPARgvX, GVout(year,2), &
         sum(P0yX(:,1))/nYears, AWENgv,GVout(year,4))
   else
    fAPARgvX=0.
  GVout(year,:) = 0.
   endif

if(isnan(fAPARgvX)) fAPARgvX = 0.

!!calculate minimum fAPAR of last 15 years to be used in ingrowth(if active) calculations
  call minFaparCalc(fAPAR,year,minFapar,fAparFactor)

!!!calculate site fAPAR and set fAPAR for preles calculations and store
   fAPARsite = fAPARtrees + fAPARgvX
   fAPARprel(:) = fAPARsite
   fAPAR(year) = fAPARtrees  !store fAPAR trees
   GVout(year,1) = fAPARgvX !store fAPAR GV
   ! if(fAPARsite>0.) then

   ! pars(4:23) = pPRELES(4:23)
     ! pars(28:30) = pPRELES(28:30)

   if(layerPRELES == 0) then

  !run preles
     call preles(weatherPRELES(year,:,:),DOY,fAPARprel,prelesOut, pars, &
    dailyPRELES((1+((year-1)*365)):(365*year),1), &  !daily GPP
    dailyPRELES((1+((year-1)*365)):(365*year),2), &  !daily ET
    dailyPRELES((1+((year-1)*365)):(365*year),3), &  !daily SW
    etmodel,CO2model)    !type of ET model

   !store ET of the ECOSYSTEM!!!!!!!!!!!!!!
     STAND_all(22,:) = prelesOut(2)    !ET
   ! STAND_all(40,:) = prelesOut(15)  !aSW
   ! STAND_all(41,:) = prelesOut(16)  !summerSW

  !store GPP

     GVout(year,3) = prelesOut(1) * fAPARgvX/fAPARsite! GV Photosynthesis in g C m-2
	 if(GVout(year,1)<0.00000001) GVout(year,:) = 0.
     STAND_all(10,:) = prelesOut(1)/1000. * fAPARtrees/fAPARsite * coeff! trees Photosynthesis in g C m-2 (converted to kg C m-2)

!initialize for next year
     pars(24) = prelesOut(3);siteInfo(4) = prelesOut(3)!SWinit
     pars(25) = prelesOut(13); siteInfo(5) = prelesOut(13) !CWinit
     pars(26) = prelesOut(4); siteInfo(6) = prelesOut(4) !SOGinit
     pars(27) = prelesOut(14); siteInfo(7) = prelesOut(14) !Sinit
   endif

   if(layerPRELES == 1) then
    fAPARlayers(1:nLayers) = fAPARtrees * coeff
  fAPARlayers(1+nLayers) = fAPARgvX
    LUElayers(1:nLayers) = LUEtrees(int(STAND_all(4,:)))
  LUElayers(1 + nLayers) = LUEgv !!fill for GV
  LUEsite = sum(fAPARlayers * LUElayers)/fAPARsite
  pars(5) = LUEsite
    call preles(weatherPRELES(year,:,:),DOY,fAPARprel,prelesOut, pars, &
    dailyPRELES((1+((year-1)*365)):(365*year),1), &  !daily GPP
    dailyPRELES((1+((year-1)*365)):(365*year),2), &  !daily ET
    dailyPRELES((1+((year-1)*365)):(365*year),3), &  !daily SW
    etmodel,CO2model)    !type of ET model

  !store ET of the ECOSYSTEM!!!!!!!!!!!!!!
    STAND_all(22,:) = prelesOut(2)    !ET
   ! STAND_all(40,:) = prelesOut(15)  !aSW
   ! STAND_all(41,:) = prelesOut(16)  !summerSW
  !store GPP
    GVout(year,3) = prelesOut(1) * fAPARlayers(1+nLayers)/fAPARsite * &
            LUElayers(1+nLayers)/LUEsite !!!GV photosynthesis

     STAND_all(10,:) = prelesOut(1)/1000. * fAPARlayers(1:nLayers)/fAPARsite * &
              LUElayers(1:nLayers)/LUEsite! trees Photosynthesis in g C m-2 (converted to kg C m-2)


!initialize for next year
     pars(24) = prelesOut(3);siteInfo(4) = prelesOut(3)!SWinit
     pars(25) = prelesOut(13); siteInfo(5) = prelesOut(13) !CWinit
     pars(26) = prelesOut(4); siteInfo(6) = prelesOut(4) !SOGinit
     pars(27) = prelesOut(14); siteInfo(7) = prelesOut(14) !Sinit

   endif


    outt(46,1,2)  = prelesOut(7) !SMI
    SMI = prelesOut(7) !SMI
	dailySW = dailyPRELES((1+((year-1)*365)):(365*year),3)

endif
!enddo !! end site loop

do ij = 1 , nLayers
 outt(30,ij,2) = 0.
 STAND=STAND_all(:,ij)
 species = int(max(1.,stand(4)))
 param = pCrobas(:,species)
 sitetype=stand(3)

   !!!reset annual litterfall
stand_all(26:29,ij) = 0.
stand(26:29) = 0.
s_fol = 0.
S_fr = 0.
S_branch = 0.
S_wood = 0.

!initialize ECMmodelling
r_RT = 0.d0
rm_aut_roots = 0.d0
litt_RT = 0.d0
exud(ij) = 0.d0

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
 par_ksi = param(38)
 par_rhof1 = 0.!param(20)
 par_Cr2 = 0.!param(24)
 par_sla0 = param(39)
 par_tsla = param(40)
 par_zb = param(41)
 par_fAa = param(45)
 par_fAb = param(46)
 par_fAc = param(47)
! do siteNo = 1, nSites  !start site loop

    !!!!update kRein and cR
   !!!!update par_kRein as a function of sitetype if parameters (param(50>-999.))) are are provided
   if(param(50)>-999.d0) call linearUpdateParam(param(50:51),stand(3),par_kRein)
   !!!!update par_cR as a function of sitetype if parameters (param(52>-999.))) are are provided
   if(param(52)>-999.d0) call linearUpdateParam(param(52:53),stand(3),par_cR)

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
  p0 = STAND(6)/1000.  ! convert    g C m-2    to    kg C m-2  !#!#
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
  p_eff_all = STAND(10)*P0yX(year, 2)/P0yX(year, 1) !!##!! 2 smoothing PHOTOSYNTHESIS
  weight = STAND(23)
  A = STAND(16)
  dHc=stand(52)
  S_fol = STAND(26)
  S_fr = STAND(27)
  S_branch = STAND(28)
  if(par_tsla .gt. 0.) then
   par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)
  else
   par_sla = par_sla
  endif

  rc = Lc / (H-1.3) !crown ratio
  B = BA / N
  ! A = rc * B

  wf_STKG = STAND(33)
  wf_treeKG = STAND(34)
  ! B = STAND(35)
  Light = STAND(36)
  hb = par_betab * Lc ** par_x
  Cw = 2. * hb

  age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - H)/par_fAc)))/par_fAa !! Root allocation of seedlings
if (N>0.) then

!!!!###here starts stand2 subroutine!!!!!!!!!!!#########
  par_alfar = modOut(year,3,ij,2) * age_factor
  ! if (sitetype <= 1.) then
   ! par_alfar = par_alfar1 * age_factor
  ! else if (sitetype==2.) then
   ! par_alfar = par_alfar2 * age_factor
  ! else if (sitetype==3.) then
   ! par_alfar = par_alfar3 * age_factor
  ! else if (sitetype==4.) then
   ! par_alfar = par_alfar4 * age_factor
  ! else
   ! par_alfar = par_alfar5 * age_factor
  ! end if

!!relate metabolic and structural parameters to site conditions

!  par_mf = par_mf0 * p0 / p0_ref
!  par_mr = par_mr0 * p0 / p0_ref
!  par_mw = par_mw0 * p0 / p0_ref


  ! par_H0 = par_H0max * (1 - exp(-par_kH * ETS/par_alfar)) !!! attempt to improve model for diameter/heigh allocation, not used currently.
  ! theta = par_thetaMax / (1. + exp(-(H - par_H0)/(par_H0*par_gamma)))   !!!! see above, get zero now
  theta = par_thetaMax

  mrFact = max(0., par_aETS * (ETS_ref-ETS)/ETS_ref) !!!new version
  if(ECMmod==1) mrFact = 0.
  par_mr = par_mr0* p0 / p0_ref + (1+par_c) * mrFact / par_vr0    !!!new version !!newX
  par_mf = par_mf0* p0 / p0_ref !!newX
  par_mw = par_mw0* p0 / p0_ref !!newX

  par_rhof0 = par_rhof1 * ETS_ref + par_rhof2 ! rho: pipe model parameter for foliage
  par_rhof = par_rhof1 * ETS + par_rhof2
  par_vf = par_vf0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref)
  par_vr = par_vr0 / (1. + par_aETS * (ETS-ETS_ref)/ETS_ref) !!new version
  par_rhor = par_alfar * par_rhof !rho: pipe model parameter for fine roots

    ! -------------------------------------
    !GPP all STAND$species   UNITS: g C  /  m2
    ! -------------------------------------
        p_eff = p_eff_all!weight * p_eff_all
    ! gpp_sp = weight * STAND(10)

    if(wf_STKG > 0.) then
!if(H < 10.) then
 !       s0 = min(0.35 * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.)
!else
        s0 = min(par_s0scale * P0 * par_k * par_sla, P_eff / wf_STKG * 10000.)
! endif
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
      age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - H)/par_fAc)))/par_fAa
            beta0 = par_beta0 * age_factor
            beta1 = (beta0 + betab + par_betas) !!newX
            beta2 = 1. - betab - par_betas     !!newX
      betaC = (beta1 + gammaC * beta2) / par_betas
      W_wsap = stand(47)
            ! wf_treeKG = par_rhof * A  !!newX
      ! wf_STKG = N * wf_treeKG  !!newX
      W_froot = stand(25)
      W_c = stand(48) !sapwood stem below Crown
      W_s = stand(49) !sapwood stem within crown
      W_branch = stand(24) !branches biomass
      W_croot = stand(32) !W_stem * (beta0 - 1.)  !coarse root biomass
      Wsh = stand(50)
      Wdb = stand(51)
      W_bh = stand(53)
      W_crh = stand(54)
      W_bs =  par_rhow * A * N * betab * Lc
      W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)  !coarse root biomass


      ! ECM modelling
      if(ECMmod==1) then !!!ECMmodelling
          normFactETS = 1. + par_aETS * (ETS-ETS_ref)/ETS_ref
        ! normFactP = p0 / p0_ref
        ! call CUEcalc(ETSmean, sitetype,par_mr0,W_froot,r_RT,rm_aut_roots,litt_RT,exud(ij),normFactP,normFactETS,P_RT,pECMmod) !!!ECMmodelling
        call CUEcalc(ETSmean, sitetype,par_mr,W_froot,r_RT,rm_aut_roots,litt_RT,exud(ij),P_RT,normFactETS,pECMmod) !!!ECMmodelling

        modOut((year+1),45,ij,1) = P_RT  !add priming to heterotrophic respiration
        Respi_m = par_mf * wf_STKG + par_mw * W_wsap + rm_aut_roots * W_froot  !!!ECMmodelling
        Cost_m  = par_mf * wf_STKG + par_mw * W_wsap + r_RT * W_froot  !!!ECMmodelling
      else !!! only consider fine root
        Respi_m = (par_mf + par_alfar*par_mr)* wf_STKG + par_mw * W_wsap
        Cost_m = Respi_m
      endif
      nppCost = max(0.,(gpp_sp - Cost_m / 10000.) / (1.+par_c))  !!newX
      npp = max(0.,(gpp_sp - Respi_m / 10000.) / (1.+par_c))  !!newX


      Respi_tot = gpp_sp - npp
         ! ! litter fall in the absence of thinning
      S_fol = S_fol + wf_STKG / par_vf  !foliage litterfall
      S_fr  = S_fr + W_froot / par_vr  !fine root litter
    S_branch = max(0.,S_branch + Wdb/Tdb)

    ! S_branch = S_branch + N * par_rhow * betab * A * (dHc + theta*Lc)

        !Height growth-----------------------
    f1 = nppCost*10000 - (wf_STKG/par_vf) - (W_froot/par_vr) - (theta * W_wsap)
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


 !-----------------------------------
        !crown rise
!         if(H - Hc > par_Cr2*100./sqrt(N)) then
!        if(2.*hb > 100./sqrt(N) ) then
        dHc = min(gammaC * dH,(H - Hc))

if(time==1)then
      dHcCum = 0.
      dHCum = 0.
endif
        dHcCum = dHcCum + dHc
        dHCum = dHCum + dH
!      else
!    dHc = 0  !CAN BE DIFFERENT FROM THE PAPER HARKONEN ET AL. 2013 CANADIAN JOURNAL, SEE THE EQUATION THERE
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
            dWsh = max(0.,par_rhow * par_z * A/Lc * dHc * Hc * N  + theta*(W_c + W_s))
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
    age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - H)/par_fAc)))/par_fAa
        beta0 = par_beta0 * age_factor
        beta1 = (beta0 + betab + par_betas) !!newX
        beta2 = 1. - betab - par_betas     !!newX
    betaC = (beta1 + gammaC * beta2) / par_betas
    W_froot = par_rhor * A * N  !!to check  !!newX
    W_wsap = par_rhow * A * N * (beta1 * H + beta2 * Hc) !!newX
    W_c = par_rhow * A * N * Hc !sapwood stem below Crown
    W_s = par_rhow * A * N * par_betas * Lc !sapwood stem within crown
    W_bs =  par_rhow * A * N * betab * Lc !branches biomass
        W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)  !coarse root sapwood biomass
    ! Wsh = Wsh + par_z * gammaC * W_c/Lc * dH + theta*(W_c + W_s)
    W_stem = W_c + W_s + Wsh
    V = W_stem / par_rhow
    W_branch = W_bs + W_bh
    W_croot = W_crs + W_crh

  age = age + step

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


!set species from thinning matrix (strart)
   species = int(thinning(countThinning,2))
   stand(4) = thinning(countThinning,2)
     !!!check if ingrowth and calculate dominant species
   if(D==0.d0 .and. H==0.d0 .and. thinning(countThinning,6)==-777.d0) then
    domSp = maxloc(STAND_all(13,:))
	layer = int(domSp(1))
	species = int(max(1.,stand_all(4,layer)))
	stand(4) = max(1.,stand_all(4,layer))
	thinning(countThinning,2) = stand(4)
	
    modOut(:,3,ij,2) = pCrobas(int(20 + stand(3)),species)
	prebasFlags(9) = ij
	
   endif
!set species from thinning matrix (end)


  ! if(year >= yearX) then
    STAND_tot = STAND
    IF (thinning(countThinning,6) < STAND_tot(13)) siteInfoDist(2) = 0

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
      species = int(max(1.,stand(4)))
  if(pCrobas(2,species)>0.) energyWood(year,ij,1) = energyWood(year,ij,2) / pCrobas(2,species)
      S_branch = max(0.,((W_branch) * (1-energyRatio) + S_branch + Wdb + &
          W_stem* (1-harvRatio)* (1-energyRatio) + &
          (0.3 * (1-energyRatio)+0.7) * W_croot *0.83))
      S_wood = S_wood + (0.3 * (1-energyRatio)+0.7) * W_croot *0.17
     else
      S_branch = max(0.,(W_branch + Wdb + W_croot*0.83 + S_branch + W_stem* (1-harvRatio)))
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

     !!!check if ingrowth and calculate the number of trees
     if(D==0.d0 .and. H==0.d0 .and. thinning(countThinning,6)==-777.d0) then
      BA = pi*((0.5d0/200.d0)**2.d0)*min((500.d0/minFapar),4000.d0)
     else
      BA = thinning(countThinning,6)
     endif

     if (thinning(countThinning,4) /= -999.) H = thinning(countThinning,4)
     if (thinning(countThinning,7) /= -999.) stand(14) = thinning(countThinning,7)
     if (thinning(countThinning,5) /= -999.) D = thinning(countThinning,5)
     Hc=stand(14)
     Lc = H - Hc !Lc
     rc = Lc / (H-1.3) !crown ratio
     Nold = N
     wf_STKG_old = wf_STKG
     W_stem_old = W_stem
     N = BA/(pi*((D/2./100.)**2.)) ! N
     Nthd = max(0.,(Nold-N)) ! number of cutted trees
     B = BA/N!(pi*((D/2/100)**2))
     if (thinning(countThinning,10) /= -999.) then
       A = thinning(countThinning,10)
     else
       A = stand(16) * B/stand(35)
     endif
     ! Update dependent variables
     hb = par_betab * Lc**par_x
     gammaC = par_cR/stand(36)
     par_rhof = par_rhof1 * ETS + par_rhof2
     par_rhor = par_alfar * par_rhof
     betab = hb/Lc
     age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - H)/par_fAc)))/par_fAa
     beta0 = par_beta0 * age_factor
     beta1 = (beta0 + betab + par_betas) !!newX
     beta2 = 1. - betab - par_betas     !!newX
     betaC = (beta1 + gammaC * beta2) / par_betas

     !!!reinitialize the stand and some variables when the thinning matrix is used to initialize the stand in the middle of the runs(start)
     if(Nold==0.) then
      flagInitWithThin = .true.
      Nold = N
      Wdb = 0.
      A = par_ksi/par_rhof * Lc**par_z
      stand(7) = Ainit
      yearX = 0.
     endif
     !!!reinitialize Nold and some variables when the thinning matrix is used to initialize the stand (end)

      if(isnan(stand(50))) stand(50) = 0
      if(isnan(stand(53))) stand(53) = 0
      if(isnan(stand(54))) stand(54) = 0
      W_bh = stand(53)
      W_crh = stand(54)
      wf_treeKG = max(par_rhof * A,0.)
      wf_STKG = max(N * wf_treeKG,0.)
      W_froot = max(par_rhor * A * N,0.)  !!to check  !!newX
      W_wsap = max(par_rhow * A * N * (beta1 * H + beta2 * Hc),0.) !!newX
      W_c = max(par_rhow * A * N * Hc,0.) !sapwood stem below Crown
      W_s = max(par_rhow * A * N * par_betas * Lc,0.) !sapwood stem within crown
      W_bs =  max(par_rhow * A * N * betab * Lc,0.) !branches biomass
      W_crs = max(par_rhow * beta0 * A * H * N,0.) !W_stem * (beta0 - 1.)  !coarse root biomass
      ! Wsh = Wsh + par_z * gammaC * W_c/Lc * dH + theta*(W_c + W_s)
      Wsh = max(stand(50) * N/Nold,0.) !!this is not proportional to the size of the standing trees
      W_stem = max(W_c + W_s + Wsh,0.)
      W_crh = max(W_crh * N/Nold,0.)
      W_bh = max(W_bh * N/Nold,0.)
      V = max(W_stem / par_rhow,0.)
      W_croot = max(W_crs + W_crh,0.)
      W_branch = max(W_bs + W_bh,0.)
      Wdb = max(Wdb * N/Nold,0.)
  !! calculate litter including residuals from thinned trees
  !energyCut
  pHarvTrees = thinning(countThinning,11)
  S_fol = max(0.,stand(26) + stand(33) - wf_STKG)
  S_fr = max(0.,stand(27) + stand(25) - W_froot)

  hW_branch = max(0.,(stand(24) - W_branch)* pHarvTrees)
  hW_croot = max(0.,(stand(32) - W_croot)* pHarvTrees)
  hW_stem = max(0.,(stand(31) - W_stem)* pHarvTrees)
  hWdb = max(0.,(stand(51) - Wdb)* pHarvTrees)
  remhW_branch = max(0.,(stand(24) - W_branch) * (1.-pHarvTrees))
  remhW_croot = max(0.,(stand(32) - W_croot) * (1.-pHarvTrees))
  remhW_stem = max(0.,(stand(31) - W_stem) * (1.-pHarvTrees))
  remhWdb = max(0.,(stand(51) - Wdb) * (1.-pHarvTrees))

  if(energyCut==1.) then
  species = int(max(1.,stand(4)))

  energyWood(year,ij,2) = max(0.,(hW_branch + hW_croot*0.3 + &
                    (hW_stem) * (1-harvRatio)) * energyRatio)
  if(pCrobas(2,species)>0.) energyWood(year,ij,1) = max(0.,energyWood(year,ij,2) / pCrobas(2,species))
  S_branch = max(0.,stand(28) + hW_branch * (1-energyRatio) + hWdb +&
           (0.3 * (1-energyRatio)+0.7) * hW_croot * 0.83 + &
           hW_stem * (1-harvRatio) * (1-energyRatio))
  S_wood = max(0.,stand(29) +(0.3 * (1-energyRatio)+0.7) * (stand(32) - W_croot) *0.17)
  else
    S_branch = max(0.,stand(28)+hW_branch+hWdb+hW_croot*0.83 + &
             hW_stem * (1-harvRatio))
    S_wood = max(0.,stand(29)  + hW_croot*0.17)
  endif

  !!! if part of the thinned trees is not harvested (e.g.,residues from disturbances) litterfall is updated
  if(pHarvTrees < 1.) then
  S_branch = S_branch + remhW_branch + remhW_croot * 0.83 + remhWdb
  S_wood = S_wood + remhW_croot*0.17 + remhW_stem
  endif


     outt(11,ij,2) = STAND_tot(11)
     outt(12,ij,2) = STAND_tot(12)
     outt(13,ij,2) = (STAND_tot(13) - BA) * pHarvTrees
     outt(14,ij,2) = STAND_tot(14)
     outt(15,ij,2) = STAND_tot(15)
     outt(16,ij,2) = STAND_tot(16)
     outt(17,ij,2) = (Nthd) * pHarvTrees
     outt(18:23,ij,2) = -999.
     outt(24,ij,2) = (STAND_tot(24) - W_branch) * pHarvTrees
     outt(25,ij,2) = (STAND_tot(25) - W_froot) * pHarvTrees
     outt(26:29,ij,2) = -999.
     outt(30,ij,2) = outt(30,ij,2) + max((STAND_tot(30) - V) * pHarvTrees,0.)
     VmortDist(ij) = stand(42) + max((STAND_tot(30) - V) * (1-pHarvTrees),0.)
     outt(31,ij,2) = outt(31,ij,2) + max((STAND_tot(31) - W_stem) * pHarvTrees,0.)
     outt(32,ij,2) = max((STAND_tot(32) - W_croot) * pHarvTrees,0.)
     outt(33,ij,2) = max((STAND_tot(33) - wf_STKG) * pHarvTrees,0.)
     outt(34,ij,2) = max((STAND_tot(34)*Nold - wf_treeKG*N)/Nthd,0.)
     outt(35,ij,2) = -999.; outt(36,ij,2)= -999.

  if(flagInitWithThin) then
     flagInitWithThin = .false.
  endif
    stand(4) = thinning(countThinning,2)
    stand(11) = H
    stand(12) = D
    stand(13) = BA
    stand(16) = A
    stand(14) = Hc
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
  ! endif
  countThinning = countThinning + 1

   End If
 End If
   STAND_all(:,ij)=STAND
 end do !!!!end loop species
 end do !!!!end loop inttimes

!Perform thinning or defoliation events for this time period using standard management routines!!!!!!!!!!!!!!!!
!!!!test for clearcut!!!!
  if(oldLayer==1) then
   ll=max((nLayers-1),1)
  else
   ll=nLayers
  endif
 domSp = maxloc(STAND_all(13,1:ll))
 layer = int(domSp(1))

if (ClCut > 0.5 .or. outdist(max(INT(year-1),1), 9) == 1.) then !outdist(,9): cc-inducing wind dist in previous year

  species = int(max(1.,stand_all(4,layer)))
  D_clearcut = clct_pars(species,1)
  A_clearcut = clct_pars(species,2)
  H_clearcut = clct_pars(species,3)
  D = stand_all(12,layer)
  age = stand_all(7,layer)
  if ((D > D_clearcut) .or. (H > H_clearcut) .or. (age > A_clearcut) &
		.or. outdist(max(INT(year-1),1), 9) == 1.) then !outdist(,9): cc-inducing wind dist in previous year

  if (outdist(max(INT(year-1),1), 9) == 1.)   outdist(year, 9) = 2. !set disturbance-induced cc flag to 2 (= 'conducted') otherwise, this triggers clearcuts over and over again...)
  ! modOut(year+1,1,2,2) = 1. !flag for clearcut
  thinClx(year,2) = 1 !flag for clearcut
   !if fertilization at thinning is active reset flagFert
  if(fertThin > 0) then
  flagFert = 0
  endif

  if(oldLayer==1) then !if oldLayer
  !!check dominant layer
   ! domSp = maxloc(STAND_all(13,1:(nLayers - 1)))
   ! layer = int(domSp(1))

   !set siteType and alfar for old layer
  modOut(:,3,nLayers,1) = modOut(year,3,species,1)
  modOut(:,3,nLayers,2) = modOut(year,3,species,2)
  outt(3,nLayers,:) = modOut(year,3,species,1:2)

  !!!!calculate percentage of trees remaining after clearcut(pDomRem)
  call random_number(randX)
   pDomRem =   max((randX*5.+5.)/100.* &  !!randomly sample between 5 and 10 %
      sum(stand_all(13,1:ll))/stand_all(13,layer),0.05)

  !update old layer
   stand_all(:,nLayers) = stand_all(:,layer)
   stand_all(42,nLayers) = 0.
   stand_all((/9,10,13,17,18,37,38,40,43,44,53,54/),nLayers) = &
  stand_all((/9,10,13,17,18,37,38,40,43,44,53,54/),layer) * (pDomRem)
   stand_all(24:34,nLayers) = stand_all(24:34,layer) * (pDomRem)
   stand_all(47:51,nLayers) = stand_all(47:51,layer) * (pDomRem)
  !update dominant layer
   stand_all((/9,10,13,17,18,37,38,40,43,44,53,54/),layer) = &
  stand_all((/9,10,13,17,18,37,38,40,43,44,53,54/),layer) * (1-pDomRem)
   stand_all(24:34,layer) = stand_all(24:34,layer) * (1-pDomRem)
   stand_all(47:51,layer) = stand_all(47:51,layer) * (1-pDomRem)

  !!!implement clearcut by layer (start) (not for old layer in oldLayer sceario)
   do ij = 1, ll

! (JH interpretation) filling harvested dimension of output with stand_all variables..
 !for some variables additively:
    outt((/9,10,13,16,17,18,24,25,30,31,32,33/),ij,2) = outt((/9,10,13,16,17,18,24,25,30,31,32,33/),ij,2) + &
                    stand_all((/9,10,13,16,17,18,24,25,30,31,32,33/),ij)
! for some by replacement. Problematic for ,42, as 42,2 is used to transfer killed V that is salvage logged to next year in case the harvest limit has already been met.
!  outt((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,42,43,44,49,50,51,52,53,54/),ij,2) = &
!    stand_all((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,42,43,44,49,50,51,52,53,54/),ij)
  outt((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,44,49,50,51,52,53,54/),ij,2) = &
      stand_all((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,44,49,50,51,52,53,54/),ij) ! excluding 42 from this replacement to fix bug where vmort (,1) is harvested in the year after the management reaction
  !energyCut
    S_fol = stand_all(33,ij) + stand_all(26,ij)
  S_fr = stand_all(25,ij) + stand_all(27,ij)
  if(energyCut==1.) then
   energyWood(year,ij,2) = energyWood(year,ij,2) + (stand_all(24,ij) + &
          stand_all(32,ij)*0.3 + stand_all(31,ij) * (1-harvRatio)) * energyRatio
  species = int(max(1.,stand_all(4,ij)))
if(pCrobas(2,species)>0.) energyWood(year,ij,1) = energyWood(year,ij,2) / pCrobas(2,species)
   S_branch = max(0.,(stand_all(28,ij) + (stand_all(24,ij)) * (1-energyRatio) + &
    stand_all(51,ij) + (0.3 * (1-energyRatio)+0.7) * stand_all(32,ij) *0.83 + &
    stand_all(31,ij)* (1-harvRatio) * (1-energyRatio)))
   S_wood = (0.3 * (1-energyRatio)+0.7) * stand_all(32,ij) *0.17 + stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts
  else
   S_branch = max(0.,(stand_all(51,ij)+stand_all(24,ij)+stand_all(28,ij)+stand_all(32,ij)* 0.83 +&
      stand_all(31,ij)* (1-harvRatio)))
   S_wood = stand_all(32,ij) *0.17 + stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts
  endif
  !energyCut
    stand_all(2,ij) = 0. !!newX
    stand_all(7,ij) = 0.
  stand_all(8,ij) = 0.
    stand_all(10:17,ij) = 0.
    stand_all(19:21,ij) = 0.
    stand_all(23:38,ij) = 0.
    stand_all(41,ij) = 0.
    stand_all(43,ij) = 0.
    stand_all(47:nVar,ij) = 0.
    stand_all(26,ij) = S_fol
    stand_all(27,ij) = S_fr
    stand_all(28,ij) = S_branch
    stand_all(29,ij) = S_wood
   enddo   !!!implement clearcut by layer (end) (not for old layer in oldLayer sceario)
  else !end if oldLayer
    !!!implement clearcut by layer (start)
   do ij = 1, nLayers

    outt((/9,10,13,16,17,18,24,25,30,31,32,33/),ij,2) = outt((/9,10,13,16,17,18,24,25,30,31,32,33/),ij,2) + &
                    stand_all((/9,10,13,16,17,18,24,25,30,31,32,33/),ij)
  ! outt((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,42,43,44,49,50,51,52,53,54/),ij,2) = & !45:47,
  !   stand_all((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,42,43,44,49,50,51,52,53,54/),ij) !45:47,
  outt((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,44,49,50,51,52,53,54/),ij,2) = & !45:47, !!! REMOVAL OF 42 her as well (see ~50 lines above)
    stand_all((/6,7,8,11,12,14,15,19,20,21,22,23,26,27,28,29,34,35,36,37,38,39,40,41,44,49,50,51,52,53,54/),ij) !45:47,

  !energyCut
    S_fol = stand_all(33,ij) + stand_all(26,ij)
    S_fr = stand_all(25,ij) + stand_all(27,ij)
  if(energyCut==1.) then
   energyWood(year,ij,2) = energyWood(year,ij,2) + (stand_all(24,ij) + &
          stand_all(32,ij)*0.3 + stand_all(31,ij) * (1-harvRatio)) * energyRatio
species = int(max(1.,stand_all(4,ij)))
if(pCrobas(2,species)>0.)   energyWood(year,ij,1) = energyWood(year,ij,2) / pCrobas(2,species)
   S_branch = max(0.,(stand_all(28,ij) + (stand_all(24,ij)) * (1-energyRatio) + &
    stand_all(51,ij) + (0.3 * (1-energyRatio)+0.7) * stand_all(32,ij) *0.83 + &
    stand_all(31,ij)* (1-harvRatio) * (1-energyRatio)))
   S_wood = (0.3 * (1-energyRatio)+0.7) * stand_all(32,ij) *0.17+ stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts

  else
   S_branch = max(0.,(stand_all(51,ij)+stand_all(24,ij)+stand_all(28,ij)+stand_all(32,ij)* 0.83 +&
      stand_all(31,ij)* (1-harvRatio)))
   S_wood = stand_all(32,ij) *0.17 + stand_all(29,ij) !(1-harvRatio) takes into account of the stem residuals after clearcuts
  endif
  !energyCut
    stand_all(2,ij) = 0. !!newX
    stand_all(7,ij) = 0.
  stand_all(8,ij) = 0.
    stand_all(10:17,ij) = 0.
    stand_all(19:21,ij) = 0.
    stand_all(23:38,ij) = 0.
    stand_all(41,ij) = 0.
    stand_all(43,ij) = 0.
    stand_all(47:nVar,ij) = 0.
    stand_all(26,ij) = S_fol
    stand_all(27,ij) = S_fr
    stand_all(28,ij) = S_branch
    stand_all(29,ij) = S_wood

   enddo !!!implement clearcut by layer (end)
  endif !!!end if oldLayer
 endif
endif

!!!!test for thinnings!!!!
 !!!!!!!for coniferous dominated stands!!!!!!
if(defaultThin == 1.) then
  if(oldLayer==1) then
   ll=max((nLayers-1),1)
  else
   ll=nLayers
  endif

 ! sitetype = siteInfo(3)
 BA_tot = sum(stand_all(13,1:ll))!+stand_all(13,2)+stand_all(13,3)
 BAr = stand_all(13,1:ll)/BA_tot
 domSp = maxloc(STAND_all(13,1:ll))
 layer = int(domSp(1))
 siteType = modOut(year,3,layer,1) !siteInfo(3)
 H = stand_all(11,layer)
 species = int(max(1.,stand_all(4,layer)))

 ! counting the dominant height of the dominant species
 Hdom = pCrobas(42,species)*exp(-1/max((H-1.3),0.001))+pCrobas(43,species)*H
 Ntot = sum(STAND_all(17,:))
  !! here we decide what thinning function to use; 3 = tapioThin, 2 = tapioFirstThin, 1 = tapioTend
 call chooseThin(species, siteType, ETSmean, Ntot, Hdom, tTapio, ftTapio, thinningType)
 ! thinx = thinningType

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

   siteInfoDist(2) = 0 !reset counter for time since thinning (wind dist model predictor)

 !!!fertilization at thinning
  if(fertThin == 3 .and. flagFert<1 .and. siteType>3. .and. siteType<6.) then
    flagFert=1

    yearsFert = max(1,min((nYears) - year,nYearsFert))
    modOut((year+1):(year+yearsFert),3,:,1) = max(1.,siteType-1.)
    call calcAlfar(sitetype,initVar(1,:),pCrobas, &
        nLayers,alfarFert,nSp,nYearsFert,npar,deltaSiteTypeFert)
    modOut((year+1):(year+yearsFert),3,:,2) = alfarFert(1:yearsFert,:,2)
  endif
!!!end fertilization at thinning

  ! modOut(year+1,1,1,2) = thinx !flag for thinning
  thinClx(year,1) = thinningType  !flag for thinning
  !if scenario = "oldLayer" do not consider the old layer in the thinnings

  do ij = 1, ll

   if(stand_all(17,ij)>0.) then
    STAND_tot = stand_all(:,ij)
    species = int(max(1.,stand_all(4,ij)))
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
    par_fAa = param(45)
    par_fAb = param(46)
    par_fAc = param(47)
    par_rhof1 = 0.!param(20)
    par_Cr2 = 0.!param(24)
    par_rhof = par_rhof1 * stand_all(5,ij) + par_rhof2
  Nold = stand_all(17,ij)

    !!!!update kRein and cR
   !!!!update par_kRein as a function of sitetype if parameters (param(50>-999.))) are are provided
   if(param(50)>-999.d0) call linearUpdateParam(param(50:51),stand(3),par_kRein)
   !!!!update par_cR as a function of sitetype if parameters (param(52>-999.))) are are provided
   if(param(52)>-999.d0) call linearUpdateParam(param(52:53),stand(3),par_cR)


  if(thinningType == 1. .or. thinningType == 2.) then
    ! N = number of trees in the current layer after thinning
    N = (stand_all(17,ij)/Ntot) * dens_thd
    H = stand_all(11,ij)
    D = stand_all(12,ij)
    BA = N*pi*(D/2./100.)**2.
  else if(thinningType == 3.) then
    BA_tot = BA_thd
    BA = BAr(ij) * BA_thd
        if(thinInt > 0.) then
      H = stand_all(11,ij) * 0.9
      D = stand_all(12,ij) * 0.9
    else
      if(par_sarShp==1.) then
        H = stand_all(11,ij) *  (1.2147-0.2086 * min(1.,(BA/ stand_all(13,ij))))
        D = stand_all(12,ij) * (1.2192 -0.2173 * min(1.,(BA/ stand_all(13,ij))))
      else
        H = stand_all(11,ij) *  (1.07386 -0.06553 * min(1.,(BA/ stand_all(13,ij))))
        D = stand_all(12,ij) * (1.1779 -0.1379 * min(1.,(BA/ stand_all(13,ij))))
      endif
    endif
    N = BA/(pi*((D/2./100.)**2.))
  endif

  stand_all(13,ij) = BA
    Nthd = max(0.,(Nold - N))
    Hc = min(stand_all(14,ij),0.9*H)
  if(siteInfo(1)==719400.) then
  endif

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
  age_factor = (1. - (1. - par_fAa)/ (1. + exp((par_fAb - h)/par_fAc)))/par_fAa
  par_alfar = modOut(year,3,ij,2) * age_factor
    ! if (sitetype <= 1.) then
     ! par_alfar = par_alfar1 * age_factor
    ! else if (sitetype==2.) then
     ! par_alfar = par_alfar2 * age_factor
    ! else if (sitetype==3.) then
     ! par_alfar = par_alfar3 * age_factor
    ! else if (sitetype==4.) then
     ! par_alfar = par_alfar4 * age_factor
    ! else
     ! par_alfar = par_alfar5 * age_factor
    ! end if

     ! Update dependent variables
   gammaC = par_cR/stand_all(36,ij)
   par_rhof = par_rhof1 * ETS + par_rhof2
    par_rhor = par_alfar * par_rhof
     betab = hb/Lc
     beta0 = par_beta0 * age_factor
     beta1 = (beta0 + betab + par_betas) !!newX
     beta2 = 1. - betab - par_betas     !!newX
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
        W_crs = par_rhow * beta0 * A * H * N !W_stem * (beta0 - 1.)  !coarse root biomass
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
  species = int(max(1.,stand_all(4,ij)))
if(pCrobas(2,species)>0.) energyWood(year,ij,1) = energyWood(year,ij,2) / pCrobas(2,species)

     S_branch = max(0.,stand_all(28,ij) + (stand_all(24,ij) - W_branch) * (1-energyRatio) +&
    stand_all(51,ij) - Wdb + &
    (0.3 * (1-energyRatio)+0.7) * (stand_all(32,ij) - W_croot) *0.83 + &
    (stand_all(31,ij)-W_stem)*(1-harvRatio)*(1-energyRatio))
   S_wood = max(0.,stand_all(29,ij)+(0.3 * (1-energyRatio)+0.7) * (stand_all(32,ij) - W_croot) *0.17)
  else
     S_branch = max(0.,stand_all(28,ij) + stand_all(24,ij) - W_branch + stand_all(51,ij) - Wdb + &
    (stand_all(32,ij) - W_croot) * 0.83+ (stand_all(31,ij) - W_stem) * (1-harvRatio))
     S_wood = max(0.,stand_all(29,ij)  + (stand_all(32,ij) - W_croot) * 0.17)
  endif
  !energyCut

    outt(11,ij,2) = STAND_tot(11)
    outt(12,ij,2) = STAND_tot(12)
    outt(13,ij,2) = outt(13,ij,2) + STAND_tot(13) - BA
    outt(14,ij,2) = STAND_tot(14)
    outt(15,ij,2) = STAND_tot(15)
    outt(16,ij,2) = outt(16,ij,2) + STAND_tot(16)
    outt(17,ij,2) = outt(17,ij,2) + Nthd
    outt(18:23,ij,2) = -999.
    outt(24,ij,2) = outt(24,ij,2) + STAND_tot(24) - W_branch
    outt(25,ij,2) = outt(25,ij,2) + STAND_tot(25) - W_froot
    outt(26:29,ij,2) = -999.
    outt(30,ij,2) = outt(30,ij,2) + STAND_tot(30) - V
    outt(31,ij,2) = outt(31,ij,2) + STAND_tot(31) - W_stem
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

 !calculate reneike and random mortality
 include 'mortalityCalc.h'
 !!model disturbances
 if(disturbance_wind) then
   include 'disturbanceCalc.h'
 endif

!add dead trees from disturbances
STAND_all(42,:) = STAND_all(42,:) + VmortDist

outt(:,:,1) = STAND_all

modOut((year+1),2,:,:) = outt(2,:,:)
modOut((year+1),4,:,:) = outt(4,:,:) !update species
modOut((year+1),7,:,:) = outt(7,:,:)
modOut((year+1),9:nVar,:,:) = outt(9:nVar,:,:)

!!!!calculate bark beetle disturbance
   include 'SBB_dist_Calc.h'

 if(oldLayer==1) then
  modOut((year+1),:,nLayers,:) = outt(:,nLayers,:)
 endif
!!!!run Yasso
 if(yassoRun==1.) then
  do ijj = 1, nLayers
   Lst(ijj) = outt(29,ijj,1)
   Lb(ijj) =  outt(28,ijj,1)
   Lf(ijj) = outt(26,ijj,1)+outt(27,ijj,1)

   species = int(max(1.,modOut((year+1),4,ijj,1)))
   call compAWENH(Lf(ijj),folAWENH(ijj,:),pAWEN(1:4,species))   !!!awen partitioning foliage
   if(GVrun==1 .and. ijj==1) then
    folAWENH(ijj,1:4) = folAWENH(ijj,1:4) + AWENgv       !!!add AWEN gv to 1st layer
   endif

   !!!ECMmodelling.
   !add W for all layer to W folAWENH(ijj,2) = folAWENH(ijj,2) + exud(ijj) !!!ECMmodelling.  exud(ijj)=0 if ECMmod= 0
   folAWENH(ijj,2) = folAWENH(ijj,2) + exud(ijj) !!!ECMmodelling

   call compAWENH(Lb(ijj),fbAWENH(ijj,:),pAWEN(5:8,species))   !!!awen partitioning branches
   call compAWENH(Lst(ijj),stAWENH(ijj,:),pAWEN(9:12,species))         !!!awen partitioning stems

   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,1,ijj),stAWENH(ijj,:),litterSize(1,species), &
  leac,soilC((year+1),:,1,ijj),0.)
   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,2,ijj),fbAWENH(ijj,:),litterSize(2,species), &
  leac,soilC((year+1),:,2,ijj),0.)
   call mod5c(pYasso,t,weatherYasso(year,:),soilC((year),:,3,ijj),folAWENH(ijj,:),litterSize(3,species), &
  leac,soilC((year+1),:,3,ijj),0.)
  enddo

  soilCtot(year+1) = sum(soilC(year+1,:,:,:))
 endif !end yassoRun if

!!!fire disturbance calculations
 ! if(fireDistFlag)
  Cpool_litter_wood =  sum(soilC((year+1),1:4,1,:)) + sum(soilC((year+1),1:4,2,:))
  Cpool_litter_green = sum(soilC((year+1),1:4,3,:)) * sum(outt(26,:,1))/sum(outt(26,:,1)+outt(27,:,1))
  livegrass = 0.!GVout(year,4)
  soil_moisture(:) = ((dailySW/pPRELES(1))-pPRELES(3))/(pPRELES(2)-pPRELES(3)) !relative extractable soil water
  Tmin = weatherPRELES(year,:,2) - 3.6
  Tmax = weatherPRELES(year,:,2) + 3.7
  FDI(:) = 0.
  call fireDist(Cpool_litter_wood,Cpool_litter_green,livegrass,soil_moisture, &
	weatherPRELES(year,:,2),NI((1+((year-1)*365)):(365*year)),weatherPRELES(year,:,4),&
 FDI,n_fire_year)
  modOut((year+1),47,:,2) = 0.
  modOut((year+1),47,1,2) = n_fire_year !maxval(FDI)
 ! endif

enddo !end year loop

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
  enddo
enddo

!compute fluxes in g C m−2 y−1
 modOut(:,44,:,1) = modOut(:,44,:,1)*1000. !*1000 coverts units to g C m−2 y−1
 modOut(:,9,:,1) = modOut(:,9,:,1)*1000.    !*1000 coverts units to g C m−2 y−1
 modOut(:,18,:,1) = modOut(:,18,:,1)*1000.    !*1000 coverts units to g C m−2 y−1

  modOut(2:(nYears+1),45,:,1) = modOut(2:(nYears+1),45,:,1) + & !! this includes priming (P_RT) calculated earlier otherwise is 0.
    modOut(1:(nYears),39,:,1)/10. - modOut(2:(nYears+1),39,:,1)/10. + &  !/10 coverts units to g C m−2 y−1
    modOut(2:(nYears+1),26,:,1)/10. + modOut(2:(nYears+1),27,:,1)/10. + &
    modOut(2:(nYears+1),28,:,1)/10. + modOut(2:(nYears+1),29,:,1)/10.
  if(GVrun==1) modOut(2:(nYears+1),45,1,1) = modOut(2:(nYears+1),45,1,1) + GVout(:,2)/10.  !/10 coverts units to g C m−2 y−1

modOut(:,46,:,1) = modOut(:,44,:,1) - modOut(:,9,:,1) - modOut(:,45,:,1)

!!!!ground vegetation Add Npp ground vegetation to the NEE first layer
!!!calculate state of GV at the last year
if(GVrun==1) then
 stand_all = modOut((nYears+1),:,:,1)
 do ij = 1, nLayers
  if(stand_all(4,ij)==0.) stand_all(4,ij)=1.
 enddo
 call Ffotos2(stand_all,nLayers,nSpec,pCrobas,&
  nVar,nPar,MeanLight,coeff,fAPARtrees)
call fAPARgv(fAPARtrees, ETSmean, siteInfo(3), lastGVout(1), lastGVout(2), &
         sum(P0yX(:,1))/nYears, AWENgv,lastGVout(4)) !reduced input output
     lastGVout(3) = prelesOut(1) * GVout(year,1)/fAPARsite!  this can be improved running the model for ground vegetation if layerPRELES==1
  if(nYears > 1) then
   GVout(1:(nYears-1),5) = GVout(2:(nYears),4)/10.d0 - GVout(1:(nYears-1),4)/10.d0 + GVout(1:(nYears-1),2)/10.d0
   GVout(nYears,5) = lastGVout(4)/10.d0 - GVout((nYears),4)/10.d0 + GVout((nYears),2)/10.d0
  ! GVout(1:(nYears-1),4) = GVout(2:(nYears),4)
  ! GVout(nYears,4) = lastGVout(4)
 else  !!!when nYears ==1 in the region multi prebas

  GVout(nYears,5) = (lastGVout(4) - GVout((nYears),4) + GVout((nYears),2))/10.
  ! GVout(nYears,4) = lastGVout(4)
  endif
  modOut(2:(nYears+1),46,1,1) = modOut(2:(nYears+1),46,1,1) + GVout(:,5)

endif
!!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
 do year = 2,(nYears +1)
  do ij = 1,nLayers
   D = modOut((year-1),12,ij,1)
   Vmort = modOut(year,42,ij,1)
    species = int(max(1.,modOut(year,4,ij,1)))
  if(Vmort>0. .and. D>pCrobas(48,species))then
     modOut(year,8,ij,1) = Vmort + modOut(year,8,ij,1)
   do i=1,(nYears+1-year)
    perVmort = exp(-exp(pCrobas(35,species) + pCrobas(36,species)*i +  &
                 pCrobas(37,species)*D + pCrobas(44,species)))
    if(perVmort > pCrobas(49,species)) then
       modOut((year+i),8,ij,1) = modOut((year+i),8,ij,1) + Vmort * perVmort
      endif
   enddo
   endif
  enddo
 enddo

 output = modOut(2:(nYears+1),:,:,:)
 output(:,5:6,:,:) = modOut(1:(nYears),5:6,:,:)
 soilCinOut = soilC(2:(nYears+1),:,:,:)
 soilCtotInOut = soilCtot(2:(nYears+1))
 output(:,1:2,1,2) = thinClx

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
