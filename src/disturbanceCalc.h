!Disturbance module

!!!! WIND RISK CALCULATIONS !!!!
! switched on by disturbanceON (in B_prebas), that in turn is switched on if there's siteInfoDist input in e.g. TransectRuns
! currently 'plugged in':
!   siteInfoDist: wspeed, tsincethin, soiltype, shallowsoil
!   dom layer spec, h
!   sitetype, tempsum
!   tsincethin counter reset by manual & tapio thinnings (tested), compensation thinnings (untested) 
!   NOTE: openedge (at least one neighbouring stand <5m) currently not implemented
!         postponed, requires dynamic spacially explicit info (e.g. via lookup table)

! PREPARING WIND RISK INPUTS
! dev: outDist(,X) 1 = dom layer #; 2 dom layer spec; 3 dom layer h, 4 sitetype, 5 ETS; 6-10 wind risks
wdistproc(1) = 1. !dominant layer
wdistproc(2) = STAND_all(4,1) ! species, layer 1
wdistproc(3) = STAND_all(11,1) !h, layer 1

!outDist(year, 10) = STAND_all(11,2)
! layer-level data: spec & h of dominant layer
IF(nLayers>1) THEN !if there's more than one layer
  do i = 2, nLayers !loop through them
     if(STAND_all(11,i) > wdistproc(3)) then !higher than previous ones
        wdistproc(1) = i ! set new dom layer dim
        wdistproc(2) = STAND_all(4,i) !set new spec
        wdistproc(3) = STAND_all(11,i) !set new h
     end if
  end do
end if

! setting wrisks to 0 (subroutine inout), to be simplified
wrisk5dd1 = 0
wrisk5dd2 = 0
wrisk5dd3 = 0
wrisk0 = 0
wrisk5 = 0
wrisk = 0





! WINDRISK SUBROUTINE
!call windrisk(siteInfoDist, spec, h, openedge, sitetype, tsum, &
!  wrisk5dd1, wrisk5dd2, wrisk5dd3, wrisk0, wrisk5, wrisk)
call windrisk(siteInfoDist, INT(wdistproc(2)), wdistproc(3), 0, STAND_all(3,1), STAND_all(5,1), INT(siteInfoDist(2)), &
  wrisk5dd1,wrisk5dd2,wrisk5dd3,wrisk0,wrisk5,wrisk)

!assigning risks

outDist(year,1) = wdistproc(2) ! dom spec
outDist(year,2) = siteInfoDist(2) !tsincethin
outDist(year,3) = wrisk !1a



!!!! END WIND RISK CALCULATIONS

!!! DRAFT FOR WIND DISTURBANCE IMPACT MODELLING
! basic idea: 3-step sampling
!    - first sample with above risk as probability if a wind disturbance occurs
!    - second: sample severity class (for now: based on 2001 post-storm inventory shares, approx 1=81, 2=14, 3=5%)
!    - then allocate share of damaged V based on list of damaged Vs of severity classes (post-storm inventory)
! current status: damaged V in outDist[,,9], needs to be tested. further implementation via thinmat or so do be discussed

! step 1: wind disturbance 0/1 based on wind risk
call random_number(rndm)
if(rndm <= wrisk) outDist(year,4) = 1 !wind disturbance occurs/ set severity class to 1
if(rndm > wrisk) outDist(year,4) = 0 !... or doesn't.

!step 2: severity class (for now with shares/probabilities from post-storm inventory )
if (outDist(year,4)==1) then 
  call random_number(rndm) ! leave sevclass at 1 or increase based on sampling
  !outDist(year,8) = 1 !set sevclass to 1 
  if(rndm <= 0.13888889) outDist(year,4) = 2 
  if(rndm <= 0.05555556) outDist(year,4) = 3 !these are now sampled even if there's no disturbance occuring; keep for dev purposes, only calculate when dist occurres in final version
endif

  ! step 3: sample from severity class-specific set of relative disturbed volumes
   if (outDist(year,4)==1) then ! sevclass 1
     call random_number(rndm)
     sevclasslength = 87 !n of plots with sc 1 in ps inventory data
     distvloc = FLOOR(sevclasslength*rndm+1) ! sample one of these
     sc1vols =  (/ 0.047246314, 0.067229849, 0.169719737, 0.318203784, 0.018104818, 0.104955687, 0.032123615, 0.092026088, &
       0.013472795, 0.166694679, 0.195763598, 0.045904633, 0.030510592, 0.257592283, 0.055838402, 0.091561805, 0.085415300, &
       0.045814558, 0.036180144, 0.019006098, 0.040064294, 0.071564091, 0.010477727, 0.019651214, 0.175753183, 0.208317710, &
       0.009252750, 0.082452182, 0.031980969, 0.087094521, 0.021563759, 0.091875943, 0.075276931, 0.057730433, 0.030528242, &
       0.113322118, 0.062922399, 0.220426425, 0.026159837, 0.033814844, 0.037818739, 0.102943860, 0.112303663, 0.095156499, &
       0.054769579, 0.135111101, 0.026199722, 0.012797435, 0.034510249, 0.041657274, 0.069087262, 0.117174984, 0.045324950, &
       0.073308835, 0.021494620, 0.034361839, 0.045795929, 0.199031553, 0.014349513, 0.035105534, 0.108747020, 0.077563323, &
       0.017255370, 0.061953846, 0.208826663, 0.429494553, 0.025722064, 0.007571658, 0.024630056, 0.314765240, 0.024407173, &
       0.027934229, 0.012025332, 0.024008892, 0.028082671, 0.043077586, 0.015088951, 0.069155659, 0.044578726, 0.037450261, &
       0.003549555, 0.030742784, 0.114136973, 0.012353190, 0.039182845, 0.061220194, 0.032197320 /)
     wdistproc(4) = sc1vols(distvloc)
 
   else if (outDist(year,4)==2) then ! sevclass 2
     call random_number(rndm)
     sevclasslength = 15 !n of plots with sc 2 in ps inventory data
   distvloc = FLOOR(sevclasslength*rndm+1) ! sample one of these
     sc2vols =  (/ 0.23033922, 0.37008936, 0.09131254, 0.25163284, 0.29842135, 0.33619803, 0.13216089, 0.43466740, 0.20183134, &
       0.08506387, 0.07892917, 0.04246911, 0.24749787, 0.02704628, 0.02121135 /)
     wdistproc(4) = sc2vols(distvloc)
 
   else if (outDist(year,4)==3) then
     call random_number(rndm)
     sevclasslength = 6 !n of plots with sc 3 in ps inventory data
     distvloc = FLOOR(sevclasslength*rndm+1) ! sample one of these
     sc3vols =  (/ 0.09466817, 0.76296213, 1.0, 0.82065451, 0.21933642, 0.33760203 /)
     wdistproc(4) = sc3vols(distvloc)
   endif


!!! END WIND IMPACT CALCULATIONS !!!!

!!! DISTRIBUTE SHARE OF VOLUME DISTURBED TO LAYERS !!!
! idea: - calculate layer-level risks for dominant layer + those with H>(domh-3m) (for now)
!       - distribute share across layers according to ratios of wrisks

! quick & dirty: calculate wind risk for all layers, delete all that have h < (hdom-3)

wriskLayers(:,:) = 0
IF(outDist(year,4) >0 .AND. nLayers>1) THEN !if there's a wind disturbance and more than one layer
  do i = 1, nLayers !loop through them
     if (STAND_all(11,i) > (wdistproc(3)-5)) then !within 5 m of highest layer (which the total wind risk is based on) !!! 
     call windrisk(siteInfoDist, INT(STAND_all(4,i)), STAND_all(11,i), 0, STAND_all(3,1), STAND_all(5,1),  INT(siteInfoDist(2)), & !calculate layer wind risk
       wrisk5dd1,wrisk5dd2,wrisk5dd3,wrisk0,wrisk5,wrisk)
       wriskLayers(i, 1) = wrisk  
     end if
  end do
end if


BA_tot = sum(STAND_all(13,:))
V_tot = sum(STAND_all(30,:))


vdam = wdistproc(4)*V_tot

if(outDist(year, 4)>0) then 
  outDist(year, 5) = vdam
  outDist(year, 6) = wdistproc(4)
endif

wriskLayers(:, 2) = STAND_all(30,:)*wriskLayers(:,1) !weighing factor for vol 'at risk': if layer-level risk and volume are equal across layers, all would receive the same amount of damage; otherwise weighed by risk AND volume share
wriskLayers(:, 3) = wriskLayers(:, 2) / sum(wriskLayers(:, 2)) !shares of disturbwd volumes
!wriskLayers(:, 4) = wriskLayers(:,2)*wriskLayers(:,3) !share of potentially affected V !! attention: doesn't add up to 1, this is too simple
! better: 
wriskLayers(:, 4) = vdam * wriskLayers(:, 3)  ! plot-level damaged volume allocated to layers 
wriskLayers(:, 5) = STAND_all(30,:)/STAND_all(13,:)!V per ba
wriskLayers(:, 6) = wriskLayers(:, 4)/wriskLayers(:,5)! convert affected vol to affected ba


do layer = 1, nLayers
  !if(wriskLayers(layer, 6) /= wriskLayers(layer, 6)) !old version
  if(wriskLayers(layer, 6) /= wriskLayers(layer, 6)) wriskLayers(layer, 6) = 0. ! NaN check (div by 0) NaN is not equal to itself...
  !if(wriskLayers(layer, 6) /= wriskLayers(layer, 6)) outDist(year, 1) = 999. !checking
end do

!outDist(year, 1:nLayers) = wriskLayers(:, 6)


!write(1,*) wriskLayers(:,1), wriskLayers(:,2), wriskLayers(:,3), wriskLayers(:,4), wriskLayers(:,5), wriskLayers(:,6) !!to write wdistdev output

!!! END DISTRIBUTE SHARE OF VOLUME DISTURBED TO LAYERS !!!



! now implementing impact in Francesco's code below (search for wdimp)
! - set if in beginning to true (everything deactivated as of now)
! - include condition to activate layer loop in case of wind disturbance : max(windrisklayer(:,1)) > 0)




! salvlog/mgmtrect module
! additional parameters in siteInfoDist; for now, due to tab issue, hardcoded in siteInfoDisttemp
! siteInfoDisttemp(1:4) = siteInfoDist 
! siteInfoDisttemp(5) = 5. !salvlogthresh
! siteInfoDisttemp(6) = 1. !salvlogshare
! siteInfoDisttemp(7) = 0.9 !pHarvTrees
! siteInfoDisttemp(8) = 10. !mgmtreactthresh
! siteInfoDisttemp(9) = 1.  !mgmtreactshare
! siteInfoDisttemp(10) = 1.!sevdistccshare



if (outDist(year,4)>0.) then !in case of disturbance 
  pHarvTrees = 0.
  ! salvage logging
  if(vdam>=siteInfoDist(5)) then
    call random_number(rndm)
    if(rndm<=siteInfoDist(6)) then 
      pHarvTrees = siteInfoDist(7)! if sampled for salvlog set pHarvTrees
      outDist(year,7) = 1. !indicate salvage logging in output
    endif
  endif

  !mgmtract/prioritisation in siteOrder
  if(vdam>=siteInfoDist(8)) then
    call random_number(rndm)
    if(rndm<=siteInfoDist(7)) then
      outDist(year,8) = 1.! if sampled for mgmtreact
      pHarvTrees = siteInfoDist(7)! force salvlog as well (very unlikely to be omitted)
      outDist(year,7) = 1.
    endif
  endif

  ! cc in severely disturbed sites (putting in action to come...)
  if((wdistproc(4)>=0.5 .OR. outDist(year,4)==3) .AND. siteInfoDist(10)>0.) then !CC if sevclass = 3 or >50% of volume disturbed
    call random_number(rndm)
    if(rndm<=siteInfoDist(10)) then
       outDist(year,9) = 1. !indicate clearcut
       outDist(year,8) = 1. !mgmtreact = T in order to include cc harvests towards meeting harvlim (and not after it's been met if lower in siteorder...)
    endif
  endif  
endif ! end salvlog/mgmtrect module

  
 if(.TRUE.) then !if XX everything is switch off for the moment !wdimp x1
 ! 
 ! !!!!!check litterfall!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if (disturbanceON) then !x2
   ! BAmort = 0.d0
   ! pMort = 0.2d0
   BAdist = wriskLayers(:,6) !disturbed layer ba/layer ba

 !   call pDistTest(ETS,1200.0d0,pMort) !!!calculate probability of fire to occur in this stand
 !   call random_number(randX)
 !   if(randX < pMort)  call intTest(pMort,perBAmort) !!!calculate the intensity of possible fire of fire to occur in this stand
 ! 
 ! ! calculate probability of the disturbance to occur and
 ! ! the intensity of the disturbance
 !   perBAmort = 0. ! deactivate Francesco's randomised mortality, seems to be very active and reduces n < 1 over rotation
 !if(perBAmort > 0.0d0 .OR. maxval(wriskLayers(:,1)) > 0) then !!! ADD CONDITION for occurence of wind disturbance wdimp x3
  ! outDist(year,1) = sum(BAdist)
   if(maxval(BAdist) > 0.) then !!! ADD CONDITION for occurence of wind disturbance wdimp x3

   !BA_tot = sum(STAND_all(13,:))
   !BAr = STAND_all(13,:)/BA_tot

     ! perBAmort = 0.1
       ! write(1,*) "disturbance", year, pMort, perBAmort
   do ij = 1 , nLayers     !loop Species xl1
    
    BAmort = BAdist(ij)
    dN=0.d0
    STAND=STAND_all(:,ij)
    species = int(stand(4))
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

     !!!!update kRein and cR   
    !!!!update par_kRein as a function of sitetype if parameters (param(50>-999.))) are are provided 
    if(param(50)>-999.d0) call linearUpdateParam(param(50:51),stand(3),par_kRein) 
    !!!!update par_cR as a function of sitetype if parameters (param(52>-999.))) are are provided 
    if(param(52)>-999.d0) call linearUpdateParam(param(52:53),stand(3),par_cR) 
!activate 
   if (year > maxYearSite) then !x4
     STAND(2) = 0. !!newX
     STAND(8:21) = 0. !#!#
     STAND(23:37) = 0. !#!#
     STAND(42:44) = 0. !#!#
     STAND(47:nVar) = 0. !#!#
   else

   ! initialize site variables
     age = STAND(7)
     H = STAND(11)
     D = STAND(12)
     BA = STAND(13)
     Hc = STAND(14)
     N = BA/(pi*((D/2/100)**2))
     B = BA/N! * par_ops2
     A = stand(16)
     Lc = H - Hc
     hb = par_betab * Lc ** par_x
     Cw = 2. * hb
     STAND(15) = Cw
     ETS = STAND(5)
     Light = STAND(36)
     V = stand(30)
     mort = stand(41)
     par_sla = par_sla + (par_sla0 - par_sla) * Exp(-ln2 * (age / par_tsla) ** 2.)

    if (N>0.) then !x5

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

   ! Mortality - use random model from siilipehto et al. 2020
      Vold = stand(30)
      Nold = stand(17)
     ! BAmort = perBAmort * BA



    if(BAmort > 0.) then !check if mortality occurs UPDATE: only activated if there is no wind disturbance wdimp
!      if(BAmort(ij) > 0. .and. maxval(wriskLayers(:,1)) == 0) then !check if mortality occurs UPDATE: only activated if there is no wind disturbance wdimp
      !dN = -Nold * (BAmort/(BA/BAr(ij)))
     dN = -Nold * (BAmort/BA)
     
    ! elseif(maxval(wriskLayers(:,1)) > 0) then !wdimp define dN based on layer-level disturbed ba
    ! !  dN = -Nold * (BAmort/(BA/BAr(ij)))
    !   dN = -Nold * (wriskLayers(ij,6)/BA) !disturbed layer ba/layer ba
    else
      dN = 0.  
     endif

   !!!update variables
       N = max(0.0, N + step*dN)

       if (dN<0. .and. Nold>0.) then !x6
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
     S_branch = max(0.,(W_branch + W_croot*0.83 + Wdb) * min(1.,-dN*step/Nold) )!rem
    !S_wood = (W_croot*0.17 + W_stem) * min(1.,-dN*step/Nold)
    S_wood = (W_croot*0.17 + W_stem) * min(1.,-dN*step/Nold)*(1-pHarvTrees)!rem
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
     STAND(26) = S_fol + STAND(26)
     STAND(27) = S_fr + STAND(27)
     STAND(28) = S_branch + STAND(28)
     STAND(29) = S_wood + STAND(29)
     STAND(31) = W_stem
     STAND(32) = W_croot
    !STAND(42) = Vold - V + STAND(42)!* min(1.,-dN*step/Nold)
    STAND(42) = (Vold - V)*(1-pHarvTrees) + STAND(42)!* min(1.,-dN*step/Nold)
     STAND(47) = W_wsap
     STAND(48) = W_c
     STAND(49) = W_s
     STAND(50) = Wsh
     STAND(53) = W_bh
     STAND(54) = W_crh
     STAND(51) = Wdb
!!!
!! allocating salvage logging to current (regionPrebas harvlimit not met when site is checked or all mgmt switched off) or next year (some mgmt allowed / harvlimit exceeded)
if(ClCut == 0. .and. defaultThin == 0.) then ! either mgmt switched off entirely or blocked due to harvest limit being met
    outt(42,ij,2) = outt(30,ij,2) + max((Vold-V)*pHarvTrees,0.)*harvRatio !salvnext save salvlogged layer-level vol here to be included in next year's harvest limit in regionPrebas (harvRatio otherwise applied when going from ,,30,,2 to ,,37,,1)
elseif(ClCut > 0. .or. defaultThin > 0.) then 
    outt(30,ij,2) = outt(30,ij,2) + max((Vold-V)*pHarvTrees,0.)
endif
    !!!
    
  endif !x6

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
  endif !x5
endif !x4
!//activate
    STAND_all(:,ij)=STAND
    end do !!!!!!!end loop layers xl1
 endif !bamort>0... x3
!  ! 
 endif !if disturbanceON x2
!  ! ! endif
endif !end if XX switch off the modules x1
 ! 

