!Disturbance module

!!!! WIND RISK CALCULATIONS !!!!
! switched on by disturbanceON (in B_prebas), that in turn is switched on if there's siteInfoDist input in e.g. TransectRuns
! currently 'plugged in':
!   siteInfoDist: wspeed, tsincethin, soiltype, shallowsoil
!   dom layer spec, h
!   sitetype, tempsum
!   NOTE: tsincethin as of now static, controlled by siteInfoDist inputs
!         counter needs to be set to 0 in all thinning events, +1 if there are none
!   NOTE: openedge (at least one neighbouring stand <5m) currently not implemented
!         postponed, requires dynamic spacially explicit info (e.g. via lookup table)

! PREPARING WIND RISK INPUTS
! dev: outDist(,X) 1 = dom layer #; 2 dom layer spec; 3 dom layer h, 4 sitetype, 5 ETS; 6-10 wind risks
outDist(year,1) = 1 !dominant layer
outDist(year,2) = STAND_all(4,1) ! species, layer 1
outDist(year,3) = STAND_all(11,1) !h, layer 1
outDist(year,4) = STAND_all(3,1) ! sitetype
outDist(year,5) = STAND_all(5,1) ! ETS/tempsum

! layer-level data: spec & h of dominant layer
IF(nLayers>1) THEN !if there's more than one layer
  do i = 2, nLayers !loop through them
     if (STAND_all(11,i) > outDist(year,3)) then !higher than previous ones
        outDist(year,1) = i ! set new dom layer dim
        outDist(year,2) = STAND_all(4,i) !set new spec
        outDist(year,3) = STAND_all(11,i) !set new h
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
call windrisk(siteInfoDist, INT(outDist(year,2)), outDist(year,3), 0, outDist(year,4), outDist(year,5), &
  wrisk5dd1,wrisk5dd2,wrisk5dd3,wrisk0,wrisk5,wrisk)

!assigning risks
outDist(year,6) = siteInfoDist(2) ! tsincethin, developing
outDist(year,7) = wrisk5dd2 !5a, damage density class 2
outDist(year,8) = wrisk5dd3 !5a, damage density class 3
outDist(year,9) = wrisk5 !5a, frequency weighted average
outDist(year,10) = wrisk !1a

!!!! END WIND RISK CALCULATIONS

!endif




if(.false.) then !if XX everything is switch off for the moment

!!!!!check litterfall!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
if (disturbanceON) then 
  ! BAmort = 0.d0
  ! pMort = 0.2d0
  perBAmort = 0.0d0
  call pDistTest(ETS,1200.0d0,pMort) !!!calculate probability of fire to occur in this stand
  call random_number(randX)
  if(randX < pMort)  call intTest(pMort,perBAmort) !!!calculate the intensity of possible fire of fire to occur in this stand
  
!! calculate probability of the disturbance to occur and
!! the intensity of the disturbance
	
  if(perBAmort > 0.0d0) then

  BA_tot = sum(STAND_all(13,:))
  BAr = STAND_all(13,:)/BA_tot

	  ! perBAmort = 0.1
      ! write(1,*) "disturbance", year, pMort, perBAmort
	 
	do ij = 1 , nLayers 		!loop Species
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

	if (year > maxYearSite) then
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

	! Mortality - use random model from siilipehto et al. 2020
		 Vold = stand(30)
		 Nold = stand(17)
		 BAmort = perBAmort * BA
		if(BAmort > 0.) then !check if mortality occurs
		 dN = -Nold * (BAmort/(BA/BAr(ij)))
		else
		 dN = 0.	
		endif

	!!!update variables
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
		S_branch = max(0.,(W_branch + W_croot*0.83 + Wdb) * min(1.,-dN*step/Nold) )
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
	  STAND(26) = S_fol + STAND(26)
	  STAND(27) = S_fr + STAND(27)
	  STAND(28) = S_branch + STAND(28)
	  STAND(29) = S_wood + STAND(29)
	  STAND(31) = W_stem
	  STAND(32) = W_croot
	  STAND(42) = Vold - V + STAND(42)!* min(1.,-dN*step/Nold)
	  STAND(47) = W_wsap
	  STAND(48) = W_c
	  STAND(49) = W_s
	  STAND(50) = Wsh
	  STAND(53) = W_bh
	  STAND(54) = W_crh
	  STAND(51) = Wdb

		  endif
	  
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
	 endif
	endif

	 STAND_all(:,ij)=STAND
	end do !!!!!!!end loop layers
 endif
	
endif !if disturbanceON
! endif
endif !end if XX switch off the modules