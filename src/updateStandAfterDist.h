!update stand after a disturbance

 ! 
 ! !!!!!check litterfall!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   if(sum(BAdist) > 0.) then !!! ADD CONDITION for occurence of wind disturbance wdimp x3

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

    outt(30,ij,2) = outt(30,ij,2) + max((Vold-V)*pHarvTrees,0.)
pHarvTrees = 0
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

