
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine bridging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine regionPrebas(siteOrder,HarvLim,minDharv,multiOut,nSites,areas,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,energyCuts,inDclct,inAclct,dailyPRELES,yassoRun,multiWood,&
		tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,GVrun,clearcuttingArea,compHarv)		!!energCuts


implicit none

integer, parameter :: nVar=54,npar=44!, nSp=3
real (kind=8), parameter :: harvRatio = 0.9, energyRatio = 0.7
integer, intent(in) :: nYears(nSites),nLayers(nSites),allSP
integer :: i,climID,ij,iz,ijj,ki,n,jj,az
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers,siteOrder(nSites,maxYears)
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5),minDharv
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP),compHarv(2)
 real (kind=8), intent(in) :: tapioPars(5,2,3,20),thdPer(nSites),limPer(nSites)
 real (kind=8), intent(in) :: tTapio(5,3,2,7), ftTapio(5,3,3,7)
 real (kind=8), intent(inout) :: siteInfo(nSites,10), areas(nSites),HarvLim(maxYears,2)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,9),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
 real (kind=8), intent(inout) :: energyCuts(nSites)	!!energCuts
 !!!ground vegetation
 integer, intent(in) :: gvRun			!!!ground vegetation
 real (kind=8), intent(inout) :: GVout(nSites,maxYears,3),clearcuttingArea(maxYears,2) !fAPAR_gv,litGV,photoGV,respGV			!!!ground vegetation
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,7,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: multiWood(nSites,maxYears,maxNlayers,2)!!energCuts
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(1,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers),wood(1,maxNlayers,2)
 real (kind=8) :: ClCutX, defaultThinX,maxState(nSites),check(maxYears), thinningX(maxThin,9)
 real (kind=8) :: energyWood, roundWood, energyCutX,thinFact	!!energCuts
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species

!!!!initialize run
! multiOut = 0.
yearX = 0.
soilC = soilCinOut
soilCtot = soilCtotInOut
multiWood = 0.
clearcuttingArea(:,2) = 0.
thinFact = compHarv(2)

    open(1,file="test1.txt")
    ! open(2,file="test2.txt")
    ! open(3,file="test3.txt")
 ! open(1,file="ftTapioREg.txt")
 ! open(2,file="tTapioReg.txt")
 ! write(1,*) ftTapio
 ! write(1,*) tTapio
 ! close(1)
 ! close(2)
!!inititialize A and biomasses
do i = 1,nSites
 do ijj = 1,nLayers(i)
	species = int(initVar(i,1,ijj))
		initVar(i,7,ijj) = pCrobas(38,species)/pCrobas(15,species) * (initVar(i,3,ijj) -&
			initVar(i,6,ijj))**pCrobas(11,species)!A = p_ksi/p_rhof * Lc^p_z
	call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),siteInfo(i,3),multiOut(i,1,:,ijj,1))
 enddo
enddo
do i = 1,nSites
 relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
enddo

do ij = 1,maxYears
 roundWood = 0.
 energyWood = 0.	!!energCuts
 do iz = 1,nSites
 	i=siteOrder(iz,ij)
! open(10,file="multiSite.txt")
 ! write(10,*) "years =",ij, "siteRun = ",iz
! close(10)
	ClCutX = ClCut(i)
	defaultThinX = defaultThin(i)
	energyCutX = energyCuts(i)		!!energCuts
	thinningX(:,:) = -999.
	az = 0

	if(ij > 1) then
	 soilC(i,ij,:,:,1:nLayers(i)) = soilC(i,(ij-1),:,:,1:nLayers(i))
	endif

!!!check if the limit has been exceeded if yes no havest (thinning or clearcut will be performed)
    ! write(1,*) clearcuttingArea(ij,:)
	if (clearcuttingArea(ij,1) > 0. .and. clearcuttingArea(ij,2) > clearcuttingArea(ij,1)) then !!!swithch off clear cuts if threshold area (clearcuttingArea(1)), has been reached
	 ClCutX = 0.
	endif
	if (HarvLim(ij,1) > 0. .and. roundWood >= HarvLim(ij,1)) then
	 ClCutX = 0.
	 defaultThinX = 0.
	endif
	if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	 energyCutX = 0.
	endif

!!!
	climID = siteInfo(i,2)
	if(ij==int(min(yearX(i),maxYears)))then
	 initClearcut(i,5) = int(min(initClearcut(i,5), initClearcut(i,5) + maxYears - yearX(i)))
	 yearX(i) = 0

	 do ijj = 1,nLayers(i)
	  species = int(multiOut(i,1,4,ijj,1))
	  initVar(i,1,ijj) = multiOut(i,1,4,ijj,1)
	  initVar(i,2,ijj) = initClearcut(i,5)
	  initVar(i,3,ijj) = initClearcut(i,1)
	  initVar(i,4,ijj) = initClearcut(i,2)
	  if(fixBAinitClarcut(i)==1) then
	   initVar(i,5,ijj) = initClearcut(i,3) * initCLcutRatio(i,ijj)
	  else
	   initVar(i,5,ijj) = initClearcut(i,3) * relBA(i,ijj)
      endif
	  initVar(i,6,ijj) = initClearcut(i,4)
	  ! initVar(i,8,ijj) = 0. !!newX
	  initVar(i,7,ijj) = max(0.,pCrobas(38,species)/pCrobas(15,species) * (initClearcut(i,1) -&
		initClearcut(i,4))**pCrobas(11,species))!A = p_ksi/p_rhof * Lc^p_z
	  do ki = 1,int(initClearcut(i,5)+1)
	   multiOut(i,int(ij-initClearcut(i,5)+ki-1),7,ijj,1) = ki !#!#
	  enddo !ki
	  call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),siteInfo(i,3),multiOut(i,(ij-1),:,ijj,1))
	 enddo !ijj
	endif

	do jj = 1, nThinning(i)
	 if(thinning(i,jj,1) == ij) then
	  az = az + 1
	  thinningX(az,:) = thinning(i,jj,:)
	  thinningX(az,1) = 1.
	 endif
	enddo
  ! if(ij==1) then
   ! write(*,*) sum(soilCinOut(i,ij,:,:,1:nLayers(i)))
  ! endif
	if(ij>2) then
		if(i==274) then 
 		 write(1,*) i,ij-1, multiOut(i,(ij-1),13,1:nLayers(i),1)
		endif
		output(1,1:7,1:nLayers(i),:) = multiOut(i,(ij-1),1:7,1:nLayers(i),:)
		output(1,9:nVar,1:nLayers(i),:) = multiOut(i,(ij-1),9:nVar,1:nLayers(i),:)
	else
		output(1,:,:,1) = multiOut(i,1,:,:,1)
		output(1,3:nVar,:,2) = multiOut(i,1,3:nVar,:,2)
	endif

	  call prebas(1,nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX(1:az,:),output(1,:,1:nLayers(i),:),az,maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,ij),P0y(climID,ij,:),&
		weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilC(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize,soilCtot(i,ij),&
		defaultThinX,ClCutX,energyCutX,inDclct(i,:),inAclct(i,:), & !!energCuts
		dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i),wood(1,1:nLayers(i),:),&
		tapioPars,thdPer(i),limPer(i),ftTapio,tTapio,GVout(i,ij,:),GVrun) !!energCuts
	
	! if clearcut occur initialize initVar and age
	if(sum(output(1,11,1:nLayers(i),1))==0 .and. yearX(i) == 0) then
	 if((maxYears-ij)<10) then
 		! if(initClearcut(i,5)<998.) then
			! Ainit = initClearcut(i,5)
		! else
			Ainit = nint(6 + 2*siteInfo(i,3) - 0.005*ETSy(climID,ij) + 2.25)
		! endif
	 else
 		! if(initClearcut(i,5)<998.) then
			! Ainit = initClearcut(i,5)
		! else
			Ainit = nint(6 + 2*siteInfo(i,3) - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
		! endif
	 endif
	 
	 clearcuttingArea(ij,2) = clearcuttingArea(ij,2) + areas(i) !calculate the clearcut area
 	! if(ij==35) then
	 ! write(1,*) clearcuttingArea(ij,2),i,iz,roundWood
	! endif
	 yearX(i) = Ainit + ij + 1
	 initClearcut(i,5) = Ainit
	 if(ij==1) then
	  relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
	 endif
	endif
	
	!!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
	do ijj = 1,nLayers(i)
	  if(output(1,8,ijj,1)>0.) then
	  multiOut(i,ij,8,ijj,1) = multiOut(i,ij,8,ijj,1) + output(1,8,ijj,1)
	  jj = int(output(1,4,ijj,1))
	    do ki = 1,(maxYears-ij)
			multiOut(i,(ki+ij),8,ijj,1) = multiOut(i,(ki+ij),8,ijj,1) + output(1,8,ijj,1) * &
				exp(-exp(pCrobas(34,jj) + pCrobas(35,jj)*ijj + pCrobas(36,jj)*output(1,12,ijj,1) + pCrobas(44,jj)))
		enddo
	  end if
	enddo

	multiWood(i,ij,1:nLayers(i),:) = wood(1,1:nLayers(i),:)
	multiOut(i,ij,1:7,1:nLayers(i),:) = output(1,1:7,1:nLayers(i),:)
	multiOut(i,ij,9:nVar,1:nLayers(i),:) = output(1,9:nVar,1:nLayers(i),:)
	! do ijj = 1,nLayers(i)
	  ! multiOut(i,ij,38,ijj,1) = sum(multiOut(i,1:ij,30,ijj,2)) + &
		! sum(multiOut(i,1:ij,42,ijj,1)) + multiOut(i,ij,30,ijj,1)
		
	  ! if(ij > 1.5) then
	! !compute gross growth
	   ! multiOut(i,ij,43,ijj,1) = multiOut(i,ij,38,ijj,1) - multiOut(i,(ij-1),38,ijj,1)
	  ! endif

	! enddo !ijj

	initVar(i,1,1:nLayers(i)) = output(1,4,1:nLayers(i),1)
	initVar(i,2,1:nLayers(i)) = output(1,7,1:nLayers(i),1)
	initVar(i,3:6,1:nLayers(i)) = output(1,11:14,1:nLayers(i),1)
	initVar(i,7,1:nLayers(i)) = output(1,16,1:nLayers(i),1)
	! initVar(i,8,1:nLayers(i)) = output(1,2,1:nLayers(i),1)  !!newX
	if(isnan(sum(output(1,37,1:nLayers(i),1)))) then
		roundWood = roundWood
		energyWood = energyWood
	else
		roundWood = roundWood + sum(output(1,37,1:nLayers(i),1))* areas(i)
		energyWood = energyWood + sum(wood(1,1:nLayers(i),1))* areas(i)   !!energCuts !!!we are looking at volumes	
	endif
 end do !iz i site loop

! write(10,*) "here3"
! write(2,*) roundWood,HarvLim(ij,1), ij
 !!! check if the harvest limit of the area has been reached otherwise clearcut the stands sorted by DBH 
 !or thin based on stand density index
if(roundWood < HarvLim(ij,1) .and. compHarv(1)>0.) then
 if(compHarv(1)==1.) then  !!!!clearcut to compensate harvest limits

  n = 0
  do while(n < nSites .and. roundWood < HarvLim(ij,1))		!!energCuts
   n = n + 1
   do i = 1, nSites
	if(ClCut(i) > 0. ) then
	 maxState(i) = maxval(multiOut(i,ij,12,1:nLayers(i),1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = int(siteInfo(siteX,2))
! write(3,*) ij,roundWood,HarvLim(ij,1),n, maxState(siteX), ClCut(siteX),areas(siteX) 
	if(maxState(siteX)>minDharv .and. ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
  ! close(10)
!!   !!clearcut!!
! write(1,*) "clearcutting", ij,maxState(siteX),minDharv
	 clearcuttingArea(ij,2) = clearcuttingArea(ij,2) + areas(siteX) !calculate the clearcut area
	   roundWood = roundWood + sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio)*areas(siteX) !!energCuts
	! write(1,*) roundWood,HarvLim(ij,1), ij,sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio),areas(siteX),n,nSites
	   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio
	   multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,31,1:nLayers(siteX),1)*harvRatio
	 multiOut(siteX,ij,2,1,2) = 2. !!!flag for clearcut compensation
     do ijj = 1, nLayers(siteX)
      multiOut(siteX,ij,6:nVar,ijj,2) = multiOut(siteX,ij,6:nVar,ijj,1)
      multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)+multiOut(siteX,ij,51,ijj,1)) * &
		(1-energyRatio) + multiOut(siteX,ij,28,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83 +&
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio)))
       multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,32,ijj,1)*0.17*(1-energyRatio)+ &
			multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  else
	   multiOut(siteX,ij,28,ijj,1) = max(0.,(multiOut(siteX,ij,24,ijj,1) + multiOut(siteX,ij,28,ijj,1) + &
		multiOut(siteX,ij,51,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83 + &
		multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)))
       multiOut(siteX,ij,29,ijj,1)=multiOut(siteX,ij,32,ijj,1)*0.17+multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  endif
!!energCuts
	  ! multiOut(siteX,ij,8,ijj,1) = 0.
	  multiOut(siteX,ij,10:17,ijj,1) = 0.
	  multiOut(siteX,ij,19:21,ijj,1) = 0.
	  multiOut(siteX,ij,2,ijj,1) = 0. !!newX
      multiOut(siteX,ij,23:36,ijj,1) = 0. !#!#
      multiOut(siteX,ij,43,ijj,1) = 0.
	  multiOut(siteX,ij,47:nVar,ijj,1) = 0.
    ! multiOut(siteX,ij,38,ijj,1) = sum(multiOut(siteX,1:ij,30,ijj,2)) + &
		! sum(multiOut(siteX,1:ij,42,ijj,1)) + multiOut(siteX,ij,30,ijj,1)
     enddo
	 if((maxYears-ij)<10) then
	  Ainit = nint(6 + 2*siteInfo(siteX,3) - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*siteInfo(siteX,3) - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 yearX(siteX) = Ainit + ij + 1
	 initClearcut(siteX,5) = Ainit
	 if(ij==1) then
	  relBA(siteX,1:nLayers(siteX)) = initVar(siteX,5,1:nLayers(siteX))/ &
		sum(initVar(siteX,5,1:nLayers(siteX)))
	 else
	  relBA(siteX,1:nLayers(siteX)) = multiOut(siteX,(ij-1),13,1:nLayers(siteX),1)/ &
		sum(multiOut(siteX,(ij-1),13,1:nLayers(siteX),1))
	 endif

    !initVar(siteX,1,1:nLayers(siteX)) = 0. !output(1,4,:,1)
     initVar(siteX,2,1:nLayers(siteX)) = 0.!output(1,7,:,1)
     initVar(siteX,3:7,1:nLayers(siteX)) = 0.!output(1,11:14,:,1)  !!newX
    endif !(maxState(i)>minDharv)
   enddo !end do while
   
 elseif(compHarv(1)==2.) then  !!!thin to compansate harvest limits
   !Perform thinning to compensate harvest levels
   !calculate SDI
  do i = 1, nSites
! !!start!! use stand density index to order the forests
   ! if(ClCut(i) > 0.) then
	! call calRein(multiOut(i,ij,:,:,1),nLayers(i),pCrobas(17,:),nVar,allSP,maxState(i))
   ! else
	! maxState(i) = 0.
   ! endif
! !!end!! use stand density index to order the forests

!!start!! use stand volume to order the forests
	if(ClCut(i) > 0. ) then
	 maxState(i) = sum(multiOut(i,ij,30,1:nLayers(i),1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
!!end!! use stand volume to order the forests
  enddo
  n = 0
  do while(n < nSites .and. roundWood < HarvLim(ij,1))		!!energCuts
   n = n + 1
   ops = maxloc(maxState)
   siteX = int(ops(1))
   maxState(siteX)=0.
	if(ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
	 !!!harvest thinFact 
	 roundWood = roundWood + sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio)* thinFact *areas(siteX) !!energCuts
     multiOut(siteX,ij,1,1,2) = 4. !!!flag for thinning compensation
	 multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,30,1:nLayers(siteX),1)*thinFact
	 multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,31,1:nLayers(siteX),1)*thinFact
     !update state of the forests
	 do ijj = 1, nLayers(siteX)
      multiOut(siteX,ij,9:10,ijj,2) = multiOut(siteX,ij,9:10,ijj,1) * thinFact
	  multiOut(siteX,ij,11:12,ijj,2) = multiOut(siteX,ij,11:12,ijj,1)
	  multiOut(siteX,ij,13,ijj,2) = multiOut(siteX,ij,13,ijj,1) * thinFact
	  multiOut(siteX,ij,14:16,ijj,2) = multiOut(siteX,ij,14:16,ijj,1)
	  multiOut(siteX,ij,17:34,ijj,2) = multiOut(siteX,ij,17:34,ijj,1) * thinFact
	  multiOut(siteX,ij,35,ijj,2) = multiOut(siteX,ij,35,ijj,1)
	  multiOut(siteX,ij,44,ijj,2) = multiOut(siteX,ij,44,ijj,1) * thinFact
	  
	  ! Litter foliage and branches
	  multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) * thinFact + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) * thinFact + multiOut(siteX,ij,27,ijj,1)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio * thinFact
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)+multiOut(siteX,ij,51,ijj,1)) * &
		(1-energyRatio)* thinFact + multiOut(siteX,ij,28,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83* thinFact +&
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio)** thinFact))
       multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,32,ijj,1)*0.17*(1-energyRatio)* thinFact+ &
			multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  else
	   multiOut(siteX,ij,28,ijj,1) = max(0.,(multiOut(siteX,ij,24,ijj,1)* thinFact + multiOut(siteX,ij,28,ijj,1) + &
		multiOut(siteX,ij,51,ijj,1)* thinFact + multiOut(siteX,ij,32,ijj,1)*0.83* thinFact + &
		multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)* thinFact))
       multiOut(siteX,ij,29,ijj,1)=multiOut(siteX,ij,32,ijj,1)*0.17* thinFact+multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  endif
!!energCuts
	  multiOut(siteX,ij,9:10,ijj,1) = multiOut(siteX,ij,9:10,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,13,ijj,1) = multiOut(siteX,ij,13,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,17:25,ijj,1) = multiOut(siteX,ij,17:25,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,30:33,ijj,1) = multiOut(siteX,ij,30:33,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,47:51,ijj,1) = multiOut(siteX,ij,47:51,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,53:nVar,ijj,1) = multiOut(siteX,ij,53:nVar,ijj,1)*(1-thinFact)
     enddo
    endif !(maxState(i)>minDharv)
   enddo !end do while
 elseif(compHarv(1)==3.) then  !!!thin to compansate harvest limits
   n = 0
  do while((n < nSites .and. roundWood < HarvLim(ij,1)) .and. &
	clearcuttingArea(ij,2) < clearcuttingArea(ij,1))		!!energCuts
   n = n + 1
   do i = 1, nSites
	if(ClCut(i) > 0. ) then
	 maxState(i) = maxval(multiOut(i,ij,12,1:nLayers(i),1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = int(siteInfo(siteX,2))
! write(3,*) ij,roundWood,HarvLim(ij,1),n, maxState(siteX), ClCut(siteX),areas(siteX) 
	if(maxState(siteX)>minDharv .and. ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
  ! close(10)
!!   !!clearcut!!
! write(1,*) "clearcutting", ij,maxState(siteX),minDharv
	 clearcuttingArea(ij,2) = clearcuttingArea(ij,2) + areas(siteX) !calculate the clearcut area
	   roundWood = roundWood + sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio)*areas(siteX) !!energCuts
	! write(1,*) roundWood,HarvLim(ij,1), ij,sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio),areas(siteX),n,nSites
	   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio
	   multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,31,1:nLayers(siteX),1)*harvRatio
	 multiOut(siteX,ij,2,1,2) = 2. !!!flag for clearcut compensation
     do ijj = 1, nLayers(siteX)
      multiOut(siteX,ij,6:nVar,ijj,2) = multiOut(siteX,ij,6:nVar,ijj,1)
      multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)+multiOut(siteX,ij,51,ijj,1)) * &
		(1-energyRatio) + multiOut(siteX,ij,28,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83 +&
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio)))
       multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,32,ijj,1)*0.17*(1-energyRatio)+ &
			multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  else
	   multiOut(siteX,ij,28,ijj,1) = max(0.,(multiOut(siteX,ij,24,ijj,1) + multiOut(siteX,ij,28,ijj,1) + &
		multiOut(siteX,ij,51,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83 + &
		multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)))
       multiOut(siteX,ij,29,ijj,1)=multiOut(siteX,ij,32,ijj,1)*0.17+multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  endif
!!energCuts
	  ! multiOut(siteX,ij,8,ijj,1) = 0.
	  multiOut(siteX,ij,10:17,ijj,1) = 0.
	  multiOut(siteX,ij,19:21,ijj,1) = 0.
	  multiOut(siteX,ij,2,ijj,1) = 0. !!newX
      multiOut(siteX,ij,23:36,ijj,1) = 0. !#!#
      multiOut(siteX,ij,43,ijj,1) = 0.
	  multiOut(siteX,ij,47:nVar,ijj,1) = 0.
    ! multiOut(siteX,ij,38,ijj,1) = sum(multiOut(siteX,1:ij,30,ijj,2)) + &
		! sum(multiOut(siteX,1:ij,42,ijj,1)) + multiOut(siteX,ij,30,ijj,1)
     enddo
	 if((maxYears-ij)<10) then
	  Ainit = nint(6 + 2*siteInfo(siteX,3) - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*siteInfo(siteX,3) - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 yearX(siteX) = Ainit + ij + 1
	 initClearcut(siteX,5) = Ainit
	 if(ij==1) then
	  relBA(siteX,1:nLayers(siteX)) = initVar(siteX,5,1:nLayers(siteX))/ &
		sum(initVar(siteX,5,1:nLayers(siteX)))
	 else
	  relBA(siteX,1:nLayers(siteX)) = multiOut(siteX,(ij-1),13,1:nLayers(siteX),1)/ &
		sum(multiOut(siteX,(ij-1),13,1:nLayers(siteX),1))
	 endif

    !initVar(siteX,1,1:nLayers(siteX)) = 0. !output(1,4,:,1)
     initVar(siteX,2,1:nLayers(siteX)) = 0.!output(1,7,:,1)
     initVar(siteX,3:7,1:nLayers(siteX)) = 0.!output(1,11:14,:,1)  !!newX
    endif !(maxState(i)>minDharv)
   enddo !end do while

  do i = 1, nSites
! !!start!! use stand density index to order the forests
   ! if(ClCut(i) > 0.) then
	! call calRein(multiOut(i,ij,:,:,1),nLayers(i),pCrobas(17,:),nVar,allSP,maxState(i))
   ! else
	! maxState(i) = 0.
   ! endif
! !!end!! use stand density index to order the forests

!!start!! use stand volume to order the forests
	if(ClCut(i) > 0. ) then
	 maxState(i) = sum(multiOut(i,ij,30,1:nLayers(i),1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
!!end!! use stand volume to order the forests

  enddo
  n = 0
  do while(n < nSites .and. roundWood < HarvLim(ij,1))		!!energCuts
   n = n + 1
   ops = maxloc(maxState)
   siteX = int(ops(1))
   maxState(siteX)=0.
	if(ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
	 !!!harvest thinFact 
	 roundWood = roundWood + sum(multiOut(siteX,ij,30,1:nLayers(siteX),1)*harvRatio)* thinFact *areas(siteX) !!energCuts
     multiOut(siteX,ij,1,1,2) = 4. !!!flag for thinning compensation
	 multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,30,1:nLayers(siteX),1)*thinFact
	 multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:nLayers(siteX),1) + &
			multiOut(siteX,ij,31,1:nLayers(siteX),1)*thinFact
     !update state of the forests
	 do ijj = 1, nLayers(siteX)
      multiOut(siteX,ij,9:10,ijj,2) = multiOut(siteX,ij,9:10,ijj,1) * thinFact
	  multiOut(siteX,ij,11:12,ijj,2) = multiOut(siteX,ij,11:12,ijj,1)
	  multiOut(siteX,ij,13,ijj,2) = multiOut(siteX,ij,13,ijj,1) * thinFact
	  multiOut(siteX,ij,14:16,ijj,2) = multiOut(siteX,ij,14:16,ijj,1)
	  multiOut(siteX,ij,17:34,ijj,2) = multiOut(siteX,ij,17:34,ijj,1) * thinFact
	  multiOut(siteX,ij,35,ijj,2) = multiOut(siteX,ij,35,ijj,1)
	  multiOut(siteX,ij,44,ijj,2) = multiOut(siteX,ij,44,ijj,1) * thinFact
	  
	  ! Litter foliage and branches
	  multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) * thinFact + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) * thinFact + multiOut(siteX,ij,27,ijj,1)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio * thinFact
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)+multiOut(siteX,ij,51,ijj,1)) * &
		(1-energyRatio)* thinFact + multiOut(siteX,ij,28,ijj,1) + multiOut(siteX,ij,32,ijj,1)*0.83* thinFact +&
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio)** thinFact))
       multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,32,ijj,1)*0.17*(1-energyRatio)* thinFact+ &
			multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  else
	   multiOut(siteX,ij,28,ijj,1) = max(0.,(multiOut(siteX,ij,24,ijj,1)* thinFact + multiOut(siteX,ij,28,ijj,1) + &
		multiOut(siteX,ij,51,ijj,1)* thinFact + multiOut(siteX,ij,32,ijj,1)*0.83* thinFact + &
		multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)* thinFact))
       multiOut(siteX,ij,29,ijj,1)=multiOut(siteX,ij,32,ijj,1)*0.17* thinFact+multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  endif
!!energCuts
	  multiOut(siteX,ij,9:10,ijj,1) = multiOut(siteX,ij,9:10,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,13,ijj,1) = multiOut(siteX,ij,13,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,17:25,ijj,1) = multiOut(siteX,ij,17:25,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,30:33,ijj,1) = multiOut(siteX,ij,30:33,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,47:51,ijj,1) = multiOut(siteX,ij,47:51,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,53:nVar,ijj,1) = multiOut(siteX,ij,53:nVar,ijj,1)*(1-thinFact)
     enddo
	 	if(multiOut(siteX,ij,1,1,1)==274.) then 
 		 write(1,*) i,ij, multiOut(i,(ij),13,1:nLayers(i),1)
		endif

    endif !(maxState(i)>minDharv)
   enddo !end do while
 
 
 endif 
endif !roundWood < HarvLim .and. HarvLim /= 0.

! write(10,*) "here4"
  !HarvLim(ij,1) = roundWood
  !HarvLim(ij,2) = energyWood
end do !end Year loop 
 do i = 1,nSites
  do ij = 1, maxYears 
    do ijj = 1,nLayers(i)
	  ! multiOut(i,ij,38,ijj,1) = sum(multiOut(i,1:ij,30,ijj,2)) + &
		! sum(multiOut(i,1:ij,42,ijj,1)) + multiOut(i,ij,30,ijj,1)
		
	  if(ij > 1.5) then
	!compute gross growth
	   multiOut(i,ij,43,ijj,1) = multiOut(i,ij,30,ijj,1) - multiOut(i,(ij-1),30,ijj,1) + & 
			multiOut(i,ij,37,ijj,1)/harvRatio + multiOut(i,ij,42,ijj,1)
	  endif

    enddo !ijj
  enddo
 enddo	
  close(1)
  ! close(2)
  ! close(3)
! write(10,*) "here5"
! close(10)
soilCinOut = soilC
soilCtotInOut = soilCtot

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


