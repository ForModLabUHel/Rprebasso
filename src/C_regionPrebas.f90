
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine bridging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine regionPrebas(siteOrder,HarvLim,minDharv,multiOut,nSites,areas,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,energyCuts,inDclct,inAclct,dailyPRELES,yassoRun,multiWood,&
		tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,GVrun,cuttingArea,compHarv,thinInt, &
		ageMitigScen, fertThin,flagFert,nYearsFert,oldLayer,mortMod)

implicit none

integer, parameter :: nVar=54,npar=47!, nSp=3
real (kind=8), parameter :: harvRatio = 0.9, energyRatio = 0.7
integer, intent(in) :: nYears(nSites),nLayers(nSites),allSP,oldLayer
integer :: i,climID,ij,iz,ijj,ki,n,jj,az
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers
integer, intent(inout) :: siteOrder(nSites,maxYears)
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5),minDharv,ageMitigScen
 integer, intent(in) :: DOY(365),etmodel,mortMod
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP)
!cuttingArea columns are clcutA target(1) simuation(2);tending target(3), sim(4);firstThin targ(5) sim(6)
 real (kind=8), intent(inout) :: compHarv(2),cuttingArea(maxYears,6) 
 real (kind=8), intent(in) :: tapioPars(5,2,3,20),thdPer(nSites),limPer(nSites)
 real (kind=8), intent(in) :: tTapio(5,3,2,7), ftTapio(5,3,3,7)
 real (kind=8), intent(inout) :: siteInfo(nSites,10), areas(nSites),HarvLim(maxYears,2)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,9),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
 real (kind=8), intent(in) :: thinInt(nSites) !site specific parameter that determines the thinning intensity; 
					!from below (thinInt>1) or above (thinInt<1);thinInt=999. uses the default value from tapio rules
 real (kind=8), intent(inout) :: energyCuts(nSites)	!!energCuts
 !!!ground vegetation
 integer, intent(in) :: gvRun			!!!ground vegetation
 real (kind=8), intent(inout) :: GVout(nSites,maxYears,5) !fAPAR_gv,litGV,photoGV,wGV			!!!ground vegetation
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
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species,layerX,domSp(1)
 real (kind=8) :: tTapioX(5,3,2,7), ftTapioX(5,3,3,7), Vmort, D,randX
 

!!! fertilization parameters
 integer, intent(inout) :: fertThin !!! flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3
 integer, intent(inout) :: flagFert(nSites) !!! flag that indicates if fertilization has already been applied along the rotation
 integer :: yearsFert !!actual number of years for fertilization (it depends if the thinning occur close to the end of the simulations)
 integer, intent(inout) :: nYearsFert !!number of years for which the fertilization is effective
 real(8) :: alfarFert(nYearsFert,maxNlayers),pDomRem, age(nSites), siteOrdX(nSites)

!!!!initialize run
! multiOut = 0.
yearX = 0.
soilC = soilCinOut
soilCtot = soilCtotInOut
multiWood = 0.
cuttingArea(:,2) = 0.
cuttingArea(:,4) = 0.
cuttingArea(:,6) = 0.
thinFact = compHarv(2)
tTapioX = tTapio
ftTapioX = ftTapio
multiOut(:,1,7,:,1) = initVar(:,2,:) !initialize age used in the mitigation scenario to select the sites to harvest
multiOut(:,1,4,:,1) = initVar(:,1,:) !initialize species 

    open(1,file="test1.txt")
    open(2,file="test2.txt")
    open(3,file="test3.txt")

! write(2,*) "compHarv",compHarv
!!inititialize A and biomasses
do i = 1,nSites
 do ijj = 1,nLayers(i)
	if(initVar(i,5,ijj) == 0.) then
		initVar(i,7,ijj) = 0. 
		multiOut(i,1,(/24,25,30,31,32,33,47,48,49,50,51,54/),ijj,1)=0.
	else		
		species = int(initVar(i,1,ijj))
		initVar(i,7,ijj) = pCrobas(38,species)/pCrobas(15,species) * (initVar(i,3,ijj) -&
			initVar(i,6,ijj))**pCrobas(11,species)!A = p_ksi/p_rhof * Lc^p_z
		call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),siteInfo(i,3),multiOut(i,1,:,ijj,1))
	endif
 enddo
 relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
enddo

do ij = 1,maxYears
 roundWood = 0.
 energyWood = 0.	!!energCuts
 

 if(ageMitigScen > 0.) then
  do i = 1,nSites
   if(oldLayer==1) then
    jj = max((nLayers(i)-1),1)
   else
    jj = nLayers(i)
   endif
   domSp = maxloc(multiOut(i,ij,13,1:jj,1))
   layerX = int(domSp(1))
   age(i) = multiOut(i,ij,7,layerX,1)
  enddo
  siteOrdX = real(siteOrder(:,ij),8)
  call changeOrder(siteOrdX,age, & 
					siteOrdX,nSites,ageMitigScen)
  siteOrder(:,ij) = int(siteOrdX)
 endif
 
 do iz = 1,nSites
 	i=siteOrder(iz,ij)
	ClCutX = ClCut(i)
	defaultThinX = defaultThin(i)
	energyCutX = energyCuts(i)		!!energCuts
	thinningX(:,:) = -999.
	az = 0

	if(ij > 1) then
	 soilC(i,ij,:,:,1:nLayers(i)) = soilC(i,(ij-1),:,:,1:nLayers(i))
	endif

!!!check if the limit has been exceeded if yes no harvest (thinning or clearcut will be performed)
    if (cuttingArea(ij,1) > 0. .and. cuttingArea(ij,2) > cuttingArea(ij,1)) then !!!swithch off clear cuts if threshold area (cuttingArea(1)), has been reached
	 ClCutX = 0.
	endif
	if (HarvLim(ij,1) > 0. .and. roundWood >= HarvLim(ij,1)) then
	 ClCutX = 0.
	 defaultThinX = 0.
	endif
	if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	 energyCutX = 0.
	endif

!!!check if the limit area for tendings has been exceeded if yes no tending havest 
	if (cuttingArea(ij,3) > 0. .and. cuttingArea(ij,4) > cuttingArea(ij,3)) then !!!swithch off tendings if threshold area (cuttingArea(3)), has been reached
	 tTapioX = tTapio * 1.e5
	else
	 tTapioX = tTapio
	endif
!!!check if the limit area for firstThin has been exceeded if yes no firstThin havest 
	if (cuttingArea(ij,5) > 0. .and. cuttingArea(ij,6) > cuttingArea(ij,5)) then !!!swithch off firstThin if threshold area (cuttingArea(5)), has been reached
	 ftTapioX = ftTapio * 1.e5
	else
	 ftTapioX = ftTapio
	endif

!!!
	climID = siteInfo(i,2)
	if(ij==int(yearX(i)))then
	!if(ij==int(min(yearX(i),maxYears)))then
	 ! initClearcut(i,5) = int(min(initClearcut(i,5), initClearcut(i,5) + maxYears - yearX(i)))
	 initClearcut(i,5) = int(initClearcut(i,5))
	 yearX(i) = 0

!if scenario = "oldLayer" do not consider the old layer
	 if(oldLayer==1) then
		jj=max((nLayers(i)-1),1)
	 else
		jj=nLayers(i)
	 endif

	 do ijj = 1,jj
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
	  call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),multiOut(i,ij,3,ijj,1),multiOut(i,(ij-1),:,ijj,1))
	 enddo !ijj
	endif

	do jj = 1, nThinning(i)
	 if(thinning(i,jj,1) == ij) then
	  az = az + 1
	  thinningX(az,:) = thinning(i,jj,:)
	  thinningX(az,1) = 1.
	 endif
	enddo

	if(ij>1) then
		if(oldLayer==1) output(1,3,:,:) = multiOut(i,(ij-1),3,:,:)
		output(1,1:7,1:nLayers(i),:) = multiOut(i,(ij-1),1:7,1:nLayers(i),:)
		output(1,9:nVar,1:nLayers(i),:) = multiOut(i,(ij-1),9:nVar,1:nLayers(i),:)
	else
		output(1,:,:,1) = multiOut(i,1,:,:,1)
		output(1,3:nVar,:,2) = multiOut(i,1,3:nVar,:,2)
	endif
 	! if(siteInfo(i,1)==411310.) write(1,*) ij,output(1,11,1:nLayers(i),1)
	! if(siteInfo(i,1)==35.) write(2,*) ij,output(1,11,1:nLayers(i),1)

	  call prebas(1,nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX(1:az,:),output(1,:,1:nLayers(i),:),az,maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,ij),P0y(climID,ij,:),&
		weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilC(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize,soilCtot(i,ij),&
		defaultThinX,ClCutX,energyCutX,inDclct(i,:),inAclct(i,:), & !!energCuts
		dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i),wood(1,1:nLayers(i),:),&
		tapioPars,thdPer(i),limPer(i),ftTapioX,tTapioX,GVout(i,ij,:),GVrun,thinInt(i), &
		fertThin,flagFert(i),nYearsFert,oldLayer,mortMod) !!energCuts
 	! if(siteInfo(i,1)==411310.) write(1,*) ij,output(1,11,1:nLayers(i),1)
	! if(siteInfo(i,1)==35.) write(2,*) ij,output(1,11,1:nLayers(i),1)
	!!!if oldLayer is active import siteType and alfar from the single site simulations simulations
	if(oldLayer==1 .and. output(1,3,nLayers(i),2)>0.) then
	 	 multiOut(i,ij:maxYears,3,nLayers(i),1) = output(1,3,nLayers(i),1)
		 multiOut(i,ij:maxYears,3,nLayers(i),2) = output(1,3,nLayers(i),2)
	endif	

if(compHarv(1)==-10. .or. compHarv(1)==2.) then
 if(siteInfo(i,1) == 454702.) write(3,*) "remaining 1", ij, output(1,11,1,1),output(1,13,1,1),&
	 output(1,37,1,1)
 if(siteInfo(i,1) == 454702.) write(3,*) "remaining 2", output(1,11,2,1), output(1,13,2,1), &
	 output(1,37,2,1)
 if(siteInfo(i,1) == 454702.) write(3,*) "remaining 3", output(1,11,3,1), output(1,13,3,1), &
	 output(1,37,3,1)
if(siteInfo(i,1) == 454702.) write(3,*) "remaining 4", output(1,11,4,1), output(1,13,4,1), &
	 output(1,37,4,1)
 if(siteInfo(i,1) == 454702.) write(3,*) "thinned 1", ij, output(1,11,1,2), output(1,13,1,2),&
	 output(1,30,1,2)
 if(siteInfo(i,1) == 454702.) write(3,*) "thinned 2", output(1,11,2,2), output(1,13,2,2), &
	 output(1,30,2,2)
 if(siteInfo(i,1) == 454702.) write(3,*) "thinned 3", output(1,11,3,2), output(1,13,3,2), &
	 output(1,30,3,2)
if(siteInfo(i,1) == 454702.) write(3,*) "thinned 4", output(1,11,4,2), output(1,13,4,2), &
	 output(1,30,4,2)
	  ! if(siteInfo(siteX,1) == 454702.) write(2,*) "thinned x", multiOut(siteX,ij,11,ijj,2), multiOut(siteX,ij,13,ijj,2), &
	  ! multiOut(siteX,ij,30,ijj,2)
endif
	
	!!!if fertilization at thinning is active,  increase siteType
	if(flagFert(i)==1 .and. fertThin>0) then 
		yearsFert = max(1,min(((nYears(i)) - ij-1),nYearsFert))
		multiOut(i,(ij+1):(ij+yearsFert),3,:,1) = siteInfo(i,3)-1.
		call calcAlfar(multiOut(i,ij,3,1:nLayers(i),:),initVar(i,1,1:nLayers(i)),pCrobas, &
				nLayers(i),alfarFert,allSP,nYearsFert,npar)
		multiOut(i,(ij+1):(ij+yearsFert),3,:,2) = alfarFert(1:yearsFert,:)
		flagFert(i)=2
	endif

	! if clearcut occur initialize initVar and age
!!calculate year of replanting after a clearcut
!if scenario = "oldLayer" do not consider the old layer
if(oldLayer==1) then
 jj=max((nLayers(i)-1),1)
else
 jj=nLayers(i)
endif
	if(sum(output(1,11,1:jj,1))==0 .and. yearX(i) == 0) then
	 if((maxYears-ij)<10) then
 		Ainit = nint(6 + 2*siteInfo(i,3) - 0.005*ETSy(climID,ij) + 2.25)
	 else
 		Ainit = nint(6 + 2*siteInfo(i,3) - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 !!!!update area of cuttings
	 cuttingArea(ij,2) = cuttingArea(ij,2) + areas(i) !calculate the clearcut area
	 yearX(i) = Ainit + ij + 1
	 initClearcut(i,5) = Ainit
	 if(ij==1) then
	  relBA(i,1:jj) = initVar(i,5,1:jj)/sum(initVar(i,5,1:jj))
	 endif
	endif
	
	multiWood(i,ij,1:nLayers(i),:) = wood(1,1:nLayers(i),:)
	multiOut(i,ij,1:2,1:nLayers(i),:) = output(1,1:2,1:nLayers(i),:)
	multiOut(i,ij,4:7,1:nLayers(i),:) = output(1,4:7,1:nLayers(i),:)
	multiOut(i,ij,9:nVar,1:nLayers(i),:) = output(1,9:nVar,1:nLayers(i),:)

	if(multiOut(i,ij,1,1,2) == 1.) then
	  cuttingArea(ij,4) = cuttingArea(ij,4) + areas(i)
	endif
	if(multiOut(i,ij,1,1,2) == 2.) then
	  cuttingArea(ij,6) = cuttingArea(ij,6) + areas(i)
	endif

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

 !!! check if the harvest limit of the area has been reached otherwise clearcut the stands sorted by DBH 
 !or thin based on stand density index
if(roundWood < HarvLim(ij,1) .and. compHarv(1)>0.5) then
! if(.false.) then
 if(compHarv(1)==1.) then  !!!!clearcut to compensate harvest limits
  n = 0
  do while(n < nSites .and. roundWood < HarvLim(ij,1))		!!energCuts
   n = n + 1
   !!!search for site with highest DBH
   do i = 1, nSites
if(oldLayer==1) then
 jj=max((nLayers(i)-1),1)
else
 jj=nLayers(i)
endif 
	if(ClCut(i) > 0. ) then
	 maxState(i) = maxval(multiOut(i,ij,12,1:jj,1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = int(siteInfo(siteX,2))
	if(maxState(siteX)>minDharv .and. ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
	 
	 !if fertilization at thinning is active reset flagFert
	 if(fertThin > 0) then
	  flagFert(siteX) = 0  
	 endif

!!if oldLayer scenario -> populate old layer with the dominant species and clearcut other layers
if(oldLayer==1) then
 jj=max((nLayers(siteX)-1),1)
 domSp = maxloc(multiOut(siteX,ij,13,1:jj,1))
 layerX = int(domSp(1))
  
  !!!!calculate percentage of trees remaining after clearcut(pDomRem)
  call random_number(randX)
   pDomRem =   (randX*5.+5.)/100.* &  !!randomly sample between 5 and 10 %
			sum(multiOut(siteX,ij,13,1:jj,1))/multiOut(siteX,ij,13,layerX,1)
   
	!update old layer
	multiOut(siteX,ij,:,nLayers(siteX),1) = multiOut(siteX,ij,:,layerX,1)
	multiOut(siteX,ij,42,nLayers(siteX),1) = 0.
	! multiOut(siteX,ij,3,nLayers(siteX),2) = multiOut(siteX,ij,3,layerX,2)
	multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),nLayers(siteX),1) = &
	 multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) * (pDomRem)
	multiOut(siteX,ij,24:34,nLayers(siteX),1) = multiOut(siteX,ij,24:34,layerX,1) * (pDomRem)
	multiOut(siteX,ij,47:51,nLayers(siteX),1) = multiOut(siteX,ij,47:51,layerX,1) * (pDomRem)
    !!!if oldLayer is active import siteType and alfar from the single site simulations simulations
	multiOut(siteX,ij:maxYears,3,nLayers(siteX),1) = multiOut(siteX,ij,3,layerX,1)
	multiOut(siteX,ij:maxYears,3,nLayers(siteX),2) = multiOut(siteX,ij,3,layerX,2)
	!update dominant layer
	multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) = &
 	 multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) * (1-pDomRem)
	multiOut(siteX,ij,24:34,layerX,1) = multiOut(siteX,ij,24:34,layerX,1) * (1-pDomRem)
	multiOut(siteX,ij,47:51,layerX,1) = multiOut(siteX,ij,47:51,layerX,1) * (1-pDomRem)

else
 jj=nLayers(siteX)
endif

!!   !!clearcut!!
	 cuttingArea(ij,2) = cuttingArea(ij,2) + areas(siteX) !calculate the clearcut area
	   roundWood = roundWood + sum(multiOut(siteX,ij,30,1:jj,1)*harvRatio)*areas(siteX) !!energCuts
	   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:jj,1) + &
			multiOut(siteX,ij,30,1:jj,1)*harvRatio
	   multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:jj,1) + &
			multiOut(siteX,ij,31,1:jj,1)*harvRatio
	 multiOut(siteX,ij,2,1,2) = 2. !!!flag for clearcut compensation
     do ijj = 1, jj
      multiOut(siteX,ij,6:23,ijj,2) = multiOut(siteX,ij,6:23,ijj,1)
      multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
      multiOut(siteX,ij,28:29,ijj,2) = multiOut(siteX,ij,28:29,ijj,1)
      multiOut(siteX,ij,35:nVar,ijj,2) = multiOut(siteX,ij,35:nVar,ijj,1)
	  !update biomasses and Volumes
	  multiOut(siteX,ij,24:25,ijj,2) = multiOut(siteX,ij,24:25,ijj,1) + &
					multiOut(siteX,ij,24:25,ijj,2)
	  multiOut(siteX,ij,30:34,ijj,2) = multiOut(siteX,ij,30:34,ijj,1) + &
					multiOut(siteX,ij,30:34,ijj,2)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)*(1-energyRatio) +   &
			multiOut(siteX,ij,51,ijj,1) + multiOut(siteX,ij,28,ijj,1)) + &
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio) + &
		 (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.83 ))
       multiOut(siteX,ij,29,ijj,1) = (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.17 + &
			multiOut(siteX,ij,29,ijj,1) 
			
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
	  relBA(siteX,1:jj) = initVar(siteX,5,1:jj)/ &
		sum(initVar(siteX,5,1:jj))
	 else
	  relBA(siteX,1:jj) = multiOut(siteX,(ij-1),13,1:jj,1)/ &
		sum(multiOut(siteX,(ij-1),13,1:jj,1))
	 endif

    !initVar(siteX,1,1:nLayers(siteX)) = 0. !output(1,4,:,1)
   	 initVar(siteX,1,1:nLayers(siteX)) = multiOut(siteX,ij,4,1:nLayers(siteX),1)
	 initVar(siteX,2,1:nLayers(siteX)) = multiOut(siteX,ij,7,1:nLayers(siteX),1)
	 initVar(siteX,3:6,1:nLayers(siteX)) = multiOut(siteX,ij,11:14,1:nLayers(siteX),1)
	 initVar(siteX,7,1:nLayers(siteX)) = multiOut(siteX,ij,16,1:nLayers(siteX),1)
     initVar(siteX,2,1:jj) = 0.!output(1,7,:,1)
     initVar(siteX,3:7,1:jj) = 0.!output(1,11:14,:,1)  !!newX
    endif !(maxState(i)>minDharv)
   enddo !end do while
   
 elseif(compHarv(1)==2.) then  !!!thin to compansate harvest limits
   !Perform thinning to compensate harvest levels
   !calculate SDI
  do i = 1, nSites
!!start!! use stand volume to order the forests
if(oldLayer==1) then
 jj=max((nLayers(i)-1),1)
else
 jj=nLayers(i)
endif 
	if(ClCut(i) > 0. ) then
	 maxState(i) = sum(multiOut(i,ij,30,1:jj,1))!!!search for site with highest volume
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
if(oldLayer==1) then
 jj=max((nLayers(siteX)-1),1)
else
 jj=nLayers(siteX)
endif 

! if(siteInfo(siteX,1) == 454702.) write(2,*) "remaining 1", ij, multiOut(siteX,ij,11,1,1), multiOut(siteX,ij,13,1,1),&
	! multiOut(siteX,ij,37,1,1)
! if(siteInfo(siteX,1) == 454702.) write(2,*) "remaining 2", multiOut(siteX,ij,11,2,1), multiOut(siteX,ij,13,2,1), &
	! multiOut(siteX,ij,37,2,1)
! if(siteInfo(siteX,1) == 454702.) write(2,*) "remaining 3", multiOut(siteX,ij,11,3,1), multiOut(siteX,ij,13,3,1), &
	! multiOut(siteX,ij,37,3,1)

	 roundWood = roundWood + sum(multiOut(siteX,ij,30,1:jj,1)*harvRatio)* thinFact *areas(siteX) !!energCuts
     multiOut(siteX,ij,1,1,2) = 4. !!!flag for thinning compensation
	 multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:jj,1) + &
			multiOut(siteX,ij,30,1:jj,1)*harvRatio*thinFact
	 multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:jj,1) + &
			multiOut(siteX,ij,31,1:jj,1)*harvRatio*thinFact
     !update state of the forests
	 do ijj = 1, jj
      multiOut(siteX,ij,9:10,ijj,2) = multiOut(siteX,ij,9:10,ijj,1) * thinFact
	  multiOut(siteX,ij,11:12,ijj,2) = multiOut(siteX,ij,11:12,ijj,1)
	  multiOut(siteX,ij,13,ijj,2) = multiOut(siteX,ij,13,ijj,1) * thinFact
	  multiOut(siteX,ij,14:16,ijj,2) = multiOut(siteX,ij,14:16,ijj,1)
	  multiOut(siteX,ij,17:23,ijj,2) = multiOut(siteX,ij,17:23,ijj,1) * thinFact
	  multiOut(siteX,ij,26:29,ijj,2) = multiOut(siteX,ij,26:29,ijj,1) * thinFact
	  
	  ! if(siteInfo(siteX,1) == 454702.) write(2,*) "thinned x", multiOut(siteX,ij,11,ijj,2), multiOut(siteX,ij,13,ijj,2), &
	  ! multiOut(siteX,ij,30,ijj,2)

	  !update biomasses and Volumes
	  multiOut(siteX,ij,24:25,ijj,2) = multiOut(siteX,ij,24:25,ijj,1) * thinFact + &
					multiOut(siteX,ij,24:25,ijj,2)
	  multiOut(siteX,ij,30:34,ijj,2) = multiOut(siteX,ij,30:34,ijj,1) * thinFact + &
					multiOut(siteX,ij,30:34,ijj,2)
	  
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
	   
	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)*(1-energyRatio) * thinFact +   &
			multiOut(siteX,ij,51,ijj,1)* thinFact + multiOut(siteX,ij,28,ijj,1)) + &
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio) * thinFact + &
		 (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.83 * thinFact))
       multiOut(siteX,ij,29,ijj,1) = (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.17 * thinFact+ &
			multiOut(siteX,ij,29,ijj,1) 
				
	  else
	   multiOut(siteX,ij,28,ijj,1) = max(0.,(multiOut(siteX,ij,24,ijj,1)* thinFact + multiOut(siteX,ij,28,ijj,1) + &
		multiOut(siteX,ij,51,ijj,1)* thinFact + multiOut(siteX,ij,32,ijj,1)*0.83* thinFact + &
		multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)* thinFact))
       multiOut(siteX,ij,29,ijj,1)=multiOut(siteX,ij,32,ijj,1)*0.17* thinFact+multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
	  endif
!!end energCuts
	  multiOut(siteX,ij,9:10,ijj,1) = multiOut(siteX,ij,9:10,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,13,ijj,1) = multiOut(siteX,ij,13,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,17:25,ijj,1) = multiOut(siteX,ij,17:25,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,30:33,ijj,1) = multiOut(siteX,ij,30:33,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,47:51,ijj,1) = multiOut(siteX,ij,47:51,ijj,1)*(1-thinFact)
	  multiOut(siteX,ij,53:nVar,ijj,1) = multiOut(siteX,ij,53:nVar,ijj,1)*(1-thinFact)

   	  initVar(siteX,1,ijj) = multiOut(siteX,ij,4,jj,1)
 	  initVar(siteX,2,ijj) = multiOut(siteX,ij,7,ijj,1)
	  initVar(siteX,3:6,ijj) = multiOut(siteX,ij,11:14,ijj,1)
	  initVar(siteX,7,ijj) = multiOut(siteX,ij,16,ijj,1)
! if(siteInfo(siteX,1) == 454702.) write(2,*) "thinned", multiOut(siteX,ij,11,ijj,2), multiOut(siteX,ij,13,ijj,2)
! if(siteInfo(siteX,1) == 454702.) write(2,*) "remaining", multiOut(siteX,ij,11,ijj,1), multiOut(siteX,ij,13,ijj,1)

     enddo !ijj layers loop
	 
 	!!!if fertilization at thinning is active,  decrease siteType
	if(flagFert(siteX)==0 .and. fertThin>0) then 
		yearsFert = max(1,min(((nYears(siteX)) - ij-1),nYearsFert))
		multiOut(siteX,(ij+1):(ij+yearsFert),3,:,1) = siteInfo(siteX,3)-1.
		call calcAlfar(multiOut(siteX,ij,3,1:nLayers(siteX),:),initVar(siteX,1,1:nLayers(siteX)),pCrobas, &
				nLayers(siteX),alfarFert,allSP,nYearsFert,npar)
		multiOut(siteX,(ij+1):(ij+yearsFert),3,:,2) = alfarFert(1:yearsFert,:)
		flagFert(siteX)=2
	endif

    endif !(maxState(i)>minDharv)
   enddo !end do while
 elseif(compHarv(1)==3.) then  !!!thin to compensate harvest limits
   n = 0
  do while((n < nSites .and. roundWood < HarvLim(ij,1)) .and. &
	cuttingArea(ij,2) < cuttingArea(ij,1))		
   n = n + 1
   do i = 1, nSites
if(oldLayer==1) then
 jj=max((nLayers(i)-1),1)
else
 jj=nLayers(i)
endif 
	if(ClCut(i) > 0. ) then
	 maxState(i) = maxval(multiOut(i,ij,12,1:jj,1))!!!search for site with highest DBH
	else
	 maxState(i) = 0.
	endif
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = int(siteInfo(siteX,2))

	if(maxState(siteX)>minDharv .and. ClCut(siteX) > 0.) then
     energyCutX = energyCuts(siteX)
	 if (HarvLim(ij,2) > 0. .and.  energyWood >= HarvLim(ij,2)) then		!!energCuts
	  energyCutX = 0.
	 endif
	 !if fertilization at thinning is active reset flagFert
	 if(fertThin > 0) then
	  flagFert(siteX) = 0  
	 endif

!!if oldLayer scenario -> populate old layer with the dominant species and clearcut other layers
if(oldLayer==1) then
 jj=max((nLayers(siteX)-1),1)
 ops = maxloc(multiOut(siteX,ij,13,1:jj,1))
 layerX = int(ops(1))
  
  !!!!calculate percentage of trees remaining after clearcut(pDomRem)
   call random_number(randX)
   pDomRem =   (randX*5.+5.)/100.* &  !!randomly sample between 5 and 10 %
   sum(multiOut(siteX,ij,13,1:jj,1))/multiOut(siteX,ij,13,layerX,1)
   

	!update old layer
	multiOut(siteX,ij,:,nLayers(i),1) = multiOut(siteX,ij,:,layerX,1)
	multiOut(siteX,ij,42,nLayers(i),1) = 0.
	multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),nLayers(i),1) = &
	 multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) * (pDomRem)
	multiOut(siteX,ij,24:34,nLayers(i),1) = multiOut(siteX,ij,24:34,layerX,1) * (pDomRem)
	multiOut(siteX,ij,47:51,nLayers(i),1) = multiOut(siteX,ij,47:51,layerX,1) * (pDomRem)
    !!!if oldLayer is active import siteType and alfar from the single site simulations simulations
	multiOut(siteX,ij:maxYears,3,nLayers(siteX),1) = multiOut(siteX,ij,3,layerX,1)
	multiOut(siteX,ij:maxYears,3,nLayers(siteX),2) = multiOut(siteX,ij,3,layerX,2)
	!update dominant layer
	multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) = &
 	 multiOut(siteX,ij,(/9,10,13,17,18,37,38,40,43,44,53,54/),layerX,1) * (1-pDomRem)
	multiOut(siteX,ij,24:34,layerX,1) = multiOut(siteX,ij,24:34,layerX,1) * (1-pDomRem)
	multiOut(siteX,ij,47:51,layerX,1) = multiOut(siteX,ij,47:51,layerX,1) * (1-pDomRem)
else
 jj=nLayers(i)
endif
!!   !!clearcut!!
	 cuttingArea(ij,2) = cuttingArea(ij,2) + areas(siteX) !calculate the clearcut area
	   roundWood = roundWood + sum(multiOut(siteX,ij,30,1:jj,1)*harvRatio)*areas(siteX) !!energCuts
	   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:jj,1) + &
			multiOut(siteX,ij,30,1:jj,1)*harvRatio
	   multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:jj,1) + &
			multiOut(siteX,ij,31,1:jj,1)*harvRatio
	 multiOut(siteX,ij,2,1,2) = 2. !!!flag for clearcut compensation
     do ijj = 1, jj
      multiOut(siteX,ij,6:23,ijj,2) = multiOut(siteX,ij,6:23,ijj,1)
      multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
      multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
      multiOut(siteX,ij,28:29,ijj,2) = multiOut(siteX,ij,28:29,ijj,1)
      multiOut(siteX,ij,35:nVar,ijj,2) = multiOut(siteX,ij,35:nVar,ijj,1)
	  !update biomasses and Volumes
	  multiOut(siteX,ij,24:25,ijj,2) = multiOut(siteX,ij,24:25,ijj,1) + &
					multiOut(siteX,ij,24:25,ijj,2)
	  multiOut(siteX,ij,30:34,ijj,2) = multiOut(siteX,ij,30:34,ijj,1) + &
					multiOut(siteX,ij,30:34,ijj,2)
!!energCuts
      if(energyCutX == 1.) then
	   multiWood(siteX,ij,ijj,2) = multiWood(siteX,ij,ijj,2) + (multiOut(siteX,ij,24,ijj,1) + &
	    multiOut(siteX,ij,32,ijj,1)*0.3 + multiOut(siteX,ij,31,ijj,1)* (1-harvRatio)) * energyRatio
	   species = int(multiOut(siteX,ij,4,ijj,1))
	   multiWood(siteX,ij,ijj,1) = multiWood(siteX,ij,ijj,2) / pCrobas(2,species)
	   energyWood = energyWood + multiWood(siteX,ij,ijj,1) * areas(siteX)   !!energCuts !!!we are looking at volumes
	   
   	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)*(1-energyRatio)  +   &
			multiOut(siteX,ij,51,ijj,1) + multiOut(siteX,ij,28,ijj,1)) + &
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio) + &
		 (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.83 ))
       multiOut(siteX,ij,29,ijj,1) = (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.17 + &
			multiOut(siteX,ij,29,ijj,1) 
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
	  relBA(siteX,1:nLayers(i)) = initVar(siteX,5,1:nLayers(i))/ &
		sum(initVar(siteX,5,1:nLayers(i)))
	 else
	  relBA(siteX,1:nLayers(i)) = multiOut(siteX,(ij-1),13,1:nLayers(i),1)/ &
		sum(multiOut(siteX,(ij-1),13,1:nLayers(i),1))
	 endif

   	 ! initVar(siteX,1,1:nLayers(siteX)) = multiOut(siteX,ij,4,1:nLayers(siteX),1)
	 ! initVar(siteX,2,1:nLayers(siteX)) = multiOut(siteX,ij,7,1:nLayers(siteX),1)
	 ! initVar(siteX,3:6,1:nLayers(siteX)) = multiOut(siteX,ij,11:14,1:nLayers(siteX),1)
	 ! initVar(siteX,7,1:nLayers(siteX)) = multiOut(siteX,ij,16,1:nLayers(siteX),1)
     initVar(siteX,2,1:jj) = 0.!output(1,7,:,1)
     initVar(siteX,3:7,1:jj) = 0.!output(1,11:14,:,1)  !!newX
    endif !(maxState(i)>minDharv)
   enddo !end do while

  do i = 1, nSites
if(oldLayer==1) then
 jj=max((nLayers(i)-1),1)
else
 jj=nLayers(i)
endif 
! !!start!! use stand density index to order the forests
   ! if(ClCut(i) > 0.) then
	! call calRein(multiOut(i,ij,:,:,1),nLayers(i),pCrobas(17,:),nVar,allSP,maxState(i))
   ! else
	! maxState(i) = 0.
   ! endif
! !!end!! use stand density index to order the forests

!!start!! use stand volume to order the forests
	if(ClCut(i) > 0. ) then
	 maxState(i) = sum(multiOut(i,ij,30,1:jj,1))!!!search for site with highest DBH
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
if(oldLayer==1) then
 jj=max((nLayers(siteX)-1),1)
else
 jj=nLayers(siteX)
endif 

	 roundWood = roundWood + sum(multiOut(siteX,ij,30,1:jj,1)*harvRatio)* thinFact *areas(siteX) !!energCuts
     multiOut(siteX,ij,1,1,2) = 4. !!!flag for thinning compensation
	 multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:jj,1) + &
			multiOut(siteX,ij,30,1:jj,1)*harvRatio*thinFact
	 multiOut(siteX,ij,38,:,1) = multiOut(siteX,ij,38,1:jj,1) + &
			multiOut(siteX,ij,31,1:jj,1)*harvRatio*thinFact
     !update state of the forests
	 do ijj = 1, jj
      multiOut(siteX,ij,9:10,ijj,2) = multiOut(siteX,ij,9:10,ijj,1) * thinFact
	  multiOut(siteX,ij,11:12,ijj,2) = multiOut(siteX,ij,11:12,ijj,1)
	  multiOut(siteX,ij,13,ijj,2) = multiOut(siteX,ij,13,ijj,1) * thinFact
	  multiOut(siteX,ij,14:16,ijj,2) = multiOut(siteX,ij,14:16,ijj,1)
	  multiOut(siteX,ij,17:23,ijj,2) = multiOut(siteX,ij,17:23,ijj,1) * thinFact
	  multiOut(siteX,ij,26:29,ijj,2) = multiOut(siteX,ij,26:29,ijj,1) * thinFact
	  
	  !update biomasses and Volumes
	  multiOut(siteX,ij,24:25,ijj,2) = multiOut(siteX,ij,24:25,ijj,1) * thinFact + &
					multiOut(siteX,ij,24:25,ijj,2)
	  multiOut(siteX,ij,30:34,ijj,2) = multiOut(siteX,ij,30:34,ijj,1) * thinFact + &
					multiOut(siteX,ij,30:34,ijj,2)

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

	   multiOut(siteX,ij,28,ijj,1) = max(0.,((multiOut(siteX,ij,24,ijj,1)*(1-energyRatio) * thinFact +   &
			multiOut(siteX,ij,51,ijj,1)* thinFact + multiOut(siteX,ij,28,ijj,1)) + &
		 multiOut(siteX,ij,31,ijj,1)* (1-harvRatio) * (1-energyRatio) * thinFact + &
		 (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.83 * thinFact))
       multiOut(siteX,ij,29,ijj,1) = (0.3 * (1-energyRatio)+0.7) * multiOut(siteX,ij,32,ijj,1) *0.17 * thinFact+ &
			multiOut(siteX,ij,29,ijj,1) 

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

	 initVar(siteX,1,ijj) = multiOut(siteX,ij,4,ijj,1)
	 initVar(siteX,2,ijj) = multiOut(siteX,ij,7,ijj,1)
	 initVar(siteX,3:6,ijj) = multiOut(siteX,ij,11:14,ijj,1)
	 initVar(siteX,7,ijj) = multiOut(siteX,ij,16,ijj,1)
     enddo !ijj layers loop
	
	! if(siteInfo(siteX,1)==411310.) write(1,*) ij,multiOut(siteX,ij,11,:,1)
	! if(siteInfo(siteX,1)==35.) write(1,*) ij,multiOut(siteX,ij,11,:,1)
	
	 	!!!if fertilization at thinning is active,  increase siteType
	if(flagFert(siteX)==0 .and. fertThin>0) then 
		yearsFert = max(1,min(((nYears(siteX)) - ij-1),nYearsFert))
		multiOut(siteX,(ij+1):(ij+yearsFert),3,:,1) = siteInfo(siteX,3)-1.
		call calcAlfar(multiOut(siteX,ij,3,1:nLayers(siteX),:),initVar(siteX,1,1:nLayers(siteX)),pCrobas, &
				nLayers(siteX),alfarFert,allSP,nYearsFert,npar)
		multiOut(siteX,(ij+1):(ij+yearsFert),3,:,2) = alfarFert(1:yearsFert,:)
		flagFert(siteX)=2
	endif

	endif !(maxState(i)>minDharv)
   enddo !end do while
 
 
 endif 
endif !roundWood < HarvLim .and. HarvLim /= 0.

  !HarvLim(ij,1) = roundWood
  !HarvLim(ij,2) = energyWood
end do !end Year loop 

 do i = 1,nSites
  do ij = 1, maxYears 
    do ijj = 1,nLayers(i)
	  ! multiOut(i,ij,38,ijj,1) = sum(multiOut(i,1:ij,30,ijj,2)) + &
		! sum(multiOut(i,1:ij,42,ijj,1)) + multiOut(i,ij,30,ijj,1)

!!!!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
	   if (ij==1) D = multiOut(i,ij,12,ijj,1)
	   if (ij>1) D = multiOut(i,(ij-1),12,ijj,1)
	   Vmort = multiOut(i,ij,42,ijj,1)
	   if(Vmort>0. .and. ij < maxYears)then
		species = int(multiOut(i,ij,4,ijj,1))
		multiOut(i,ij,8,ijj,1) = Vmort + multiOut(i,ij,8,ijj,1)
		do ki=1,(maxYears-ij)
		 multiOut(i,(ij+ki),8,ijj,1) = multiOut(i,(ij+ki),8,ijj,1) + Vmort * & 
		   exp(-exp(pCrobas(35,species) + pCrobas(36,species)*ki + &
					 pCrobas(37,species)*D + pCrobas(44,species)))
		enddo
	   endif
		
	! !!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
	! do ijj = 1,nLayers(i)
	  ! if(output(1,8,ijj,1)>0.) then
	  ! multiOut(i,ij,8,ijj,1) = multiOut(i,ij,8,ijj,1) + output(1,8,ijj,1)
	  ! jj = int(output(1,4,ijj,1))
	    ! do ki = 1,(maxYears-ij)
			! multiOut(i,(ki+ij),8,ijj,1) = multiOut(i,(ki+ij),8,ijj,1) + output(1,8,ijj,1) * &
				! exp(-exp(pCrobas(34,jj) + pCrobas(35,jj)*ki + pCrobas(36,jj)*output(1,12,ijj,1) + pCrobas(44,jj)))
		! enddo
	  ! end if
	! enddo

	  if(ij > 1.5) then
	!compute gross growth
	   multiOut(i,ij,43,ijj,1) = multiOut(i,ij,30,ijj,1) - multiOut(i,(ij-1),30,ijj,1) + & 
			multiOut(i,ij,37,ijj,1)/harvRatio + multiOut(i,ij,42,ijj,1)
	  endif

    enddo !ijj
  enddo
 enddo	
  close(1)
  close(2)
  close(3)
soilCinOut = soilC
soilCtotInOut = soilCtot

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
