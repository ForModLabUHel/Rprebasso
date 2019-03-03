
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine bridging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine regionPrebas(siteOrder,HarvLim,minDharv,multiOut,nSites,areas,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,inDclct,inAclct,dailyPRELES,yassoRun,prebasVersion)

implicit none

integer, parameter :: nVar=46,npar=37!, nSp=3
integer, intent(in) :: nYears(nSites),nLayers(nSites),allSP
integer :: i,climID,ij,iz,ijj,ki,n,jj,az
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers,siteOrder(nSites,maxYears)
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5),HarvLim(maxYears),minDharv
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP)
 real (kind=8), intent(inout) :: siteInfo(nSites,7), areas(nSites)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,8),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites),prebasVersion(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,6,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(1,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers)
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears), thinningX(maxThin,8)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1)

!!!!initialize run
multiOut = 0.
yearX = 0.
soilC = soilCinOut
soilCtot = soilCtotInOut

do i = 1,nSites
 relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
enddo

do ij = 1,maxYears
 HarvArea = 0.
 do iz = 1,nSites
 	i=siteOrder(iz,ij)
! open(10,file="multiSite.txt")
 ! write(10,*) "years =",ij, "siteRun = ",iz
! close(10)
	ClCutX = ClCut(i)
	defaultThinX = defaultThin(i)
	thinningX(:,:) = -999.
	az = 0

	if(ij > 1) then
	 soilC(i,ij,:,:,1:nLayers(i)) = soilC(i,(ij-1),:,:,1:nLayers(i))
	endif

!!!check if the limit has been exceeded if yes no havest (thinning or clearcut will be performed)
	if (HarvLim(ij) > 0. .and. HarvArea >= HarvLim(ij)) then
	 ClCutX = 0.
	 defaultThinX = 0.
	endif
!!!
	climID = siteInfo(i,2)
	if(ij==int(min(yearX(i),maxYears)))then
	 initClearcut(i,5) = int(min(initClearcut(i,5), initClearcut(i,5) + maxYears - yearX(i)))
	 yearX(i) = 0

	 do ijj = 1,nLayers(i)
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
	  do ki = 1,int(initClearcut(i,5)+1)
	   multiOut(i,int(ij-initClearcut(i,5)+ki-1),7,ijj,1) = ki !#!#
	  enddo !ki
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
  
  
	if(prebasVersion(i)==0.) then
	  call prebas_v0(1,nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX(1:az,:),output(1,:,1:nLayers(i),:),az,maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,ij),P0y(climID,ij,:),&
		weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilC(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize,soilCtot(i,ij),&
		defaultThinX,ClCutX,inDclct(i,:),inAclct(i,:),dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i))
	elseif(prebasVersion(i)==1.) then
	  call prebas_v1(1,nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX(1:az,:),output(1,:,1:nLayers(i),:),az,maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,ij),P0y(climID,ij,:),&
		weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilC(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize,soilCtot(i,ij),&
		defaultThinX,ClCutX,inDclct(i,:),inAclct(i,:),dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i))
	endif
	
	! if clearcut occur initialize initVar and age
	if(sum(output(1,11,1:nLayers(i),1))==0 .and. yearX(i) == 0) then
	 if((maxYears-ij)<15) then
	  Ainit = nint(6 + 2*3.5 - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*3.5 - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 yearX(i) = Ainit + ij + 1
	 initClearcut(i,5) = Ainit
	 if(ij==1) then
	  relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
	 endif
	endif
	
	! write(10,*) "here1"
	  !!!calculate deadWood using Gompetz function (Makinen et al. 2006)!!!!
	do ijj = 1,nLayers(i)
	  if(output(1,8,ijj,1)>0.) then
	  multiOut(i,ij,8,ijj,1) = multiOut(i,ij,8,ijj,1) + output(1,8,ijj,1)
	  jj = int(output(1,4,ijj,1))
	    do ki = 1,(maxYears-ij)
			multiOut(i,(ki+ij),8,ijj,1) = multiOut(i,(ki+ij),8,ijj,1) + output(1,8,ijj,1) * &
				exp(-exp(pCrobas(34,jj) + pCrobas(35,jj)*ijj + pCrobas(36,jj)*output(1,12,ijj,1) + 0.))
		enddo
	  end if
	enddo
	
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
! write(10,*) "here2"

	initVar(i,1,1:nLayers(i)) = output(1,4,1:nLayers(i),1)
	initVar(i,2,1:nLayers(i)) = output(1,7,1:nLayers(i),1)
	initVar(i,3:6,1:nLayers(i)) = output(1,11:14,1:nLayers(i),1)
	HarvArea = HarvArea + 300!sum(output(1,37,1:nLayers(i),1))*0.5!areas(i)
 end do !iz i


! write(10,*) "here3"


 !!! check if the harvest limit of the area has been reached otherwise clearcut the stands sorted by basal area
 if (HarvArea < HarvLim(ij)) then
  n = 0
  do while(n < nSites .and. HarvArea < HarvLim(ij))
   n = n + 1
   do i = 1, nSites
	maxState(i) = maxval(multiOut(i,ij,12,1:nLayers(i),1))!!!search for site with highest DBH
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = int(siteInfo(siteX,2))
if(maxState(siteX)>minDharv .and. ClCut(siteX) > 0.) then
  ! close(10)
!!   !!clearcut!!
   HarvArea = HarvArea + 300!sum(multiOut(siteX,ij,30,1:nLayers(siteX),1))*areas(i)
   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,1:nLayers(siteX),1) + multiOut(siteX,ij,30,1:nLayers(siteX),1)
   do ijj = 1, nLayers(siteX)
    multiOut(siteX,ij,6:nVar,ijj,2) = multiOut(siteX,ij,6:nVar,ijj,1)
    multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
    multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
    multiOut(siteX,ij,28,ijj,1) = multiOut(siteX,ij,24,ijj,1) + multiOut(siteX,ij,28,ijj,1)
    multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,31,ijj,1)* 0.1 + &
	multiOut(siteX,ij,32,ijj,1) + multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
    multiOut(siteX,ij,8:21,ijj,1) = 0.
    multiOut(siteX,ij,23:36,ijj,1) = 0. !#!#
    multiOut(siteX,ij,43:44,ijj,1) = 0.
    multiOut(siteX,ij,38,ijj,1) = sum(multiOut(siteX,1:ij,30,ijj,2)) + &
		sum(multiOut(siteX,1:ij,42,ijj,1)) + multiOut(siteX,ij,30,ijj,1)
   enddo
	 if((maxYears-ij)<10) then
	  Ainit = nint(6 + 2*3.5 - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*3.5 - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
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
  initVar(siteX,3:6,1:nLayers(siteX)) = 0.!output(1,11:14,:,1)
endif !(maxState(i)>minDharv)
  enddo !end do while
 endif !HarvArea < HarvLim .and. HarvLim /= 0.
! write(10,*) "here4"
end do
do i = 1,nSites
  do ij = 1, maxYears 
	do ijj = 1,nLayers(i)
	  ! multiOut(i,ij,38,ijj,1) = sum(multiOut(i,1:ij,30,ijj,2)) + &
		! sum(multiOut(i,1:ij,42,ijj,1)) + multiOut(i,ij,30,ijj,1)
		
	  if(ij > 1.5) then
	!compute gross growth
	   multiOut(i,ij,43,ijj,1) = multiOut(i,ij,30,ijj,1) - multiOut(i,(ij-1),30,ijj,1) + & 
			multiOut(i,ij,37,ijj,1) + multiOut(i,ij,42,ijj,1)
	  endif

	enddo !ijj
  enddo
enddo	
! close(10)
! write(10,*) "here5"
! close(10)
soilCinOut = soilC
soilCtotInOut = soilCtot

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


