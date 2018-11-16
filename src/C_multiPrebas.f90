 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut, fixBAinitClearcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,inDclct,inAclct,dailyPRELES,yassoRun,prebasVersion)

implicit none

integer, parameter :: nVar=46,npar=37!
integer, intent(in) :: nYears(nSites),nLayers(nSites),allSP
integer :: i,climID
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP)
 real (kind=8), intent(inout) :: siteInfo(nSites,7),fixBAinitClearcut(nSites),initCLcutRatio(nSites,maxNlayers)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,8),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites),prebasVersion(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(in) :: initVar(nSites,6,maxNlayers),P0y(nClimID,maxYears),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(maxYears,nVar,maxNlayers,2)
 integer :: maxYearSite = 100000000

 multiOut = 0.
 do i = 1,nSites
 ! write(*,*) i
	climID = siteInfo(i,2)
	if(prebasVersion(i)==0.) then
	  call prebas_v0(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinning(i,1:nThinning(i),:),output(1:nYears(i),:,1:nLayers(i),:),nThinning(i),maxYearSite,&
		fAPAR(i,1:nYears(i)),initClearcut(i,:),&
		fixBAinitClearcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),P0y(climID,1:nYears(i)),&
		weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES,etmodel, &
		soilCinOut(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
		litterSize,soilCtotInOut(i,1:nYears(i)),&
		defaultThin(i),ClCut(i),inDclct(i,:),inAclct(i,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i))
	elseif(prebasVersion(i)==1.) then
	  call prebas_v1(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinning(i,1:nThinning(i),:),output(1:nYears(i),:,1:nLayers(i),:),nThinning(i),maxYearSite,&
		fAPAR(i,1:nYears(i)),initClearcut(i,:),&
		fixBAinitClearcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),P0y(climID,1:nYears(i)),&
		weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES,etmodel, &
		soilCinOut(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
		litterSize,soilCtotInOut(i,1:nYears(i)),&
		defaultThin(i),ClCut(i),inDclct(i,:),inAclct(i,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i))
	endif
	multiOut(i,:,:,1:nLayers(i),:) = output(:,:,1:nLayers(i),:)
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


