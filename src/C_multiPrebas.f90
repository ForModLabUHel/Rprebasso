 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,inDclct,inAclct,dailyPRELES,yassoRun,prebasVersion)

implicit none

integer, parameter :: nVar=46,npar=38!, nSp=3
integer, intent(in) :: nYears(nSites),nLayers(nSites),allSP

 integer :: i,climID,ij,iz,ijj,ki,n,jj,az
integer, intent(in) :: nSites, maxYears,maxThin,nClimID,maxNlayers
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP)
 real (kind=8), intent(inout) :: siteInfo(nSites,7)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,8),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites),prebasVersion(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,7,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(maxYears,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers)
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears), thinningX(maxThin,8)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species

!!!!initialize run
multiOut = 0.
output = 0.
yearX = 0.
soilC = soilCinOut
soilCtot = soilCtotInOut
do i = 1,nSites
 do ijj = 1,nLayers(i)
	species = int(multiOut(i,1,4,ijj,1))
		initVar(i,7,ijj) = pCrobas(38,species)/pCrobas(15,species) * (initVar(i,3,ijj) -&
			initVar(i,6,ijj))**pCrobas(11,species)!A = p_ksi/p_rhof * Lc^p_z
 enddo
enddo

 do i = 1,nSites
 ! write(*,*) i
	climID = siteInfo(i,2)
	thinningX = thinning(i,:,:)
	! nYears(i) = nYears(i)
	  call prebas_v0(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX,output(1:nYears(i),:,1:nLayers(i),:),maxThin,maxYearSite,fAPAR(i,1:nYears(i)),initClearcut(i,:),&
		fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),P0y(climID,1:nYears(i),:),&
		weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES,etmodel, &
		soilC(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
		litterSize,soilCtot(i,1:nYears(i)),&
		defaultThinX,ClCutX,inDclct(i,:),inAclct(i,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i))

		multiOut(i,1:nYears(i),:,1:nLayers(i),:) = output(1:nYears(i),:,1:nLayers(i),:)
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
