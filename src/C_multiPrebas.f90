 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,energyCuts,inDclct,inAclct,dailyPRELES,yassoRun,multiEnergyWood, &
		tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,GVrun) !!energCut

implicit none

integer, parameter :: nVar=54,npar=43!, nSp=3
integer, intent(in) :: nSites, maxYears,maxThin,nClimID,maxNlayers,allSP
integer, intent(in) :: nYears(nSites),nLayers(nSites)

 integer :: i,climID,ij,iz,ijj,ki,n,jj,az
 real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP),tapioPars(5,2,3,20)
 real (kind=8), intent(in) :: tTapio(5,3,2,7), ftTapio(5,3,3,7)
 real (kind=8), intent(inout) :: siteInfo(nSites,10),thdPer(nSites),limPer(nSites)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,9),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP),energyCuts(nSites)	!!energCuts
 !!!ground vegetation
 integer, intent(in) :: gvRun			!!!ground vegetation
 real (kind=8), intent(inout) :: GVout(nSites,maxYears,3) !fAPAR_gv,litGV,photoGV,respGV			!!!ground vegetation
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,7,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: multiEnergyWood(nSites,maxYears,maxNlayers,2)!!energCuts
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(maxYears,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers)
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears), thinningX(maxThin,9)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species

!!!!initialize run
! multiOut = 0.
! open(1,file="test1.txt")
 ! open(1,file="ftTapioMsite.txt")
 ! open(2,file="tTapioMsite.txt")
 ! write(1,*) ftTapio
 ! write(2,*) tTapio
 ! close(1)
 ! close(2)

output = 0.
yearX = 0.
multiEnergyWood = 0.
soilC = soilCinOut
soilCtot = soilCtotInOut
do i = 1,nSites
 do ijj = 1,nLayers(i)
	species = int(initVar(i,1,ijj))
		initVar(i,7,ijj) = max(0.,pCrobas(38,species)/pCrobas(15,species) * (initVar(i,3,ijj) -&
			initVar(i,6,ijj))**pCrobas(11,species))!A = p_ksi/p_rhof * Lc^p_z
	call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),siteInfo(i,3),multiOut(i,1,:,ijj,1))
 enddo
enddo

do i = 1,nSites
 ! write(*,*) i
 output(1,:,:,:) = multiOut(i,1,:,:,:)
 ! write(1,*) i

	climID = siteInfo(i,2)
	defaultThinX = defaultThin(i)
	ClCutX = ClCut(i)
	thinningX = thinning(i,:,:)
	! nYears(i) = nYears(i)
	  call prebas(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinningX(1:nThinning(i),:),output(1:nYears(i),:,1:nLayers(i),:),nThinning(i),maxYearSite,fAPAR(i,1:nYears(i)), &
		initClearcut(i,:),fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),&
		P0y(climID,1:nYears(i),:),weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES,etmodel, &
		soilC(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
		litterSize,soilCtot(i,1:nYears(i)),defaultThinX,&
		ClCutX,energyCuts(i),inDclct(i,:),inAclct(i,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i),&
		multiEnergyWood(i,1:nYears(i),1:nLayers(i),:),tapioPars,thdPer(i),limPer(i),ftTapio,tTapio,&
		GVout(i,1:nYears(i),:),GVrun) !energyCut)
		
		multiOut(i,1:nYears(i),:,1:nLayers(i),:) = output(1:nYears(i),:,1:nLayers(i),:)
end do
! close(1)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
