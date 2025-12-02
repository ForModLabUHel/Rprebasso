 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,maxYears,maxThin, &
    nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
    nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
    weatherPRELES,DOY,pPRELES, soilC,pYasso,&
    pAWEN,weatherYasso,litterSize,soilCtot, &
    defaultThin,ClCut,energyCuts,clct_pars,dailyPRELES,yassoRun,multiEnergyWood, &
    tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,thinInt, &
    flagFert,nYearsFert,mortMod,pECMmod,& !protect removed btw nYearsFert and mortModX, neither in prebas subroutine nor multiPrebas() R function
    layerPRELES,LUEtrees,LUEgv, siteInfoDist, outDist, prebasFlags, &
	latitude, TsumSBBs)


implicit none

integer, parameter :: nVar=54,npar=53!, nSp=3
integer, intent(in) :: nSites, maxYears,maxThin,nClimID,maxNlayers,allSP
integer, intent(in) :: nYears(nSites),nLayers(nSites) !protect removed; neither in prebas subroutine nor multiPrebas() R function

 integer :: i,climID,ij,iz,ijj,ki,n,jj,az
 real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),layerPRELES !, ECMmod fvec
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP),tapioPars(5,2,3,20),pECMmod(12)
 real (kind=8), intent(inout) :: tTapio(5,allSP,2,7), ftTapio(5,allSP,3,7),mortMod(2)
 real (kind=8), intent(inout) :: siteInfo(nSites,11),thdPer(nSites),limPer(nSites)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,11),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: LUEtrees(allSP),LUEgv
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)  !initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: clct_pars(nSites,allSP,3),energyCuts(nSites)  !!energCuts
 real (kind=8), intent(in) :: thinInt(nSites) !site specific parameter that determines the thinning intensity; 
          !from below (thinInt>1) or above (thinInt<1);thinInt=999. uses the default value from tapio rules


! logical :: disturbanceON !!!this could be site specific but to block dist. in some sites you can work on the inputs
real (kind=8), intent(inout) :: siteInfoDist(nSites,10), outDist(nSites,maxYears,10) !inputs(siteInfoDist) & outputs(outDist) of disturbance modules

 !!! fertilization parameters
 !integer, intent(inout) :: fertThin !!! flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 fvec
 integer, intent(inout) :: flagFert !!! flag that indicates if fertilization has already been applied along the rotation
 integer, intent(inout) :: nYearsFert !!number of years for which the fertilization is effective

!!!ground vegetation
 !integer, intent(in) :: gvRun      !!!ground vegetation fvec
 real (kind=8), intent(inout) :: GVout(nSites,maxYears,5) !fAPAR_gv,litGV,photoGV,wGV      !!!ground vegetation
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,7,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2),latitude(nSites), TsumSBBs(nSites,4)
 real (kind=8), intent(inout) :: multiEnergyWood(nSites,maxYears,maxNlayers,2)!!energCuts
 real (kind=8), intent(inout) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 ! real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(maxYears,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers),mortModX
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears), thinningX(maxThin,11)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species

 integer :: etmodel,CO2model, gvRun, fertThin, ECMmod, oldLayer !not direct inputs anymore, but in prebasFlags fvec !wdimpl pflags
 integer, intent(inout) :: prebasFlags(10)

!!! 'un-vectorise' flags, fvec
etmodel = prebasFlags(1)
gvRun = prebasFlags(2)
fertThin = prebasFlags(3)
oldLayer = prebasFlags(4)
ECMmod = prebasFlags(5)
CO2model = prebasFlags(7)
! if(prebasFlags(6)==0) disturbanceON = .FALSE.
! if(prebasFlags(6)==1) disturbanceON = .TRUE.

!outDist(:,:,:) = 99!prebasFlags(6)
!outDist(1,10) = siteInfoDist(1,1)
!!!!initialize run
! multiOut = 0.
! open(1,file="test1.txt")
! open(2,file="test2.txt")

output = 0.
yearX = 0.
multiEnergyWood = 0.
!soilC = soilCinOut
!soilCtot = soilCtotInOut
do i = 1,nSites
 do ijj = 1,nLayers(i)
  species = int(initVar(i,1,ijj))
!    initVar(i,7,ijj) = max(0.,pCrobas(38,species)/pCrobas(15,species) * (initVar(i,3,ijj) -&
!      initVar(i,6,ijj))**pCrobas(11,species))!A = p_ksi/p_rhof * Lc^p_z
  call initBiomasses(pCrobas(:,species),initVar(i,:,ijj),siteInfo(i,3),multiOut(i,1,:,ijj,1),nVar,npar)
 enddo
enddo

do i = 1,nSites
 prebasFlags(8) = int(multiOut(i,1,7,1,2))
 multiOut(i,1,7,1,2) = 0.
 output(:,:,:,:) = multiOut(i,:,:,:,:)

  climID = siteInfo(i,2)
  defaultThinX = defaultThin(i)
  ClCutX = ClCut(i)
  
  !!!##set mortality model for managed and unmanaged forests
  mortModX = mortMod(1) !!mortality model to be used in the managed forests
  if(ClCut(i) < 0.5 .and. defaultThin(i) < 0.5) mortModX = mortMod(2) !!mortality model to be used in the unmanaged forests
  
  thinningX = thinning(i,:,:)
  ! nYears(i) = nYears(i)
    call prebas(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
    thinningX(1:nThinning(i),:),output(1:nYears(i),:,1:nLayers(i),:),nThinning(i),maxYearSite,fAPAR(i,1:nYears(i)), &
    initClearcut(i,:),fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),&
    P0y(climID,1:nYears(i),:),weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES, &
    soilC(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
    litterSize,soilCtot(i,1:nYears(i)),defaultThinX,&
    ClCutX,energyCuts(i),clct_pars(i,:,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i),&
    multiEnergyWood(i,1:nYears(i),1:nLayers(i),:),tapioPars,thdPer(i),limPer(i),ftTapio,tTapio,&
    GVout(i,1:nYears(i),:),thinInt(i), &
    flagFert,nYearsFert,mortModX,pECMmod,layerPRELES,LUEtrees,LUEgv, & !protect removed btw nYearsFert and mortModX, neither in prebas subroutine nor multiPrebas() R function
    siteInfoDist(i,:), outDist(i,1:nYears(i),:), prebasFlags,latitude(i), TsumSBBs(i,:))
    
    ! pre flag vectorisatio:
    ! call prebas(nYears(i),nLayers(i),allSP,siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
    ! thinningX(1:nThinning(i),:),output(1:nYears(i),:,1:nLayers(i),:),nThinning(i),maxYearSite,fAPAR(i,1:nYears(i)), &
    ! initClearcut(i,:),fixBAinitClarcut(i),initCLcutRatio(i,1:nLayers(i)),ETSy(climID,1:nYears(i)),&
    ! P0y(climID,1:nYears(i),:),weatherPRELES(climID,1:nYears(i),:,:),DOY,pPRELES,etmodel, &
    ! soilC(i,1:nYears(i),:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,1:nYears(i),:),&
    ! litterSize,soilCtot(i,1:nYears(i)),defaultThinX,&
    ! ClCutX,energyCuts(i),inDclct(i,:),inAclct(i,:),dailyPRELES(i,1:(nYears(i)*365),:),yassoRun(i),&
    ! multiEnergyWood(i,1:nYears(i),1:nLayers(i),:),tapioPars,thdPer(i),limPer(i),ftTapio,tTapio,&
    ! GVout(i,1:nYears(i),:),GVrun,thinInt(i), &
    ! fertThin,flagFert,nYearsFert,protect,mortModX,ECMmod,pECMmod,layerPRELES,LUEtrees,LUEgv, &
    ! disturbanceON, siteInfoDist(i,:), outDist(i,1:nYears(i),:))
    
    
    
    multiOut(i,1:nYears(i),:,1:nLayers(i),:) = output(1:nYears(i),:,1:nLayers(i),:)
end do
 ! close(1)
 ! close(2)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
