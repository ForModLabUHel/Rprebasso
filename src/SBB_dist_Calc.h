!!!! Bark beetle CALCULATIONS !!!!
  call TsumSBBfun(latitude,weatherPRELES(year,:,2),TsumSBBs(4)) 
  if(TsumSBBs(1)<=-998.) then   !!!initialize the first three years this will be done only in the first year of the simulations if the inputs are not provvided
   TsumSBBs(1:3) = TsumSBBs(4)
  endif
  if(year > 1) then
   SMIt0 = modOut(year,46,1,2)
  else
   if(SMIt0 < -998.d0) SMIt0 = SMI !!if year 1 SMIt0 is the same of first year  
  endif
  call spruceVars(outt((/4,7,13/),:,1),nLayers,(/2,10/),2,spruceStandVars,rBAspruce)
  call riskBB(pBB,TsumSBBs,spruceStandVars(1),spruceStandVars(3),spruceStandVars(2),(SMI+SMIt0)/2.)
  !update output
  !modOut((year+1),45,:,2) = 0.
  !modOut((year+1),45,1,2) = pBB(1)
  outt(45,:,2) = 0.
  outt(45,1,2) = pBB(1)
  TsumSBBs(1:3) = TsumSBBs(2:4)  
  
!calculate intensity
  if(spruceStandVars(3)>0.) then
   call bb_imp_mod(SMIt0,spruceStandVars(1),intenSpruce)
 !old version (start)
   ! SHI = (spruceStandVars(1)/spruceStandVars(3))*(1.0-SMIt0)/0.2093014 !(spruceStandVars(1)/spruceStandVars(3)) = Baspruce fraction
   ! intenSpruce = 1.d0/(1.d0+exp(3.9725-2.9673*SHI))
 !old version (end)
   if((spruceStandVars(1)/spruceStandVars(3)) < 0.05) intenSpruce = 0. ! If no spruce, damage intensity set zero
  endif
  outt(48,1,2) = intenSpruce !* spruceStandVars(1)

if(disturbance_bb) then !!!!mortality caused by bark beetle is switched off for now
!!!sample to see if bb attack occurs
  call random_number(rndm)
  if(rndm <= pBB(1)) then !if BB attack
!   outt(48,:,2) = intenSpruce * rBAspruce
   pHarvTrees = 0. !!!for now is 0 then we need to modify this if salvage logging, see jonathan code
   call baBBdist_bylay_fun(outt((/4,11,13/),:,1),nLayers,(/2,10/),2,BAdist,((intenSpruce) * spruceStandVars(1)))
   outt(43,:,2)=BAdist! BAdist = intenSpruce * rBAspruce * STAND_all(13,:)
   
   include 'updateStandAfterDist.h'

  endif !endif BB attack
endif
  
!!reupdate output after bark beetle dist
  outt(:,:,1) = STAND_all
  modOut((year+1),2,:,:) = outt(2,:,:)
  modOut((year+1),4,:,:) = outt(4,:,:) !update species
  modOut((year+1),7,:,:) = outt(7,:,:)
  modOut((year+1),9:nVar,:,:) = outt(9:nVar,:,:)
