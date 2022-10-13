###
#' calBioIndices function 
#' @Description 
#'  Computed some habitat suitability indices 
#'
#' @param modOut The Prebass output 
#'
#' @return the habitat suitability indices
#' @export
#'
#' @examples
calBioIndices <- function(modOut){
  ### test if is is a multiSite version or a single site version
  if(class(modOut) %in% c("multiPrebas","regionPrebas")){
    ##calculate inputs to the indices
    Vpine = apply(modOut$multiOut,1:2,calVspec,speciesID=1)
    Vspruce = apply(modOut$multiOut,1:2,calVspec,speciesID=2)
    Vdec = apply(modOut$multiOut,1:2,calVspec,speciesID=3)
    VdecDead = apply(modOut$multiOut,1:2,calVdeadSpec,speciesID=3)
    Vtot = apply(modOut$multiOut[,,30,,1],1:2,sum,na.rm=T)
    VtotDead= apply(modOut$multiOut[,,8,,1],1:2,sum,na.rm=T)
    Ppine = Vpine * 100 / Vtot; Ppine[which(is.na(Ppine))]=0
    Pspruce = Vspruce / Vtot *100; Pspruce[which(is.na(Pspruce))]=0
    Pdec = Vdec / Vtot *100; Pdec[which(is.na(Pdec))]=0
    ageoldest = apply(modOut$multiOut[,,7,,1],1:2,max,na.rm=T)
    Nstems = apply(modOut$multiOut[,,17,,1],1:2,sum,na.rm=T)
    BAtot = apply(modOut$multiOut[,,13,,1],1:2,sum,na.rm=T)
    BAmortTot = (Vtot/VtotDead)*BAtot; BAmortTot[which(is.na(BAmortTot))]=0
    BAdec = apply(modOut$multiOut,1:2,calBAspec,speciesID=3)
    BAdecMort = apply(modOut$multiOut,1:2,calBADspec,speciesID=3); BAdecMort[which(is.na(BAdecMort))]=0
    ###Computing the indices
    HSIcaperRun <- matrix(mapply(HSIcaper,Vpine,Vspruce,Nstems),modOut$nSites,modOut$nYears)
    HSIhgRun <- matrix(mapply(HSIhg,ageoldest,Pdec,Pspruce),modOut$nSites,modOut$nYears)
    #HSIttwoRun <- matrix(mapply(HSIttwo,BAtot,Vtot),modOut$nSites,modOut$nYears)
    HSIttwoRun <- matrix(mapply(HSIttwo,BAmortTot,Vtot),modOut$nSites,modOut$nYears)
    #HSIlswoRun <- matrix(mapply(HSIlswo,BAdec,ageoldest),modOut$nSites,modOut$nYears)
    HSIlswoRun <- matrix(mapply(HSIlswo,BAdecMort,ageoldest),modOut$nSites,modOut$nYears)
    HSIlttRun <- matrix(mapply(HSIltt,ageoldest,BAtot,Pdec),modOut$nSites,modOut$nYears)
    HSIfsRun <- matrix(mapply(HSIfs,Vspruce,Pspruce,Vdec),modOut$nSites,modOut$nYears)
    #resource availability is not still implemented
    
  }else if(class(modOut) == "prebas"){
    ##calculate inputs to the indices
    Vpine = apply(PREBASout$output,1,calVspec,speciesID=1)
    Vspruce = apply(PREBASout$output,1,calVspec,speciesID=2)
    Vdec = apply(PREBASout$output,1,calVspec,speciesID=3)
    VdecDead = apply(modOut$multiOut,1:2,calVdeadSpec,speciesID=3)
    Vtot = apply(PREBASout$output[,30,,1],1,sum,na.rm=T)
    VtotDead= apply(modOut$multiOut[,,8,,1],1:2,sum,na.rm=T)
    Ppine = Vpine / Vtot *100 
    Pspruce = Vspruce / Vtot *100
    Pdec = Vdec / Vtot *100
    ageoldest = apply(PREBASout$output[,7,,1],1,max,na.rm=T)
    Ppine[which(is.na(Ppine))] <- 0
    Pspruce[which(is.na(Pspruce))] <- 0
    Pdec[which(is.na(Pdec))] <- 0
    Nstems = apply(PREBASout$output[,17,,1],1,sum,na.rm=T)
    BAtot = apply(PREBASout$output[,13,,1],1,sum,na.rm=T)
    BAmortTot = (Vtot/VtotDead)*BAtot; BAmortTot[which(is.na(BAmortTot))]=0
    BAdec = apply(PREBASout$output,1,calBAspec,speciesID=3)
    BAdecMort = apply(modOut$multiOut,1:2,calBADspec,speciesID=3); BAdecMort[which(is.na(BAdecMort))]=0
    HSIcaperRun <- mapply(HSIcaper,Vpine,Vspruce,Nstems)
    HSIhgRun <- mapply(HSIhg,ageoldest,Pdec,Pspruce)
    #HSIttwoRun <- mapply(HSIttwo,BAtot,Vtot)
    HSIttwoRun <- matrix(mapply(HSIttwo,BAmortTot,Vtot),modOut$nSites,modOut$nYears)
    #HSIlswoRun <- mapply(HSIlswo,BAdec,ageoldest)
    HSIlswoRun <- matrix(mapply(HSIlswo,BAdecMort,ageoldest),modOut$nSites,modOut$nYears)
    HSIlttRun <- mapply(HSIltt,ageoldest,BAtot,Pdec)
    HSIfsRun <- mapply(HSIfs,Vspruce,Pspruce,Vdec)
    #resource availability is not still implemented
  }
  return(list(HSIcaper=HSIcaperRun,HSIhg=HSIhgRun,HSIttwo=HSIttwoRun,HSIlswo=HSIlswoRun,HSIltt=HSIlttRun,HSIfs=HSIfsRun))
}
## TODO: Add a second way, by option in the funcion calBioIndices(option), 0 = default (BA deciduous recently died), 1= total BA recently died
