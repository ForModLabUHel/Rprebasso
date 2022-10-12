#List of bioindicators in this file:
# 1.HSI capercaillie 2. HSI hazel grouse 3. HSI Three-toed woodpecker 
# 4. HSI Lesser-spotted woodpecker 5. HSI Long-tailed tit 
# 6. HSI Siberian flying squirrel 7. Resource availability (not yet available)
#List of other function in this file
# 1. calVspec 2. calBAspec 3. calVdeadSpec 4. calBADspec

#' Habitat suitability for the Capercaillie
#' @Description 
#'  Biodiversity indicator for habitat suitability of Capercaillie  (Tetrao uralensis) 
#'
#' @param Vpine volume (m3/ha) of Pine
#' @param Vspruce volume (m3/ha) of Spruce
#' @param Nstems total number of trees per hectare
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIcaper <- function(Vpine,Vspruce,Nstems){
  #wpine
  if (Vpine <= 60) Wpine = 0
  else if (Vpine <= 80) Wpine = 0.05 * Vpine -3
  else if (Vpine > 80) Wpine = 1
  #wspruce
  if (Vspruce <= 5) Wspruce= 0
  else if (Vspruce <= 10) Wspruce = 0.2 * Vspruce - 1
  else if (Vspruce <= 20) Wspruce = 1
  else if (Vspruce <= 30) Wspruce = -0.1 * Vspruce + 3
  else if (Vspruce > 30) Wspruce = 0
  #Nstems
  if (Nstems <= 500) Wdensity= 0
  else if (Nstems <= 600) Wdensity = 0.01 * Nstems - 5
  else if (Nstems <= 800) Wdensity = 1
  else if (Nstems <= 1000) Wdensity = -0.005 * Nstems + 5
  else if (Nstems > 1000) Wdensity = 0
  #eq
  HSIcap <- Wpine * Wspruce * Wdensity
  return(HSIcap)
}

#' Habitat suitability for the Hazel grouse
#' @Description 
#'  Biodiversity indicator for habitat suitability of The Hazel grouse (Bonasa bonasia)
#'
#' @param ageoldest Age of the oldest layer 
#' @param Pdec Proportion of deciduous tress (percentage) of the total tree volume
#' @param Pspruce Proportion of spruce (percentage) of the total tree volume
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIhg <- function(ageoldest,Pdec,Pspruce){
  #wage
  if (ageoldest <= 20) Wage = 0
  else if (ageoldest <= 30) Wage = 0.1*ageoldest-2
  else if (ageoldest <= 60) Wage = 1
  else if (ageoldest > 60) Wage = -0.012*ageoldest+1.72
  #wdec
  if (is.nan(Pdec)) Wdec = 0
  else if (Pdec <= 5) Wdec = 0
  else if (Pdec <= 20) Wdec = 0.066*Pdec-0.33
  else if (Pdec <= 40) Wdec = 1
  else if (Pdec <= 60) Wdec = -0.05*Pdec+3
  else if (Pdec >60) Wdec = 0
  #wspruce
  if (is.nan(Pspruce)) Wspruce = 0
  else if (Pspruce <= 20) Wspruce = 0
  else if (Pspruce <= 25) Wspruce = 0.2*Pspruce-4
  else if (Pspruce >25) Wspruce = 1
  #eq
  HSIhgint = Wage*Wdec*Wspruce
  return(HSIhgint)
}

#' Habitat suitability for the Three-toed woodpecker
#' @Description 
#'  Biodiversity indicator for habitat suitability of Three-toed woodpecker (Picoides tridactylus)
#'
#' @param BArecd Total basal area (m2/ha) of recently died trees 
#' @param Vtotal Total timber volume of living trees
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIttwo <- function(BArecd,Vtotal) {
  #wdw
  Wdw = 1 / (1 + exp(-(3.55*BArecd-4.46)))
  #wage
  if (Vtotal < 60) Wvol = 0
  else if (Vtotal <= 200) Wvol = Vtotal/200
  else if (Vtotal > 200) Wvol = 1
  #eq
  HSIttwoint = Wdw*Wvol
  return(HSIttwoint)
}

#' Habitat suitability for the Lesser-spotted woodpecker
#' @Description 
#'  Biodiversity indicator for habitat suitability of Lesser-spotted woodpecker (Dendrocopos minor) 
#'
#' @param BArecddec Total basal area (m2/ha) of recently died deciduous trees
#' @param ageoldest Age of the oldest layer
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIlswo <- function(BArecddec,ageoldest) {
  #wdw
  Wdw = 1 / (1 + exp(-(6.32*BArecddec-2.96)))
  #wage
  if (ageoldest < 60) Wage = 0
  else if (ageoldest <= 200) Wage = ageoldest/200
  else if (ageoldest > 200) Wage = 1
  #eq
  HSIlswoint = Wdw*Wage
  return(HSIlswoint)
}

#' Habitat suitability for the Long-tailed tit
#' @Description 
#'  Biodiversity indicator for habitat suitability of Long-tailed tit (Aegithalos caudatus)   
#'
#' @param ageoldest Age of the oldest layer
#' @param BAtot Total basal area (m2/ha) 
#' @param Pdec Proportion of deciduous tress (percentage) of the total tree volume
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIltt <- function(ageoldest,BAtot,Pdec){
  #age
  if (ageoldest < 30) Wage = 0
  else if (ageoldest <60) Wage = 0.033*ageoldest-1
  else if (ageoldest >= 60) Wage = 1
  #BA
  if (BAtot <= 10) Wba = 0
  else if (BAtot <= 15) Wba = 0.2*BAtot-2
  else if (BAtot > 15) Wba = 1 
  #wdec
  if (is.nan(Pdec)) Wdec = 0
  else if (Pdec <= 20) Wdec = 0
  else if (Pdec <= 60) Wdec = 0.025*Pdec-0.5
  else if (Pdec >60) Wdec = 1
  #eq
  HSIlttint = Wage*Wdec*Wba
  return(HSIlttint)
}

#' Habitat suitability for the Siberian flying squirrel
#' @Description 
#'  Biodiversity indicator for habitat suitability of Long-tailed tit (Aegithalos caudatus)   
#'
#' @param Vspruce volume (m3/ha) of Spruce
#' @param Pspruce Proportion of spruce (percentage) of the total tree volume
#' @param Vdec volume (m3/ha) of deciduous trees
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIfs <- function(Vspruce,Pspruce,Vdec) {
  #wsprucevol
  if (Vspruce <= 140) Wsprucevol = 0
  else if (Vspruce <= 175) Wsprucevol = 0.028*Vspruce-4
  else if (Vspruce > 175) Wsprucevol = 1
  #wsprucep
  if (is.nan(Pspruce)) Wsprucep = 0
  else if (Pspruce <= 50) Wsprucep = 0
  else if (Pspruce <= 60) Wsprucep = 0.1*Pspruce-5
  else if (Pspruce > 60) Wsprucep = 1
  #wdec
  if (Vdec <= 12) Wdec = 0
  else if (Vdec <= 15) Wdec = 0.33*Vdec-4
  else if (Vdec > 15) Wdec = 1
  #eq
  HSIfsint = Wsprucevol*Wsprucep*Wdec
  return(HSIfsint)
}##Add poc from Hurme et al., 2007 ??

#HSIresav <- function(sunny,BA,diameter,s,r) {
#  if (sunny == 1) {
#    Fmicroclima <- -0.07*sqrt(BA)+1
#  } else {
#    Fmicroclima <- 0.15*sqrt(BA)
#  }
#  Fdiameter <- 1.014 / (1+exp(-1*((diameter-13.621)/2.397)^5.056))
#  HSIresavint = Fmicroclima * Fdiameter * s * r
#}


#' calVspec function 
#' @Description 
#'  Computed the Volume of a given specie of tree by it SpecieID
#'
#' @param SpecieID ID number of the layer for the specie
#' @param Prebasoutput route or variable to the data needed
#'
#' @return the volume for the specie given
#' @export
#'
#' @examples
calVspec <- function(prebout,speciesID){
  speciesLoc <- which(prebout[4,,1]==speciesID)
  if (!is.null(speciesLoc)){
    Vtot <- sum(prebout[30,speciesLoc,1],na.rm=T)
  }else {
    Vtot = 0
  }
  return(Vtot)
}

#' calVdeadSpec function 
#' @Description 
#'  Computed the Dead wood volume of a given specie of tree by it SpecieID
#'
#' @param SpecieID ID number of the layer for the specie
#' @param Prebasoutput route or variable to the data needed
#'
#' @return the volume for the specie given
#' @export
#'
#' @examples
calVdeadSpec <- function(prebout,speciesID){
  speciesLoc <- which(prebout[4,,1]==speciesID)
  if (!is.null(speciesLoc)){
    Vtot <- sum(prebout[8,speciesLoc,1],na.rm=T)
  }else {
    Vtot = 0
  }
  return(Vtot)
}

#' calBAspec function 
#' @Description 
#'  Computed the Basal Area of a given specie of tree by it SpecieID
#'
#' @param SpecieID ID number of the layer for the specie
#' @param Prebasoutput route or variable to the data needed
#'
#' @return the Basal area for the Specie given
#' @export
#'
#' @examples
calBAspec <- function(prebout,speciesID){
  speciesLoc <- which(prebout[4,,1]==speciesID)
  if (!is.null(speciesLoc)){
    BAspec <- sum(prebout[13,speciesLoc,1],na.rm=T)
  }else {
    BAspec = 0
  }
  return(BAspec)
}

#' calBADspec function 
#' @Description 
#'  Computed the Basal Area of died trees of a given specie of tree by it SpecieID
#'
#' @param SpecieID ID number of the layer for the specie
#' @param Prebasoutput route or variable to the data needed
#'
#' @return the Basal area for the Specie given
#' @export
#'
#' @examples
calBADspec <- function(prebout,speciesID){
  speciesLoc <- which(prebout[4,,1]==speciesID)
  if (!is.null(speciesLoc)){
    BAspec <- sum(prebout[13,speciesLoc,1],na.rm=T)
    Vspec <- sum(prebout[30,speciesLoc,1],na.rm=T)
    Vdeadspec <- sum(prebout[8,speciesLoc,1],na.rm=T)
  }else {
    BAspec = 0
    Vspec = 1
    Vdeadspec = 1
  }
  BAspecMort = (Vdeadspec/Vspec)*BAspec
  return(BAspecMort)
}
