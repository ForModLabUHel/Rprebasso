#' Habitat suitability for the capercaillie
#' @Description 
#'  Biodiversity indicator for habitat suitability of Capercaillie
#' 
#' @param Vpine volume (m3/ha) of Pine 
#' @param Vspruce volume (m3/ha) of Spruce 
#' @param Nstems total number of trees per hectar
#'
#' @return the habitat suitability index
#' @export
#'
#' @examples
HSIcaper <- function(Vpine,Vspruce,Nstems){
  
  if(Vpine <= 60) Wpine= 0
  if(Vpine > 60 & Vpine <= 80) Wpine = 0.05 * Vpine -3
  if(Vpine > 80) Wpine = 1

  if(Vspruce <= 5) Wspruce= 0
  if(Vspruce > 5 & Vspruce <= 10) Wspruce= 0.2 * Vspruce - 1
  if(Vspruce > 10 & Vspruce <= 20) Wspruce= 1
  if(Vspruce > 20 & Vspruce <= 30) Wspruce= -0.1 * Vspruce + 3
  if(Vspruce > 30) Wspruce = 0
  
  if(Nstems <= 500) Wdensity= 0
  if(Nstems > 500 & Nstems <= 600) Wdensity= 0.01 * Nstems - 5
  if(Nstems > 600 & Nstems <= 800) Wdensity= 1
  if(Nstems > 800 & Nstems <= 1000) Wdensity= -0.005 * Nstems + 5
  if(Nstems > 1000) Wdensity = 0
  
  HSIcap <- Wpine * Wspruce * Wdensity
  return(HSIcap)
}