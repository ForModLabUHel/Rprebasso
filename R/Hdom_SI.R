#' Estimate dominant height based on average height
#'
#' @param H.av Average height (m)
#' @param pCROBAS The matrix of CROBAS parameters
#' @param SpeciesNo The column number of the given species in pCROB
#'
#' @return The dominant height (m)
#' @export
#'
#' @examples Hdom.fun(H.av=2:50,SpeciesNo=1)
Hdom.fun<-function(H.av=NA,SpeciesNo=NA,pCROBAS=pCROB){
  Hdominant<-pCROBAS[42,SpeciesNo]*exp(-1/pmax((H.av-1.3),0.001))+pCROBAS[43,SpeciesNo]*H.av
  return(Hdominant)
}
