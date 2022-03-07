#' SW
#'
#' @param a  precipitation
#' @param b  rate
#' @param c  previous water storage
#'
#' @return the soil water
#' @export
#'
#' @examples
xx<-function(a,b,c){
  y<-a*b+c
  return(y)
  }