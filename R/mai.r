#' compMAI
#'
#' @param modRun multiSite PREBAS output
#' @param 
#' #'
#' @return mean annual increment of the stand (standMAI) and bylayer (layerMAI)
#' 
#' @export
#'
#' @examples prebOut <- TransectRun(species="Mixed"); compMAI(prebOut)
compMAI <- function(modRun){
  if(class(modRun)=="multiPrebas"){
    
    nSites <- modRun$nSites
    nLayers <- modRun$maxNlayers
    standMAI <- rep(NA,nSites)
    layerMAI <- matrix(NA,nSites,nLayers)
    ##calculate standing Volume for each site and each year (sum of layers)
    standV <- apply(modRun$multiOut[,,30,,1],1:2,sum)
    
    ##calculate harvested Volume for each site and each year (sum of layers)
    harvV <- apply(modRun$multiOut[,,37,,1],1:2,sum)
    
    ###identify the year of Clcut: standV is 0 and harvested Volume is > 0
    yearClCt <- which(standV==0. & harvV>0,arr.ind=T)
    
    for(i in 1:nSites){
      yearClCtX <- yearClCut[which(yearClCt[,1]==i),2][1]
      ageClCtX <- modRun$multiOut[i,(yearClCtX-1),7,,1] + 1
      maiX <- modRun$multiOut[i,yearClCtX,37,,1]/ageClCtX
      layerMAI[i,] <- maiX
      standMAI[i] <- sum(maiX,na.rm = T)  
    }  
  }
return(list(standMAI=standMAI,layerMAI=layerMAI))  
}
# 
# 
# ##function to compute mai from model output
# mai <- function(xx){
#   end <- which(xx[,2,1]==0)[1]-1
#   rot <- xx[end,1,1]
#   mai=(xx[end,2,1] + sum(xx[2:end,2,2]))/rot
#   return(as.double(mai))
# }
# 
# ##function to compute mai from model output + predictive uncertainty
# maiPU <- function(xx,errP){
#   end <- which(xx[,2,1]==0)[1]-1
#   rot <- xx[end,1,1]
#   totV <- xx[end,2,1] + sum(xx[2:end,2,2])
#   errV <- rnorm(1,0,(errP[1] + errP[2]*totV))
#   mai=(totV + errV)/rot
#   return(as.double(mai))
# }
# 
# 
# 
