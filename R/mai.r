
##function to compute mai from model output
mai <- function(xx){
  end <- which(xx[,2,1]==0)[1]-1
  rot <- xx[end,1,1]
  mai=(xx[end,2,1] + sum(xx[2:end,2,2]))/rot
  return(as.double(mai))
}

##function to compute mai from model output + predictive uncertainty
maiPU <- function(xx,errP){
  end <- which(xx[,2,1]==0)[1]-1
  rot <- xx[end,1,1]
  totV <- xx[end,2,1] + sum(xx[2:end,2,2])
  errV <- rnorm(1,0,(errP[1] + errP[2]*totV))
  mai=(totV + errV)/rot
  return(as.double(mai))
}



