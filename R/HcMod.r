# Function to Compute Hc based on ksi parameter
ksiHcMod <- function(initVar){
  h <- initVar[3] - 1.3
  b <- pi * initVar[4]^2 / 40000
  p_ksi <- pCROBAS[38,initVar[1]]
  p_rhof <- pCROBAS[15,initVar[1]]
  p_z <- pCROBAS[11,initVar[1]]
  Lc <- h* ((p_rhof * b)/(p_ksi * h^p_z))^(1/(p_z-1))
  Hc <- max(0.,(h-Lc))
  return(Hc)
}

###function to replace HC NAs in initial variable initVar
findHcNAs <- function(initVar,pHcMod,hcFactor=1.){
  if(is.vector(initVar)){
    if(is.na(initVar[6])){
        Ntot <- initVar[5]/(pi*((initVar[4]/2/100)**2))
        inModHc <- c(pHcMod[,initVar[1]],initVar[3],
                     initVar[4],initVar[2],initVar[5],initVar[5],Ntot,initVar[1])
        initVar[6] <- model.Hc(inModHc) * hcFactor
    }
  } else if(any(is.na(initVar[6,]))){
    initVar[1,][which(initVar[1,]==0)] <- 1 ###deals with 0 species ID
    HcNAs <- which(is.na(initVar[6,]))
    BAtot <- sum(initVar[5,],na.rm = T)
    Ntot <- sum(initVar[5,]/(pi*((initVar[4,]/2/100)**2)),na.rm = T)
    if(length(HcNAs)==1){
        inModHc <- c(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                     initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],
                     BAtot,Ntot,initVar[1,HcNAs])
        initVar[6,HcNAs] <- model.Hc(inModHc) * hcFactor
    }else{
        inModHc <- rbind(pHcMod[,initVar[1,HcNAs]],initVar[3,HcNAs],
                         initVar[4,HcNAs],initVar[2,HcNAs],initVar[5,HcNAs],
                         BAtot,Ntot,initVar[1,HcNAs])
        initVar[6,HcNAs] <- apply(inModHc,2,model.Hc) * hcFactor
    }
  }
  return(initVar)
}

model.Hc <- function(inputs){ 
  spID = inputs[13]
  Hc = HcModDef[[spID]](inputs)
  return(Hc)
}


HcModDef <- list()
###default HcModel for pisy, piab, beal and pipi
HcModDef[[1]] <- HcModDef[[2]] <-HcModDef[[3]] <- HcModDef[[5]]<- function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
    lnHc_sim <- pValues[1]+pValues[2]*log(H)+pValues[3]*D/H+
    pValues[4]*log(age)+ pValues[5]*log(BA_sp)+
    pValues[6]*(BA_sp/BA_tot)
  Hc_sim <- exp(lnHc_sim)
  return(pmax(Hc_sim,0.,na.rm = T)) 
} 

###default HCmodel for fasy
HcModDef[[4]] <-function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
  N_tot = inputs[12]
  Hc_sim <-  H/(1+exp(pValues[1]-pValues[1]*N_tot/1000-
                        pValues[1]*log(BA_tot)+pValues[1]*D/H-pValues[1]*log(H)))
  return(pmax(Hc_sim,0.,na.rm = T)) 
} 

###default HCmodel for eugl
HcModDef[[6]] <-function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
  N_tot = inputs[12]
  Hc_sim <- H*exp(pValues[1]+pValues[2]*N_tot/1000+
              pValues[3]*log(BA_tot)+pValues[4]*D/H +pValues[5]*H)
  return(pmax(Hc_sim,0.,na.rm = T)) 
} 

###default HCmodel for rops
HcModDef[[7]] <-function(inputs){ 
  pValues=inputs[1:6]
  H=inputs[7]
  D=inputs[8]
  age=inputs[9]
  BA_sp=inputs[10]
  BA_tot=inputs[11]
  Hc_sim <- H/(1+exp(pValues[1]-pValues[2]*log(BA_tot)+pValues[3]*D/H))
  return(pmax(Hc_sim,0.,na.rm = T)) 
} 



# ##Height of the crown base model
# model.Hc <- function(inputs){ 
#   pValues=inputs[1:6]
#   H=inputs[7]
#   D=inputs[8]
#   age=inputs[9]
#   BA_sp=inputs[10]
#   BA_tot=inputs[11]
#   N = 
#   lnHc_sim <- pValues[1]+pValues[2]*log(H)+pValues[3]*D/H+
#     pValues[4]*log(age)+ pValues[5]*log(BA_sp)+
#     pValues[6]*(BA_sp/BA_tot)
#   Hc_sim <- exp(lnHc_sim)
#   return(pmax(Hc_sim,0.,na.rm = T)) 
# } 