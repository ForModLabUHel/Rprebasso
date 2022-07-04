devtools::install("C:/Users/minunno/Documents/github/Rprebasso/", build_vignettes = FALSE)
devtools::install_github("ForModLabUHel/Rprebasso", ref="master")

library(Rprebasso)
library(ggplot2)
library(ggpubr)

ETS <- 1000
st = 3 #site type
r_r <- 0.53
W_RT <- seq(1,10000,by=200)

#initialize output variables
CN <- 0 #initial C:N
rhoM <- 0 #ratio ECM biomass to fine root biomass
r_RT=rep(0,length(W_RT))
rm_ECM=rep(0,length(W_RT))
litt_RT=rep(0,length(W_RT))
exud=rep(0,length(W_RT))

  
CNratio <- .Fortran("CNratio",CN=as.double(CN),
                    ETS=as.double(ETS),st=as.double(st))
rhoM <- .Fortran("rhoMcalc",rhoM=as.double(rhoM),
                 CN=as.double(CNratio$CN))

for(i in 1:length(W_RT)){
  CUEcalc <- .Fortran("CUEcalc",
                      ETS=as.double(ETS), 
                      siteType=as.double(st),
                      r_r=as.double(r_r),
                      W_RT=as.double(W_RT[i]),
                      r_RT=as.double(0),
                      rm_ECM=as.double(0),
                      litt_RT=as.double(0),
                      exud=as.double(0))
  r_RT[i] <- CUEcalc$r_RT
  rm_ECM[i] <- CUEcalc$rm_ECM
  litt_RT[i] <- CUEcalc$litt_RT
  exud[i] <- CUEcalc$exud
}

p1 = ggplot() + geom_point(aes(x=W_RT,y=r_RT))
p2 = ggplot() + geom_point(aes(x=W_RT,y=rm_ECM))
p3 = ggplot() + geom_point(aes(x=W_RT,y=litt_RT))
p4 = ggplot() + geom_point(aes(x=W_RT,y=exud))

ggarrange(p1,p2,p3,p4)


C:\HYAPP\rtools40v2
Sys.setenv(PATH = paste("C:HYAPP/rtools40v2/usr/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/HYAPP/rtools40v2/mingw$(WIN)/bin/")
                      C:\HYAPP\rtools40v2\mingw64\bin