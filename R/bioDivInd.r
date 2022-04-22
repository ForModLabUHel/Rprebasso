library(Rprebasso)
ciao1 <- TransectRun(species = "Pine",mortMod=1)
ciao2 <- TransectRun(species = "Pine",mortMod=2)
ciao3 <- TransectRun(species = "Pine",mortMod=3)

siteX=2
varX=17
spX=1

ylim=range(ciao1$multiOut[siteX,,varX,spX,1],
           ciao2$multiOut[siteX,,varX,spX,1],
           ciao3$multiOut[siteX,,varX,spX,1])
plot(ciao1$multiOut[siteX,,varX,spX,1],ylim=ylim)
points(ciao2$multiOut[siteX,,varX,spX,1],col=2,pch=20)
points(ciao3$multiOut[siteX,,varX,spX,1],col=3,pch=20)



x1=300
y1=0.18
x2=1300
y2=0.05
b<- (y2/y1)^(1/(x2-x1))
a=y1/(b^x1)
x=1500
a*b^x

test <- list()
N=5000
ba=100
pMort=0.
perBAmort=0.

test$ba <- ba
for(i in 1:5){
  test <- .Fortran("randMort",N=N,ba=test$ba,
                   pMort=pMort,perBAmort=perBAmort)
  print(test$ba)
}
test#$ba/ba

library(data.table)
load("C:/Users/checcomi/Documents/research/IBC-carbon/test/data.all_maakunta_5.rdata")
data.all[,totCover:=(pine + spruce+decid)]
age = data.all$age
ba=data.all$ba
d=data.all$dbh
N=ba/(pi*(d/200)^2)
rBApine = data.all$pine/data.all$totCover
rBAbrd = data.all$decid/data.all$totCover
slope = 0.048
pSize = 314.16
pMort=rep(0,nrow(data.all))
perBAmort= BAmort= rep(0,nrow(data.all))
step=5
siteX <- union(union(which(is.na(N)),which(ba==0)),which(d==0))
sites <- (1:nrow(data.all))[-siteX]
for(i in sites){
  ciao <- .Fortran("randMort",age=as.double(age[i]),
                   d=as.double(d[i]),
                   ba=as.double(ba[i]),
                   N=as.double(N[i]),
                   rBApine=as.double(rBApine[i]),
                   rBAbrd=as.double(rBAbrd[i]),
                   slope=as.double(slope),
                   pSize=as.double(pSize),
                   pMort=as.double(0),
                   perBAmort=as.double(0),
                   step=as.double(step),
                   BAmort=as.double(0))
  pMort[i] <- ciao$pMort
  perBAmort[i] <- ciao$perBAmort
  BAmort[i] <- ciao$BAmort
}

plot(age[sites],pMort[sites],pch='.')
lines(spline(age[sites], pMort[sites],n=20),col=2)
plot(d[sites],pMort[sites],pch='.')
lines(spline(d[sites], pMort[sites],n=20),col=2)
plot(N[sites],pMort[sites],pch='.',xlim=c(0,2000))
lines(spline(N[sites], pMort[sites],n=102000),col=2)
plot(ba[sites],pMort[sites],pch='.')
lines(spline(ba[sites], pMort[sites],n=20),col=2)
abline(h=c(0.1,.2,.3))

plot(ba[sites],perBAmort[sites],pch='.')
plot(N[sites],perBAmort[sites],pch='.',xlim=c(0,2000))
plot(d[sites],perBAmort[sites],pch='.')


plot(ba[sites],BAmort[sites],pch='.')
abline(h=0.7)
plot(N[sites],BAmort[sites],pch='.',xlim=c(0,2000))
plot(d[sites],BAmort[sites],pch='.')







###average condition
age = 80 #84.1
ba=21.7 #25.2
d=23.6 #23.3
N=685.8 #826.4
rBApine = 0.5 #0.38
rBAbrd = 0.2 #0.2
slope = 0.048 #0.05
pSize = 314.16
perBAmort=0
pMort=0.
step=5
# sites <- (1:nrow(data.all))[-which(is.na(N))]
# for(i in sites[(1:10000)]){
  ciao <- .Fortran("randMort",age=as.double(age),
                   d=as.double(d),
                   ba=as.double(ba),
                   N=as.double(N),
                   rBApine=as.double(rBApine),
                   rBAbrd=as.double(rBAbrd),
                   slope=as.double(slope),
                   pSize=as.double(pSize),
                   pMort=as.double(0),
                   perBAmort=as.double(0),
                   step=as.double(step))
  ciao$pMort
  ciao$perBAmort
   ciao$perBAmort * ba
  