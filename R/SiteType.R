#' Estimate site type based on age, height, and climate.
#'
#' @param Age Age of dominant trees or mean trees (yr). 
#' @param Hdom Dominant height (m).
#' @param Hmean Average height (m). If dominant height is not given, the dominant height will be calculated based on average height.
#' @param Species Species, ranging from 1 to 10, according to the column of pCROB.
#' @param PAR A matrix of daily PPFD, nrow=5, ncol=365*years.
#' @param TAir A matrix of daily air temperature, nrow=5, ncol=365*years.
#' @param VPD A matrix of daily VPD, nrow=5, ncol=365*years.
#' @param Precip A matrix of daily precipitation, nrow=5, ncol=365*years.
#' @param CO2 A matrix of daily CO2 concentration, nrow=5, ncol=365*years.
#' @param pCROBAS Parameters of CROBAS.
#' @param Age.1.5m Age of seedlings when the height reach 1.5m.
#'
#' @return Site type
#' @export
#'
#' @examples SiteType(Age=c(50,70),Hdom=c(25,30),Species=3,
#'      PAR = PAR,TAir= TAir,VPD= VPD,Precip= Precip,CO2= CO2,
#'      pCROBAS=pCROB,Age.1.5m=5)
#'
SiteType<-function(Age=NA,Hdom=NA,Hmean=NA,Species=NA,PAR=NA,TAir=NA,VPD=NA,Precip=NA,CO2=NA,
                   pCROBAS=pCROB,Age.1.5m=5){
  nyears<-max(Age) # Years of the simulation
  
  # Use measurements of average heights if dominant height is not available
  if(any(is.na(Hdom)) & all(!is.na(Hmean))){
    H=Hmean
  }else{
    H=(Hdom-pCROB['aHdom',Species])/pCROB['bHdom',Species]
  }
  
  nSites<-5
  siteInfo<- data.frame(siteID=c(1:nSites),
                        climID=c(1:nSites),
                        siteType=1:5,
                        SWinit=rep(200,nSites),
                        CWinit=rep(0,nSites),
                        SOGinit=rep(0,nSites),
                        Sinit =rep(10,nSites),
                        nLayers =rep(1,nSites),
                        nSpecies=rep(1,nSites),
                        Dsoil=rep(413,nSites),
                        FC=rep(0.450,nSites),
                        WP=rep(0.118,nSites) )
  
  multiInitVar <- aperm(array(c(Species, Age.1.5m, initSeedling.def), dim = c(7, 5, 1)), c(2, 1, 3))
  
  pPRELES <- switch(
    Species,
    pPREL, # "pisy"
    pPREL, # "piab" 
    pPREL, # "beal"
    pPRELESfasy, #  "fasy" 
    pPRELESpipi, #  "pipi" 
    pPRELESeugl, # "eugl"
    pPRELES.Df.DBF, # "rops"
    pPRELES.Df.DBF, # "popu"
    pPRELESeugl, # "eugrur" 
    pPRELESpiabDE # "piab(DE)"
    )

  initPrebas <- InitMultiSite(nYearsMS = rep(nyears,5),
                              siteInfo=as.matrix(siteInfo),
                              pCROBAS = pCROBAS,
                              pPRELES=pPRELES,
                              defaultThin=0.,
                              ClCut = 0.,
                              multiInitVar = multiInitVar,
                              PAR = PAR,
                              TAir= TAir,
                              VPD= VPD,
                              Precip= Precip,
                              CO2= CO2)
  
  output<- multiPrebas(initPrebas)$multiOut
  
  pl<-list()
  for (site in 1:5) {
    pl[[site]]<-data.frame(A=output[site,,7,1,1],
                           H=output[site,,11,1,1],
                           site.type=site
    )
  }
  pldata<-rbind(pl[[1]],pl[[2]],pl[[3]],pl[[4]],pl[[5]])
  
  Site.Type<-rep(NA,length(Age))
  for (i in 1:length(Site.Type)) {
    age.i<-which(pldata$A==Age[i])
    Site.Type[i]<-pldata$site.type[age.i][which.min(abs(H[i]-pldata$H[age.i]))]
  }
  
  return(Site.Type)
}
