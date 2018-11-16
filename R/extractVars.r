extractVars <- function(x,variables=NA,sites=NA,speciesIDs=NA,
                        speciesNam = NA,years=NA){

  varNam <- getVarNam()
  if (anyNA(speciesNam)) speciesNam <- as.character(paste("sp",1:nSp))
  if (any(variables == "all") | anyNA(variables)) variables <-c(5,6,8:18,22,24:34,37:46)
  variables <- c(7,variables)

  if(inherits(x,"prebas")){
    if(anyNA(speciesIDs)) speciesIDs <- 1:dim(x$output)[3]
    nSp <- length(speciesIDs)
    if (anyNA(years)) years <- 1:(dim(x$output)[1])

    out <- x$output[years,variables,speciesIDs,1]
    dimnames(out) <- list(NULL,varNam,speciesNam)

  }


  if(inherits(x,"multiPrebas")){
    if(anyNA(speciesIDs)) speciesIDs <- 1:dim(x$multiOut)[4]
    nSp <- length(speciesIDs)
    if (anyNA(years)) years <- 1:(dim(x$output)[2])
    if (anyNA(sites)) sites <- 1:(dim(x$multiOut)[1])

    out <- x$multiOut[sites,years,variables,speciesIDs,1]
    dimnames(out) <- list(x$multiOut[sites,1,1,1,1],NULL,varNam[variables],speciesNam[speciesIDs])

  }
return(out)
}
