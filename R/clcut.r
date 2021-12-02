ClCutD_Pine <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_pine){   #TEST DELETE THIS COMMENT LAuri
    if(siteType<=3 & ETSmean>=ETSthres) inDclct <- pClcut[1,1]
    if(siteType==4 & ETSmean>=ETSthres) inDclct <- pClcut[2,1]
    if(siteType>=5 & ETSmean>=ETSthres) inDclct <- pClcut[3,1]
    if(siteType<=3 & ETSmean<ETSthres) inDclct <- pClcut[1,3]
    if(siteType==4 & ETSmean<ETSthres) inDclct <- pClcut[2,3]
    if(siteType>=5 & ETSmean<ETSthres) inDclct <- pClcut[3,3]
    return(as.double(inDclct))
}

ClCutD_Spruce <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_spruce){
  if(siteType<=2 & ETSmean>=ETSthres) inDclct <- pClcut[1,1]
  if(siteType>=3 & ETSmean>=ETSthres) inDclct <- pClcut[2,1]
  if(siteType<=2 & ETSmean<ETSthres) inDclct <- pClcut[1,3]
  if(siteType>=3 & ETSmean<ETSthres) inDclct <- pClcut[2,3]
  return(as.double(inDclct))
}

ClCutD_Birch <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_birch){
  if(siteType<=2 & ETSmean>=ETSthres) inDclct <- pClcut[1,1]
  if(siteType>=3 & ETSmean>=ETSthres) inDclct <- pClcut[2,1]
  if(siteType<=2 & ETSmean<ETSthres) inDclct <- pClcut[1,3]
  if(siteType>=3 & ETSmean<ETSthres) inDclct <- pClcut[2,3]
  return(as.double(inDclct))
}

ClCutA_Pine <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_pine){
    if(siteType<=3 & ETSmean>=ETSthres) inAclct <- pClcut[1,2]
    if(siteType==4 & ETSmean>=ETSthres) inAclct <- pClcut[2,2]
    if(siteType>=5 & ETSmean>=ETSthres) inAclct <- pClcut[3,2]
    if(siteType<=3 & ETSmean<ETSthres) inAclct <- pClcut[1,4]
    if(siteType==4 & ETSmean<ETSthres) inAclct <- pClcut[2,4]
    if(siteType>=5 & ETSmean<ETSthres) inAclct <- pClcut[3,4]
    return(as.double(inAclct))
}

ClCutA_Spruce <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_spruce){
  if(siteType<=2 & ETSmean>=ETSthres) inAclct <- pClcut[1,2]
  if(siteType>=3 & ETSmean>=ETSthres) inAclct <- pClcut[2,2]
  if(siteType<=2 & ETSmean<ETSthres) inAclct <- pClcut[1,4]
  if(siteType>=3 & ETSmean<ETSthres) inAclct <- pClcut[2,4]
  return(as.double(inAclct))
}

ClCutA_Birch <- function(ETSmean,ETSthres,siteType,pClcut=ClCut_birch){
  if(siteType<=2 & ETSmean>=ETSthres) inAclct <- pClcut[1,2]
  if(siteType>=3 & ETSmean>=ETSthres) inAclct <- pClcut[2,2]
  if(siteType<=2 & ETSmean<ETSthres) inAclct <- pClcut[1,4]
  if(siteType>=3 & ETSmean<ETSthres) inAclct <- pClcut[2,4]
  return(as.double(inAclct))
}


