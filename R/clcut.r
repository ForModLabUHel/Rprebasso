ClCutD_Pine <- function(ETSmean,ETSthres,siteType){   #TEST DELETE THIS COMMENT LAuri
    if(siteType<=3 & ETSmean>=ETSthres) inDclct <- ClCut_pine[1,1]
    if(siteType==4 & ETSmean>=ETSthres) inDclct <- ClCut_pine[2,1]
    if(siteType>=5 & ETSmean>=ETSthres) inDclct <- ClCut_pine[3,1]
    if(siteType<=3 & ETSmean<ETSthres) inDclct <- ClCut_pine[1,3]
    if(siteType==4 & ETSmean<ETSthres) inDclct <- ClCut_pine[2,3]
    if(siteType>=5 & ETSmean<ETSthres) inDclct <- ClCut_pine[3,3]
    return(as.double(inDclct))
}

ClCutD_Spruce <- function(ETSmean,ETSthres,siteType){
  if(siteType<=2 & ETSmean>=ETSthres) inDclct <- ClCut_spruce[1,1]
  if(siteType>=3 & ETSmean>=ETSthres) inDclct <- ClCut_spruce[2,1]
  if(siteType<=2 & ETSmean<ETSthres) inDclct <- ClCut_spruce[1,3]
  if(siteType>=3 & ETSmean<ETSthres) inDclct <- ClCut_spruce[2,3]
  return(as.double(inDclct))
}

ClCutD_Birch <- function(ETSmean,ETSthres,siteType){
  if(siteType<=2 & ETSmean>=ETSthres) inDclct <- ClCut_birch[1,1]
  if(siteType>=3 & ETSmean>=ETSthres) inDclct <- ClCut_birch[2,1]
  if(siteType<=2 & ETSmean<ETSthres) inDclct <- ClCut_birch[1,3]
  if(siteType>=3 & ETSmean<ETSthres) inDclct <- ClCut_birch[2,3]
  return(as.double(inDclct))
}

ClCutA_Pine <- function(ETSmean,ETSthres,siteType){
    if(siteType<=3 & ETSmean>=ETSthres) inAclct <- ClCut_pine[1,2]
    if(siteType==4 & ETSmean>=ETSthres) inAclct <- ClCut_pine[2,2]
    if(siteType>=5 & ETSmean>=ETSthres) inAclct <- ClCut_pine[3,2]
    if(siteType<=3 & ETSmean<ETSthres) inAclct <- ClCut_pine[1,4]
    if(siteType==4 & ETSmean<ETSthres) inAclct <- ClCut_pine[2,4]
    if(siteType>=5 & ETSmean<ETSthres) inAclct <- ClCut_pine[3,4]
    return(as.double(inAclct))
}

ClCutA_Spruce <- function(ETSmean,ETSthres,siteType){
  if(siteType<=2 & ETSmean>=ETSthres) inAclct <- ClCut_spruce[1,2]
  if(siteType>=3 & ETSmean>=ETSthres) inAclct <- ClCut_spruce[2,2]
  if(siteType<=2 & ETSmean<ETSthres) inAclct <- ClCut_spruce[1,4]
  if(siteType>=3 & ETSmean<ETSthres) inAclct <- ClCut_spruce[2,4]
  return(as.double(inAclct))
}

ClCutA_Birch <- function(ETSmean,ETSthres,siteType){
  if(siteType<=2 & ETSmean>=ETSthres) inAclct <- ClCut_birch[1,2]
  if(siteType>=3 & ETSmean>=ETSthres) inAclct <- ClCut_birch[2,2]
  if(siteType<=2 & ETSmean<ETSthres) inAclct <- ClCut_birch[1,4]
  if(siteType>=3 & ETSmean<ETSthres) inAclct <- ClCut_birch[2,4]
  return(as.double(inAclct))
}


