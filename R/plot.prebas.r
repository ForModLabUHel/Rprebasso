plot.prebas <- function(x,variableIDs=NA,siteIDs=NA,layerIDs=NA,leg=T,
                        layerNam = NA,obsData=NA){

  varNam <- getVarNam()
  if(all(is.na(obsData))) obsData <- matrix(NA,2,10)
  if (any(variableIDs == "all") | anyNA(variableIDs)) variableIDs <-c(5,6,8:18,22,24:34,37:46)

  if(inherits(x,"prebas")){
    if(anyNA(layerIDs)) layerIDs <- 1:dim(x$output)[3]
    nLayers <- length(layerIDs)
    if (anyNA(layerNam)) layerNam <- as.character(paste("layer",1:nLayers))

    count <- 0
    if (length(variableIDs) > 1) par(mfrow=c(2,3)) else par(mfrow=c(1,1))
    for(vars in variableIDs){
      count <- count + 1
      plot(x$output[,vars,layerIDs[1],1],type='l',xaxt='n',
       main=varNam[vars],ylab = "units",xlab="age (y)",col=layerIDs[1],
       ylim=c(min(x$output[,vars,,1]),max(x$output[,vars,,1])))
      selObs <- which(obsData[,4]==vars & obsData[,2]==layerIDs[1])
      if(length(selObs)>0) points(obsData[selObs,3],obsData[(selObs),5],col=layerIDs[1])

      if (nLayers>1) for(ij in layerIDs[2:nLayers]){
        lines(x$output[,vars,ij,1],col=ij)
        selObs <- which(obsData[,4]==vars & obsData[,2]==ij)
        if(length(selObs)>0) points(obsData[selObs,3],obsData[(selObs),5],col=ij)
      }# points(data)
      axis(1, at=seq(1,(dim(x$output)[1]),length.out=6), labels=x$output[seq(1,(dim(x$output)[1]),length.out=6),7,1,1])
      if (leg==TRUE) legend("topleft",c(layerNam[layerIDs]),lty=1, col=layerIDs)
      if (count %% 6 == 0 & vars!=tail(variableIDs,n=1)) pause()
    }
  }


  if(inherits(x,"multiPrebas")){
    if(anyNA(layerIDs)) layerIDs <- 1:dim(x$multiOut)[4]
    nLayers <- length(layerIDs)
    if (anyNA(siteIDs)) siteIDs <- 1:dim(x$multiOut)[1]

    for(iz in siteIDs){
      if (anyNA(layerNam[iz])) layerNam <- as.character(paste("layer",1:x$nLayers[iz]))
      count <- 0
      if (length(variableIDs) > 1) par(mfrow=c(2,3)) else par(mfrow=c(1,1))
      for(vars in variableIDs){

        plot(x$multiOut[iz,,vars,layerIDs[1],1],type='l',xaxt='n',
             main=varNam[vars],ylab = "units",xlab="age (y)",col=layerIDs[1],
             ylim=c(min(x$multiOut[iz,,vars,,1]),max(x$multiOut[iz,,vars,,1])))
        if (nLayers>1) for(ij in layerIDs[2:nLayers]) lines(x$multiOut[iz,,vars,ij,1],col=ij)
        axis(1, at=seq(1,(dim(x$multiOut)[2]),length.out=6), labels=x$multiOut[iz,seq(1,(dim(x$multiOut)[2]),length.out=6),7,1,1])
        if (leg==TRUE) legend("topleft",c(layerNam[layerIDs]),lty=1,col=layerIDs)
        if(count%%6==0) title(paste('Site:', x$multiOut[iz,1,1,1,1]), line = -18, outer = TRUE,cex.main=2)
        count <- count + 1
        if (count %% 6 == 0 & vars!=tail(variableIDs,n=1)) pause()
      }
    if (length(siteIDs)>1 & iz != tail(siteIDs,n=1)) pause()}
  }
}



# title_plot <- "NRMSE for the PGE Pine dataset"
# PSPpine_nrmse <- ggplot(data=pspPine_nrmse, aes(x=variableIDs, y=NRMSE, fill=calibration)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   ggtitle(title_plot) + scale_fill_manual(values=gray.colors(2)) + theme_classic() +
#   # ylim(0,10) +
#   theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),
#         axis.text.x=element_text(size=rel(1.7)),axis.title.y=element_text(size=15),
#         axis.text.y=element_text(size=rel(1.5)))
#
#
# dd<-data.frame(matrix(c(x$output[,7,1,1],x$output[,vars,1,1]),
#                       length(x$output[,vars,1,1]),2))
#
#   pg <- ggplot(data = dd, aes(x=X1, y=X2)) + geom_point()
#   # geom_point(data = ddx, aes(x=observed, y=simulated), color = 1,shape=19) +
#   # annotate("text", x=min(c(dd$observed,dd$simulated)), y= max(c(dd$observed,dd$simulated)), label = legendSP, hjust=0) +
#   # ggtitle(paste("PGEcal:",varNam[vars])) +
#   # geom_abline(slope=1, intercept=0) +
#   # scale_x_continuous(limits = c(min(c(dd$observed,dd$simulated)), max(c(simulated,observed)))) +
#   # scale_y_continuous(limits = c(min(c(dd$observed,dd$simulated)), max(c(simulated,observed)))) +
#   # theme(plot.title = element_text(hjust = 0.5))
#
#
#   grid.newpage()
#   pushViewport(viewport(layout = grid.layout(6, 6, heights = unit(c(0.8, 8,0.8,8), "null"))))
#   grid.text("NFI calibration", gp=gpar(cex=1.3),
#             vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
#   print(pg, vp = viewport(layout.pos.row = 2, layout.pos.col = 6))
#
#
#   print(NFIplots[[2]], vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
#   print(NFIplots[[1]], vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
#   grid.text("PGE calibration", gp=gpar(cex=1.3),
#             vp = viewport(layout.pos.row = 3, layout.pos.col = 1:3))
#   print(PGEplots[[3]], vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
#   print(PGEplots[[2]], vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
#   print(PGEplots[[1]], vp = viewport(layout.pos.row = 4, layout.pos.col = 3))
#
