# # Empirical wind risk model
# # source: Suvanto et al 2019
# # Demo/testing script
# # jh, 1/24
# 
# UNCOMMENT ALL (if active package building fails)
#
# # Test inputs & var explanation
# spec <- 1 # 1 pine, 2 spruce, 3 other
# tsincethin <- 1# time since last thinning in years; ref: 0:5)
# openedge <- 0 # 0 = no open edge, 1 = open edge; ref: 0
# soiltype <- 0 # 0 = mineral, coarse; 1 = mineral, fine; 3 = organic; ref: 0 
# shallowsoil <- 0 # 0 = soil depth >30cm, 1 = <30cm; ref: 0
# sitetype <- 4 # prebas site type 1:5 (converted to site fertility with 1:3 as fertile, 4:5 as infertile (ref)
# h <- 35 # in m, avg: 16.4
# wspeed <- 12.2 # m/s (10a max), avg: 12.2
# tsum <- 1187 # effective temperature sum in degree days (note: 100 dd in Suvanto 2019)
# 
# # single 'site' demo
# ftest <- .Fortran("windrisk",
#                        spec=as.integer(spec),
#                        h=as.double(h),
#                        tsincethin=as.integer(tsincethin),
#                        wspeed=as.double(wspeed),
#                        openedge=as.integer(openedge),
#                        soiltype=as.integer(soiltype),
#                        shallowsoil=as.integer(shallowsoil),
#                        sitetype=as.integer(sitetype),
#                        tsum=as.double(tsum),
#                        wrisk5dd1=as.double(0), # 5a risk for damage density class 1 (0-2)
#                        wrisk5dd2=as.double(0), # 5a risk for damage density class 1 (2-3)
#                        wrisk5dd3=as.double(0), # 5a risk for damage density class 3 (<3)
#                        wrisk0=as.double(0), # pre-logit transformation value 
#                        wrisk5=as.double(0), # 5a weighted average of all damage density classes
#                        wrisk=as.double(0)) # annual risk
# 
# ftest
# 
# 
# # calculate risk for species-specific risk as a function of height
# # with other vars set to reference
# library(data.table)
# library(ggplot2)
# htest <- data.table(h=rep(1:40,3), spec=rep(1:3, each=40))
# 
# for(specx in c(1:3)){
#    for(hx in c(1:40)){
#     ftest <- .Fortran("windrisk",
#                     spec=as.integer(specx),
#                     h=as.double(hx),
#                     tsincethin=as.integer(tsincethin),
#                     wspeed=as.double(wspeed),
#                     openedge=as.integer(openedge),
#                     soiltype=as.integer(soiltype),
#                     shallowsoil=as.integer(shallowsoil),
#                     sitetype=as.integer(sitetype),
#                     tsum=as.double(tsum),
#                     wrisk5dd1=as.double(0),
#                     wrisk5dd2=as.double(0),
#                     wrisk5dd3=as.double(0),
#                     wrisk0=as.double(0),
#                     wrisk5=as.double(0),
#                     wrisk=as.double(0))
# 
# htest[h==hx & spec==specx, wrisk:=ftest$wrisk5]
#  }
# }
# 
# htest
# ggplot(data=htest[h<=35,], aes(x=h, y=wrisk, col=as.factor(spec)))+
#   geom_line()+
#   ggtitle("5a wind disturbance risk")
# 
# ggplot(data=htest[h<=35,], aes(x=h, y=wrisk/5, col=as.factor(spec)))+
#   geom_line()+
#   ggtitle("Annual wind disturbance risk")