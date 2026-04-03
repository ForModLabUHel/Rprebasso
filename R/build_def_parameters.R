# ####code to create parameters objects  ##start
# 
# library(data.table)
# library(readxl)
# 
# ####default pcrobas parameters
# pCROB <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pCROBAS"))
# pCROBASnames <- pCROB$pNames
# pCROB <- pCROB[,2:(dim(pCROB)[2])]
# speciesNam <- names(pCROB)
# nparsAll <- length(speciesNam)
# rownames(pCROB) <- pCROBASnames
# 
# ###new parameters from ritika's calibration for the boreal species
# pCROBAS_Ritika_new <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pCROBAS_Ritika"))
# pCROBAS_Ritika <- pCROB
# pCROBAS_Ritika[,1:3] <- pCROBAS_Ritika_new[,2:4]
# 
# #####default clearcut parameters
# pars_clcut <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pClearcut"))
# inDclct_def <- rep(NA,(dim(pCROB)[2]))
# inAclct_def <- pars_clcut$inAclct_def
# 
# ####read preles parameters
# pPREL <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pPRELES_matrix"))
# pPREL <- as.matrix(pPREL[,2:ncol(pPREL)])
# # pPRELES_tab <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pPRELES_matrix"))
# # pPREL <- pPRELES_tab$pPREL
# # pPRELESeugl <-  pPRELES_tab$eucalyptus
# # pPRELESpiabDE <- pPRELES_tab$pPRELESpiabDE
# # pPRELESfasy <- pPRELES_tab$pPRELESfasy
# # pPRELES.Df.DBF <- pPRELES_tab$pPRELES.Df.DBF
# # pPRELESpipi <- pPRELES_tab$pPRELESpipi
# # pPRELES_Ritika <- pPRELES_tab$pPRELES_Ritika
# # pPRELES_catalonia <- pPRELES_tab$Catalonia
# 
# ####read preles parameters
# pPeattp_tab <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pPeattp"))
# pPeattp_def <- cbind(pPeattp_tab$peattype_1,pPeattp_tab$peattype_2)
# 
# ####light use efficiency parameters by species
# pLUEtrees <- as.numeric(pPREL[11,])
# # pLUEtrees <- rep(NA, dim(pCROB)[2])
# # pLUEtrees <- c(rep(pPREL[5],3),pPRELESfasy[5],pPRELESpipi[5],
# #                pPRELESeugl[5],pPRELES.Df.DBF[5],pPRELES.Df.DBF[5],
# #                pPRELESeugl[5],pPRELESpiabDE[5],pPRELES.Df.DBF[5],pPRELESfasy[5],
# #                rep(pPRELES_tab$Catalonia[5],10))
# names(pLUEtrees) <- colnames(pCROB)
# pLUEgv <- pPREL[5]
# 
# ####read parsAWEN
# parsAWEN <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pAWEN"))
# 
# ####read Hcmodel parameters
# pHcM <- data.table(read_excel("data/PREBAS_parameters.xlsx", sheet = "pHcM"))
# 
# ###convert to matrices
# pCROB = as.matrix(pCROB)
# pCROBAS_Ritika = as.matrix(pCROBAS_Ritika)
# parsAWEN = as.matrix(parsAWEN)
# pHcM = as.matrix(pHcM)
# 
# save(pCROB,pCROBAS_Ritika,
#      speciesNam, nparsAll,
#      inDclct_def,inAclct_def,pPeattp_def,
#      pPREL,#pPRELESeugl,pPRELESpiabDE,pPRELESfasy,pPRELES.Df.DBF,
#      # pPRELESpipi,pPRELES_Ritika,pPRELES_catalonia,
#      pLUEtrees,pLUEgv,
#      parsAWEN,
#      pHcM,file="data/parameters.rda")
# 
# ####code to create parameters objects  ##end