library(eulerr)
################################################################################
################################################################################
FFLs <-  unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2]))
DMDs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2]))
FBLs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2]))
None <- unique(c(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,2],
                 read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,2],
                 read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,2]))
None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))

ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )


vd1 <- euler(ee, shape="circle")
plot(vd1, quantities=TRUE)


# Non plast
FFLs <-  c(subset(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2])
DMDs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2])
FBLs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2])
None <- c(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,2],
          read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,2],
          read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,2])
None <- unique(None)
None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))

ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )

vd2 <- euler(ee, shape="circle")
plot(vd2, quantities=TRUE)


pdf(paste0("figures/Loop_Venn_plast",".pdf"), width=4, height=3)
par(mgp=c(2.5, 1.2, 0), mar = c(0,2, 0,0), mfrow=c(2,1))

plot(vd1, quantities = list(type = c("percent"), round=2), main = "\nE. coli plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))

plot(vd2, quantities = list(type = c("percent"), round=2), main = "E. coli non-plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))
dev.off()

#SIMULATIONS####################################################################
FFLpl <- read.csv("scripts/data/full_netw_005_Pl_nbrFFL.csv", sep = ",")
FFLpl$Gene_ID <- paste0(FFLpl$X, "Pl")
DMDpl <- read.csv("scripts/data/full_netw_005_Pl_nbrDMD.csv", sep = ",")
DMDpl$Gene_ID <- paste0(DMDpl$X, "Pl")
FBLpl <- read.csv("scripts/data/full_netw_005_Pl_nbrFBL.csv", sep = ",")
FBLpl$Gene_ID <- paste0(FBLpl$X, "Pl")

FFLs <- c(subset(FFLpl, Loop_number!=0)[,5])
DMDs <- c(subset(DMDpl, Loop_number!=0)[,5])
FBLs <- c(subset(FBLpl, FBL_number!=0)[,5])
colnames(FBLpl) <- colnames(FFLpl)
None <- subset(rbind(FFLpl, DMDpl, FBLpl), !(Gene_ID %in% c(FFLs, DMDs, FBLs)))[,5]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))


ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )


vd1 <- euler(ee, shape="circle")
plot(vd1, quantities=TRUE)

####
FFLnp <- read.csv("scripts/data/full_netw_005_NP_nbrFFL.csv", sep = ",")
FFLnp$Gene_ID <- paste0(FFLnp$X, "Np")
DMDnp <- read.csv("scripts/data/full_netw_005_NP_nbrDMD.csv", sep = ",")
DMDnp$Gene_ID <- paste0(DMDnp$X, "Np")
FBLnp <- read.csv("scripts/data/full_netw_005_NP_nbrFBL.csv", sep = ",")
FBLnp$Gene_ID <- paste0(FBLnp$X, "Np")

FFLs <- c(subset(FFLnp, Loop_number!=0)[,5])
DMDs <- c(subset(DMDnp, Loop_number!=0)[,5])
FBLs <- c(subset(FBLnp, FBL_number!=0)[,5])
colnames(FBLnp) <- colnames(FFLnp)
None <- subset(rbind(FFLnp, DMDnp, FBLnp), !(Gene_ID %in% c(FFLs, DMDs, FBLs)))[,5]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))


ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )


vd2 <- euler(ee, shape="circle")
plot(vd2, quantities=TRUE)



pdf(paste0("figures/Loop_Venn_plast_simul",".pdf"), width=4, height=3)
par(mgp=c(2.5, 1.2, 0), mar = c(0,2, 0,0), mfrow=c(2,1))

plot(vd1, quantities = list(type = c("percent"), round=2), main = "\n\nSimulated plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))
plot(vd2, quantities = list(type = c("percent"), round=2), main = "Simulated non-plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))
dev.off()


#DRIFT####################################################################
FFLpl <- read.csv("scripts/data/full_netw_005_drift_Pl_nbrFFL.csv", sep = ",")
FFLpl$Gene_ID <- paste0(FFLpl$X, "Pl")
DMDpl <- read.csv("scripts/data/full_netw_005_drift_Pl_nbrDMD.csv", sep = ",")
DMDpl$Gene_ID <- paste0(DMDpl$X, "Pl")
FBLpl <- read.csv("scripts/data/full_netw_005_drift_Pl_nbrFBL.csv", sep = ",")
FBLpl$Gene_ID <- paste0(FBLpl$X, "Pl")

FFLs <- c(subset(FFLpl, Loop_number!=0)[,5])
DMDs <- c(subset(DMDpl, Loop_number!=0)[,5])
FBLs <- c(subset(FBLpl, FBL_number!=0)[,5])
colnames(FBLpl) <- colnames(FFLpl)
None <- subset(rbind(FFLpl, DMDpl, FBLpl), !(Gene_ID %in% c(FFLs, DMDs, FBLs)))[,5]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))


ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )


vd1 <- euler(ee, shape="circle")
plot(vd1, quantities=TRUE)

####
FFLnp <- read.csv("scripts/data/full_netw_005_drift_NP_nbrFFL.csv", sep = ",")
FFLnp$Gene_ID <- paste0(FFLnp$X, "Np")
DMDnp <- read.csv("scripts/data/full_netw_005_drift_NP_nbrDMD.csv", sep = ",")
DMDnp$Gene_ID <- paste0(DMDnp$X, "Np")
FBLnp <- read.csv("scripts/data/full_netw_005_NP_nbrFBL.csv", sep = ",")
FBLnp$Gene_ID <- paste0(FBLnp$X, "Np")

FFLs <- c(subset(FFLnp, Loop_number!=0)[,5])
DMDs <- c(subset(DMDnp, Loop_number!=0)[,5])
FBLs <- c(subset(FBLnp, FBL_number!=0)[,5])
colnames(FBLnp) <- colnames(FFLnp)
None <- subset(rbind(FFLnp, DMDnp, FBLnp), !(Gene_ID %in% c(FFLs, DMDs, FBLs)))[,5]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))


ee <- c("FFL"=length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]),
        "DMD"=length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]),
        "FBL"=length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), 
        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "None"=length(None)  )


vd2 <- euler(ee, shape="circle")
plot(vd2, quantities=TRUE)



pdf(paste0("figures/Loop_Venn_plast_drift",".pdf"), width=4, height=3)
par(mgp=c(2.5, 1.2, 0), mar = c(0,2, 0,0), mfrow=c(2,1))

plot(vd1, quantities = list(type = c("percent"), round=2), main = "\n\nDrift plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))
plot(vd2, quantities = list(type = c("percent"), round=2), main = "Drift non-plastic genes", lty = 0, fills = c("powderblue", "khaki1", "sienna1", "grey"))
dev.off()

################################################################################
################################################################################
# ################################################################################
# library(venneuler)
# library(plotrix)
# ################################################################################
# #E coli plast and non plast
# ################################################################################
# 
# pdf(paste0("figures/Loop_Venn_plast",".pdf"), width=5, height=4)
# par(mgp=c(2.5, 1.2, 0), mar = c(0,2, 0,0), mfrow=c(1,2))
# FFLs <-  unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2]))
# DMDs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2]))
# FBLs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2]))
# None <- unique(c(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,2],
#                  read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,2],
#                  read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,2]))
# None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
# Total <- length(unique(c(FFLs, DMDs, FBLs, None)))
# MyVenn <- venneuler(c(FFL=length(FFLs),
#                       DMD=length(DMDs),
#                       FBL=length(FBLs), 
#                       "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
#                       "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
#                       "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
# 
# MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n\n"),paste0("\n\n\n\nDMD\n"),"FBL                        ")
# 
# plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"), add=TRUE, xpd=TRUE, main="\n\n\nPlastic genes")
# text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8 )
# text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/Total, 1) ,"%"), cex=0.8)
# text(0.36,0.59, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8)
# text(0.36,0.41, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
#                        round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8) 
# 
# text(0.58,0.72, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8)
# text(0.58,0.22, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8) 
# text(0.28,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8) 
# text(0.23,0.58, "FBL", xpd=TRUE) 
# 
# 
# barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", add=TRUE)
# draw.circle(x=0.55, y=-0.10, radius=0.21, col="gray", border = 0)
# text(0.55,0.05, "No loop", xpd=TRUE)
# text(0.55,-0.02, paste0(length(None), "\n", round(length( None)*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )
# ##########
# 
# FFLs <-  c(subset(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2])
# DMDs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2])
# FBLs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2])
# None <- c(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,2],
#           read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,2],
#           read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,2])
# None <- unique(None)
# None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
# Total <- length(unique(c(FFLs, DMDs, FBLs, None)))
# MyVenn2 <- venneuler(c(FFL=length(FFLs),
#                        DMD=length(DMDs),
#                        FBL=length(FBLs), 
#                        "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
#                        "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
#                        "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
# MyVenn2$labels <- c(paste0("FFL\n\n\n\n\n\n\n"),paste0("\n\n\n\nDMD\n"),"FBL                        ")
# plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"), xpd=TRUE, main="\n\n\nNon plastic genes")
# text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8 )
# text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/Total, 1) ,"%"), cex=0.8)
# text(0.36,0.59, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
#                        round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8)
# text(0.36,0.41, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
#                        round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8) 
# 
# text(0.58,0.72, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8)
# text(0.58,0.22, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8) 
# text(0.28,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )
# text(0.23,0.58, "FBL", xpd=TRUE) 
# 
# barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", add=TRUE)
# draw.circle(x=0.55, y=-0.3, radius=0.4, col="gray", border = 0)
# text(0.55,0.05, "No loop", xpd=TRUE)
# text(0.55,-0.02, paste0(length( None), "\n", round(length( None)*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )
# 
# dev.off()
