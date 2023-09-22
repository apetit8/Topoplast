################################################################################
library(venneuler)
library(plotrix)

FFLs <-  c(subset(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2])
DMDs <- c(subset(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2])
FBLs <- c(subset(read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2])

pdf(paste0("figures/Loop_Venn",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2,0, 0,0))
MyVenn <- venneuler(c(FFL=length(FFLs),
                      DMD=length(DMDs),
                      FBL=length(FBLs), 
                      "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n\n\n"),paste0("\n\n\n\n\nDMD\n"),"FBL                        ")
plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"))
text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%") )
text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.35,0.57, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.35,0.42, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 

text(0.58,0.70, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.58,0.22, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 
text(0.28,0.49, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 
dev.off()



################################################################################
################################################################################
library(venneuler)

##And now for simulations ?#####################################################
#Genes are in the same order ; Plasticity+row = ID !

FFLpl <- read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")
FFLpl$Gene_ID <- paste0(FFLpl$X, "Pl")
DMDpl <- read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")
DMDpl$Gene_ID <- paste0(DMDpl$X, "Pl")
FBLpl <- read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ",")
FBLpl$Gene_ID <- paste0(FBLpl$X, "Pl")
#
FFLnp <- read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")
FFLnp$Gene_ID <- paste0(FFLnp$X, "Np")
DMDnp <- read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")
DMDnp$Gene_ID <- paste0(DMDnp$X, "Np")
FBLnp <- read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ",")
FBLnp$Gene_ID <- paste0(FBLnp$X, "Np")

FFLs <- c(subset(FFLpl, Loop_number!=0)[,5], subset(FFLnp, Loop_number!=0)[,5])
DMDs <- c(subset(DMDpl, Loop_number!=0)[,5], subset(DMDnp, Loop_number!=0)[,5])
FBLs <- c(subset(FBLpl, FBL_number!=0)[,5], subset(FBLnp, FBL_number!=0)[,5])

pdf(paste0("figures/Loop_Venn_simulations",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2,0, 0,0))
MyVenn <- venneuler(c(FFL=length(FFLs),
                      DMD=length(DMDs),
                      FBL=length(FBLs), 
                      "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n\n\n"),paste0("\n\n\n\n\nDMD\n"),"FBL                              ")
plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"))
text(0.63,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%") )
text(0.45,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.38,0.6, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.40,0.40, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 

text(0.56,0.7, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%"))
text(0.57,0.24, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 
text(0.3,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/length(unique(c(FFLs, DMDs, FBLs))), 1) ,"%")) 
dev.off()


################################################################################
#E coli plast and non plast
################################################################################

pdf(paste0("figures/Loop_Venn_plast",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(0,2, 0,0), mfrow=c(1,2))
FFLs <-  unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2]))
DMDs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2]))
FBLs <- unique(c(subset(read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2]))
None <- unique(c(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,2],
          read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,2],
          read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,2]))
None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))
MyVenn <- venneuler(c(FFL=length(FFLs),
                      DMD=length(DMDs),
                      FBL=length(FBLs), 
                      "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))

MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n\n"),paste0("\n\n\n\nDMD\n"),"FBL                        ")

plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"), add=TRUE, xpd=TRUE, main="\n\n\nPlastic genes")
text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8 )
text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/Total, 1) ,"%"), cex=0.8)
text(0.36,0.59, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8)
text(0.36,0.41, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8) 

text(0.58,0.72, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8)
text(0.58,0.22, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8) 
text(0.28,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8) 
text(0.23,0.58, "FBL", xpd=TRUE) 


barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", add=TRUE)
draw.circle(x=0.55, y=-0.10, radius=0.21, col="gray", border = 0)
text(0.55,0.05, "No loop", xpd=TRUE)
text(0.55,-0.02, paste0(length(None), "\n", round(length( None)*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )
##########

FFLs <-  c(subset(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2])
DMDs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2])
FBLs <- c(subset(read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2])
None <- c(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,2],
          read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,2],
          read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,2])
None <- unique(None)
None <- None[!(None %in% FFLs | None %in% DMDs | None %in% FBLs)]
Total <- length(unique(c(FFLs, DMDs, FBLs, None)))
MyVenn2 <- venneuler(c(FFL=length(FFLs),
                       DMD=length(DMDs),
                       FBL=length(FBLs), 
                       "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                       "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                       "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
MyVenn2$labels <- c(paste0("FFL\n\n\n\n\n\n\n"),paste0("\n\n\n\nDMD\n"),"FBL                        ")
plot(MyVenn, col=c("orange", "palegreen3", "hotpink2"), xpd=TRUE, main="\n\n\nNon plastic genes")
text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8 )
text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/Total, 1) ,"%"), cex=0.8)
text(0.36,0.59, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8)
text(0.36,0.41, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/Total, 1) ,"%"), cex=0.8) 

text(0.58,0.72, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8)
text(0.58,0.22, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/Total, 1) ,"%"), cex=0.8) 
text(0.28,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )
text(0.23,0.58, "FBL", xpd=TRUE) 

barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", add=TRUE)
draw.circle(x=0.55, y=-0.3, radius=0.4, col="gray", border = 0)
text(0.55,0.05, "No loop", xpd=TRUE)
text(0.55,-0.02, paste0(length( None), "\n", round(length( None)*100/Total, 1) ,"%"), xpd=TRUE, cex=0.8 )

dev.off()

