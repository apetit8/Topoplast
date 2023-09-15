################################################################################
library(venneuler)

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
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%") )
text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%"))
text(0.35,0.57, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%"))
text(0.35,0.42, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%")) 

text(0.58,0.22, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%"))
text(0.58,0.70, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%")) 
text(0.28,0.49, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/sum(length(FFLs),length(DMDs),length(FBLs)), 1) ,"%")) 
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

MyVenn <- venneuler(c(FFL=length(FFLs),
                      DMD=length(DMDs),
                      FBL=length(FBLs), 
                      "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n"),paste0("\n\n\n\n\n\n\nDMD\n"),"FBL                                      ")
plot(MyVenn)
text(0.58,0.52, paste0(length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%") )
text(0.40,0.50, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 2) ,"%"))
text(0.35,0.55, paste0(length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n",
                       round((length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%"))
text(0.35,0.43, paste0((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))), "\n",
                       round((length(which(FBLs %in% DMDs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 

text(0.58,0.21, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%"))
text(0.59,0.78, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 
text(0.3,0.5, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 



################################################################################
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


x <- list(
  FFL = c(subset(FFLpl, Loop_number!=0)[,5], subset(FFLnp, Loop_number!=0)[,5]),
  DMD = c(subset(DMDpl, Loop_number!=0)[,5], subset(DMDnp, Loop_number!=0)[,5]), 
  FBL = c(subset(FBLpl, FBL_number!=0)[,5], subset(FBLnp, FBL_number!=0)[,5])
)

pdf(paste0("figures/Loop_Venn_simulations",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
ggvenn(
  x, fill_alpha = 0.7, stroke_color = "white",
  fill_color = c("khaki1", "darkseagreen2", "darkolivegreen2"),
  stroke_size = 0.5, set_name_size = 6, text_size = 4.5
)
dev.off()

################################################################################
library(venneuler)

FFLs <- c(subset(FFLpl, Loop_number!=0)[,5], subset(FFLnp, Loop_number!=0)[,5])
DMDs <- c(subset(DMDpl, Loop_number!=0)[,5], subset(DMDnp, Loop_number!=0)[,5])
FBLs <- c(subset(FBLpl, FBL_number!=0)[,5], subset(FBLnp, FBL_number!=0)[,5])

MyVenn <- venneuler(c(FFL=length(FFLs),
                      DMD=length(DMDs),
                      FBL=length(FBLs), 
                      "FFL&DMD"=length(which(FFLs %in% DMDs)) - length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),"FFL&FBL"=length(which(FFLs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "DMD&FBL"=length(which(DMDs %in% FBLs))- length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))),
                      "FFL&DMD&FBL"= length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))  ))
MyVenn$labels <- c(paste0("FFL\n\n\n\n\n\n"),paste0("\n\n\n\n\n\n\nDMD\n"),"FBL                                      ")
plot(MyVenn)
text(0.68,0.52, paste0(length(which(FFLs %in% DMDs)), "\n", round(length(which(FFLs %in% DMDs))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%") )
text(0.48,0.51, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n", round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 2) ,"%"))
text(0.4,0.63, paste0(length(which(FFLs %in% FBLs)), "\n", round(length(which(FFLs %in% FBLs))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%"))
text(0.4,0.35, paste0(length(which(FBLs %in% DMDs)), "\n", round(length(which(FBLs %in% DMDs))*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 
text(0.48,0.51, paste0(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs))), "\n", round(length(which(which(FFLs %in% DMDs) %in% which(FFLs %in% FBLs)))*100/sum(length(FFLs),length(DMDs),length(FBLs)), 2) ,"%"))



text(0.56,0.25, paste0(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)]), "\n", round(length( FFLs[!(FFLs %in% DMDs  | FFLs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%"))
text(0.63,0.70, paste0(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)]), "\n", round(length( DMDs[!(DMDs %in% FFLs  | DMDs %in% FBLs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 
text(0.24,0.45, paste0(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)]), "\n", round(length( FBLs[!(FBLs %in% FFLs  | FBLs %in% DMDs)])*100/sum(length(FFLs),length(DMDs),length(FBLs))) ,"%")) 







