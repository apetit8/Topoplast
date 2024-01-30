source("scripts/functions/functions.R")
#Plot combining both empirical data and theoretical data
################################################################################
#Whisker plot###############################
nbrloop_thpl <- as.data.frame( cbind(read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")[,2],
                      read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")[,2],
                      read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thpl$Type <- "4thpl"
nbrloop_thpl$Sum <- rowSums2(as.matrix(nbrloop_thpl[,1:3]))

nbrloop_thnp <- as.data.frame( cbind(read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thnp$Type <- "3thnp"
nbrloop_thnp$Sum <- rowSums2(as.matrix(nbrloop_thnp[,1:3]))
  
nbrloop_empl <- as.data.frame(cbind(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,c(2,3)]))
nbrloop_empl$Type <- "2empl"
nbrloop_empl$Sum <- rowSums2(as.matrix(nbrloop_empl[,c(2,4,6)]))

nbrloop_emnp <- as.data.frame(cbind(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,c(2,3)]))
nbrloop_emnp$Type <- "1emnp"
nbrloop_emnp$Sum <- rowSums2(as.matrix(nbrloop_emnp[,c(2,4,6)]))

df <- rbind(nbrloop_emnp[,7:8], nbrloop_empl[,7:8], nbrloop_thnp[,4:5], nbrloop_thpl[,4:5] )
df <- df[(df$Sum !=0),]


pdf(paste0("figures/Distrib_loop_nbr",".pdf"), width=5, height=4)
par(mar = c(3.5,3.5, 1,1))
boxplot(df$Sum ~ df$Type, log = "y", col = NA, border = NA, axes = FALSE,
        yaxt="", space=c(0.3,0.1,0.4,0.1,0.4), ylab = "", xlab = "")
title(ylab = "Number of motif", line=2)
polygon(x=c(0.35, 0.35, 3, 3),  y=c(1000, 1e-5, 1e-5, 1000),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(2.5, 2.5, 4.65, 4.65),  y=c(1000, 1e-5, 1e-5, 1000),  col="cornsilk", border=NA, xpd=TRUE)
boxplot(df$Sum ~ df$Type,  log = "y", add=TRUE, col=c("lavender", "lightskyblue"), 
        names=c("Non plastic","Plastic","Non plastic","Plastic"),
        space=c(0.3,0.1,0.4,0.1,0.4), frame=FALSE)
text(1.5, 0.21, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1)
text(3.5,0.21, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1)
text(1.5,500, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
dev.off()

mean(subset(df, Type=="3thnp")$Sum) - mean(subset(df, Type=="4thpl")$Sum)
mean(subset(df, Type=="1emnp")$Sum) - mean(subset(df, Type=="2empl")$Sum)

################################################################################
#Drift









