#Fig supp total number of regulations
#Take E_coli matrix, subset plastic genes, get row sum and compare to non plastic gene row sum
source("scripts/functions/functions.R")
library(purrr)
#########################################
pdfname <- "figures/fig_supp_reg_nbr"
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd_2024_01_29.csv") #List of regulations from Ecocyc
ec_genes <- read.csv("e_coli/ncbi_dataset_K-12_annotation.csv", sep ="\t") #List of E coli genes from Ecocyc

#Transcriptions factors
TF_genes <- unique(ec_cyc[,1]) # Here, TF = all regulating genes from Ecocyc reg data
#Adding non-regulated genes in regulation data
nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,2] %in% ec_cyc[,1]),]) #Regulators
nonreg_genes <- nonreg_genes[!(nonreg_genes[,2] %in% ec_cyc[,2]),] #Regulatees
#1639 genes that are non-regulated
nr <- as.data.frame(cbind(nonreg_genes[,2],nonreg_genes[,2], rep(0, nrow(nonreg_genes)))) #Same file format as ec_cyc
setnames(nr, 1:3, colnames(ec_cyc))
ec_reg_full <- rbind(ec_cyc, nr) #Reg file with every genes
##########
g <- graph.data.frame(ec_reg_full, directed=TRUE)
E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3'))) #t() to have regulators in columns and regulees in rows

E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
E_coli_mat[is.na(E_coli_mat)] <- 0 #fill NA to 0, mandatory for later analyses
################################################################################
##E. coli
plastgenes <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,2]
npgenes <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,2]
E_coli_mat_plast <- E_coli_mat[plastgenes, ]
E_coli_mat_np <- E_coli_mat[npgenes, ]
#
plastmat <- E_coli_mat_plast
plastmat[(plastmat == -1)] <- 1
reg_number_plast <- as.data.frame(rowSums2(plastmat))
reg_number_plast$Type <- "2empl"
colnames(reg_number_plast) <- c("nbr","Type")
#
plastmat <- E_coli_mat_np
plastmat[(plastmat == -1)] <- 1
reg_number_np <- as.data.frame(rowSums2(plastmat))
reg_number_np$Type <- "1emnp"
colnames(reg_number_np) <- c("nbr","Type")

##Simulations
topo_plastic <- readRDS("scripts/data/list_plastic_topo_full_netw.Rds")
topo_nnplast <- readRDS("scripts/data/list_nnplast_topo_full_netw.Rds")

THreg_sum_plast <- as.data.frame(unlist(lapply(topo_plastic, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_plast$Type <- "4thpl"
colnames(THreg_sum_plast) <- c("nbr","Type")
#
THreg_sum_np <- as.data.frame(unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_np$Type <- "3thnp"
colnames(THreg_sum_np) <- c("nbr","Type")

df1 <- rbind(reg_number_plast, reg_number_np, THreg_sum_np, THreg_sum_plast )

mean(reg_number_plast[,1])
mean(reg_number_np[,1])
mean(THreg_sum_np[,1])
mean(THreg_sum_plast[,1])

################################################################################
#Number of motif per gene
################################################################################
#Whisker plot###############################
nbrloop_thpl2 <- as.data.frame( cbind(read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thpl2$Type <- "4thpl"
nbrloop_thpl2$Sum <- rowSums2(as.matrix(nbrloop_thpl2[,1:3]))

nbrloop_thnp2 <- as.data.frame( cbind(read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")[,2],
                                     read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thnp2$Type <- "3thnp"
nbrloop_thnp2$Sum <- rowSums2(as.matrix(nbrloop_thnp2[,1:3]))

nbrloop_empl2 <- as.data.frame(cbind(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,c(2,3)]))
nbrloop_empl2$Type <- "2empl"
nbrloop_empl2$Sum <- rowSums2(as.matrix(nbrloop_empl2[,c(2,4,6)]))

nbrloop_emnp2 <- as.data.frame(cbind(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,c(2,3)],
                                    read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,c(2,3)]))
nbrloop_emnp2$Type <- "1emnp"
nbrloop_emnp2$Sum <- rowSums2(as.matrix(nbrloop_emnp2[,c(2,4,6)]))

df <- rbind(nbrloop_emnp2[,7:8], nbrloop_empl2[,7:8], nbrloop_thnp2[,4:5], nbrloop_thpl2[,4:5] )
df <- df[(df$Sum !=0),]


reg_number_plast <- rowSums2(abs(E_coli_mat_plast))
E_coli_mat_np <- E_coli_mat[npgenes, ]
reg_number_np <- rowSums2(abs(E_coli_mat_np))

randomTf <- as.data.frame(cbind(read.csv("scripts/data/plast_genes_E_coli_random3_nffl.csv", sep = ",")[,c(2,3)],
                                     read.csv("scripts/data/plast_genes_E_coli_random3_nDMD.csv", sep = ",")[,c(2,3)],
                                     read.csv("scripts/data/plast_genes_E_coli_random3_nFBL.csv", sep = ",")[,c(2,3)]))
randomTf$Type <- "RandomTF"
randomTf$Sum <- rowSums2(as.matrix(randomTf[,c(2,4,6)]))
randomTf <- randomTf[order(as.character(randomTf[,1]), method = c("radix")),c(1,8)]
reg_number_ranTf <- reg_number_plast[order(names(reg_number_plast), method = c("radix"))] #Put genes of the different dataset in the same order

#Without the 0 for a log axis
RTF <- as.data.frame(cbind(reg_number_ranTf, randomTf[,2]))
PL <- as.data.frame(cbind(reg_number_plast, nbrloop_empl2[,8]))
NP <- as.data.frame(cbind(reg_number_np, nbrloop_emnp2[,8]))

#FAIRE PAREIL AVEC SIMU

randomTf <- as.data.frame( cbind(read.csv("scripts/data/full_netw_randomTF_Pl_nbrFFL.csv", sep = ",")[,2],
                                 read.csv("scripts/data/full_netw_randomTF_Pl_nbrDMD.csv", sep = ",")[,2],
                                 read.csv("scripts/data/full_netw_randomTF_Pl_nbrFBL.csv", sep = ",")[,2]))
randomTf$Type <- "RandomTF"
randomTf$Sum <- rowSums2(as.matrix(nbrloop_thpl2[,1:3]))

PL2 <- as.data.frame(cbind(THreg_sum_plast[,1], nbrloop_thpl2[,5]))
NP2 <- as.data.frame(cbind(THreg_sum_np[,1], nbrloop_thnp2[,5]))



#####################
RTF <- subset(RTF, V2 != 0)
PL <- subset(PL, V2 != 0)
NP <- subset(NP, V2 != 0)
RTF <- subset(RTF, V2 != 0)
PL2 <- subset(PL2, V2 != 0)
NP2 <- subset(NP2, V2 != 0)

#To display the same number of plastic genes in E coli and simu
PL2 <- PL2[sample(nrow(PL2), nrow(PL)),]

#Drift##########
topo_plasticD <- readRDS("scripts/data/list_plastic_topo_full_netw_drift.Rds")
topo_nnplastD <- readRDS("scripts/data/list_nnplast_topo_full_netw_drift.Rds")

THreg_sum_plastD <- as.data.frame(unlist(lapply(topo_plasticD, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_plastD$Type <- "6thplD"
colnames(THreg_sum_plastD) <- c("nbr","Type")
#
THreg_sum_npD <- as.data.frame(unlist(lapply(topo_nnplastD, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_npD$Type <- "5thnpD"
colnames(THreg_sum_npD) <- c("nbr","Type")

df1 <- rbind(df1, THreg_sum_npD, THreg_sum_plastD )

nbrloop_thplD <- as.data.frame( cbind(read.csv("scripts/data/full_netw_drift_Pl_nbrFFL.csv", sep = ",")[,2],
                                      read.csv("scripts/data/full_netw_drift_Pl_nbrDMD.csv", sep = ",")[,2],
                                      read.csv("scripts/data/full_netw_drift_Pl_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thplD$Type <- "4thpl"
nbrloop_thplD$Sum <- rowSums2(as.matrix(nbrloop_thplD[,1:3]))

nbrloop_thnpD <- as.data.frame( cbind(read.csv("scripts/data/full_netw_drift_NP_nbrFFL.csv", sep = ",")[,2],
                                      read.csv("scripts/data/full_netw_drift_NP_nbrDMD.csv", sep = ",")[,2],
                                      read.csv("scripts/data/full_netw_drift_NP_nbrFBL.csv", sep = ",")[,2]))
nbrloop_thnpD$Type <- "3thnp"
nbrloop_thnpD$Sum <- rowSums2(as.matrix(nbrloop_thnpD[,1:3]))

PLD <- as.data.frame(cbind(THreg_sum_plastD[,1], nbrloop_thplD[,5]))
NPD <- as.data.frame(cbind(THreg_sum_npD[,1], nbrloop_thnpD[,5]))


#FIGURES########################################################################
pdf(paste0("figures/Fig1_Reg_nbr",".pdf"), width=4, height=3.5)
par(mar = c(4.8,3.5, 2,0))
boxplot(df1$nbr ~ df1$Type,  main="", col = NA, border = NA, axes = FALSE,
        yaxt="", at=c(1,2, 3.5,4.5,6,7), ylab = "", xlab = "")
title(ylab = "Number of regulators", line=2)
polygon(x=c(0.2, 0.2, 2.75, 2.75),  y=c(50, -20, -20, 50),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(2.75, 2.75, 120, 120),  y=c(50, -20, -20, 50),  col="cornsilk", border=NA, xpd=TRUE)
boxplot(df1$nbr ~ df1$Type, add=TRUE, col=c("salmon", "grey"),
        at=c(1,2, 3.5,4.5,6,7), frame=FALSE, xaxt="n", pch=19, outcol=alpha("black", 0.1))
text(x = c(1,2, 3.5,4.5,6,7), y = par("usr")[3] - 0.3,
     labels = c("Non-Plastic", "Plastic", "Non-Plastic", "Plastic", "Non-Plastic", "Plastic"), xpd = NA, srt = 90, cex = 1, adj=1.0)
text(1.5, 32.5, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1)
text(5.5, 32.5, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1)
text(4, 27, substitute(paste("  Under\nSelection")), xpd=TRUE, cex=1)
text(6.5, 29, substitute(paste("Drift")), xpd=TRUE, cex=1)
text(1.5,25, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
dev.off()

t.test(reg_number_np, reg_number_plast)

pdf(paste0("figures/Motif_Reg_both_log",".pdf"), width=6, height=3.5)
#E. coli
layout(matrix(c(1,2,3), 1,3, byrow = TRUE))
#Plot 1
par(mar=c(3, 3, 2, 0.3), mgp = c(1.75, 0.75, 0), las=0)
plot(PLD[,1], PLD[,2], log = "y", col=alpha("darkblue", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19, frame.plot = FALSE)
polygon(x=c(-15, -15, 35, 35),  y=c(50000, 0.01, 0.01, 50000),  col="honeydew2", border=NA, xpd=TRUE)
axis(side = 1, cex=1.4)
axis(side = 2, cex=1.4)
points( PL[,1], PL[,2], col=alpha("black", 0.5), pch=19)
points( NP[,1], NP[,2], col=alpha("tomato", 0.5), pch=19)
title(ylab="Number of Loops", line=2, cex.lab=1.4)
text(15, 11000, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1.4)
#Plot 2
par(mar=c(3, 0.5, 2, 0), mgp = c(1.75, 0.75, 0), las=0)
plot(PLD[,1], PLD[,2], log = "y", col=alpha("black", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19, yaxt='n', frame.plot = FALSE)
polygon(x=c(-5, -5, 35, 35),  y=c(50000, 0.01, 0.01, 50000),  col="cornsilk", border=NA, xpd=TRUE)
axis(side = 1)
#axis(side = 2)
title(ylab="Number of Loops", xlab="Number of regulators", line=2, cex.lab=1.4)
points(PL2[,1], PL2[,2], col=alpha("black", 0.5), pch=19)
points(NP2[,1], NP2[,2], col=alpha("tomato", 0.5), pch=19)
text(31.5, 11000, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.4)
legend(12, 5, legend=c("Plastic genes", "Non plastic genes"), col=c("black", "tomato"), bty=1, cex=1, bg="white", pch = c(19))
text(15, 6500, substitute(paste("Under Selection")), xpd=TRUE, cex=1.4)
#Plot 3
plot(PLD[,1], PLD[,2], log = "y", col=alpha("black", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19, yaxt='n', frame.plot = FALSE)
polygon(x=c(-15, -15, 32, 32),  y=c(50000, 0.01, 0.01, 50000),  col="cornsilk", border=NA, xpd=TRUE)
axis(side = 1)
points(PLD[,1], PLD[,2], col=alpha("black", 0.5), pch=19)
points(NPD[,1], NPD[,2], col=alpha("tomato", 0.5), pch=19)
text(-2, 11000, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.4)
text(15, 6500, substitute(paste("Drift")), xpd=TRUE, cex=1.4)
dev.off()



pdf(paste0("figures/Motif_Reg_log_drift",".pdf"), width=9, height=3.5)
# par(mar = c(3.5,3.5, 1,1))
#E. coli
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
#par(mar=c(3, 3, 1.5, 0.3), mgp = c(1.75, 0.75, 0), las=0)
#Plot 1
par(mar = c(4.8,3.5, 2,0))
boxplot(dfD$nbr ~ dfD$Type,  main="", col = NA, border = NA, axes = FALSE,
        yaxt="", space=c(0.3,0.1,0.4,0.1,0.4), ylab = "", xlab = "")
title(ylab = "Number of regulation toward gene", line=2)
polygon(x=c(-5, -5, 2.5, 2.5),  y=c(50, -20, -20, 50),  col="cornsilk", border=NA, xpd=TRUE)
boxplot(dfD$nbr ~ dfD$Type, add=TRUE, col=c("salmon", "grey"),
        space=c(0.3,0.1,0.4,0.1,0.4), frame=FALSE, xaxt="n")
text(x = c(1,2,3,4), y = par("usr")[3] - 0.3,
     labels = c("Non-Plastic", "Plastic", "Non-Plastic", "Plastic"), xpd = NA, srt = 90, cex = 1, adj=1.0)
text(1.5, 29, substitute(paste(bold("Simulations: Drift"))), xpd=TRUE, cex=1)
#Plot 2
par(mar=c(3, 3, 1.5, 0.3), mgp = c(1.75, 0.75, 0), las=0)
plot(PLD[,1], PLD[,2], log = "y", col=alpha("black", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19, yaxt='n')
polygon(x=c(-15, -15, 32, 32),  y=c(50000, 0.01, 0.01, 50000),  col="cornsilk", border=NA, xpd=TRUE)
axis(side = 1)
axis(side = 2)
title(ylab="Number of Loop", xlab="Number of regulation", line=2)
points(PLD[,1], PLD[,2], col=alpha("black", 0.5), pch=19)
points(NPD[,1], NPD[,2], col=alpha("tomato", 0.5), pch=19)
#legend(15, 5, legend=c("Plast genes", "Non plast genes"), col=c("blue", "green"), lty=1, cex=0.8)
text(15, 9000, substitute(paste(bold("Simulations: Drift"))), xpd=TRUE, cex=1)
legend(12, 5, legend=c("Plastic genes", "Non plastic genes"),
       col=c("black", "tomato"), bty=1, cex=0.8, bg="white", pch = c(19))
dev.off()




##Non-log:
# pdf(paste0("figures/Motif_Reg_both",".pdf"), width=6, height=4)
# #E. coli
# layout(matrix(c(1,2), 1, 2, byrow = TRUE))
# #Plot 1
# par(mar=c(3, 3, 1.5, 0.3), mgp = c(1.75, 0.75, 0), las=0)
# plot(PL2[,1], PL2[,2], col=alpha("darkblue", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19)
# polygon(x=c(-15, -15, 35, 35),  y=c(3000, -200, -200, 3000),  col="honeydew2", border=NA, xpd=TRUE)
# axis(side = 1)
# axis(side = 2)
# points( PL[,1], PL[,2], col=alpha("black", 0.5), pch=19)
# points( NP[,1], NP[,2], col=alpha("tomato", 0.5), pch=19)
# title(ylab="Number of Loop", xlab="Number of regulation", line=2)
# # abline(lm( log10(PL[,2]) ~ PL[,1]), col="black")
# # abline(lm( log10(NP[,2]) ~ NP[,1]), col="tomato")
# text(15, 1500, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1)
# #Simu
# par(mar=c(3, 0.5, 1.5, 0), mgp = c(1.75, 0.75, 0), las=0)
# plot(PL2[,1], PL2[,2], col=alpha("black", 0.5), ylab="", xlab="", xlim=c(0,30), pch=19, yaxt='n')
# polygon(x=c(-5, -5, 32, 32),  y=c(3000, -200, -200, 3000),  col="cornsilk", border=NA, xpd=TRUE)
# axis(side = 1)
# #axis(side = 2)
# title(ylab="Number of Loop", xlab="Number of regulation", line=2)
# points(PL2[,1], PL2[,2], col=alpha("black", 0.3), pch=19)
# points(NP2[,1], NP2[,2], col=alpha("tomato", 0.5), pch=19)
# # abline(lm(log10(PL2[,2]) ~ PL2[,1]), col="black")
# # abline(lm(log10(NP2[,2]) ~ NP2[,1]), col="tomato")
# #legend(15, 5, legend=c("Plast genes", "Non plast genes"), col=c("blue", "green"), lty=1, cex=0.8)
# text(15, 1500, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1)
# legend(12, 700, legend=c("Plastic genes", "Non plastic genes"),
#        col=c("black", "tomato"), bty=1, cex=0.8, bg="white", pch = c(19))
# dev.off()


# plot(THreg_sum_plast[,1], randomTf[,5], log = "y", col="darkred", ylab="", xlab="", xlim=c(0,30))
# title(ylab="Number of Loop", xlab="Number of regulation", line=2)
# points(THreg_sum_plast[,1], nbrloop_thpl2[,5], col="blue")
# points(THreg_sum_np[,1], nbrloop_thnp2[,5], col="green")
# abline(lm(log10(nbrloop_thpl2[,5]) ~ THreg_sum_plast[,1]), col="blue")
# abline(lm(log10(nbrloop_thnp2[,5]) ~ THreg_sum_np[,1]), col="green")
# legend(1, 300, legend=c("Plast genes: reg from random TF", "Plast genes", "Non plast genes"),
#        col=c("darkred", "blue", "green"), lty=1, cex=0.8)
# 
# 
# 
# #DO FFL AND DIAMONDS SEPARATELY
# 
# 
# randomTf <- as.data.frame(cbind(read.csv("scripts/data/plast_genes_E_coli_random3_nffl.csv", sep = ",")[,c(2,3)]))
# randomTf <- randomTf[order(as.character(randomTf[,1]), method = c("radix")),c(1,2)]
# reg_number_ranTf <- reg_number_plast[order(names(reg_number_plast), method = c("radix"))] #Put genes of the different dataset in the same order
# 
# 
# #Only FFL, should use all loops
# par(mar = c(3.5,3.5, 1,1))
# plot(reg_number_ranTf, randomTf[,2], col="darkred", ylab="", xlab="", ylim=c(0,100), xlim=c(0,30))
# title(ylab="Number of Loop", xlab="Number of regulation", line=2)
# points(reg_number_plast, read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3], col="blue")
# points(reg_number_np, read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3], col="green")
# abline(lm( randomTf[,2] ~ reg_number_ranTf), col="darkred")
# abline(lm(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3]  ~ reg_number_plast), col="blue")
# abline(lm( read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3] ~ reg_number_np), col="green")
# legend(1, 300, legend=c("Plast genes: reg from random TF", "Plast genes", "Non plast genes"),
#        col=c("darkred", "blue", "green"), lty=1, cex=0.8)
# 
# 
# 
# randomTf <- as.data.frame(cbind(read.csv("scripts/data/plast_genes_E_coli_random3_nDMD.csv", sep = ",")[,c(2,3)]))
# randomTf <- randomTf[order(as.character(randomTf[,1]), method = c("radix")),c(1,2)]
# reg_number_ranTf <- reg_number_plast[order(names(reg_number_plast), method = c("radix"))] #Put genes of the different dataset in the same order
# #Only DMD
# par(mar = c(3.5,3.5, 1,1))
# plot(reg_number_ranTf, randomTf[,2], col="darkred", ylab="", xlab="", ylim=c(0,100), xlim=c(0,30))
# title(ylab="Number of Loop", xlab="Number of regulation", line=2)
# points(reg_number_plast, read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3], col="blue")
# points(reg_number_np, read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3], col="green")
# abline(lm( randomTf[,2] ~ reg_number_ranTf), col="darkred")
# abline(lm(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3]  ~ reg_number_plast), col="blue")
# abline(lm( read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3] ~ reg_number_np), col="green")
# legend(1, 300, legend=c("Plast genes: reg from random TF", "Plast genes", "Non plast genes"),
#        col=c("darkred", "blue", "green"), lty=1, cex=0.8)
# 
# 
# 
# 
# pdf(paste0("figures/Distrib_loop_nbr",".pdf"), width=5, height=4)
# par(mar = c(3.5,3.5, 1,1))
# boxplot(df$Sum ~ df$Type, log = "y", col = NA, border = NA, axes = FALSE,
#         yaxt="", space=c(0.3,0.1,0.4,0.1,0.4), ylab = "", xlab = "")
# title(ylab = "Number of motif", line=2)
# polygon(x=c(0.35, 0.35, 3, 3),  y=c(1000, 1e-5, 1e-5, 1000),  col="honeydew2", border=NA, xpd=TRUE)
# polygon(x=c(2.5, 2.5, 4.65, 4.65),  y=c(1000, 1e-5, 1e-5, 1000),  col="cornsilk", border=NA, xpd=TRUE)
# boxplot(df$Sum ~ df$Type,  log = "y", add=TRUE, col=c("lavender", "lightskyblue"), 
#         names=c("Non plastic","Plastic","Non plastic","Plastic"),
#         space=c(0.3,0.1,0.4,0.1,0.4), frame=FALSE)
# text(1.5, 0.21, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1)
# text(3.5,0.21, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1)
# text(1.5,500, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
# dev.off()
# 
# mean(subset(df, Type=="3thnp")$Sum) - mean(subset(df, Type=="4thpl")$Sum)
# mean(subset(df, Type=="1emnp")$Sum) - mean(subset(df, Type=="2empl")$Sum)