source("scripts/functions/functions.R")
pdfname <- "figures/fig_full_netw"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
##FFL motifs Data###############################################################
control_np <- read.csv("scripts/data/nonplast_E_coli_random2_FFL.csv", sep = ",")
control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_FFL.csv", sep = ",")
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_FFL.csv", sep = ",")[,1:11]
drift <- read.csv("scripts/data/full_netw_drift_FFL.csv", sep = ",")[,1:11]
s_shuffled <- read.csv("scripts/data/full_netw_shuffled_FFL.csv", sep = ",")[,1:11]

df_ffl <- as.data.frame(rbind(colSums(control_np[,5:12])*100/sum(control_np[,3]), colSums(control_pl[,5:12])*100/sum(control_pl[,3]),
                              colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                              theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2],
                              s_shuffled[2,4:11]*100/s_shuffled[2,2],s_shuffled[1,4:11]*100/s_shuffled[1,2],
                        drift[2,4:11]*100/drift[2,2],drift[1,4:11]*100/drift[1,2]))
rownames(df_ffl) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nShuffled","Plastic\nShuffled", "Non-plastic\nDrift","Plastic\nDrift")
df_ffl <- df_ffl[,order(as.character(colnames(df_ffl)), method = c("radix"))]


#Stat test to know if plastic and non plastic are different from each other
rstatix::row_wise_fisher_test(t(df_ffl[c(3,4),]))

# ##Each motif topology
# barplot(t(df_ffl), main="A) FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
#         space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))

################################################################################
##DMD motifs Data###############################################################
control_np <- read.csv("scripts/data/nonplast_E_coli_random2_diamond.csv", sep = ",")
control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_diamond.csv", sep = ",")
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_DMD.csv", sep = ",")[,1:13]
drift <- read.csv("scripts/data/full_netw_drift_DMD.csv", sep = ",")[,1:13]
s_shuffled <- read.csv("scripts/data/full_netw_shuffled_DMD.csv", sep = ",")[,1:13]

df1 <- as.data.frame(rbind(colSums(control_np[,5:14])*100/sum(control_np[,3]), colSums(control_pl[,5:14])*100/sum(control_pl[,3]),
                              colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                              theory[2,4:13]*100/theory[2,2],theory[1,4:13]*100/theory[1,2],
                              s_shuffled[2,4:13]*100/s_shuffled[2,2],s_shuffled[1,4:13]*100/s_shuffled[1,2],
                              drift[2,4:13]*100/drift[2,2],drift[1,4:13]*100/drift[1,2]))
rownames(df1) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nShuffled","Plastic\nShuffled", "Non-plastic\nDrift","Plastic\nDrift")
df1 <- df1[,order(as.character(colnames(df1)), method = c("radix"))]


#Stat test to know if plastic and non plastic are different from each other
rstatix::row_wise_fisher_test(t(df1[c(3,4),]))

# #Each motif topology
# barplot(t(df1[1:10]), main="B) Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
#         space=c(0.3,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))

##FIGURE########################################################################
pdf(paste0("figures/Fig3_motifs_control",".pdf"), width=5, height=8)
par(mar = c(8,3, 2,1), mfrow=c(2,1))

barplot(t(df_ffl), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="FFL motifs", width = c(0.6, 0.6,1.1,1.1,1.1,1.1,0.6,0.6,0.6,0.6))
#polygon(x=c(7.6, 7.6, 2.6, 2.6),  y=c(103, -14, -14, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(0.1, 0.1, 5, 5),  y=c(103, -70, -70, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(4.28, 4.28, 10.5, 10.5),  y=c(103, -70, -70, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.3, -60, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.2, -60, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df_ffl), main="FFL motif distribution", col=c("blue4","blue","dodgerblue3","skyblue3","lightskyblue","deepskyblue","cadetblue2","cyan3"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1), add=TRUE,
        srt=45, xaxt = "n", width = c(0.6, 0.6,1.1,1.1,1.1,1.1,0.6,0.6,0.6,0.6))
title(ylab = "%", line=2.1)
text(x = c(0.5,1.2,2.4,3.5,5,6.1,7.3,8,8.9,9.6), y = par("usr")[3] - 0.45,
     labels = c("Non-Plastic", "Plastic", "Non-Plastic", "Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic"),
     xpd = NA, srt = 90, cex = 1, adj=1.05)
text(0.9, 105, substitute(paste("Shuffled")), xpd=TRUE, font=3, cex=0.85)
text(7.7, 105, substitute(paste("Shuffled")), xpd=TRUE, font=3, cex=0.85)
text(9.2, 105, substitute(paste("Drift")), xpd=TRUE, font=3, cex=0.85)
text(10.19, 9, substitute(paste("C1")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 22, substitute(paste("C2")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 33, substitute(paste("C3")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 43, substitute(paste("C4")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 59, substitute(paste("I1")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 69, substitute(paste("I2")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 78, substitute(paste("I3")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 91, substitute(paste("I4")), xpd=TRUE, font=3, col="black", cex=0.7)

barplot(t(df1[1:10]), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="DMD motifs", width = c(0.6, 0.6,1.1,1.1,1.1,1.1,0.6,0.6,0.6,0.6))
polygon(x=c(0.1, 0.1, 5, 5),  y=c(103, -70, -70, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(4.28, 4.28, 10.5, 10.5),  y=c(103, -70, -70, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.3, -60, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.2, -60, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("#E65100","#FF6C00","#FB8C00","#FFB300","#FFCC80","#FDD835","#FFEE58","#FFF59D","#EEFF91","#FFFDE7"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), xaxt="n", add=TRUE, width = c(0.6, 0.6,1.1,1.1,1.1,1.1,0.6,0.6,0.6,0.6))
title(ylab = "%", line=2.1)
text(x = c(0.5,1.2,2.4,3.5,5,6.1,7.3,8,8.9,9.6), y = par("usr")[3] - 0.45,
     labels = c("Non-Plastic", "Plastic", "Non-Plastic", "Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic"),
     xpd = NA, srt = 90, cex = 1, adj=1.05)
text(0.9, 105, substitute(paste("Shuffled")), xpd=TRUE, font=3, cex=0.85)
text(7.7, 105, substitute(paste("Shuffled")), xpd=TRUE, font=3, cex=0.85)
text(9.2, 105, substitute(paste("Drift")), xpd=TRUE, font=3, cex=0.85)
text(10.19, 8, substitute(paste("MM1")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 21, substitute(paste("MM2")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 30, substitute(paste("MN")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 41, substitute(paste("MP")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 54, substitute(paste("NM")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 62, substitute(paste("NN")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 68, substitute(paste("NP")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 78, substitute(paste("PM")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 87, substitute(paste("PN")), xpd=TRUE, font=3, col="black", cex=0.7)
text(10.19, 95, substitute(paste("PP")), xpd=TRUE, font=3, col="black", cex=0.7)

dev.off()

# 
# pdf(paste0("figures/Fig2_motifs_control_test",".pdf"), width=5.5, height=8)
# par(mar = c(8,3, 2,1), mfrow=c(2,1))
# 
# barplot(t(df_ffl), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="FFL motif frequency")
# #polygon(x=c(7.6, 7.6, 2.6, 2.6),  y=c(103, -14, -14, 103),  col="honeydew2", border=NA, xpd=TRUE)
# polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -70, -70, 103),  col="honeydew2", border=NA, xpd=TRUE)
# polygon(x=c(5.1, 5.1, 12.59, 12.59),  y=c(103, -70, -70, 103),  col="cornsilk", border=NA, xpd=TRUE)
# text(2.5, -60, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
# text(9, -60, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
# barplot(t(df_ffl), main="FFL motif distribution", col=c("blue4","blue","dodgerblue3","skyblue3","lightskyblue","deepskyblue","cadetblue2","cyan3"),
#         space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1), add=TRUE,
#         srt=45, xaxt = "n")
# title(ylab = "Frequency %", line=2.1)
# text(x = c(0.9,2,3.4,4.4,5.9,6.9,8.3,9.4,10.8,11.9), y = par("usr")[3] - 0.45,
#      labels = c("Non-Plastic", "Plastic", "Non-Plastic", "Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic", "Non-plastic","Plastic"),
#      xpd = NA, srt = 80, cex = 1, adj=1.05)
# text(1.3, 105, substitute(paste(bold("Control"))), xpd=TRUE, font=3, cex=0.85)
# text(8.6, 105, substitute(paste(bold("Control"))), xpd=TRUE, font=3, cex=0.85)
# text(11.3, 105, substitute(paste(bold("Drift"))), xpd=TRUE, font=3, cex=0.85)
# text(11.75, 10, "C1", xpd=TRUE, font=3, col="white")
# text(11.75, 25, "C2", xpd=TRUE, font=3, col="white")
# text(11.75, 39, "C3", xpd=TRUE, font=3, col="white")
# text(11.75, 47, "C4", xpd=TRUE, font=3, col="white")
# text(11.75, 55, "I1", xpd=TRUE, font=3, col="black")
# text(11.75, 65, "I2", xpd=TRUE, font=3, col="black")
# text(11.75, 80, "I3", xpd=TRUE, font=3, col="black")
# text(11.75, 95, "I4", xpd=TRUE, font=3, col="black")
# 
# barplot(t(df1[1:10]), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="DMD motif frequency")
# polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -70, -70, 103),  col="honeydew2", border=NA, xpd=TRUE)
# polygon(x=c(5.1, 5.1, 12.59, 12.59),  y=c(103, -70, -70, 103),  col="cornsilk", border=NA, xpd=TRUE)
# text(2.5, -60, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
# text(9, -60, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
# barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("chocolate", "tomato2", "indianred1","lightsalmon1","orange","darkgoldenrod1","gold","yellow2","khaki1","cornsilk"),
#         space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), xaxt="n", add=TRUE)
# title(ylab = "Frequency %", line=2.1)
# text(x = c(0.9,2,3.4,4.4,5.9,6.9,8.3,9.4,10.8,11.9), y = par("usr")[3] - 0.45,
#      labels = c("Non-Plastic\nControl", "Plastic\nControl", "Non-Plastic", "Plastic", "Non-plastic","Plastic", "Non-plastic\nControl","Plastic\nControl", "Non-plastic\nDrift","Plastic\nDrift"),
#      xpd = NA, srt = 80, cex = 1, adj=1.05)
# text(11.75, 5, "MM1", xpd=TRUE, font=3, col="black")
# text(11.75, 19, "MM2", xpd=TRUE, font=3, col="black")
# text(11.75, 29, "MN", xpd=TRUE, font=3, col="black")
# text(11.75, 40, "MP", xpd=TRUE, font=3, col="black")
# text(11.75, 52, "NM", xpd=TRUE, font=3, col="black")
# text(12.6, 60, "NN", xpd=TRUE, font=3, col="black")
# text(11.75, 65, "NP", xpd=TRUE, font=3, col="black")
# text(11.75, 78, "PM", xpd=TRUE, font=3, col="black")
# text(12.6, 87, "PN", xpd=TRUE, font=3, col="black")
# text(11.75, 95, "PP", xpd=TRUE, font=3, col="black")
# 
# dev.off()

################################################################################
#PCA
#Empiric
non_plast <- cbind(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,c(5:12)],  read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")[,c(5:14)],  read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")[,c(7:10)])
all_plast <- cbind(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")[,c(5:14)], read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")[,c(7:10)])
drift_np <- cbind(read.csv("scripts/data/nonplast_E_coli_random2_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/nonplast_E_coli_random2_diamond.csv", sep = ",")[,c(5:14)], read.csv("scripts/data/nonplast_E_coli_random2_FBL.csv", sep = ",")[,c(7:10)])
drift_pl <- cbind(read.csv("scripts/data/plast_genes_E_coli_random2_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/plast_genes_E_coli_random2_diamond.csv", sep = ",")[,c(5:14)], read.csv("scripts/data/plast_genes_E_coli_random2_FBL.csv", sep = ",")[,c(7:10)])

#FFL
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/nonplast_E_coli_random2_nffl.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/plast_genes_E_coli_random2_nffl.csv", sep = ",")[,3]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
drift_np[,1:8] <- drift_np[,1:8]*nbrloop_emdrnp
drift_pl[,1:8] <- drift_pl[,1:8]*nbrloop_emdrpl
#DMD
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/nonplast_E_coli_random2_nDMD.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/plast_genes_E_coli_random2_nDMD.csv", sep = ",")[,3]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
drift_np[,9:18] <- drift_np[,9:18]*nbrloop_emdrnp
drift_pl[,9:18] <- drift_pl[,9:18]*nbrloop_emdrpl
#FBL
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/nonplast_E_coli_random2_nFBL.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/plast_genes_E_coli_random2_nFBL.csv", sep = ",")[,3]
non_plast[,19:22] <- non_plast[,19:22]*nbrloop_emnp
all_plast[,19:22] <- all_plast[,19:22]*nbrloop_empl
drift_np[,19:22] <- drift_np[,19:22]*nbrloop_emdrnp
drift_pl[,19:22] <- drift_pl[,19:22]*nbrloop_emdrpl

all_plast <- all_plast[which(rowSums2(as.matrix(all_plast))!=0),]
non_plast <- non_plast[which(rowSums2(as.matrix(non_plast))!=0),]
drift_np <- drift_np[which(rowSums2(as.matrix(drift_np))!=0),]
drift_pl <- drift_pl[which(rowSums2(as.matrix(drift_pl))!=0),]

df.ACPe <- as.data.frame(rbind(colMeans2(as.matrix(all_plast[,1:22]*100/rowSums2(as.matrix(all_plast[,1:22])))),
                               colMeans2(as.matrix(non_plast[,1:22]*100/rowSums2(as.matrix(non_plast[,1:22])))),
                               colMeans2(as.matrix(drift_np[,1:22]*100/rowSums2(as.matrix(drift_np[,1:22])))),
                               colMeans2(as.matrix(drift_pl[,1:22]*100/rowSums2(as.matrix(drift_pl[,1:22]))))))
row.names(df.ACPe) <- c("EcoPl", "EcoNP", "ShuffledNP","ShuffledPl")
df.ACPe$Genes <- c("Plastic","Non plastic","Shuffled Non Plastic","Shuffled Plastic")




##Simulations #Add shuffle
non_plast <- cbind(read.csv("scripts/data/full_netw_control_FFL.csv", sep = ",")[,c(4:11)],  read.csv("scripts/data/full_netw_control_DMD.csv", sep = ",")[,c(4:13)],  read.csv("scripts/data/full_netw_control_FBL.csv", sep = ",")[,c(6:9)])
all_plast <- cbind(read.csv("scripts/data/full_netw_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_plast_DMD.csv", sep = ",")[,c(4:13)], read.csv("scripts/data/full_netw_plast_FBL.csv", sep = ",")[,c(6:9)])
drift_np <- cbind(read.csv("scripts/data/full_netw_drift_control_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_drift_control_DMD.csv", sep = ",")[,c(4:13)], read.csv("scripts/data/full_netw_drift_control_FBL.csv", sep = ",")[,c(6:9)])
drift_pl <- cbind(read.csv("scripts/data/full_netw_drift_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_drift_plast_DMD.csv", sep = ",")[,c(4:13)], read.csv("scripts/data/full_netw_drift_plast_FBL.csv", sep = ",")[,c(6:9)])
shuff_np <- cbind(read.csv("scripts/data/full_netw_shuffled_control_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_shuffled_control_DMD.csv", sep = ",")[,c(4:13)], read.csv("scripts/data/full_netw_shuffled_control_FBL.csv", sep = ",")[,c(6:9)])
shuff_pl <- cbind(read.csv("scripts/data/full_netw_shuffled_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_shuffled_plast_DMD.csv", sep = ",")[,c(4:13)], read.csv("scripts/data/full_netw_shuffled_plast_FBL.csv", sep = ",")[,c(6:9)])

#FFL
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")[,2]
nbrloop_emdrnp <- read.csv("scripts/data/full_netw_drift_NP_nbrFFL.csv", sep = ",")[,2]
nbrloop_emdrpl <- read.csv("scripts/data/full_netw_drift_Pl_nbrFFL.csv", sep = ",")[,2]
nbrloop_emshnp <- read.csv("scripts/data/full_netw_shuffled_NP_nbrFFL.csv", sep = ",")[,2]
nbrloop_emshpl <- read.csv("scripts/data/full_netw_shuffled_Pl_nbrFFL.csv", sep = ",")[,2]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
drift_np[,1:8] <- drift_np[,1:8]*nbrloop_emdrnp
drift_pl[,1:8] <- drift_pl[,1:8]*nbrloop_emdrpl
shuff_np[,1:8] <- shuff_np[,1:8]*nbrloop_emshnp
shuff_pl[,1:8] <- shuff_pl[,1:8]*nbrloop_emshpl
#DMD
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")[,2]
nbrloop_emdrnp <- read.csv("scripts/data/full_netw_drift_NP_nbrDMD.csv", sep = ",")[,2]
nbrloop_emdrpl <- read.csv("scripts/data/full_netw_drift_Pl_nbrDMD.csv", sep = ",")[,2]
nbrloop_emshnp <- read.csv("scripts/data/full_netw_shuffled_NP_nbrDMD.csv", sep = ",")[,2]
nbrloop_emshpl <- read.csv("scripts/data/full_netw_shuffled_Pl_nbrDMD.csv", sep = ",")[,2]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
drift_np[,9:18] <- drift_np[,9:18]*nbrloop_emdrnp
drift_pl[,9:18] <- drift_pl[,9:18]*nbrloop_emdrpl
shuff_np[,9:18] <- shuff_np[,9:18]*nbrloop_emshnp
shuff_pl[,9:18] <- shuff_pl[,9:18]*nbrloop_emshpl
#FBL
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/full_netw_drift_NP_nbrFBL.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/full_netw_drift_Pl_nbrFBL.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/full_netw_shuffled_NP_nbrFBL.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/full_netw_shuffled_Pl_nbrFBL.csv", sep = ",")[,3]
non_plast[,19:22] <- non_plast[,19:22]*nbrloop_emnp
all_plast[,19:22] <- all_plast[,19:22]*nbrloop_empl
drift_np[,19:22] <- drift_np[,19:22]*nbrloop_emdrnp
drift_pl[,19:22] <- drift_pl[,19:22]*nbrloop_emdrpl
shuff_np[,19:22] <- shuff_np[,19:22]*nbrloop_emshnp
shuff_pl[,19:22] <- shuff_pl[,19:22]*nbrloop_emshpl
all_plast <- all_plast[which(rowSums2(as.matrix(all_plast))!=0),]
non_plast <- non_plast[which(rowSums2(as.matrix(non_plast))!=0),]
drift_np <- drift_np[which(rowSums2(as.matrix(drift_np))!=0),]
drift_pl <- drift_pl[which(rowSums2(as.matrix(drift_pl))!=0),]
shuff_np <- shuff_np[which(rowSums2(as.matrix(shuff_np))!=0),]
shuff_pl <- shuff_pl[which(rowSums2(as.matrix(shuff_pl))!=0),]

#What is needed : frequency by row (sum of every roww = 1) ; and mean of the frequency.

df.ACPt <- as.data.frame(rbind(colMeans2(as.matrix(all_plast[,1:22]*100/rowSums2(as.matrix(all_plast[,1:22])))),
                 colMeans2(as.matrix(non_plast[,1:22]*100/rowSums2(as.matrix(non_plast[,1:22])))),
                 colMeans2(as.matrix(drift_np[,1:22]*100/rowSums2(as.matrix(drift_np[,1:22])))),
                 colMeans2(as.matrix(drift_pl[,1:22]*100/rowSums2(as.matrix(drift_pl[,1:22])))),
                 colMeans2(as.matrix(shuff_np[,1:22]*100/rowSums2(as.matrix(shuff_np[,1:22])))),
                 colMeans2(as.matrix(shuff_pl[,1:22]*100/rowSums2(as.matrix(shuff_pl[,1:22]))))))
row.names(df.ACPt) <- c("SimuPl", "SimuNP", "DriftNP","DriftPl","ShuffNP","ShuffPL")
df.ACPt$Genes <- c("Plastic","Non plastic","Non plastic", "Plastic","Shuffled Non Plastic","Shuffled Plastic")



########################
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(ggfortify)
library(ggplot2)

df.ACP <- rbind(df.ACPe,df.ACPt)
df.ACP$Datatype <- c(rep("1E.coli",4),rep("2Simulation",2),rep("3Drift",2),rep("2Simulation",2))
colnames(df.ACP)[19:22] <- c("FBL2","FBL3","FBL4","FBL5")
pca_rest <- prcomp(df.ACP[,1:22], scale. = FALSE)



gg <- autoplot(pca_rest, data = df.ACP, colour = 'Genes', shape="Datatype", size = 4, label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'lightgrey',alpha=1, cex=.8,
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "darkgrey")+  
  ggtitle("PCA on motif proportions")+theme_bw() +
  scale_color_manual(values = c("red", "black","rosybrown1","darkgrey"), guide = guide_legend(title.position = "left", nrow = 2))+
  theme(legend.position = c(0.43,-0.3), legend.box = "vertical", legend.direction = "horizontal", plot.margin = unit(c(0.1,0.1,3,0.1), "cm"))+
  scale_shape_manual(labels = c("E.coli", "Simulation","Drift"), values = c(15,17,24), guide = guide_legend(title.position = "left"))



pdf(paste0("figures/Figure3B_ACP",".pdf"), width=4, height=5)
gg
dev.off()




