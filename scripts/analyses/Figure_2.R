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
drift <- read.csv("scripts/data/Stab/full_drift_20000_FFL.csv", sep = ",")[,1:11]

df_ffl <- as.data.frame(rbind(colSums(control_np[,5:12])*100/sum(control_np[,3]), colSums(control_pl[,5:12])*100/sum(control_pl[,3]),
                              colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                              theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2],
                        drift[2,4:11]*100/drift[2,2],drift[1,4:11]*100/drift[1,2]))
rownames(df_ffl) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nDrift","Plastic\nDrift")

#Stat test to know if plastic and non plastic are different from each other
#empiric
df_anova <- rbind(subset(non_plast, FFL==1)[,5:12], subset(all_plast, FFL==1)[,5:12])
df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))

#colnames(df_ffl)
# t.test(subset(non_plast, FFL==1)$I2,  subset(all_plast, FFL==1)$I2)$p.value
legendtext <- c(paste0("C3 ", ifelse(anova(lm(C3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C1 ", ifelse(anova(lm(C1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C2 ", ifelse(anova(lm(C2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C4 ", ifelse(anova(lm(C4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I2 ", ifelse(anova(lm(I2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I4 ", ifelse(anova(lm(I4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I3 ", ifelse(anova(lm(I3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I1 ", ifelse(anova(lm(I1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))

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
drift <- read.csv("scripts/data/Stab/full_drift_20000_DMD.csv", sep = ",")[,1:13]

df1 <- as.data.frame(rbind(colSums(control_np[,5:14])*100/sum(control_np[,3]), colSums(control_pl[,5:14])*100/sum(control_pl[,3]),
                              colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                              theory[2,4:13]*100/theory[2,2],theory[1,4:13]*100/theory[1,2],
                              drift[2,4:13]*100/drift[2,2],drift[1,4:13]*100/drift[1,2]))
rownames(df1) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nDrift","Plastic\nDrift")


#Stat test to know if plastic and non plastic are different from each other
df_anova <- rbind(subset(non_plast, FFL==1)[,5:14], subset(all_plast, FFL==1)[,5:14])
df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))

# colnames(df1)
legendtext1 <- c(paste0("PP", ifelse(anova(lm(PP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PM", ifelse(anova(lm(PM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PN", ifelse(anova(lm(PN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NP", ifelse(anova(lm(NP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NM", ifelse(anova(lm(NM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NN", ifelse(anova(lm(NN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MP", ifelse(anova(lm(MP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM2", ifelse(anova(lm(MM2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM1", ifelse(anova(lm(MM1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MN", ifelse(anova(lm(MN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))


# #Each motif topology
# barplot(t(df1[1:10]), main="B) Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
#         space=c(0.3,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))

##FIGURE########################################################################

pdf(paste0("figures/Fig2_motifs_control",".pdf"), width=11, height=10)
par(mar = c(8,3, 2,1), mfrow=c(2,1))

barplot(t(df_ffl), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="FFL motif frequency")
#polygon(x=c(7.6, 7.6, 2.6, 2.6),  y=c(103, -14, -14, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -20, -20, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(5.1, 5.1, 10, 10),  y=c(103, -20, -20, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.5, -15, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.5, -15, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df_ffl), main="FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.3, 1.22)), add=TRUE
        , names.arg=c("Non-Plastic\nControl", "Plastic\nControl", "Non-Plastic\n", "Plastic\n", "Non-plastic\n","Plastic\n", "Non-plastic\nControl","Plastic\nControl"))
title(ylab = "Frequency %", line=2.1)

barplot(t(df1[1:10]), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="DMD motif frequency")
polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -20, -20, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(5.1, 5.1, 10, 10),  y=c(103, -20, -20, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.5, -15, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.5, -15, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
                 space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(0.27, 1.22)), add=TRUE,
        names.arg=c("Non-Plastic\nControl", "Plastic\nControl", "Non-Plastic\n", "Plastic\n", "Non-plastic\n","Plastic\n", "Non-plastic\nControl","Plastic\nControl"))
title(ylab = "Frequency %", line=2.1)

dev.off()

################################################################################
#PCA
#Empiric
non_plast <- cbind(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,c(5:12)],  read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")[,c(5:14)])
all_plast <- cbind(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")[,c(5:14)])
drift_np <- cbind(read.csv("scripts/data/nonplast_E_coli_random2_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/nonplast_E_coli_random2_diamond.csv", sep = ",")[,c(5:14)])
drift_pl <- cbind(read.csv("scripts/data/plast_genes_E_coli_random2_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/plast_genes_E_coli_random2_diamond.csv", sep = ",")[,c(5:14)])

#FFL
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/nonplast_E_coli_random2_nffl.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/plast_genes_E_coli_random2_nffl.csv", sep = ",")[,3]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
drift_np[,1:8] <- drift_np[,1:8]*nbrloop_emdrnp
drift_pl[,1:8] <- drift_pl[,1:8]*nbrloop_empl
#DMD
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emdrnp <- read.csv("scripts/data/nonplast_E_coli_random2_nDMD.csv", sep = ",")[,3]
nbrloop_emdrpl <- read.csv("scripts/data/plast_genes_E_coli_random2_nDMD.csv", sep = ",")[,3]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
drift_np[,9:18] <- drift_np[,9:18]*nbrloop_emdrnp
drift_pl[,9:18] <- drift_pl[,9:18]*nbrloop_empl

all_plast <- all_plast[which(rowSums2(as.matrix(all_plast))!=0),]
non_plast <- non_plast[which(rowSums2(as.matrix(non_plast))!=0),]
drift_np <- drift_np[which(rowSums2(as.matrix(drift_np))!=0),]
drift_pl <- drift_pl[which(rowSums2(as.matrix(drift_pl))!=0),]

df.ACPe <- as.data.frame(rbind(colMeans2(as.matrix(all_plast[,1:18]*100/rowSums2(as.matrix(all_plast[,1:18])))),
                               colMeans2(as.matrix(non_plast[,1:18]*100/rowSums2(as.matrix(non_plast[,1:18])))),
                               colMeans2(as.matrix(drift_np[,1:18]*100/rowSums2(as.matrix(drift_np[,1:18])))),
                               colMeans2(as.matrix(drift_pl[,1:18]*100/rowSums2(as.matrix(drift_pl[,1:18]))))))
row.names(df.ACPe) <- c("EcoPl", "EcoNP", "ControlNP","ContolPl")
df.ACPe$Genes <- c("E.coli Plastic","E.coli Non plastic","E.coli control Non plastic", "E. coli control Plastic")




##Simulations
non_plast <- cbind(read.csv("scripts/data/full_netw_control_FFL.csv", sep = ",")[,c(4:11)],  read.csv("scripts/data/full_netw_control_DMD.csv", sep = ",")[,c(4:13)])
all_plast <- cbind(read.csv("scripts/data/full_netw_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_plast_DMD.csv", sep = ",")[,c(4:13)])
drift_np <- cbind(read.csv("scripts/data/Stab/full_drift_20000_control_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/Stab/full_drift_20000_control_DMD.csv", sep = ",")[,c(4:13)])
drift_pl <- cbind(read.csv("scripts/data/Stab/full_drift_20000_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/Stab/full_drift_20000_plast_DMD.csv", sep = ",")[,c(4:13)])

#FFL
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")[,2]
nbrloop_emdrnp <- read.csv("scripts/data/Stab/full_drift_20000_NP_nbrFFL.csv", sep = ",")[,2]
nbrloop_emdrpl <- read.csv("scripts/data/Stab/full_drift_20000_Pl_nbrFFL.csv", sep = ",")[,2]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
drift_np[,1:8] <- drift_np[,1:8]*nbrloop_emdrnp
drift_pl[,1:8] <- drift_pl[,1:8]*nbrloop_empl
#DMD
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")[,2]
nbrloop_emdrnp <- read.csv("scripts/data/Stab/full_drift_20000_NP_nbrDMD.csv", sep = ",")[,2]
nbrloop_emdrpl <- read.csv("scripts/data/Stab/full_drift_20000_Pl_nbrDMD.csv", sep = ",")[,2]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
drift_np[,9:18] <- drift_np[,9:18]*nbrloop_emdrnp
drift_pl[,9:18] <- drift_pl[,9:18]*nbrloop_empl

all_plast <- all_plast[which(rowSums2(as.matrix(all_plast))!=0),]
non_plast <- non_plast[which(rowSums2(as.matrix(non_plast))!=0),]
drift_np <- drift_np[which(rowSums2(as.matrix(drift_np))!=0),]
drift_pl <- drift_pl[which(rowSums2(as.matrix(drift_pl))!=0),]

#What is needed : frequency by row (sum of every roww = 1) ; and mean of the frequency.

df.ACPt <- as.data.frame(rbind(colMeans2(as.matrix(all_plast[,1:18]*100/rowSums2(as.matrix(all_plast[,1:18])))),
                 colMeans2(as.matrix(non_plast[,1:18]*100/rowSums2(as.matrix(non_plast[,1:18])))),
                 colMeans2(as.matrix(drift_np[,1:18]*100/rowSums2(as.matrix(drift_np[,1:18])))),
                 colMeans2(as.matrix(drift_pl[,1:18]*100/rowSums2(as.matrix(drift_pl[,1:18]))))))
row.names(df.ACPt) <- c("SimuPl", "SimuNP", "DriftNP","DriftPl")
df.ACPt$Genes <- c("Simu Plastic","Simu Non plastic","Drift Non plastic", "Drift Plastic")

########################""
df.ACP <- rbind(df.ACPe,df.ACPt)
pca_rest <- prcomp(df.ACP[,1:18], scale. = FALSE)

g2 <- autoplot(pca_rest, data = df.ACP, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey',alpha=0, 
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("ACP on proportions")+theme_bw() +scale_color_manual(values = c("orange","cyan3","pink","green","darkorange","darkblue","darkred","darkgreen"))







