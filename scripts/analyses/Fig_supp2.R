source("scripts/functions/functions.R")
pdfname <- "figures/fig_full_netw"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
##FFL motifs####################################################################
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

##Each motif topology
barplot(t(df_ffl), main="A) FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))
#dev.off()

################################################################################
##DMD motifs####################################################################
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


#Each motif topology
barplot(t(df1[1:10]), main="B) Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
        space=c(0.3,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))




pdf(paste0("figures/Figsupp2_motifs_control",".pdf"), width=10, height=10)
par(mar = c(8,2, 2,1), mfrow=c(2,1))

barplot(t(df_ffl), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="FFL motif distribution")
text(4, 9, expression(hat(beta) == (X^t * X)^{-1} * X^t * y))
#polygon(x=c(7.6, 7.6, 2.6, 2.6),  y=c(103, -14, -14, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -20, -20, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(5.1, 5.1, 10, 10),  y=c(103, -20, -20, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.5, -15, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.5, -15, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df_ffl), main="FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.3, 1.22)), add=TRUE
        , names.arg=c("Non-Plastic\nControl", "Plastic\nControl", "Non-Plastic\n", "Plastic\n", "Non-plastic\n","Plastic\n", "Non-plastic\nControl","Plastic\nControl"))

barplot(t(df1[1:10]), col = NA, border = NA, axes = FALSE, ylim = c(0,100), xaxt = "n", main="DMD motif distribution")
polygon(x=c(0.1, 0.1, 5.5, 5.5),  y=c(103, -20, -20, 103),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(5.1, 5.1, 10, 10),  y=c(103, -20, -20, 103),  col="cornsilk", border=NA, xpd=TRUE)
text(2.5, -15, substitute(paste(bold("E. coli"))), xpd=TRUE, font=3)
text(7.5, -15, substitute(paste(bold("Simulations"))), xpd=TRUE, font=3)
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
                 space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(0.27, 1.22)), add=TRUE,
        names.arg=c("Non-Plastic\nControl", "Plastic\nControl", "Non-Plastic\n", "Plastic\n", "Non-plastic\n","Plastic\n", "Non-plastic\nControl","Plastic\nControl"))

dev.off()

################################################################################
# E. coli and simu in different plots
################################################################################
##FFL motifs####################################################################
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
rownames(df_ffl) <- c("E. coli\nrandom np", "E. coli\nrandom pl","Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nDrift","Plastic\nDrift")

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


##DMD motifs####
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
rownames(df1) <- c("E. coli\nrandom np", "E. coli\nrandom pl", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nDrift","Plastic\nDrift")


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





pdf(paste0("figures/Fig2_e_coli",".pdf"), width=10, height=4)
par(mar = c(8,2, 2,1), mfrow=c(1,1))
par(mar = c(2,2, 2,25))
barplot(as.matrix(c(0,0,0,0)), col = NA, border = NA, axes = FALSE, )
polygon(x=c(100, 100, 0, 0),
        y=c(100, 0, 0, 100),
        col="honeydew2", border=NA)
barplot(t(df_ffl[c(3,4,1,2),]), main="E. coli motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1), add=TRUE )
barplot(t(df1[c(3,4,1,2),]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
        space=c(5.5,0.1,0.4,0.1), add=TRUE)
dev.off()



pdf(paste0("figures/Fig2_simu",".pdf"), width=10, height=5)
par(mar = c(8,2, 2,25))
barplot(t(df_ffl[c(5,6,7,8),]), main="Simulation motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.2, 1.15)))
barplot(t(df1[c(5,6,7,8),]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
        space=c(5.5,0.1,0.4,0.1), legend.text = legendtext1, args.legend = list(ncol=4, x = "topright", inset = c(-1, 1.15)), add=TRUE)

dev.off()
 






