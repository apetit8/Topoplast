source("scripts/functions/functions.R")
pdfname <- "figures/fig_full_netw"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
##FFL motifs####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_FFL.csv", sep = ",")[,1:11]
theory03 <- read.csv("scripts/data/full_a0-3_FFL.csv", sep = ",")[,1:11]
theory05 <- read.csv("scripts/data/full_a0-5_FFL.csv", sep = ",")[,1:11]
theory07 <- read.csv("scripts/data/full_a0-7_FFL.csv", sep = ",")[,1:11]

df_ffl <- as.data.frame(rbind(colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                              theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2],
                              theory03[2,4:11]*100/theory03[2,2],theory03[1,4:11]*100/theory03[1,2],
                              theory05[2,4:11]*100/theory05[2,2],theory05[1,4:11]*100/theory05[1,2],
                        theory07[2,4:11]*100/theory07[2,2],theory07[1,4:11]*100/theory07[1,2]  ))
rownames(df_ffl) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nTheory 0.3","Plastic\nTheory 0.3", "Non-plastic\nTheory 0.5","Plastic\nTheory 0.5", "Non-plastic\nTheory 0.7","Plastic\nTheory 0.7")

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

# pdf(paste0(pdfname,"_FFL_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df_ffl), main="FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1, 0.4, 0.1, 0.4, 0.1, 0.4, 0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))
# dev.off()
################################################################################
##DMD motifs####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_DMD.csv", sep = ",")[,1:13]
theory03 <- read.csv("scripts/data/full_a0-3_DMD.csv", sep = ",")[,1:13]
theory05 <- read.csv("scripts/data/full_a0-5_DMD.csv", sep = ",")[,1:13]
theory07 <- read.csv("scripts/data/full_a0-7_DMD.csv", sep = ",")[,1:13]

df1 <- as.data.frame(rbind(colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                           theory[2,4:13]*100/theory[2,2], theory[1,4:13]*100/theory[1,2],
                           theory03[2,4:13]*100/theory03[2,2], theory03[1,4:13]*100/theory03[1,2],
                           theory05[2,4:13]*100/theory05[2,2], theory05[1,4:13]*100/theory05[1,2],
                           theory07[2,4:13]*100/theory07[2,2], theory07[1,4:13]*100/theory07[1,2] ))
rownames(df1) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nTheory 0.3","Plastic\nTheory 0.3", "Non-plastic\nTheory 0.5","Plastic\nTheory 0.5", "Non-plastic\nTheory 0.7","Plastic\nTheory 0.7")

#Stat test to know if plastic and non plastic are different from each other
df_anova <- rbind(subset(non_plast, FFL==1)[,5:14], subset(all_plast, FFL==1)[,5:14])
df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))

# colnames(df1)
legendtext <- c(paste0("PP", ifelse(anova(lm(PP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PM", ifelse(anova(lm(PM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PN", ifelse(anova(lm(PN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NP", ifelse(anova(lm(NP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NM", ifelse(anova(lm(NM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NN", ifelse(anova(lm(NN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MP", ifelse(anova(lm(MP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM2", ifelse(anova(lm(MM2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM1", ifelse(anova(lm(MM1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MN", ifelse(anova(lm(MN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))

# pdf(paste0(pdfname,"_DMD_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
        space=c(0.3,0.1,0.4,0.1, 0.4, 0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))
# dev.off()
