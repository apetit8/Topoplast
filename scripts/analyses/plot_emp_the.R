source("scripts/functions/functions.R")
pdfname <- "figures/fig_VS"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
####FFL prop
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]


df_ffl <- as.data.frame(rbind(theory[1,2:11], theory[2,2:11],
                           colSums(all_plast[,3:12])*100/nrow(all_plast), colSums(non_plast[,3:12])*100/nrow(non_plast) ))
rownames(df_ffl) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")

pdf(paste0(pdfname,"_prop_FFL",".pdf"), width=5, height=5)
par(mar = c(7,2, 2,1))
barplot(t(df_ffl[,1:2]), main = "FFL motif enrichment", col=c("gold", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "at least 1 FFL","No FFL"), 
        args.legend = list(x = "topright", inset = c(0.3, 1.2)))
dev.off()

##FFL motifs
df_ffl <- as.data.frame(rbind(theory[1,4:11]*100/theory[1,2],theory[2,4:11]*100/theory[2,2], colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                           colSums(non_plast[,5:12])*100/sum(non_plast[,3]) ))
rownames(df_ffl) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


#Stat test to know if plastic and non plastic are different from each other
#empiric
df_anova <- rbind(subset(non_plast, FFL==1)[,5:12], subset(all_plast, FFL==1)[,5:12])
df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))
#theory
theoryP <- read.csv("scripts/data/10g_plast_FFL.csv", sep = ",")[,1:11]
theoryNP <- read.csv("scripts/data/10g_control_FFL.csv", sep = ",")[,1:11]
df_anova2 <- rbind(subset(theoryNP, FFL==1)[,2:11], subset(theoryP, FFL==1)[,2:11])
df_anova2$Plasticity <- c(rep("No", nrow(subset(theoryNP, FFL==1))),rep("Yes", nrow(subset(theoryP, FFL==1))))

#colnames(df_ffl)
legendtext <- c(paste0("C3 ",ifelse(anova(lm(C3 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(C3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C1 ",ifelse(anova(lm(C1 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(C1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C2 ",ifelse(anova(lm(C2 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(C2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("C4 ",ifelse(anova(lm(C4 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(C4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I2 ",ifelse(anova(lm(I2 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(I2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I4 ",ifelse(anova(lm(I4 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(I4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I3 ",ifelse(anova(lm(I3 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(I3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("I1 ",ifelse(anova(lm(I1 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(I1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))


pdf(paste0(pdfname,"_FFL_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df_ffl), main="FFL motif distribution", col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))
dev.off()


###Same plot but for random empirical groups
df <- (rbind(non_plast, all_plast))

df1 <- df[sample(nrow(df), nrow(all_plast)), ]
df2 <- subset(df, !(V1 %in% df1$V1))

dfrandom <- as.data.frame(rbind(colSums(df1[,5:12])*100/nrow(df1), colSums(df2[,5:12])*100/nrow(df2) ))
rownames(dfrandom) <- c("Random 1\nEmpiric","Random2\nEmpiric")

pdf(paste0(pdfname,"_random_FFL_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(dfrandom), main="FFL motif distribution", col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1), legend.text = colnames(dfrandom), args.legend = list(ncol=4, x = "topright", inset = c(0.2, 1.2)))
dev.off()


################################################################################
#Do the same thing but for diamond
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/10g_diamond.csv", sep = ",")[c(1,5),1:13]

df1 <- as.data.frame(rbind(theory[1,2:13], theory[2,2:13], colSums(all_plast[,3:14])*100/nrow(all_plast),
                           colSums(non_plast[,3:14])*100/nrow(non_plast) ))
rownames(df1) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


pdf(paste0(pdfname,"_prop_DMD",".pdf"), width=5, height=5)
par(mar = c(7,2, 2,1))
barplot(t(df1[,1:2]), main = "Diamond motif enrichment", col=c("palegreen", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "at least 1 DMD","No DMD"), 
        args.legend = list(x = "topright", inset = c(0.3, 1.2)))
dev.off()

####Motifs
df1 <- as.data.frame(rbind(theory[1,4:13]*100/theory[1,2], theory[2,4:13]*100/theory[2,2], 
                           colSums(all_plast[,5:14])*100/sum(all_plast[,3]), colSums(non_plast[,5:14])*100/sum(non_plast[,3]) ))
rownames(df1) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")

#Stat test to know if plastic and non plastic are different from each other
df_anova <- rbind(subset(non_plast, FFL==1)[,5:14], subset(all_plast, FFL==1)[,5:14])
df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))
#theory
theoryP <- read.csv("scripts/data/10g_plast_diamond.csv", sep = ",")[,1:13]
theoryNP <- read.csv("scripts/data/10g_control_diamond.csv", sep = ",")[,1:13]
df_anova2 <- rbind(subset(theoryNP, FFL==1)[,2:13], subset(theoryP, FFL==1)[,2:13])
df_anova2$Plasticity <- c(rep("No", nrow(subset(theoryNP, FFL==1))),rep("Yes", nrow(subset(theoryP, FFL==1))))

# colnames(df1)
legendtext <- c(paste0("PP",ifelse(anova(lm(PP ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""), ifelse(anova(lm(PP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PM",ifelse(anova(lm(PM ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(PM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("PN",ifelse(anova(lm(PN ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(PN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NP",ifelse(anova(lm(NP ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(NP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NM",ifelse(anova(lm(NM ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(NM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("NN",ifelse(anova(lm(NN ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(NN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MP",ifelse(anova(lm(MP ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(MP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM2",ifelse(anova(lm(MM2 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(MM2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MM1",ifelse(anova(lm(MM1 ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(MM1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                paste0("MN",ifelse(anova(lm(MN ~ Plasticity, df_anova2))[1,5] <= pval, "(¨)", ""),ifelse(anova(lm(MN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))

pdf(paste0(pdfname,"_DMD_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1","plum1","lightsalmon1","indianred1","lightpink","chocolate1"),
        space=c(0.3,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))
dev.off()

################################################################################
#FBL
################################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FBL.csv", sep = ",")[c(1,5),]

df <- as.data.frame(rbind(theory[1,2:5], theory[2,2:5],
                           colSums(all_plast[,3:6])*100/nrow(all_plast), colSums(non_plast[,3:6])*100/nrow(non_plast) ))
rownames(df) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


pdf(paste0(pdfname,"_FBL",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df[,2:4]), main="Diamond motif distribution", col=c("grey","olivedrab1","palegreen"),
        space=c(0.3,0.1,0.4,0.1), legend.text = colnames(df[2:4]), args.legend = list(ncol=3, x = "topright", inset = c(0, 1.13)))
dev.off()

# #Stat test to know if FBL are different between non plastic and plastic genes
# df <- rbind(non_plast, all_plast)
# df$Plasticity <- c(rep("No", nrow(non_plast)),rep("Yes", nrow(all_plast)))
# anova(lm(FBL ~ Plasticity, df)) #Significant difference for FBL
# 
# anova(lm(Activating ~ Plasticity, subset(df, FBL==1) ))  #Almost significative ? Non significative enough ? for sign of FBL
# 



