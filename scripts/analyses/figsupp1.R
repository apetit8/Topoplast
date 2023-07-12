#Supplementary figures
source("scripts/functions/functions.R")
pdfname <- "figures/fig_supp"
################################################################################
#Different parameters: 20000 generations ; 30 genes ; with self-regulation
################################################################################
####FFL prop####################################################################
theory1 <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]
theory2 <- read.csv("scripts/data/10g_selfreg_FFL.csv", sep = ",")[c(1,5),1:11]
theory3 <- read.csv("scripts/data/10g_20k_FFL.csv", sep = ",")[c(1,5),1:11]
theory4 <- read.csv("scripts/data/30g_FFL.csv", sep = ",")[c(1,5),1:11]


df_ffl <- as.data.frame(rbind(theory1[1,2:11], theory1[2,2:11], theory2[1,2:11], theory2[2,2:11], theory3[1,2:11], theory3[2,2:11], theory4[1,2:11], theory4[2,2:11]))
rownames(df_ffl) <- c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g")
#
theory1 <- read.csv("scripts/data/10g_diamond.csv", sep = ",")[c(1,5),1:11]
theory2 <- read.csv("scripts/data/10g_selfreg_diamond.csv", sep = ",")[c(1,5),1:11]
theory3 <- read.csv("scripts/data/10g_20k_diamond.csv", sep = ",")[c(1,5),1:11]
theory4 <- read.csv("scripts/data/30g_diamond.csv", sep = ",")[c(1,5),1:11]

df_dmd <- as.data.frame(rbind(theory1[1,2:11], theory1[2,2:11], theory2[1,2:11], theory2[2,2:11], theory3[1,2:11], theory3[2,2:11], theory4[1,2:11], theory4[2,2:11]))
rownames(df_dmd) <- c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g")

#
theory1 <- read.csv("scripts/data/10g_FBL.csv", sep = ",")[c(1,5),1:11]
theory2 <- read.csv("scripts/data/10g_selfreg_FBL.csv", sep = ",")[c(1,5),1:11]
theory3 <- read.csv("scripts/data/10g_20k_FBL.csv", sep = ",")[c(1,5),1:11]
theory4 <- read.csv("scripts/data/30g_FBL.csv", sep = ",")[c(1,5),1:11]

df_fbl <- as.data.frame(rbind(theory1[1,2:11], theory1[2,2:11], theory2[1,2:11], theory2[2,2:11], theory3[1,2:11], theory3[2,2:11], theory4[1,2:11], theory4[2,2:11]))
rownames(df_fbl) <- c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g")

  
pdf(paste0(pdfname,"_1",".pdf"), width=10, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df_ffl[,1:2]), main = "FFL motif enrichment", col=c("gold", "grey"), space=c(0.0,0.1,0.35,0.1,0.35,0.1,0.35,0.1), legend.text = c( "at least 1 FFL","No FFL"), 
      args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
barplot(t(df_dmd[,1:2]), main = "DMD motif enrichment", col=c("palegreen", "grey"), space=c(0.0,0.1,0.35,0.1,0.35,0.1,0.35,0.1), legend.text = c( "at least 1 DMD","No DMD"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
barplot(t(df_fbl[,1:2]), main = "FBL motif enrichment", col=c("olivedrab1", "grey"), space=c(0.0,0.1,0.35,0.1,0.35,0.1,0.35,0.1), legend.text = c( "at least 1 FBL","No FBL"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
dev.off()


################################################################################
#Whisker plot###############################
nbrloop_1 <- read.csv("scripts/data/10g_Pl_nbrffl.csv", sep = ",")
nbrloop_1$Type <- "g1"
nbrloop_2 <- read.csv("scripts/data/10g_NP_nbrffl.csv", sep = ",")
nbrloop_2$Type <- "g2"
nbrloop_3 <- read.csv("scripts/data/10g_selfreg_Pl_nbrffl.csv", sep = ",")
nbrloop_3$Type <- "g3"
nbrloop_4 <- read.csv("scripts/data/10g_selfreg_NP_nbrffl.csv", sep = ",")
nbrloop_4$Type <- "g4"
nbrloop_5 <- read.csv("scripts/data/10g_20k_Pl_nbrffl.csv", sep = ",")
nbrloop_5$Type <- "g5"
nbrloop_6 <- read.csv("scripts/data/10g_20k_NP_nbrffl.csv", sep = ",")
nbrloop_6$Type <- "g6"
nbrloop_7 <- read.csv("scripts/data/30g_Pl_nbrffl.csv", sep = ",")
nbrloop_7$Type <- "g7"
nbrloop_8 <- read.csv("scripts/data/30g_NP_nbrffl.csv", sep = ",")
nbrloop_8$Type <- "g8"

df_ffl <- rbind(nbrloop_1, nbrloop_2, nbrloop_3, nbrloop_4,nbrloop_5, nbrloop_6, nbrloop_7, nbrloop_8)
df_ffl <- subset(df_ffl, Loop_number != 0)
#
nbrloop_1 <- read.csv("scripts/data/10g_Pl_nbrDMD.csv", sep = ",")
nbrloop_1$Type <- "g1"
nbrloop_2 <- read.csv("scripts/data/10g_NP_nbrDMD.csv", sep = ",")
nbrloop_2$Type <- "g2"
nbrloop_3 <- read.csv("scripts/data/10g_selfreg_Pl_nbrDMD.csv", sep = ",")
nbrloop_3$Type <- "g3"
nbrloop_4 <- read.csv("scripts/data/10g_selfreg_NP_nbrDMD.csv", sep = ",")
nbrloop_4$Type <- "g4"
nbrloop_5 <- read.csv("scripts/data/10g_20k_Pl_nbrDMD.csv", sep = ",")
nbrloop_5$Type <- "g5"
nbrloop_6 <- read.csv("scripts/data/10g_20k_NP_nbrDMD.csv", sep = ",")
nbrloop_6$Type <- "g6"
nbrloop_7 <- read.csv("scripts/data/30g_Pl_nbrDMD.csv", sep = ",")
nbrloop_7$Type <- "g7"
nbrloop_8 <- read.csv("scripts/data/30g_NP_nbrDMD.csv", sep = ",")
nbrloop_8$Type <- "g8"

df_DMD <- rbind(nbrloop_1, nbrloop_2, nbrloop_3, nbrloop_4,nbrloop_5, nbrloop_6, nbrloop_7, nbrloop_8)
df_DMD <- subset(df_DMD, Loop_number != 0)
#
nbrloop_1 <- read.csv("scripts/data/10g_Pl_nbrFBL.csv", sep = ",")
nbrloop_1$Type <- "g1"
nbrloop_2 <- read.csv("scripts/data/10g_NP_nbrFBL.csv", sep = ",")
nbrloop_2$Type <- "g2"
nbrloop_3 <- read.csv("scripts/data/10g_selfreg_Pl_nbrFBL.csv", sep = ",")
nbrloop_3$Type <- "g3"
nbrloop_4 <- read.csv("scripts/data/10g_selfreg_NP_nbrFBL.csv", sep = ",")
nbrloop_4$Type <- "g4"
nbrloop_5 <- read.csv("scripts/data/10g_20k_Pl_nbrFBL.csv", sep = ",")
nbrloop_5$Type <- "g5"
nbrloop_6 <- read.csv("scripts/data/10g_20k_NP_nbrFBL.csv", sep = ",")
nbrloop_6$Type <- "g6"
nbrloop_7 <- read.csv("scripts/data/30g_Pl_nbrFBL.csv", sep = ",")
nbrloop_7$Type <- "g7"
nbrloop_8 <- read.csv("scripts/data/30g_NP_nbrFBL.csv", sep = ",")
nbrloop_8$Type <- "g8"

df_FBL <- rbind(nbrloop_1, nbrloop_2, nbrloop_3, nbrloop_4,nbrloop_5, nbrloop_6, nbrloop_7, nbrloop_8)
df_FBL <- subset(df_FBL, FBL_number != 0)




pdf(paste0(pdfname,"_2",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(log(df_ffl$Loop_number) ~ df_ffl$Type,
        at = c(1, 1.9, 3, 3.9, 5, 5.9, 7, 7.9), ylab = "Log of FFL number per gene", xlab = "",
        names = c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g"),
        tck=-0.1, las = 1, col = c("orange","red"), border = "brown")
boxplot(log(df_DMD$Loop_number) ~ df_DMD$Type,
        at = c(1, 1.9, 3, 3.9, 5, 5.9, 7, 7.9), ylab = "Log of DMD number per gene", xlab = "",
        names = c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g"),
        tck=-0.1, las = 1, col = c("orange","red"), border = "brown")
boxplot(log(df_FBL$FBL_number) ~ df_FBL$Type,
        at = c(1, 1.9, 3, 3.9, 5, 5.9, 7, 7.9), ylab = "Log of FBL number per gene", xlab = "",
        names = c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g"),
        tck=-0.1, las = 1, col = c("orange","red"), border = "brown")
dev.off()



################################################################################
####FFL prop####################################################################
theory1 <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]
theory2 <- read.csv("scripts/data/10g_selfreg_FFL.csv", sep = ",")[c(1,5),1:11]
theory3 <- read.csv("scripts/data/10g_20k_FFL.csv", sep = ",")[c(1,5),1:11]
theory4 <- read.csv("scripts/data/30g_FFL.csv", sep = ",")[c(1,5),1:11]

df_ffl <- as.data.frame(rbind(theory1[1,4:11]*100/theory1[1,2], theory1[2,4:11]*100/theory1[2,2], theory2[1,4:11]*100/theory2[1,2], theory2[2,4:11]*100/theory2[2,2],
                              theory3[1,4:11]*100/theory3[1,2], theory3[2,4:11]*100/theory3[2,2], theory4[1,4:11]*100/theory4[1,2], theory4[2,4:11]*100/theory4[2,2]))
rownames(df_ffl) <- c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g")
#
theory1 <- read.csv("scripts/data/10g_diamond.csv", sep = ",")[c(1,5),1:13]
theory2 <- read.csv("scripts/data/10g_selfreg_diamond.csv", sep = ",")[c(1,5),1:13]
theory3 <- read.csv("scripts/data/10g_20k_diamond.csv", sep = ",")[c(1,5),1:13]
theory4 <- read.csv("scripts/data/30g_diamond.csv", sep = ",")[c(1,5),1:13]

df_dmd <- as.data.frame(rbind(theory1[1,4:13]*100/theory1[1,2], theory1[2,4:13]*100/theory1[2,2], theory2[1,4:13]*100/theory2[1,2], theory2[2,4:13]*100/theory2[2,2],
                              theory3[1,4:13]*100/theory3[1,2], theory3[2,4:13]*100/theory3[2,2], theory4[1,4:13]*100/theory4[1,2], theory4[2,4:13]*100/theory4[2,2]))
rownames(df_dmd) <- c("Plastic\nDefault", "Non-plastic\nDefault","Plastic\nSelfreg","Non-Plastic\nSelfreg","Plastic\n20k","Non-Plastic\n20k","Plastic\n30g","Non-Plastic\n30g")
#

pdf(paste0(pdfname,"_3",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
barplot(t(df_ffl), main="FFL motif distribution", col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = colnames(df_ffl), args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))
barplot(t(df_dmd), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1","plum1","lightsalmon1","indianred1","lightpink","chocolate1"),
        space=c(0.3,0.1,0.4,0.1,0.4,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))
dev.off()










