source("scripts/functions/functions.R")
pdfname <- "figures/fig_VS"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
####FFL prop####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]


df_ffl <- as.data.frame(rbind(theory[1,2:11], theory[2,2:11],
                           colSums(all_plast[,3:12])*100/nrow(all_plast), colSums(non_plast[,3:12])*100/nrow(non_plast) ))
rownames(df_ffl) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")

pdf(paste0(pdfname,"_prop_FFL",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df_ffl[,1:2]), main = "FFL motif enrichment", col=c("gold", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "at least 1 FFL","No FFL"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
dev.off()

#Whisker plot###############################
nbrloop_thpl <- read.csv("scripts/data/10g_Pl_nbrffl.csv", sep = ",")
nbrloop_thpl$Type <- "1thpl"
nbrloop_thnp <- read.csv("scripts/data/10g_NP_nbrffl.csv", sep = ",")
nbrloop_thnp$Type <- "2thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "3empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "4emnp"

df <- rbind(nbrloop_thpl, nbrloop_thnp, nbrloop_empl, nbrloop_emnp)
df <- subset(df, Loop_number != 0)

pdf(paste0(pdfname,"_wisker",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(df$Loop_number ~ df$Type,
        at = c(1,1.9,3,3.9), ylab = "Number of FBL per gene", xlab = "",
        names = c("Plastic\nprediction", "Non plastic\nprediction", "Plastic\nempiric", "Non plastic\nempiric"),
        tck=-0.1,
        las = 1,
        col = c("orange","red"),
        border = "brown"
)
dev.off()

################################################################################
####DMD prop####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/10g_diamond.csv", sep = ",")[c(1,5),1:13]

df1 <- as.data.frame(rbind(theory[1,2:13], theory[2,2:13], colSums(all_plast[,3:14])*100/nrow(all_plast),
                           colSums(non_plast[,3:14])*100/nrow(non_plast) ))
rownames(df1) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


pdf(paste0(pdfname,"_prop_DMD",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df1[,1:2]), main = "Diamond motif enrichment", col=c("palegreen", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "at least 1 DMD","No DMD"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
dev.off()


#Whisker plot#############################
nbrloop_thpl <- read.csv("scripts/data/10g_Pl_nbrDMD.csv", sep = ",")
nbrloop_thpl$Type <- "1thpl"
nbrloop_thnp <- read.csv("scripts/data/10g_NP_nbrDMD.csv", sep = ",")
nbrloop_thnp$Type <- "2thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "3empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "4emnp"

df <- rbind(nbrloop_thpl, nbrloop_thnp, nbrloop_empl, nbrloop_emnp)
df <- subset(df, Loop_number != 0)

pdf(paste0(pdfname,"_wisker_DMD",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(df$Loop_number ~ df$Type,
        at = c(1,1.9,3,3.9), ylab = "Number of DMD per gene", xlab = "",
        names = c("Plastic\nprediction", "Non plastic\nprediction", "Plastic\nempiric", "Non plastic\nempiric"),
        tck=-0.1,
        las = 1,
        col = c("orange","red"),
        border = "brown"
)
dev.off()

################################################################################
####FBL prop####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FBL.csv", sep = ",")[c(1,5),]

df <- as.data.frame(rbind(theory[1,2:5], theory[2,2:5],
                           colSums(all_plast[,3:6])*100/nrow(all_plast), colSums(non_plast[,3:6])*100/nrow(non_plast) ))
rownames(df) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


pdf(paste0(pdfname,"_prop_FBL",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
#Each motif topology
barplot(t(df[,1:2]), main="FBL motif enrichment", col=c("olivedrab1", "grey"),
        space=c(0.3,0.1,0.4,0.1), legend.text = c("At least 1 FBL", "No FBL"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)))
dev.off()

#Whisker plot###########################
nbrloop_thpl <- read.csv("scripts/data/10g_Pl_nbrFBL.csv", sep = ",")
nbrloop_thpl$Type <- "1thpl"
nbrloop_thnp <- read.csv("scripts/data/10g_NP_nbrFBL.csv", sep = ",")
nbrloop_thnp$Type <- "2thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "3empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "4emnp"

df <- rbind(nbrloop_thpl, nbrloop_thnp, nbrloop_empl, nbrloop_emnp)
df <- subset(df, FBL_number != 0)

pdf(paste0(pdfname,"_wisker_FBL",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(df$FBL_number ~ df$Type, ylim=c(0,1000),
        at = c(1,1.9,3,3.9), ylab = "Number of FBL per gene", xlab = "",
        names = c("Plastic\nprediction", "Non plastic\nprediction", "Plastic\nempiric", "Non plastic\nempiric"),
        tck=-0.1,
        las = 1,
        col = c("orange","red"),
        border = "brown"
)
dev.off()

# #Stat test to know if FBL are different between non plastic and plastic genes
# df <- rbind(non_plast, all_plast)
# df$Plasticity <- c(rep("No", nrow(non_plast)),rep("Yes", nrow(all_plast)))
# anova(lm(FBL ~ Plasticity, df)) #Significant difference for FBL
# 
# anova(lm(Activating ~ Plasticity, subset(df, FBL==1) ))  #Almost significative ? Non significative enough ? for sign of FBL
# 



