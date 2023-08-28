source("scripts/functions/functions.R")
pdfname <- "figures/fig_full_VS"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
#Function to rbind tables, modified from https://www.pogol.net/rbind-fill-for-1-dimensional-tables-in-r
rbind1dtable=function(tab1,tab2,tab3,tab4,fail=0){
  sapply(unique(c(names(tab1),names(tab2),names(tab3),names(tab4))),function(n){
    t1=tab1[which(names(tab1)==n)]
    t2=tab2[which(names(tab2)==n)]
    t3=tab3[which(names(tab3)==n)]
    t4=tab4[which(names(tab4)==n)]
    c(ifelse(length(t1)==0,fail,t1),ifelse(length(t2)==0,fail,t2),ifelse(length(t3)==0,fail,t3),ifelse(length(t4)==0,fail,t4))
  })
}
####FFL prop####################################################################
non_plast_ffl <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast_ffl <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/full_20k_0-01_FFL.csv", sep = ",")[,1:11]


df_ffl <- as.data.frame(rbind( colSums(non_plast_ffl[,3:12])*100/nrow(non_plast_ffl), colSums(all_plast_ffl[,3:12])*100/nrow(all_plast_ffl),
                               theory[2,2:11], theory[1,2:11]))
rownames(df_ffl) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")

pdf(paste0(pdfname,"_prop_FFL",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df_ffl[,1:2]), main = "FFL motif enrichment", col=c("khaki1", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "At least 1 FFL","No FFL"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
text(1.05, 103, paste0( ifelse(t.test(non_plast_ffl[,3], all_plast_ffl[,3], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off() #gold

#Whisker plot###############################
nbrloop_thpl <- read.csv("scripts/data/full_20k_0-01_Pl_nbrFFL.csv", sep = ",")
nbrloop_thpl$Type <- "4thpl"
nbrloop_thnp <- read.csv("scripts/data/full_20k_0-01_NP_nbrFFL.csv", sep = ",")
nbrloop_thnp$Type <- "3thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "2empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "1emnp"

df <- rbind(nbrloop_emnp, nbrloop_empl, nbrloop_thnp, nbrloop_thpl )
df <- subset(df, Loop_number != 0)
table(df[,2])
rbind.fill(table(nbrloop_thpl[,2]), table(nbrloop_thnp[,2]))

pdf(paste0(pdfname,"_wisker_FFL",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(log(df$Loop_number) ~ df$Type,
        at = c(1,1.9,3,3.9), ylab = "Log of FFL number per gene", xlab = "",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        tck=-0.1,
        las = 1,
        col = c("khaki1"),
        border = "black"
)
text(1.4, 4, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

################################################################################
####DMD prop####################################################################
non_plast_dmd <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast_dmd <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/full_20k_0-01_DMD.csv", sep = ",")[,1:13]

df1 <- as.data.frame(rbind(colSums(non_plast_dmd[,3:14])*100/nrow(non_plast_dmd), colSums(all_plast_dmd[,3:14])*100/nrow(all_plast_dmd),
                           theory[2,2:13], theory[1,2:13] ))
rownames(df1) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")


pdf(paste0(pdfname,"_prop_DMD",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df1[,1:2]), main = "Diamond motif enrichment", col=c("darkseagreen2", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "At least 1 DMD","No DMD"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
text(1.05, 103, paste0( ifelse(t.test(non_plast_dmd[,3], all_plast_dmd[,3], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()


#Whisker plot#############################
nbrloop_thpl <- read.csv("scripts/data/full_20k_0-01_Pl_nbrDMD.csv", sep = ",")
nbrloop_thpl$Type <- "4thpl"
nbrloop_thnp <- read.csv("scripts/data/full_20k_0-01_NP_nbrDMD.csv", sep = ",")
nbrloop_thnp$Type <- "3thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "2empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "1emnp"

df <- rbind(nbrloop_emnp, nbrloop_empl, nbrloop_thnp, nbrloop_thpl )
df <- subset(df, Loop_number != 0)

pdf(paste0(pdfname,"_wisker_DMD",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(log(df$Loop_number) ~ df$Type,
        at = c(1,1.9,3,3.9), ylab = "Log of DMD number per gene", xlab = "",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        tck=-0.1,
        las = 1,
        col = c("darkseagreen2"),
        border = "black"
)
text(1.4, 5.2, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

################################################################################
####FBL prop####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")
theory <- read.csv("scripts/data/full_20k_0-01_FBL.csv", sep = ",")[,]

df <- as.data.frame(rbind(colSums(non_plast[,3:6])*100/nrow(non_plast), colSums(all_plast[,3:6])*100/nrow(all_plast),
                          theory[2,2:5], theory[1,2:5] ))
rownames(df) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")


pdf(paste0(pdfname,"_prop_FBL",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
#Each motif topology
barplot(t(df[,1:2]), main="FBL motif enrichment", col=c("darkolivegreen2", "grey"),
        space=c(0.3,0.1,0.4,0.1), legend.text = c("At least 1 FBL", "No FBL"), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
text(1.35, 103, paste0( ifelse(t.test(non_plast[,3], all_plast[,3], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

#Whisker plot###########################
nbrloop_thpl <- read.csv("scripts/data/full_20k_0-01_Pl_nbrFBL.csv", sep = ",")
nbrloop_thpl$Type <- "4thpl"
nbrloop_thnp <- read.csv("scripts/data/full_20k_0-01_NP_nbrFBL.csv", sep = ",")
nbrloop_thnp$Type <- "3thnp"
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_empl$Type <- "2empl"
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp$Type <- "1emnp"

df <- rbind(nbrloop_emnp, nbrloop_empl, nbrloop_thnp, nbrloop_thpl )
df <- subset(df, FBL_number != 0)

pdf(paste0(pdfname,"_wisker_FBL",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(log(df$FBL_number) ~ df$Type,
        at = c(1,1.9,3,3.9), ylab = "Log of FBL per gene", xlab = "",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        tck=-0.1,
        las = 1,
        col = c("darkolivegreen2"),
        border = "black"
)
text(1.4, 7, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

# #Stat test to know if FBL are different between non plastic and plastic genes
# df <- rbind(non_plast, all_plast)
# df$Plasticity <- c(rep("No", nrow(non_plast)),rep("Yes", nrow(all_plast)))
# anova(lm(FBL ~ Plasticity, df)) #Significant difference for FBL
# 
# anova(lm(Activating ~ Plasticity, subset(df, FBL==1) ))  #Almost significative ? Non significative enough ? for sign of FBL
# 


################################################################################
#Number of gene with no loop####################################################
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_size10.csv", sep = ",")
length(unique(rbind(subset(all_plast_ffl, FFL==1)[,1:4],subset(all_plast_dmd, FFL==1)[,1:4],subset(all_plast, FFL==1)[,1:4])[,2]))*100/nrow(all_plast_ffl)
# unique(rbind(subset(all_plast_ffl, FFL==0)[,1:4],subset(all_plast_dmd, FFL==0)[,1:4])[,2])
#

100-length(unique(subset(all_plast, FFL==1)[,2]))*100/nrow(all_plast)


