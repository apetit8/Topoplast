source("scripts/functions/functions.R")
pdfname <- "figures/fig_VS"
#Plot combining both empirical data and theoretical data
################################################################################
####FFL prop
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]


df1 <- as.data.frame(rbind(theory[1,2:11], theory[2,2:11],
                           colSums(all_plast[,3:12])*100/nrow(all_plast), colSums(non_plast[,3:12])*100/nrow(non_plast) ))
rownames(df1) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")

pdf(paste0(pdfname,"_prop_FFL",".pdf"), width=5, height=5)
par(mar = c(7,2, 2,1))
barplot(t(df1[,1:2]), main = "FFL motif enrichment", col=c("gold", "grey"), space=c(0.0,0.1,0.35,0.1), legend.text = c( "at least 1 FFL","No FFL"), 
        args.legend = list(x = "topright", inset = c(0.3, 1.2)))
dev.off()

##FFL motifs
df1 <- as.data.frame(rbind(theory[1,4:11]*100/theory[1,2],theory[2,4:11]*100/theory[2,2], colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                           colSums(non_plast[,5:12])*100/sum(non_plast[,3]) ))
rownames(df1) <- c("Plastic\nPrediction", "Non-plastic\nPrediction","Plastic\nEmpiric","Non-Plastic\nEmpiric")


pdf(paste0(pdfname,"_FFL_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df1), main="FFL motif distribution", col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
        space=c(0.3,0.1,0.4,0.1), legend.text = colnames(df1), args.legend = list(ncol=4, x = "topright", inset = c(0.2, 1.2)))
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


pdf(paste0(pdfname,"_DMD_motifs",".pdf"), width=5, height=6)
par(mar = c(8,2, 2,1))
#Each motif topology
barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1","plum1","lightsalmon1","indianred1","lightpink","chocolate1"),
        space=c(0.3,0.1,0.4,0.1), legend.text = colnames(df1[1:10]), args.legend = list(ncol=3, x = "topright", inset = c(0, 1.13)))
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
barplot(t(df[2:4]), main="Diamond motif distribution", col=c("grey","olivedrab1","palegreen"),
        space=c(0.3,0.1,0.4,0.1), legend.text = colnames(df[2:4]), args.legend = list(ncol=3, x = "topright", inset = c(0, 1.13)))
dev.off()


