#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")
pdfname <- "figures/fig_Ecoli"

##FFL
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")

df1 <- as.data.frame(rbind(colSums(non_plast[,3:12])*100/nrow(non_plast),
                           colSums(all_plast[,3:12])*100/nrow(all_plast)
))
rownames(df1) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df1$coherent <- rowSums2(as.matrix(df1[,c(3:6)])) #Count of 
df1$incoherent <- rowSums2(as.matrix(df1[,c(7:10)]))
df1$homogenous <- rowSums2(as.matrix(df1[,c(3,4,7,8)]))
df1$heterogenous <- rowSums2(as.matrix(df1[,c(5,6,9,10)]))


pdf(paste0(pdfname,"_prop_FFL",".pdf"), width=5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df1[,1:2]), col=c("gold", "grey"), main="E coli")
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("gold", "grey"),
       legend=c( "at least 1 FFL","No FFL") )
dev.off()


################################################################################
##FFL with the X gene being plastic as well
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL_from_allplast.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_from_allplast.csv", sep = ",")

#FFL FROM percent
df <- as.data.frame(rbind(colSums(non_plast[,3:13])*100/(sum(non_plast[,3])-sum(non_plast[,13])),
                          colSums(all_plast[,3:13])*100/(sum(all_plast[,3])-sum(all_plast[,13]))
))
rownames(df) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))

pdf(paste0(pdfname,"_percent_FFL_from",".pdf"), width=5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c( "C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()


##
df <- as.data.frame(rbind(colSums(non_plast[,3:13]),
                          colSums(all_plast[,3:13])
))
rownames(df) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))

pdf(paste0(pdfname,"_abs_FFL_from",".pdf"), width=5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c( "C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()

#FFL FROM NP percent
################################################################################
##FFL with the X gene being plastic as well
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL_from_NP.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_from_NP.csv", sep = ",")

df <- as.data.frame(rbind(colSums(non_plast[,3:13])*100/(sum(non_plast[,3])-sum(non_plast[,13])),
                          colSums(all_plast[,3:13])*100/(sum(all_plast[,3])-sum(all_plast[,13]))
))
rownames(df1) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))

pdf(paste0(pdfname,"_percent_FFL_from_NP",".pdf"), width=5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c("C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()

#
df <- as.data.frame(rbind(colSums(non_plast[,3:13]),
                          colSums(all_plast[,3:13])
))
rownames(df1) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))

pdf(paste0(pdfname,"_abs_FFL_from_NP",".pdf"), width=5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c("C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()

#########################
##Diamond
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")

df <- as.data.frame(rbind(colSums(colSums(non_plast[,3:14])*100/nrow(non_plast),
                          colSums(all_plast[,3:14])*100/nrow(all_plast))
))
rownames(df1) <- c("Non-plastic\ngenes", "Plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3,5,6,8,10,11)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(4,7,9,12)]))
df$homogenous_pos <- rowSums2(as.matrix(df[,c(3,6,9)]))
df$homogenous_neg <- rowSums2(as.matrix(df[,c(5,8,12)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(4,7,10,11)]))


pdf(paste0(pdfname,"_diamond_from",".pdf"), width=6, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,13,14)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,15:17)]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()

