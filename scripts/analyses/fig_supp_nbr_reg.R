#Fig supp total number of regulations
#Take E_coli matrix, subset plastic genes, get row sum and compare to non plastic gene row sum
source("scripts/functions/functions.R")
library(purrr)
#########################################
pdfname <- "figures/fig_supp_reg_nbr"
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv") #List of regulations from Ecocyc
ec_genes <- read.csv("e_coli/ncbi_dataset_K-12_annotation.csv", sep ="\t") #List of E coli genes from Ecocyc

#Transcriptions factors
TF_genes <- unique(ec_cyc[,1]) # Here, TF = all regulating genes from Ecocyc reg data
#Adding non-regulated genes in regulation data
nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,2] %in% ec_cyc[,1]),]) #Regulators
nonreg_genes <- nonreg_genes[!(nonreg_genes[,2] %in% ec_cyc[,2]),] #Regulatees
#1639 genes that are non-regulated
nr <- as.data.frame(cbind(nonreg_genes[,2],nonreg_genes[,2], rep(0, nrow(nonreg_genes)))) #Same file format as ec_cyc
setnames(nr, 1:3, colnames(ec_cyc))
ec_reg_full <- rbind(ec_cyc, nr) #Reg file with every genes
##########
g <- graph.data.frame(ec_reg_full, directed=TRUE)
E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3'))) #t() to have regulators in columns and regulees in rows

E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
E_coli_mat[is.na(E_coli_mat)] <- 0 #fill NA to 0, mandatory for later analyses
################################################################################
#E. coli
plastgenes <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,2]
npgenes <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,2]
E_coli_mat_plast <- E_coli_mat[plastgenes, ]
reg_number_plast <- rowSums2(abs(E_coli_mat_plast))
E_coli_mat_np <- E_coli_mat[npgenes, ]
reg_number_np <- rowSums2(abs(E_coli_mat_np))

#Simulations
topo_plastic <- readRDS("scripts/data/list_plastic_topo.Rds")
topo_nnplast <- readRDS("scripts/data/list_nnplast_topo.Rds")
  
THreg_sum_plast <- unlist(lapply(topo_plastic, function(i){rowSums2(abs(i))} ))
THreg_sum_np <- unlist(lapply(topo_nnplast, function(i){rowSums2(abs(i))} ))

#
pdf(paste0(pdfname,".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(cbind(reg_number_np, reg_number_plast, THreg_sum_np, THreg_sum_plast), at = c(1,1.9,3,3.9), ylab = "Number of regulation per gene",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"))
text(1.35, 25, paste0( ifelse(t.test(reg_number_plast, reg_number_np, var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

################################################################################
#Positive and neg regulations
#E. coli
plastmat <- E_coli_mat_plast
plastmat[(plastmat == -1)] <- 0
reg_number_pos_plast <- rowSums2(plastmat)
plastmat <- E_coli_mat_plast
plastmat[(plastmat == 1)] <- 0
reg_number_neg_plast <- rowSums2(abs(plastmat))
#
plastmat <- E_coli_mat_np
plastmat[(plastmat == -1)] <- 0
reg_number_pos_np <- rowSums2(plastmat)
plastmat <- E_coli_mat_np
plastmat[(plastmat == 1)] <- 0
reg_number_neg_np <- rowSums2(abs(plastmat))


#Simulations
THreg_sum_pos_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == -1)] <- 0
return(rowSums2(abs(i)))} ))
THreg_sum_neg_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i)))} ))
#
THreg_sum_pos_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 0
return(rowSums2(abs(i)))} ))
THreg_sum_neg_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i)))} ))

df <- data.frame(Pos=c(sum(reg_number_pos_np)*100/sum(reg_number_np),  sum(reg_number_pos_plast)*100/sum(reg_number_plast),  sum(THreg_sum_pos_np)*100/sum(THreg_sum_np), sum(THreg_sum_pos_plast)*100/sum(THreg_sum_plast) ),
                 Neg=c(sum(reg_number_neg_np)*100/sum(reg_number_np), sum(reg_number_neg_plast)*100/sum(reg_number_plast), sum(THreg_sum_neg_np)*100/sum(THreg_sum_np), sum(THreg_sum_neg_plast)*100/sum(THreg_sum_plast) ))

pdf(paste0(pdfname,"_sign",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
barplot(t(df), main = "Regulation type", names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        space=c(0.0,0.1,0.35,0.1), legend.text = c( "Activating","Inhibiting"), args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
text(1.05, 103, paste0( ifelse(t.test(reg_number_neg_plast, reg_number_neg_np, var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()


