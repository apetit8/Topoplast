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
E_coli_mat_np <- E_coli_mat[npgenes, ]

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
return(rowSums2(abs(i))[2])} )) # [2] : the target of the network is in position #2
THreg_sum_neg_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i))[2])} ))
#
THreg_sum_pos_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 0
return(rowSums2(abs(i))[2])} ))
THreg_sum_neg_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i))[2])} ))

pdf(paste0(pdfname,"_pos",".pdf"), width=5, height=4)
layout(matrix(c(1), 1, 1, byrow = TRUE))
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(cbind(reg_number_pos_np, reg_number_pos_plast, THreg_sum_pos_np, THreg_sum_pos_plast), at = c(1,1.9,3,3.9), ylab = "Number of activation per gene",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        col = c("orange","cyan3","orange","cyan3"))
text(1.35, 12, paste0( ifelse(t.test(reg_number_plast, reg_number_np, var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

pdf(paste0(pdfname,"_neg",".pdf"), width=5, height=4)
layout(matrix(c(1), 1, 1, byrow = TRUE))
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
boxplot(cbind(reg_number_neg_np, reg_number_neg_plast, THreg_sum_neg_np, THreg_sum_neg_plast), at = c(1,1.9,3,3.9), ylab = "Number of inhibition per gene",
        names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"),
        col = c("orange","cyan3","orange","cyan3"))
text(1.35, 15, paste0( ifelse(t.test(reg_number_plast, reg_number_np, var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
dev.off()

################################################################################
# 
reg_number_plast <- rowSums2(abs(E_coli_mat_plast))
E_coli_mat_np <- E_coli_mat[npgenes, ]
reg_number_np <- rowSums2(abs(E_coli_mat_np))




plot(reg_number_plast, read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3])
points(reg_number_np, read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3], col="red")

abline(lm( read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3] ~ reg_number_plast))
abline(lm( read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3] ~ reg_number_np))

# 
# #Simulations
# topo_plastic <- readRDS("scripts/data/list_plastic_topo_full_netw.Rds")
# topo_nnplast <- readRDS("scripts/data/list_nnplast_topo_full_netw.Rds")
#   
# THreg_sum_plast <- unlist(lapply(topo_plastic, function(i){rowSums2(abs(i))[2]} )) # [2] : the target of the network is in position #2
# THreg_sum_np <- unlist(lapply(topo_nnplast, function(i){rowSums2(abs(i))[2]} ))
# #
# length(reg_number_np) <- max(length(reg_number_np), length(reg_number_plast), length(THreg_sum_np), length(THreg_sum_plast))
# length(reg_number_plast) <- max(length(reg_number_np), length(reg_number_plast), length(THreg_sum_np), length(THreg_sum_plast))
# length(THreg_sum_plast) <- max(length(reg_number_np), length(reg_number_plast), length(THreg_sum_np), length(THreg_sum_plast))
# length(THreg_sum_np) <- max(length(reg_number_np), length(reg_number_plast), length(THreg_sum_np), length(THreg_sum_plast))
# #
# pdf(paste0(pdfname,".pdf"), width=5, height=4)
# par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
# boxplot(cbind(reg_number_np, reg_number_plast, as.array(THreg_sum_np), as.array(THreg_sum_plast)), at = c(1,1.9,3,3.9), ylab = "Number of regulation per gene",
#         names = c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory"))
# text(1.35, 25, paste0( ifelse(t.test(reg_number_plast, reg_number_np, var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
# dev.off()





