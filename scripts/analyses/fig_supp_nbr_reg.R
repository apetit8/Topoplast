#Fig supp total number of regulations
#Take E_coli matrix, subset plastic genes, get row sum and compare to non plastic gene row sum
source("scripts/functions/functions.R")
library(purrr)
#########################################
pdfname <- "figures/fig_supp_reg_nbr"
gen <- 10000 
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
#Empiric
plastgenes <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,2]
npgenes <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,2]
E_coli_mat_plast <- E_coli_mat[plastgenes, ]
reg_number_plast <- rowSums2(abs(E_coli_mat_plast))
E_coli_mat_np <- E_coli_mat[npgenes, ]
reg_number_np <- rowSums2(abs(E_coli_mat_np))

#Prediction
sims.dirs1 <- list.dirs("simul/10g_selfreg", recursive = TRUE) #c(list.dirs("simul/10g", recursive = TRUE), list.dirs("simul/10g_a0.15", recursive = TRUE), list.dirs("simul/10g_rep", recursive = TRUE))
df.10 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10$envir <- str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Anticorrelated"), target=target,
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))
topo.corr10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Correlated"), target=target,
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))
topo.sel10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Control_sel"), target=target,
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

list_mat <- c(lapply(topo.anticor10, function(i){abs(i)} ),lapply(topo.corr10, function(i){abs(i)}))
list_mat <- list_mat %>% map(~.x[2, ])
THreg_sum_plast <- sapply(list_mat, sum)

list_mat <- lapply(topo.sel10, function(i){abs(i)} )
list_mat <- list_mat %>% map(~.x[2, ])
THreg_sum_np <- sapply(list_mat, sum)

#####
pdf(paste0(pdfname,"_reg_nbr",".pdf"), width=5, height=4)
par(mar = c(5,2, 2,1))
boxplot(cbind(THreg_sum_plast, THreg_sum_np, reg_number_plast, reg_number_np), main = "Total number of regulation", at = c(1,1.9,3,3.9), names = c("Plastic\nprediction", "Non plastic\nprediction", "Plastic\nempiric", "Non plastic\nempiric"))
dev.off()


