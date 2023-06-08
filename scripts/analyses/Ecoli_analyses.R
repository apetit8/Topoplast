source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
memory.limit(size=2500)
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
##########
igraph_options(return.vs.es=F)
#########################################
#FROM = FALSE ; meaning that the output will be the the FFLS coming from and to plastic genes.
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_FFL"
##
source("scripts/analyses/E_coli_FFL.R")
#########################################
#FROM = FALSE ; Diamond motifs
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_diamond"
##
source("scripts/analyses/E_coli_FFL.R")
#########################################
#FROM = TRUE ; meaning that the output will be the the FFLS coming from and to plastic genes.
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 1
from <- TRUE
all_plast_genes <- data.frame()
csvname <- "E_coli_FFL_from"
##
source("scripts/analyses/E_coli_FFL.R")
from <- all_plast_genes$V1
#Genes that were not find as responding to the environment in our data corpus
nonplast_FFL <- e_coli_prep_analyses(nonplast_genes, g, E_coli_mat, fun="FFL", edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(nonplast_FFL, paste0("scripts/data/nonplast_",csvname,".csv"))
print("Non plastic genes done!")
#########################################
#FROM = Non plastic gene names
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 1
from <- nonplast_genes
all_plast_genes <- data.frame()
csvname <- "E_coli_FFL_from_NP"
##
source("scripts/analyses/E_coli_FFL.R")





