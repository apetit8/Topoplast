source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd_2024_01_29.csv") #List of regulations from Ecocyc
ec_genes <- read.csv("e_coli/ncbi_dataset_K-12_annotation.csv", sep ="\t") #List of E coli genes from Ecocyc

#Keep only regulations from genes in the annotation
ec_cyc <- subset(ec_cyc, V1 %in% ec_genes[,2] & V2 %in% ec_genes[,2])

#Transcriptions factors
TF_genes <- unique(ec_cyc[,1]) # Here, TF = all regulating genes from Ecocyc reg data
#Adding non-regulated genes in regulation data
nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,2] %in% ec_cyc[,1]),]) #Genes that are not Regulators
nonreg_genes <- nonreg_genes[!(nonreg_genes[,2] %in% ec_cyc[,2]),] #Genes that are not Regulatees
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
#FROM = FALSE
#########################################
#Analyses of plastic genes from different sources
fun <- "FFL"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_FFL"
##
source("scripts/analyses/E_coli.R")
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
source("scripts/analyses/E_coli.R")
#########################################
#FROM = FALSE ; count FFL
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_nffl"
##
source("scripts/analyses/E_coli.R")

#########################################
#FROM = FALSE ; count DMD
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_nDMD"
##
source("scripts/analyses/E_coli.R")

#########################################
#FBL
#########################################
#Analyses of plastic genes from different sources
fun <- "FBL"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_FBL")
##
source("scripts/analyses/E_coli.R")

#########################################
#FBL count
#########################################
#Analyses of plastic genes from different sources
fun <- "FBLcount"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_nFBL")
##
source("scripts/analyses/E_coli.R")

#########################################

################################################################################
phgenes1 <-e_coli_gene_name(read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)[,1])
phgenes2 <- e_coli_gene_name(read.table("e_coli/Plast_genes/pH_2/plastic_genes.txt", sep ="\t", header=FALSE)[,1])
stress_genes <- e_coli_gene_name(as.data.frame(unique(c(subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S5.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S6.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S7.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S8.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S9.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5])))[,1])
medgrowth_genes <- e_coli_gene_name(read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S1.csv", sep ="\t", header=TRUE)[,1])
mg_c_genes <- e_coli_gene_name(as.data.frame(unique(subset(read.csv("e_coli/Plast_genes/C_Mg_Na/41598_2017_BFsrep45303_MOESM60_ESM.csv"), dataType=="mrna")[,2]))[,1])
ox_genes <- e_coli_gene_name(read.csv("e_coli/Plast_genes/Oxydative_stress/table_5_genes_updated_names.txt", sep ="\t", header=FALSE)[,1])
aero_genes <- e_coli_gene_name(read.table("e_coli/Plast_genes/Aero_liquid/Table_2_3_4.txt", sep ="\t", header=FALSE)[,1])
juice_genes <- e_coli_gene_name(read.table("e_coli/Plast_genes/apple_juice/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)[,1])
temptr_genes <- e_coli_gene_name(read.table("e_coli/Plast_genes/human_temp/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)[,1])
temptr_genes2 <- e_coli_gene_name(read.table("e_coli/Plast_genes/Temperature/plastic_genes.txt", sep ="\t", header=FALSE)[,1])
stringent_genes <- e_coli_gene_name(as.data.frame(unique(read.table("e_coli/Plast_genes/Stringent_response/All_genes.txt", sep ="\t", header=FALSE)[,1]))[,1])

length(as.character(as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                                  ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2)))[,1]))
length(as.character(subset( as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                                                   ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2))), Freq >=2)[,1]))

all_plast_genes <- c(phgenes1, phgenes2, stress_genes, medgrowth_genes , mg_c_genes, ox_genes, aero_genes, juice_genes, temptr_genes, temptr_genes2, stringent_genes)
#Number of plastic genes in annotation
all_plast_genes1 <- e_coli_gene_name(all_plast_genes)
length(unique(all_plast_genes1))

###
Defplastic_genes <- as.character(subset( as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                                               ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2))), Freq >=2)[,1])
length(Defplastic_genes)
#sum(E_coli_mat[(rownames(E_coli_mat) %in% unique(all_plast_genes1)) , ]==1)/(sum(E_coli_mat[(rownames(E_coli_mat) %in% unique(all_plast_genes1)) , ]==1)+sum(E_coli_mat[(rownames(E_coli_mat) %in% unique(all_plast_genes1)) , ]==-1))

################################################################################
################################################################################
#CONTROL shuffled matrix
################################################################################
##########
ec_reg_full <- rbind(ec_cyc, nr) #TO SHUFFLE

for (i in 1:nrow(ec_reg_full)) {
  if(ec_reg_full[i,3] !=0 ) ec_reg_full[i,3] <- sample(subset(ec_reg_full, V3 != 0)[,3], 1)
}

g <- graph.data.frame(ec_reg_full, directed=TRUE)

##############
#########################################
#FROM = FALSE
#########################################
#Analyses of plastic genes from different sources
fun <- "FFL"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random2_FFL"
##
source("scripts/analyses/E_coli.R")
#########################################
#FROM = FALSE ; Diamond motifs
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random2_diamond"
##
source("scripts/analyses/E_coli.R")
#########################################
#FROM = FALSE ; count FFL
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random2_nffl"
##
source("scripts/analyses/E_coli.R")

#########################################
#FROM = FALSE ; count DMD
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random2_nDMD"
##
source("scripts/analyses/E_coli.R")

#########################################
#FBL
#########################################
#Analyses of plastic genes from different sources
fun <- "FBL"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_random2_FBL")
##
source("scripts/analyses/E_coli.R")

#########################################
#FBL count
#########################################
#Analyses of plastic genes from different sources
fun <- "FBLcount"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_random2_nFBL")
##
source("scripts/analyses/E_coli.R")

################################################################################
################################################################################
################################################################################
#CONTROL random TF
################################################################################
##########
ec_reg_full <- rbind(ec_cyc, nr) #TO SHUFFLE
g <- graph.data.frame(ec_reg_full, directed=TRUE)

##############
#########################################
#FROM = FALSE
#########################################
#Analyses of plastic genes from different sources
fun <- "FFL"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random3_FFL"
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))
#########################################
#FROM = FALSE ; Diamond motifs
#########################################
#Analyses of plastic genes from different sources
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random3_diamond"
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))
#########################################
#FROM = FALSE ; count FFL
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 1
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random3_nffl"
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))
#########################################
#FROM = FALSE ; count DMD
#########################################
#Analyses of plastic genes from different sources
fun <- "FFLcount"
edges1 <- 2
edges2 <- 2
from <- FALSE
all_plast_genes <- data.frame()
csvname <- "E_coli_random3_nDMD"
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))
#########################################
#FBL
#########################################
#Analyses of plastic genes from different sources
fun <- "FBL"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_random3_FBL")
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))
#########################################
#FBL count
#########################################
#Analyses of plastic genes from different sources
fun <- "FBLcount"
edges1 <- c(2:5) #More FBL with edges1=3 than=2 ; why ?
edges2 <- 0
from <- FALSE
all_plast_genes <- data.frame()
csvname <- paste0("E_coli_random3_nFBL")
##
ran.genes <- e_coli_prep_analyses_randomized(Defplastic_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, from = from)
write.csv(ran.genes, paste0("scripts/data/plast_genes_",csvname,".csv"))








