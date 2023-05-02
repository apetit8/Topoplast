source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv") #List of regulations from Ecocyc
ec_genes <- read.csv("e_coli/All-genes-of-E.-coli-K-12-substr.-MG1655---GO.csv", sep ="\t") #List of E coli genes from Ecocyc

#Adding non-regulated genes in regulation data
nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,1] %in% ec_cyc[,1]),])
nonreg_genes <- nonreg_genes[!(nonreg_genes[,1] %in% ec_cyc[,2]),] #1704 genes that are non-regulated
nr <- as.data.frame(cbind(nonreg_genes[,1],nonreg_genes[,1], rep(0, length(nonreg_genes)))) #Same file format as ec_cyc
setnames(nr, 1:3, colnames(ec_cyc))
ec_reg_full <- rbind(ec_cyc, nr) #Reg file with every genes
##########
g <- graph.data.frame(ec_reg_full, directed=TRUE)
E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3'))) #t() to have regulators in columns and regulees in rows

E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
E_coli_mat[is.na(E_coli_mat)] <- 0 #fill NA to 0, mandatory for later analyses
#########################################
#
#
#########################################
#Analyses of plastic genes from different sources
cutoff.max <- 5
cutoff.min <- 1
#########################################
#
##PH
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002

phgenes <- read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)
phgenes <- phgenes[(phgenes[,1] %in% colnames(E_coli_mat)),] 

phloops <- mclapply(phgenes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 4)
write.csv(rbindlist(phloops), "scripts/data/ph_plast_h_ffloops.csv")

########################################
##Medium of growth
#Feugeas et al., 2016 ; DOI : 10.1093/molbev/msw105

medgrowth_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S1.csv", sep ="\t", header=TRUE)
medgrowth_genes <- medgrowth_genes[(medgrowth_genes[,1] %in% colnames(E_coli_mat)),1] 

non_envir_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S3.csv", sep ="\t", header=TRUE)
non_envir_genes <- non_envir_genes[(non_envir_genes[,1] %in% colnames(E_coli_mat)),1] 

plast_medgrowthloops <- mclapply(medgrowth_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_medgrowthloops), "scripts/data/plast_medium_growth_h_ffloops.csv")

#Non plastic genes (expression depending on strains and not on tested environment)
np_medgrowthloops <- mclapply(non_envir_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(np_medgrowthloops), "scripts/data/np_medium_growth_h_ffloops.csv")


########################################
#10.1038/srep45303
#List of plastic genes for 3 different conditions: carbone source, Mg stress, Na+ stress
mg_c_genes <- as.data.frame(unique(subset(read.csv("e_coli/Plast_genes/C_Mg_Na/41598_2017_BFsrep45303_MOESM60_ESM.csv"), dataType=="mrna")[,2]))
mg_c_genes <- mg_c_genes[(mg_c_genes[,1] %in% colnames(E_coli_mat)),] 

plast_mg_c <- mclapply(mg_c_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_mg_c), "scripts/data/plast_mg_c_h_ffloops.csv")

########################################
#10.1128/AEM.00914-09
#List of plastic genes oxydative stress
ox_genes <- read.csv("e_coli/Plast_genes/Oxydative_stress/table_5_genes_updated_names.txt", sep ="\t", header=FALSE)
ox_genes <- ox_genes[(ox_genes[,1] %in% colnames(E_coli_mat)),] 

plast_ox <- mclapply(ox_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_ox), "scripts/data/plast_mg_c_h_ffloops.csv")

########################################
#10.1007/s00253-018-9083-5
#List of DEGS between aerozolisation and liquid suspension
aero_genes <- read.table("e_coli/Plast_genes/Aero_liquid/Table_2_3_4.txt", sep ="\t", header=FALSE)
aero_genes <- aero_genes[(aero_genes[,1] %in% colnames(E_coli_mat)),] 

plast_aero <- mclapply(aero_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_aero), "scripts/data/plast_aero_h_ffloops.csv")

########################################
#10.1128/AEM.02841-08
#Apple Juice
juice_genes <- read.table("e_coli/Plast_genes/apple_juice/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
juice_genes <- juice_genes[(juice_genes[,1] %in% colnames(E_coli_mat)),] 

plast_juice <- mclapply(juice_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_juice), "scripts/data/plast_juice_h_ffloops.csv")

########################################
#10.1128/JB.01929-06
#Human body temperature
temptr_genes <- read.table("e_coli/Plast_genes/human_temp/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
temptr_genes <- temptr_genes[(temptr_genes[,1] %in% colnames(E_coli_mat)),] 

plast_temptr <- mclapply(temptr_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_temptr), "scripts/data/plast_temptr_h_ffloops.csv")

########################################
#10.1038/s41598-022-12463-3
#Compendium of stressor
stringent_genes <- as.data.frame(unique(read.table("e_coli/Plast_genes/Stringent_response/All_genes.txt", sep ="\t", header=FALSE)[,1]))
stringent_genes <- stringent_genes[(stringent_genes[,1] %in% colnames(E_coli_mat)),] 

plast_stringent <- mclapply(stringent_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_stringent), "scripts/data/plast_stringent_h_ffloops.csv")

########################################
#10.1128/JB.01092-07
#stringent response
stress_genes <- as.data.frame(unique(c(subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S5.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S6.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S7.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S8.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S9.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5])))
stress_genes <- stress_genes[(stress_genes[,1] %in% colnames(E_coli_mat)),] 

plast_stress <- mclapply(stress_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(plast_stress), "scripts/data/plast_stress_h_ffloops.csv")


########################################
#Every plastic genes
plast_genes <- unique(c(mg_c_genes,medgrowth_genes,phgenes$V1,ox_genes,aero_genes,juice_genes,temptr_genes))

plast_ffloops <- mclapply(plast_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 100)
write.csv(rbindlist(nonplast_ffloops), "scripts/data/plast_genes_h_ffloops.csv")

########################################
#Genes that were not find as responding to the environment in our data corpus
non_plast_genes <- as.data.frame(colnames(E_coli_mat))
nonplast_genes <- non_plast_genes[!(non_plast_genes[,1] %in% c(mg_c_genes,medgrowth_genes,phgenes$V1,ox_genes,aero_genes,juice_genes,temptr_genes)),]

nonplast_h_ffloops <- mclapply(nonplast_genes, function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 50)
write.csv(rbindlist(nonplast_h_ffloops), "scripts/data/nonplast_h_ffloops.csv")


########################################
#Control : every genes

all_ffloops <- mclapply(colnames(E_coli_mat), function(gene) {
  cc <- homog.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
}, mc.cores = 100)
write.csv(rbindlist(nonplast_ffloops), "scripts/data/all_genes_h_ffloops.csv")








