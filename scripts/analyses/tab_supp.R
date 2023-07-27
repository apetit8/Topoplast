source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
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
##PH ; Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002
phgenes1 <- read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)
phgenes1 <- phgenes1[(phgenes1[,1] %in% colnames(E_coli_mat)),] 
########################################
#PH ; Maurer et al., 2004 ; 
phgenes2 <- read.table("e_coli/Plast_genes/pH_2/plastic_genes.txt", sep ="\t", header=FALSE)
phgenes2 <- phgenes2[(phgenes2[,1] %in% colnames(E_coli_mat)),] 
########################################
#10.1128/JB.01092-07 Bhatia et al., 2022
#Compendium of stressor
stress_genes <- as.data.frame(unique(c(subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S5.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S6.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S7.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S8.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S9.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5])))
stress_genes <- stress_genes[(stress_genes[,1] %in% colnames(E_coli_mat)),] 
########################################
##Medium of growth
#Feugeas et al., 2016 ; DOI : 10.1093/molbev/msw105
medgrowth_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S1.csv", sep ="\t", header=TRUE)
medgrowth_genes <- medgrowth_genes[(medgrowth_genes[,1] %in% colnames(E_coli_mat)),1] 
########################################
#10.1038/srep45303 Caglar et al., 2017
#List of plastic genes for 3 different conditions: carbone source, Mg stress, Na+ stress
mg_c_genes <- as.data.frame(unique(subset(read.csv("e_coli/Plast_genes/C_Mg_Na/41598_2017_BFsrep45303_MOESM60_ESM.csv"), dataType=="mrna")[,2]))
mg_c_genes <- mg_c_genes[(mg_c_genes[,1] %in% colnames(E_coli_mat)),1]
########################################
#10.1128/AEM.00914-09 Wang et al., 2009
#List of plastic genes oxydative stress
ox_genes <- read.csv("e_coli/Plast_genes/Oxydative_stress/table_5_genes_updated_names.txt", sep ="\t", header=FALSE)
ox_genes <- ox_genes[(ox_genes[,1] %in% colnames(E_coli_mat)),] 
########################################
#10.1007/s00253-018-9083-5 Ng et al., 2018
#List of DEGS between aerozolisation and liquid suspension
aero_genes <- read.table("e_coli/Plast_genes/Aero_liquid/Table_2_3_4.txt", sep ="\t", header=FALSE)
aero_genes <- aero_genes[(aero_genes[,1] %in% colnames(E_coli_mat)),] 
########################################
#10.1128/AEM.02841-08 Bergholz et al., 2009
#Apple Juice
juice_genes <- read.table("e_coli/Plast_genes/apple_juice/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
juice_genes <- juice_genes[(juice_genes[,1] %in% colnames(E_coli_mat)),] 
########################################
#10.1128/JB.01929-06 White-Ziegler et al., 2007
#Human body temperature
temptr_genes <- read.table("e_coli/Plast_genes/human_temp/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
temptr_genes <- temptr_genes[(temptr_genes[,1] %in% colnames(E_coli_mat)),] 
########################################
#https://doi.org/10.1038/s41598-020-74606-8 Kim et al., 2020 ; Temperature
temptr_genes2 <- read.table("e_coli/Plast_genes/Temperature/plastic_genes.txt", sep ="\t", header=FALSE)
temptr_genes2 <- temptr_genes2[(temptr_genes2[,1] %in% colnames(E_coli_mat)),]
########################################
#10.1128/JB.01092-07 Durfee et al., 2008
#stringent response
stringent_genes <- as.data.frame(unique(read.table("e_coli/Plast_genes/Stringent_response/All_genes.txt", sep ="\t", header=FALSE)[,1]))
stringent_genes <- stringent_genes[(stringent_genes[,1] %in% colnames(E_coli_mat)),] 
########################################

#PROBLEM: %in% also count it if the string of character in inside a larger string. Not what I want at all ; 

study_mat <- matrix(data=NA, nrow=11, ncol=11)
colnames(study_mat) <- c("Tucker2002","Maurer2004","Bhatia2022","Feugas2016","Caglar2017","Wang2009","Ng2018","Bergholz2009","White-Ziegler2007","Kim2020","Durfee2008")
rownames(study_mat) <- colnames(study_mat)
# c(phgenes1,phgenes2,stress_genes,medgrowth_genes,mg_c_genes,ox_genes,aero_genes,juice_genes,temptr_genes,temptr_genes2,stringent_genes)

study_mat[,1] <- c( length(which(phgenes1 %in% phgenes1)), length(which(phgenes1 %in% phgenes2)),
                    length(which(phgenes1 %in% stress_genes)), length(which(phgenes1 %in% medgrowth_genes)),
                    length(which(phgenes1 %in% mg_c_genes)), length(which(phgenes1 %in% ox_genes)),
                    length(which(phgenes1 %in% aero_genes)), length(which(phgenes1 %in% juice_genes)),
                    length(which(phgenes1 %in% temptr_genes)), length(which(phgenes1 %in% temptr_genes2)),
                    length(which(phgenes1 %in% stringent_genes)))

study_mat[,2] <- c( length(which(phgenes2 %in% phgenes1)), length(which(phgenes2 %in% phgenes2)),
                    length(which(phgenes2 %in% stress_genes)), length(which(phgenes2 %in% medgrowth_genes)),
                    length(which(phgenes2 %in% mg_c_genes)), length(which(phgenes2 %in% ox_genes)),
                    length(which(phgenes2 %in% aero_genes)), length(which(phgenes2 %in% juice_genes)),
                    length(which(phgenes2 %in% temptr_genes)), length(which(phgenes2 %in% temptr_genes2)),
                    length(which(phgenes2 %in% stringent_genes)))

study_mat[,3] <- c( length(which(stress_genes %in% phgenes1)), length(which(stress_genes %in% phgenes2)),
                    length(which(stress_genes %in% stress_genes)), length(which(stress_genes %in% medgrowth_genes)),
                    length(which(stress_genes %in% mg_c_genes)), length(which(stress_genes %in% ox_genes)),
                    length(which(stress_genes %in% aero_genes)), length(which(stress_genes %in% juice_genes)),
                    length(which(stress_genes %in% temptr_genes)), length(which(stress_genes %in% temptr_genes2)),
                    length(which(stress_genes %in% stringent_genes)))

study_mat[,4] <- c( length(which(medgrowth_genes %in% phgenes1)), length(which(medgrowth_genes %in% phgenes2)),
                    length(which(medgrowth_genes %in% stress_genes)), length(which(medgrowth_genes %in% medgrowth_genes)),
                    length(which(medgrowth_genes %in% mg_c_genes)), length(which(medgrowth_genes %in% ox_genes)),
                    length(which(medgrowth_genes %in% aero_genes)), length(which(medgrowth_genes %in% juice_genes)),
                    length(which(medgrowth_genes %in% temptr_genes)), length(which(medgrowth_genes %in% temptr_genes2)),
                    length(which(medgrowth_genes %in% stringent_genes)))

study_mat[,5] <- c( length(which(mg_c_genes %in% phgenes1)), length(which(mg_c_genes %in% phgenes2)),
                    length(which(mg_c_genes %in% stress_genes)), length(which(mg_c_genes %in% medgrowth_genes)),
                    length(which(mg_c_genes %in% mg_c_genes)), length(which(mg_c_genes %in% ox_genes)),
                    length(which(mg_c_genes %in% aero_genes)), length(which(mg_c_genes %in% juice_genes)),
                    length(which(mg_c_genes %in% temptr_genes)), length(which(mg_c_genes %in% temptr_genes2)),
                    length(which(mg_c_genes %in% stringent_genes)))


study_mat[,6] <- c( length(which(ox_genes %in% phgenes1)), length(which(ox_genes %in% phgenes2)),
                    length(which(ox_genes %in% stress_genes)), length(which(ox_genes %in% medgrowth_genes)),
                    length(which(ox_genes %in% mg_c_genes)), length(which(ox_genes %in% ox_genes)),
                    length(which(ox_genes %in% aero_genes)), length(which(ox_genes %in% juice_genes)),
                    length(which(ox_genes %in% temptr_genes)), length(which(ox_genes %in% temptr_genes2)),
                    length(which(ox_genes %in% stringent_genes)))

study_mat[,7] <- c( length(which(aero_genes %in% phgenes1)), length(which(aero_genes %in% phgenes2)),
                    length(which(aero_genes %in% stress_genes)), length(which(aero_genes %in% medgrowth_genes)),
                    length(which(aero_genes %in% mg_c_genes)), length(which(aero_genes %in% ox_genes)),
                    length(which(aero_genes %in% aero_genes)), length(which(aero_genes %in% juice_genes)),
                    length(which(aero_genes %in% temptr_genes)), length(which(aero_genes %in% temptr_genes2)),
                    length(which(aero_genes %in% stringent_genes)))

study_mat[,8] <- c( length(which(juice_genes %in% phgenes1)), length(which(juice_genes %in% phgenes2)),
                    length(which(juice_genes %in% stress_genes)), length(which(juice_genes %in% medgrowth_genes)),
                    length(which(juice_genes %in% mg_c_genes)), length(which(juice_genes %in% ox_genes)),
                    length(which(juice_genes %in% aero_genes)), length(which(juice_genes %in% juice_genes)),
                    length(which(juice_genes %in% temptr_genes)), length(which(juice_genes %in% temptr_genes2)),
                    length(which(juice_genes %in% stringent_genes)))

study_mat[,9] <- c( length(which(temptr_genes %in% phgenes1)), length(which(temptr_genes %in% phgenes2)),
                    length(which(temptr_genes %in% stress_genes)), length(which(temptr_genes %in% medgrowth_genes)),
                    length(which(temptr_genes %in% mg_c_genes)), length(which(temptr_genes %in% ox_genes)),
                    length(which(temptr_genes %in% aero_genes)), length(which(temptr_genes %in% juice_genes)),
                    length(which(temptr_genes %in% temptr_genes)), length(which(temptr_genes %in% temptr_genes2)),
                    length(which(temptr_genes %in% stringent_genes)))

study_mat[,10] <- c( length(which(temptr_genes2 %in% phgenes1)), length(which(temptr_genes2 %in% phgenes2)),
                    length(which(temptr_genes2 %in% stress_genes)), length(which(temptr_genes2 %in% medgrowth_genes)),
                    length(which(temptr_genes2 %in% mg_c_genes)), length(which(temptr_genes2 %in% ox_genes)),
                    length(which(temptr_genes2 %in% aero_genes)), length(which(temptr_genes2 %in% juice_genes)),
                    length(which(temptr_genes2 %in% temptr_genes)), length(which(temptr_genes2 %in% temptr_genes2)),
                    length(which(temptr_genes2 %in% stringent_genes)))

study_mat[,11] <- c( length(which(stringent_genes %in% phgenes1)), length(which(stringent_genes %in% phgenes2)),
                    length(which(stringent_genes %in% stress_genes)), length(which(stringent_genes %in% medgrowth_genes)),
                    length(which(stringent_genes %in% mg_c_genes)), length(which(stringent_genes %in% ox_genes)),
                    length(which(stringent_genes %in% aero_genes)), length(which(stringent_genes %in% juice_genes)),
                    length(which(stringent_genes %in% temptr_genes)), length(which(stringent_genes %in% temptr_genes2)),
                    length(which(stringent_genes %in% stringent_genes)))



















