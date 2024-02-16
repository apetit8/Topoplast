source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd_2024_01_29.csv") #List of regulations from Ecocyc
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

########################################

#PROBLEM: %in% also count it if the string of character in inside a larger string. Not what I want at all ; 

study_mat <- matrix(data=NA, nrow=11, ncol=12)
colnames(study_mat) <- c("Tested environment","Tucker2002","Maurer2004","Bhatia2022","Feugas2016","Caglar2017","Wang2009","Ng2018","Bergholz2009","White-Ziegler2007","Kim2020","Durfee2008")
rownames(study_mat) <- colnames(study_mat)[2:12]
# c(phgenes1,phgenes2,stress_genes,medgrowth_genes,mg_c_genes,ox_genes,aero_genes,juice_genes,temptr_genes,temptr_genes2,stringent_genes)

study_mat[,2] <- c( sum(phgenes1 %in% phgenes1), sum(phgenes1 %in% phgenes2),
                    sum(phgenes1 %in% stress_genes), sum(phgenes1 %in% medgrowth_genes),
                    sum(phgenes1 %in% mg_c_genes), sum(phgenes1 %in% ox_genes),
                    sum(phgenes1 %in% aero_genes), sum(phgenes1 %in% juice_genes),
                    sum(phgenes1 %in% temptr_genes), sum(phgenes1 %in% temptr_genes2),
                    sum(phgenes1 %in% stringent_genes))

study_mat[,3] <- c( sum(phgenes2 %in% phgenes1), sum(phgenes2 %in% phgenes2),
                    sum(phgenes2 %in% stress_genes), sum(phgenes2 %in% medgrowth_genes),
                    sum(phgenes2 %in% mg_c_genes), sum(phgenes2 %in% ox_genes),
                    sum(phgenes2 %in% aero_genes), sum(phgenes2 %in% juice_genes),
                    sum(phgenes2 %in% temptr_genes), sum(phgenes2 %in% temptr_genes2),
                    sum(phgenes2 %in% stringent_genes))

study_mat[,4] <- c( sum(stress_genes %in% phgenes1), sum(stress_genes %in% phgenes2),
                    sum(stress_genes %in% stress_genes), sum(stress_genes %in% medgrowth_genes),
                    sum(stress_genes %in% mg_c_genes), sum(stress_genes %in% ox_genes),
                    sum(stress_genes %in% aero_genes), sum(stress_genes %in% juice_genes),
                    sum(stress_genes %in% temptr_genes), sum(stress_genes %in% temptr_genes2),
                    sum(stress_genes %in% stringent_genes))

study_mat[,5] <- c( sum(medgrowth_genes %in% phgenes1), sum(medgrowth_genes %in% phgenes2),
                    sum(medgrowth_genes %in% stress_genes), sum(medgrowth_genes %in% medgrowth_genes),
                    sum(medgrowth_genes %in% mg_c_genes), sum(medgrowth_genes %in% ox_genes),
                    sum(medgrowth_genes %in% aero_genes), sum(medgrowth_genes %in% juice_genes),
                    sum(medgrowth_genes %in% temptr_genes), sum(medgrowth_genes %in% temptr_genes2),
                    sum(medgrowth_genes %in% stringent_genes))

study_mat[,6] <- c( sum(mg_c_genes %in% phgenes1), sum(mg_c_genes %in% phgenes2),
                    sum(mg_c_genes %in% stress_genes), sum(mg_c_genes %in% medgrowth_genes),
                    sum(mg_c_genes %in% mg_c_genes), sum(mg_c_genes %in% ox_genes),
                    sum(mg_c_genes %in% aero_genes), sum(mg_c_genes %in% juice_genes),
                    sum(mg_c_genes %in% temptr_genes), sum(mg_c_genes %in% temptr_genes2),
                    sum(mg_c_genes %in% stringent_genes))


study_mat[,7] <- c( sum(ox_genes %in% phgenes1), sum(ox_genes %in% phgenes2),
                    sum(ox_genes %in% stress_genes), sum(ox_genes %in% medgrowth_genes),
                    sum(ox_genes %in% mg_c_genes), sum(ox_genes %in% ox_genes),
                    sum(ox_genes %in% aero_genes), sum(ox_genes %in% juice_genes),
                    sum(ox_genes %in% temptr_genes), sum(ox_genes %in% temptr_genes2),
                    sum(ox_genes %in% stringent_genes))

study_mat[,8] <- c( sum(aero_genes %in% phgenes1), sum(aero_genes %in% phgenes2),
                    sum(aero_genes %in% stress_genes), sum(aero_genes %in% medgrowth_genes),
                    sum(aero_genes %in% mg_c_genes), sum(aero_genes %in% ox_genes),
                    sum(aero_genes %in% aero_genes), sum(aero_genes %in% juice_genes),
                    sum(aero_genes %in% temptr_genes), sum(aero_genes %in% temptr_genes2),
                    sum(aero_genes %in% stringent_genes))

study_mat[,9] <- c( sum(juice_genes %in% phgenes1), sum(juice_genes %in% phgenes2),
                    sum(juice_genes %in% stress_genes), sum(juice_genes %in% medgrowth_genes),
                    sum(juice_genes %in% mg_c_genes), sum(juice_genes %in% ox_genes),
                    sum(juice_genes %in% aero_genes), sum(juice_genes %in% juice_genes),
                    sum(juice_genes %in% temptr_genes), sum(juice_genes %in% temptr_genes2),
                    sum(juice_genes %in% stringent_genes))

study_mat[,10] <- c( sum(temptr_genes %in% phgenes1), sum(temptr_genes %in% phgenes2),
                    sum(temptr_genes %in% stress_genes), sum(temptr_genes %in% medgrowth_genes),
                    sum(temptr_genes %in% mg_c_genes), sum(temptr_genes %in% ox_genes),
                    sum(temptr_genes %in% aero_genes), sum(temptr_genes %in% juice_genes),
                    sum(temptr_genes %in% temptr_genes), sum(temptr_genes %in% temptr_genes2),
                    sum(temptr_genes %in% stringent_genes))

study_mat[,11] <- c( sum(temptr_genes2 %in% phgenes1), sum(temptr_genes2 %in% phgenes2),
                    sum(temptr_genes2 %in% stress_genes), sum(temptr_genes2 %in% medgrowth_genes),
                    sum(temptr_genes2 %in% mg_c_genes), sum(temptr_genes2 %in% ox_genes),
                    sum(temptr_genes2 %in% aero_genes), sum(temptr_genes2 %in% juice_genes),
                    sum(temptr_genes2 %in% temptr_genes), sum(temptr_genes2 %in% temptr_genes2),
                    sum(temptr_genes2 %in% stringent_genes))

study_mat[,12] <- c( sum(stringent_genes %in% phgenes1), sum(stringent_genes %in% phgenes2),
                    sum(stringent_genes %in% stress_genes), sum(stringent_genes %in% medgrowth_genes),
                    sum(stringent_genes %in% mg_c_genes), sum(stringent_genes %in% ox_genes),
                    sum(stringent_genes %in% aero_genes), sum(stringent_genes %in% juice_genes),
                    sum(stringent_genes %in% temptr_genes), sum(stringent_genes %in% temptr_genes2),
                    sum(stringent_genes %in% stringent_genes))
study_mat[,1] <- c("pH","pH","Compendium of stressors","Medium of growth","Sources of C and Mg", "Oxidative stress","Aerosolization","Exposure to apple juice","Temperature change","Temperature change","Stringent response")

study_mat
write.csv(study_mat, "figures/Supplemental_tab1.csv")


