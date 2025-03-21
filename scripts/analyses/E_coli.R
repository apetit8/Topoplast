# source("scripts/functions/functions.R")
# source("scripts/functions/detectloops.R")
# memory.limit(size=2500)
# #########################################
# #Genetic data
# ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv") #List of regulations from Ecocyc
# ec_genes <- read.csv("e_coli/ncbi_dataset_K-12_annotation.csv", sep ="\t") #List of E coli genes from Ecocyc
# 
# #Transcriptions factors
# TF_genes <- unique(ec_cyc[,1]) # Here, TF = all regulating genes from Ecocyc reg data
# #Adding non-regulated genes in regulation data
# nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,2] %in% ec_cyc[,1]),]) #Regulators
# nonreg_genes <- nonreg_genes[!(nonreg_genes[,2] %in% ec_cyc[,2]),] #Regulatees
# nr <- as.data.frame(cbind(nonreg_genes[,2],nonreg_genes[,2], rep(0, nrow(nonreg_genes)))) #Same file format as ec_cyc
# setnames(nr, 1:3, colnames(ec_cyc))
# ec_reg_full <- rbind(ec_cyc, nr) #Reg file with every genes
# ##########
# g <- graph.data.frame(ec_reg_full, directed=TRUE)
# E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3'))) #t() to have regulators in columns and regulees in rows
# 
# E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
# E_coli_mat[is.na(E_coli_mat)] <- 0 #fill NA to 0, mandatory for later analyses
# ##########
# igraph_options(return.vs.es=F)
# #########################################
# #DEBUG
# # ff <- as.matrix(E_coli_mat[1:200,1:200])
# # Rprof()
# # cc <- FFL.coherence(list(ff), cutoff.max = 4, cutoff.min = 1, target = 35)
# # Rprof(NULL)
# # summaryRprof()
# #########################################
# #Analyses of plastic genes from different sources
# edges1 <- 2
# edges2 <- 1
# from=FALSE
# all_plast_genes <- data.frame()
# #########################################
#
##PH
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002
#Maurer et al., 2004 ; 

phgenes1 <- read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)
phgenes1 <- e_coli_gene_name(phgenes1[,1])

phloops <- e_coli_prep_analyses(phgenes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)

write.csv(rbind(phloops), paste0("scripts/data/ph1_plast_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(phloops))

print("Ph1 done!")
########################################
#Maurer et al., 2004 ; 
phgenes2 <- read.table("e_coli/Plast_genes/pH_2/plastic_genes.txt", sep ="\t", header=FALSE)
phgenes2 <- unique(e_coli_gene_name(phgenes2[,1]))
phgenes <- phgenes2[!(phgenes2 %in% all_plast_genes$V1)] 

phloops <- e_coli_prep_analyses(phgenes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(phloops), paste0("scripts/data/ph2_plast_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(phloops))

print("Ph2 done!")
########################################
#10.1128/JB.01092-07
#Compendium of stressor
stress_genes <- as.data.frame(unique(c(subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S5.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S6.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S7.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S8.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5],
                                       subset(read.csv("e_coli/Plast_genes/Mult_stress/Table_S9.csv", sep ="\t", header=TRUE), Differential.expression != "ns")[,5])))
stress_genes <- e_coli_gene_name(stress_genes[,1])
stress_genes1 <- stress_genes[!(stress_genes %in% all_plast_genes$V1)] 

plast_stress <- e_coli_prep_analyses(stress_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_stress, subset(all_plast_genes, V1 %in% stress_genes)), paste0("scripts/data/plast_stress_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_stress))

print("Stressors done!")
########################################
##Medium of growth
#Feugeas et al., 2016 ; DOI : 10.1093/molbev/msw105
medgrowth_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S1.csv", sep ="\t", header=TRUE)
medgrowth_genes <- e_coli_gene_name(medgrowth_genes[,1])
medgrowth_genes1 <- medgrowth_genes[!(medgrowth_genes %in% all_plast_genes$V1)]
#
plast_medgrowthloops <- e_coli_prep_analyses(medgrowth_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)

write.csv(rbind(plast_medgrowthloops, subset(all_plast_genes, V1 %in% medgrowth_genes)), paste0("scripts/data/plast_medium_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_medgrowthloops))

print("Medgrowth plastic genes done!")

########################################
#10.1038/srep45303
#List of plastic genes for 3 different conditions: carbone source, Mg stress, Na+ stress
mg_c_genes <- as.data.frame(unique(subset(read.csv("e_coli/Plast_genes/C_Mg_Na/41598_2017_BFsrep45303_MOESM60_ESM.csv"), dataType=="mrna")[,2]))
mg_c_genes <- e_coli_gene_name(mg_c_genes[,1]) 
mg_c_genes1 <- mg_c_genes[!(mg_c_genes %in% all_plast_genes$V1)] 

plast_mg_c  <- e_coli_prep_analyses(mg_c_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
  
write.csv(rbind(plast_mg_c , subset(all_plast_genes, V1 %in% mg_c_genes)), paste0("scripts/data/mg_c_growth_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_mg_c))

print("Mg C done!")
########################################
#10.1128/AEM.00914-09
#List of plastic genes oxydative stress
ox_genes <- read.csv("e_coli/Plast_genes/Oxydative_stress/table_5_genes_updated_names.txt", sep ="\t", header=FALSE)
#ox_genes <- ox_genes[(ox_genes[,1] %in% colnames(E_coli_mat)),] 
ox_genes <- e_coli_gene_name(ox_genes[,1])
ox_genes1 <- ox_genes[!(ox_genes %in% all_plast_genes$V1)] 

plast_ox <- e_coli_prep_analyses(ox_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_ox, subset(all_plast_genes, V1 %in% ox_genes)), paste0("scripts/data/plast_ox_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_ox))

print("Ox done!")
########################################
#10.1007/s00253-018-9083-5
#List of DEGS between aerozolisation and liquid suspension
aero_genes <- read.table("e_coli/Plast_genes/Aero_liquid/Table_2_3_4.txt", sep ="\t", header=FALSE)
aero_genes <- e_coli_gene_name(aero_genes[,1])
aero_genes1 <- aero_genes[!(aero_genes %in% all_plast_genes$V1)] 

plast_aero <- e_coli_prep_analyses(aero_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_aero, subset(all_plast_genes, V1 %in% aero_genes)), paste0("scripts/data/plast_aero_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_aero))

print("Aero done!")
########################################
#10.1128/AEM.02841-08
#Apple Juice
juice_genes <- read.table("e_coli/Plast_genes/apple_juice/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
juice_genes <- e_coli_gene_name(juice_genes[,1]) 
juice_genes1 <- juice_genes[!(juice_genes %in% all_plast_genes$V1)] 


plast_juice <- e_coli_prep_analyses(juice_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_juice, subset(all_plast_genes, V1 %in% juice_genes)), paste0("scripts/data/plast_juice_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_juice))

print("Juice done!")
########################################
#10.1128/JB.01929-06
#Human body temperature
temptr_genes <- read.table("e_coli/Plast_genes/human_temp/Table_1_genes_updated_names.txt", sep ="\t", header=FALSE)
temptr_genes <- e_coli_gene_name(temptr_genes[,1])
temptr_genes1 <- temptr_genes[!(temptr_genes %in% all_plast_genes$V1)] 

plast_temptr <- e_coli_prep_analyses(temptr_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_temptr, subset(all_plast_genes, V1 %in% temptr_genes)), paste0("scripts/data/plast_temptr1_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_temptr))

print("Temperature done!")
########################################
#
#Temperature
temptr_genes2 <- read.table("e_coli/Plast_genes/Temperature/plastic_genes.txt", sep ="\t", header=FALSE)
temptr_genes2 <- e_coli_gene_name(temptr_genes2[,1])

temptr_genes1 <- temptr_genes2[!(temptr_genes2 %in% all_plast_genes$V1)] 

plast_temptr <- e_coli_prep_analyses(temptr_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_temptr, subset(all_plast_genes, V1 %in% temptr_genes2)), paste0("scripts/data/plast_temptr2_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_temptr))

print("Temperature done!")


########################################
#10.1128/JB.01092-07
#stringent response
stringent_genes <- as.data.frame(unique(read.table("e_coli/Plast_genes/Stringent_response/All_genes.txt", sep ="\t", header=FALSE)[,1]))
stringent_genes <- e_coli_gene_name(stringent_genes[,1]) 
stringent_genes1 <- stringent_genes[!(stringent_genes %in% all_plast_genes$V1)] 

plast_stringent <- e_coli_prep_analyses(stringent_genes1, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(rbind(plast_stringent, subset(all_plast_genes, V1 %in% stringent_genes)), paste0("scripts/data/plast_stringent_",csvname,".csv"))
all_plast_genes <- rbind(all_plast_genes, as.data.frame(plast_stringent))

print("Stringent response done!")
########################################
#Plastic genes : only genes occuring in at least 2 studies are kept
###
Subplastic_genes <- as.character(subset( as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                                               ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2))), Freq >=2)[,1])

write.csv(subset(all_plast_genes, V1 %in% Subplastic_genes), paste0("scripts/data/plast_genes_",csvname,".csv"))

print("All plastic genes done!")
########################################
#Genes that were not find as responding to the environment in our data corpus
nonplast_genes <- as.data.frame(colnames(E_coli_mat))
nonplast_genes <- nonplast_genes[!(nonplast_genes[,1] %in% all_plast_genes$V1),]

nonplast_FFL <- e_coli_prep_analyses(nonplast_genes, g, E_coli_mat, fun=fun, edges1=edges1, edges2=edges2, cores=1, from = from)
write.csv(nonplast_FFL, paste0("scripts/data/nonplast_",csvname,".csv"))

print("Non plastic genes done!")

########################################
print(paste0("Plast genes: ", nrow(all_plast_genes), "; Plast genes occuring more than twice: ", nrow(subset(all_plast_genes, V1 %in% Subplastic_genes)), "; Non-plastic genes : ", nrow(nonplast_FFL)))

length(unique(all_plast_genes[,1]))
