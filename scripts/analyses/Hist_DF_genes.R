#Nbr de gène x nbr de fois reporté plastique
#Fig supp
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


all_plast_genes <- c(phgenes1, phgenes2, stress_genes, medgrowth_genes , mg_c_genes, ox_genes, aero_genes, juice_genes, temptr_genes, temptr_genes2, stringent_genes)
#Number of plastic genes in annotation
all_plast_genes1 <- e_coli_gene_name(all_plast_genes)
length(unique(all_plast_genes1))

###
DF_genes <- as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                                               ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2)))


pdf(paste0("figures/Hist_DF_genes",".pdf"), width=3, height=4)
hist(DF_genes$Freq, breaks = seq(0.5,7.5, 1), xlab="Times being reported as DF", ylab="Number of gene", main="")
dev.off()

