source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(igraph)
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
#########################################
#Genetic data
ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd_2024_01_29.csv") #List of regulations from Ecocyc
ec_cyc <- subset(ec_cyc, V1!=V2)
##########
freq_genes <- as.data.frame(table(c(stringent_genes, temptr_genes, temptr_genes2, juice_genes, aero_genes,
                                    ox_genes,mg_c_genes,medgrowth_genes,stress_genes,phgenes1,phgenes2)))
for (i in 1:nrow(ec_cyc)) {
  if(ec_cyc[i,1] %in% freq_genes[,1]) ec_cyc[i,4] <- freq_genes[which(freq_genes[,1]==ec_cyc[i,1]),2]
  else ec_cyc[i,4] <- 0
}
##########
TF_genes <- unique(ec_cyc[,1]) # Here, TF = all regulating genes from Ecocyc reg data
ec_cyc <- subset(ec_cyc, V2 %in% TF_genes)
g1 <- graph.data.frame(ec_cyc, directed=TRUE)
##########
igraph_options(return.vs.es=F)
##########


cairo_pdf("figures/Ecoli_TFs.pdf", width=20, height=20)
plot(g1, layout=layout_nicely, edge.color=ifelse(E(g1)$V3 == 1, "black",ifelse(E(g1)$V3 == -1, "red","white")),
     vertex.size=5, main="Ecoli TFs regulatory network", edge.curved=TRUE,
     vertex.color=ifelse(E(g1)$V4 >= 3, "steelblue",ifelse(E(g1)$V4 >= 2, "lightblue",ifelse(E(g1)$V4 == 1, "white","grey"))))
legend(0.5, 1, legend=c("Reported as DE\nat least 3 times\n", "Reported as DE\n2 times\n", "Reported as DE\n1 time\n", "Never reported\nas DE\n"),
       col=c( "steelblue", "lightblue", "black","grey"), bty=1, cex=2, bg="white", pch = c(19,19,21,19))
legend(0.5, 0.4, legend=c("Activation", "Inhibition"),
       col=c("black", "tomato"), lty=1, cex=2, bg="white")
dev.off()

################################################################################
ec_cyc <- subset(ec_cyc, V2 %in% TF_genes)
g1 <- graph.data.frame(ec_cyc, directed=FALSE)
##########
igraph_options(return.vs.es=F)
##########
##SLideshow defense
cairo_pdf("figures/Ecoli_TFs_COUVERTURE3.pdf", width=30, height=20)
plot(g1, layout=layout.fruchterman.reingold, edge.color=ifelse(E(g1)$V3 == 1, "lightseagreen",ifelse(E(g1)$V3 == -1, "darkseagreen","white")),
     edge.curved=TRUE, vertex.frame.color="steelblue", #vertex.label=NA,
     vertex.color=ifelse(E(g1)$V4 >= 3, "aquamarine3",ifelse(E(g1)$V4 >= 2, "darkolivegreen2",ifelse(E(g1)$V4 == 1, "white","grey"))),
     vertex.size=ifelse(E(g1)$V4 >= 3, 4,ifelse(E(g1)$V4 >= 2, 4,ifelse(E(g1)$V4 == 1, 4,4))))
dev.off()

