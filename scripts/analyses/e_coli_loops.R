source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(igraph)

ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv")
ec_genes <- read.csv("e_coli/Genes-in-E-coli.csv")

nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,1] %in% ec_cyc[,1]),])
nonreg_genes <- nonreg_genes[!(nonreg_genes[,1] %in% ec_cyc[,2]),] #1544 genes that are non-regulated

nr <- as.data.frame(cbind(nonreg_genes,nonreg_genes, rep(0, length(nonreg_genes))))
setnames(nr, 1:3, colnames(ec_cyc))

ec_reg_full <- rbind(ec_cyc, nr)

g <- graph.data.frame(ec_reg_full, directed=TRUE)
E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3')))


##PH
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002

phgenes <- c("cfa","yeaQ","ompC","yfbE","yfbF","slp","hdeB","hdeD","osmY","cbpA",
             "hdeA","yahO","yccJ","asr","ydiZ","yebV","yhiM","wrbA","ygfR","yifC",
             "ybaS","ycaC","gadB","gadA","yhiE","yhiX","yiaG","dps")
dist <- 4 #Max distance for loop

phloops <- mclapply(phgenes, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  cc <- c.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, randomFF=FALSE,
          target = which(colnames(E_coli_mat2)==gene))
  return(cc)
},mc.cores = 4)

write.csv(rbindlist(phloops), "ph_plast_loops.csv")


##Medium of growth
#Feugeas et al., 2016 ; DOI : 10.1093/molbev/msw105

medgrowth_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S1.csv", sep ="\t", header=TRUE)[,1]
non_envir_genes <- read.csv("e_coli/Plast_genes/Growth_envir/Supp_tables_S3.csv", sep ="\t", header=TRUE)[,1]


plast_medgrowthloops <- mclapply(medgrowth_genes, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  cc <- c.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, randomFF=FALSE,
                target = which(colnames(E_coli_mat2)==gene))
  return(cc)
}, mc.cores = 4)
write.csv(rbindlist(plast_medgrowthloops), "plast_medium_growth_loops.csv")

#Non plastic genes (expression depending on strains and not on tested environment)
np_medgrowthloops <- mclapply(non_envir_genes, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  cc <- c.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, randomFF=FALSE,
                target = which(colnames(E_coli_mat2)==gene))
  return(cc)
}, mc.cores = 4)
write.csv(rbindlist(np_medgrowthloops), "np_medium_growth_loops.csv")





