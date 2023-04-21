source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(igraph)

ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv")
ec_genes <- read.csv("e_coli/Genes-in-E-coli.csv")

# nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,1] %in% ec_cyc[,1]),])
# nonreg_genes <- nonreg_genes[!(nonreg_genes[,1] %in% ec_cyc[,2]),] #1544 genes that are non-regulated
# 
# nr <- as.data.frame(cbind(nonreg_genes,nonreg_genes, rep(0, length(nonreg_genes))))
# setnames(nr, 1:3, colnames(ec_cyc))
# 
# ec_reg_full <- rbind(ec_cyc, nr)

# g <- graph.data.frame(ec_reg_full, directed=TRUE)
g <- graph.data.frame(ec_cyc, directed=TRUE)

E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3')))


#Loop analyses for list of plastic genes and non-plastic genes

#W in the right order ?

#Take only subnetw connected to target, otherwise too heavy to compute !!

gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, 3, nodes = which(colnames(E_coli_mat)=="cfa"), mode = 'all'))))

E_coli_mat2 <- t((get.adjacency(gg,sparse=FALSE, attr='V3')))

feedforward.to(gg, to=which(colnames(E_coli_mat2)=="cfa"), cutoff.max=16, cutoff.min=1)


###
# I have to add every genes in the matrix even though they are not regulated by/regulating anything !!!
# How to do that ?

################################################################################
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002

phgenes <- c("cfa","yeaQ","ompC","yfbE","yfbF","slp","hdeB","hdeD","osmY","cbpA",
             "hdeA","yahO","yccJ","asr","ydiZ","yebV","yhiM","wrbA","ygfR","yifC","ybaS","ycaC","gadB","gadA","yhiE","yhiX","yiaG","dps") #"CFA" not in regulonDB ?
phkgenes <- c("gadC")

which(colnames(E_coli_mat)=="cfa")

test <- c.count(list(E_coli_mat), cutoff.max = 15, cutoff.min = 1, randomFF=TRUE,
                target = which(colnames(E_coli_mat)=="cfa"))

for (i in phgenes) {
  c.count(list(E_coli_mat), cutoff.max = 20, cutoff.min = 1, randomFF=FALSE,
          target = which(colnames(E_coli_mat)==i))
}






# #RegulonDB database
# library("regutools")
# 
# regulondb_conn <- connect_database()
# 
# e_coli_regulondb <-  regulondb(database_conn = regulondb_conn,
#     organism = "chr", database_version = "11.1", genome_version = "11.1")
# 
# get_regulatory_network(e_coli_regulondb, cytograph = TRUE)
# 
# #Getting adjency matrix from regulonDB edges list
# library(igraph)
# el1 <- read.table("e_coli/reg_11.adj")
# el <- el1
# el[,1] <- el1[,1]
# el[,2] <- el1[,3]
# el[,3] <- el1[,2]
# g <- graph.data.frame(el, directed=TRUE)
