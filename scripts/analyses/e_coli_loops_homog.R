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

E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
E_coli_mat[is.na(E_coli_mat)] <- 0
#####################################

##PH
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002

phgenes <- read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)
dist <- 4 #Max distance for loop

phloops <- mclapply(phgenes$V1, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  E_coli_mat2 <- matrix(as.numeric(E_coli_mat2), ncol = ncol(E_coli_mat2), dimnames = dimnames(E_coli_mat2)) #convert to numeric matrix
  E_coli_mat2[is.na(E_coli_mat2)] <- 0
  cc <- homog.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, randomFF=FALSE,
                target = which(colnames(E_coli_mat2)==gene))
  return(cc)
}, mc.cores = 4)

write.csv(rbindlist(phloops), "ph_plast_loops_homog.csv")




###Loop number distibution
phloops <- mclapply(phgenes$V1, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  E_coli_mat2 <- matrix(as.numeric(E_coli_mat2), ncol = ncol(E_coli_mat2), dimnames = dimnames(E_coli_mat2)) #convert to numeric matrix
  E_coli_mat2[is.na(E_coli_mat2)] <- 0
  cc <- loops_n.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, target = which(colnames(E_coli_mat2)==gene))
  return(cc)
}, mc.cores = 4)

df <- rbindlist(phloops)


hist(t(df[,1]))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFl", "Incoherent FFL","No FFL"))


##########################################
#10.1038/srep45303
#List of plastic genes for 3 different conditions: carbone source, Mg stress, Na+ stress
mg_carbon <- unique(subset(read.csv("e_coli/Plast_genes/C_Mg_Na/41598_2017_BFsrep45303_MOESM60_ESM.csv"), dataType=="mrna")[,2])








