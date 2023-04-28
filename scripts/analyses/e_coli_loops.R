source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(igraph)

ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv")
ec_genes <- read.csv("e_coli/All-genes-of-E.-coli-K-12-substr.-MG1655---GO.csv", sep ="\t")

nonreg_genes <- as.data.frame(ec_genes[!(ec_genes[,1] %in% ec_cyc[,1]),])
nonreg_genes <- nonreg_genes[!(nonreg_genes[,1] %in% ec_cyc[,2]),] #1704 genes that are non-regulated

nr <- as.data.frame(cbind(nonreg_genes[,1],nonreg_genes[,1], rep(0, length(nonreg_genes))))
setnames(nr, 1:3, colnames(ec_cyc))

ec_reg_full <- rbind(ec_cyc, nr)

g <- graph.data.frame(ec_reg_full, directed=TRUE)
E_coli_mat <- t((get.adjacency(g,sparse=FALSE, attr='V3')))

E_coli_mat <- matrix(as.numeric(E_coli_mat), ncol = ncol(E_coli_mat), dimnames = dimnames(E_coli_mat)) #convert to numeric matrix
E_coli_mat[is.na(E_coli_mat)] <- 0


##PH
#Tucker et al., 2002 ; DOI : 10.1128/JB.184.23.6551-6558.2002

phgenes <- read.table("e_coli/Plast_genes/Ph/plastic_genes.txt", sep ="\t", header=FALSE)
dist <- 4 #Max distance for loop

phloops <- mclapply(phgenes$V1, function(gene) {
  gg <- induced.subgraph(g, vids = as.vector(unlist(neighborhood(g, dist, nodes = which(colnames(E_coli_mat)==gene), mode = 'all'))))
  E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
  E_coli_mat2 <- matrix(as.numeric(E_coli_mat2), ncol = ncol(E_coli_mat2), dimnames = dimnames(E_coli_mat2)) #convert to numeric matrix
  E_coli_mat2[is.na(E_coli_mat2)] <- 0
  cc <- c.count(list(E_coli_mat2), cutoff.max = 4, cutoff.min = 1, randomFF=FALSE,
          target = which(colnames(E_coli_mat2)==gene))
  return(cc)
}, mc.cores = 4)

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







np_medgrowthloops <- read.csv("np_medium_growth_loops.csv")[,2:5]
plast_medgrowthloops <- read.csv("plast_medium_growth_loops.csv")[,2:5]

df <- rbind(colSums(np_medgrowthloops)/nrow(np_medgrowthloops), colSums(plast_medgrowthloops)/nrow(plast_medgrowthloops), colSums(rbindlist(phloops))/nrow(rbindlist(phloops)))
rownames(df) <- c("np_medgrowth", "plast_medgrowth", "pH")

layout(matrix(c(1:1), 1, 1, byrow = TRUE))

barplot(t(df[,1:3]), col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFl", "Incoherent FFL","No FFL"))


pdf("figures/FFL_distrib_e_coli.pdf", width=4, height=4)
  layout(matrix(c(1), 1, 1, byrow = TRUE))
  par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
  barplot(t(df[,1:3]), col=c(7,3,"grey"))
  legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
         legend=c("Coherent FFL", "Incoherent FFL", "No FFL"))
dev.off()





