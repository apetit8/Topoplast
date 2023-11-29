#Fig supp total number of regulations
#Take E_coli matrix, subset plastic genes, get row sum and compare to non plastic gene row sum
source("scripts/functions/functions.R")
library(purrr)
#########################################
pdfname <- "figures/fig_supp_reg_nbr"
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
################################################################################
#E. coli
plastgenes <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,2]
npgenes <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,2]
E_coli_mat_plast <- E_coli_mat[plastgenes, ]
E_coli_mat_np <- E_coli_mat[npgenes, ]

################################################################################
#Positive and neg regulations
#E. coli
plastmat <- E_coli_mat_plast
plastmat[(plastmat == -1)] <- 0
reg_number_pos_plast <- rowSums2(plastmat)
plastmat <- E_coli_mat_plast
plastmat[(plastmat == 1)] <- 0
reg_number_neg_plast <- rowSums2(abs(plastmat))
#
plastmat <- E_coli_mat_np
plastmat[(plastmat == -1)] <- 0
reg_number_pos_np <- rowSums2(plastmat)
plastmat <- E_coli_mat_np
plastmat[(plastmat == 1)] <- 0
reg_number_neg_np <- rowSums2(abs(plastmat))

#Simulations

#Simulations
topo_plastic <- readRDS("scripts/data/list_plastic_topo_full_netw.Rds")
topo_nnplast <- readRDS("scripts/data/list_nnplast_topo_full_netw.Rds")

THreg_sum_pos_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == -1)] <- 0
return(rowSums2(abs(i))[2])} )) # [2] : the target of the network is in position #2
THreg_sum_neg_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i))[2])} ))
#
THreg_sum_pos_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 0
return(rowSums2(abs(i))[2])} ))
THreg_sum_neg_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == 1)] <- 0
return(rowSums2(abs(i))[2])} ))



pdf(paste0("figures/Fig_histo_pos",".pdf"), width=5, height=4)
par(mar = c(4,2,1,0.5), mfrow=c(1,4), mgp=c(1,0.5,0))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="honeydew2", border=NA, xpd=TRUE)
barplot(as.vector(table(reg_number_pos_np))*100/sum(as.vector(table(reg_number_pos_np))), col="lavender", names.arg=names(table(reg_number_pos_np)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE)
legend(-25, 15.5, fill=c("lavender", "lightskyblue"), legend = c("Non plastic gene","Plastic gene"), bg='white', xpd=TRUE, cex=1)
text(119.5,-1.8, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1.6)
text(90,12, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="honeydew2", border=NA, xpd=TRUE)
barplot(as.vector(table(reg_number_pos_plast))*100/sum(as.vector(table(reg_number_pos_plast))), col="lightskyblue", names.arg=names(table(reg_number_pos_plast)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(-27,-1.8, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1.6)
text(119.5,17, substitute(paste(bold("Number of activation per gene"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="cornsilk", border=NA, xpd=TRUE)
barplot(as.vector(table( THreg_sum_pos_np))*100/sum(as.vector(table( THreg_sum_pos_np))), col="lavender", names.arg=names(table( THreg_sum_pos_np)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(119.5,-1.8, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.6)
text(-27,17, substitute(paste(bold("Number of activation per gene"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="cornsilk", border=NA, xpd=TRUE)
barplot(as.vector(table(THreg_sum_pos_plast))*100/sum(as.vector(table(THreg_sum_pos_plast))), col="lightskyblue", names.arg=names(table(THreg_sum_pos_plast)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(-27,-1.8, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.6)
text(-28,17, substitute(paste(bold("e"))), xpd=TRUE, cex=1.6)
dev.off()


pdf(paste0("figures/Fig_histo_neg",".pdf"), width=5, height=4)
par(mar = c(4,2,1,0.5), mfrow=c(1,4), mgp=c(1,0.5,0))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="honeydew2", border=NA, xpd=TRUE)
barplot(as.vector(table(reg_number_neg_np))*100/sum(as.vector(table(reg_number_neg_np))), col="lavender", names.arg=names(table(reg_number_neg_np)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE)
legend(-25, 15.5, fill=c("lavender", "lightskyblue"), legend = c("Non plastic gene","Plastic gene"), bg='white', xpd=TRUE, cex=1)
text(119.5,-1.8, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1.6)
text(90,12, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="honeydew2", border=NA, xpd=TRUE)
barplot(as.vector(table(reg_number_neg_plast))*100/sum(as.vector(table(reg_number_neg_plast))), col="lightskyblue", names.arg=names(table(reg_number_neg_plast)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(-27,-1.8, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1.6)
text(119.5,17, substitute(paste(bold("Number of inhibition per gene"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="cornsilk", border=NA, xpd=TRUE)
barplot(as.vector(table( THreg_sum_neg_np))*100/sum(as.vector(table( THreg_sum_neg_np))), col="lavender", names.arg=names(table( THreg_sum_neg_np)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(119.5,-1.8, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.6)
text(-27,17, substitute(paste(bold("Number of inhibition per gene"))), xpd=TRUE, cex=1.6)
#
par(mar = c(4,2,1,0.5))
barplot(t(as.data.frame(c(0,0,0))), col = NA, border = NA, axes = FALSE, ylim = c(0,17), xlim = c(0,100), xaxt = "n")
polygon(x=c(-40, -40, 120, 120),  y=c(16, -20, -20, 16),  col="cornsilk", border=NA, xpd=TRUE)
barplot(as.vector(table(THreg_sum_neg_plast))*100/sum(as.vector(table(THreg_sum_neg_plast))), col="lightskyblue", names.arg=names(table(THreg_sum_neg_plast)), axes = TRUE, space = 0, horiz=TRUE, xlab= "Counts", add=TRUE, xpd=TRUE)
text(-27,-1.8, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1.6)
text(-28,17, substitute(paste(bold("e"))), xpd=TRUE, cex=1.6)
dev.off()

################################################################################
################################################################################
#All regulations
#E. coli
plastmat <- E_coli_mat_plast
plastmat[(plastmat == -1)] <- 1
reg_number_plast <- as.data.frame(rowSums2(plastmat))
reg_number_plast$Type <- "2empl"
colnames(reg_number_plast) <- c("nbr","Type")
#
plastmat <- E_coli_mat_np
plastmat[(plastmat == -1)] <- 1
reg_number_np <- as.data.frame(rowSums2(plastmat))
reg_number_np$Type <- "1emnp"
colnames(reg_number_np) <- c("nbr","Type")

#Simulations
topo_plastic <- readRDS("scripts/data/list_plastic_topo_full_netw.Rds")
topo_nnplast <- readRDS("scripts/data/list_nnplast_topo_full_netw.Rds")

THreg_sum_plast <- as.data.frame(unlist(lapply(topo_plastic, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_plast$Type <- "4thpl"
colnames(THreg_sum_plast) <- c("nbr","Type")

#
THreg_sum_np <- as.data.frame(unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 1
return(rowSums2(abs(i))[2])} )))
THreg_sum_np$Type <- "3thnp"
colnames(THreg_sum_np) <- c("nbr","Type")

df <- rbind(reg_number_plast, reg_number_np, THreg_sum_np, THreg_sum_plast )

pdf(paste0("figures/Reg_nbr",".pdf"), width=5, height=4)
par(mar = c(3.5,3.5, 1,1))
boxplot(df$nbr ~ df$Type,  main="", col = NA, border = NA, axes = FALSE,
        yaxt="", space=c(0.3,0.1,0.4,0.1,0.4), ylab = "", xlab = "")
title(ylab = "Number of regulation toward gene", line=2)
polygon(x=c(0.35, 0.35, 2.5, 2.5),  y=c(28, -20, -20, 28),  col="honeydew2", border=NA, xpd=TRUE)
polygon(x=c(2.5, 2.5, 120, 120),  y=c(28, -20, -20, 28),  col="cornsilk", border=NA, xpd=TRUE)
boxplot(df$nbr ~ df$Type, add=TRUE, col=c("lavender", "lightskyblue"), 
        names=c("Non plastic","Plastic","Non plastic","Plastic"),
        space=c(0.3,0.1,0.4,0.1,0.4), frame=FALSE)
text(1.5, -6.5, substitute(paste(bold("E. coli"))), xpd=TRUE, cex=1)
text(3.5, -6.5, substitute(paste(bold("Simulations"))), xpd=TRUE, cex=1)
text(1.5,25, substitute(paste(bold("***"))), xpd=TRUE, cex=1.6)
dev.off()

mean(reg_number_plast[,1])
mean(reg_number_np[,1])
mean(THreg_sum_np[,1])
mean(THreg_sum_plast[,1])
