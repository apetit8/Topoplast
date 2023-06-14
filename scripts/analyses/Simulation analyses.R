source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/10g", recursive = TRUE)
genes <- 10
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
filename <- "10g"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
gen <- 10000 #max(df.10$Gen)
#####################
#DATA
##################
df.10 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10$envir <- str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,3]
df.10$anc_id <- str_split(str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]
#############
topo.anticor10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Anticorrelated"), target=target,
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Correlated"), target=target,
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.no_sel10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Control_no_sel"), target=target,
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Control_sel"), target=target,
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

# ifelse(!dir.exists(file.path("scripts/data/essential_motif")), dir.create("scripts/data/essential_motif"), FALSE)
# lapply(topo.anticor10, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.anticor10.csv'  , append= T, sep=',' ))
# lapply(topo.anticor10, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.corr10.csv'  , append= T, sep=',' ))
# lapply(topo.anticor10, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.no_sel10.csv'  , append= T, sep=',' ))
# lapply(topo.anticor10, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.sel10.csv'  , append= T, sep=',' ))
################
##### FFL
Anticor10 <- FFL.type2(topo.anticor10, edges1 = 2, edges2 = 1)
Corr10 <- FFL.type2(topo.corr10, edges1 = 2, edges2 = 1)
No_sel10 <- FFL.type2(topo.no_sel10, edges1 = 2, edges2 = 1) #
Sel10 <- FFL.type2(topo.sel10, edges1 = 2, edges2 = 1)

df <- as.data.frame(rbind(colSums(rbind(Anticor10,Corr10))/2*100/nrow(Anticor10), colSums(Anticor10)*100/nrow(Anticor10), colSums(Corr10)*100/nrow(Anticor10), colSums(No_sel10)*100/nrow(No_sel10), colSums(Sel10)*100/nrow(Sel10)))
rownames(df) <- c("Plastic\ngenes","Anticor","Cor","No_sel","Non-plastic\ngenes")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))
write.csv(df, paste0("scripts/data/",filename,"_FFL",".csv"))

pdf(paste0(filename,"prop_FFL",".pdf"), width=3.5, height=4)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[c(5,1),1:2]), col=c("gold", "grey"), main="Theoretical Prediction")
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("gold", "grey"),
       legend=c( "at least 1 FFL","No FFL") )
dev.off()

pdf(paste0("figures/",filename,"_FFL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:10]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos") )
#Coherence
barplot(t(df[,c(2,12,13)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,14,15)]), col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

##### FFL from plastic gene
Anticor10 <- FFL.type2(topo.anticor10, edges1 = 2, edges2 = 1, from = c(1))
Corr10 <- FFL.type2(topo.corr10, edges1 = 2, edges2 = 1, from = c(1))
No_sel10 <- FFL.type2(topo.no_sel10, edges1 = 2, edges2 = 1, from = c(1)) #
Sel10 <- FFL.type2(topo.sel10, edges1 = 2, edges2 = 1, from = c(1))

df <- as.data.frame(rbind(colSums(rbind(Anticor10,Corr10))*100/sum(rbind(Anticor10,Corr10)[,1]),
                          colSums(Anticor10)*100/sum(Anticor10[,1]),
                          colSums(Corr10)*100/sum(Corr10[,1]),
                          colSums(No_sel10)*100/sum(No_sel10[,1]),
                          colSums(Sel10)*100/sum(Sel10[,1])))
rownames(df) <- c("All_plast","Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))
write.csv(df, paste0("scripts/data/",filename,"_FFL_from",".csv"))


pdf(paste0("figures/",filename,"_FFL_from",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:11]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c("C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(11,12,13)]), col=c("lightslategrey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("lightslategrey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(11,14,15)]), col=c("lightslategrey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("lightslategrey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

##### FFL from NP genes
Anticor10 <- FFL.type2(topo.anticor10, edges1 = 2, edges2 = 1, from = c(3:10))
Corr10 <- FFL.type2(topo.corr10, edges1 = 2, edges2 = 1, from = c(3:10))
No_sel10 <- FFL.type2(topo.no_sel10, edges1 = 2, edges2 = 1, from = c(3:10)) #
Sel10 <- FFL.type2(topo.sel10, edges1 = 2, edges2 = 1, from = c(3:10))

df <- as.data.frame(rbind(colSums(rbind(Anticor10,Corr10))*100/sum(rbind(Anticor10,Corr10)[,1]),
                          colSums(Anticor10)*100/sum(Anticor10[,1]),
                          colSums(Corr10)*100/sum(Corr10[,1]),
                          colSums(No_sel10)*100/sum(No_sel10[,1]),
                          colSums(Sel10)*100/sum(Sel10[,1])))
rownames(df) <- c("All_plast","Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))
write.csv(df, paste0("scripts/data/",filename,"_FFL_from_NP",".csv"))

pdf(paste0("figures/",filename,"_FFL_from_NP",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:11]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c("C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(11,12,13)]), col=c("lightslategrey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("lightslategrey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(11,14,15)]), col=c("lightslategrey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("lightslategrey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

############Diamond
Anticor10 <- FFL.type2(topo.anticor10, edges1 = 2, edges2 = 2)
Corr10 <- FFL.type2(topo.corr10, edges1 = 2, edges2 = 2)
No_sel10 <- FFL.type2(topo.no_sel10, edges1 = 2, edges2 = 2) #
Sel10 <- FFL.type2(topo.sel10, edges1 = 2, edges2 = 2)

df <- as.data.frame(rbind(colSums(rbind(Anticor10,Corr10))/2*100/nrow(Anticor10), colSums(Anticor10)*100/nrow(Anticor10), colSums(Corr10)*100/nrow(Anticor10), colSums(No_sel10)*100/nrow(Anticor10), colSums(Sel10)*100/nrow(Sel10)))
rownames(df) <- c("All_plast","Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3,5,6,8,10,11)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(4,7,9,12)]))
df$homogenous_pos <- rowSums2(as.matrix(df[,c(3,6,9)]))
df$homogenous_neg <- rowSums2(as.matrix(df[,c(5,8,12)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(4,7,10,11)]))
write.csv(df, paste0("scripts/data/",filename,"_diamond",".csv"))

pdf(paste0("figures/", filename,"_diamond",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12]), col=c("grey","olivedrab1","palegreen","mediumseagreen","plum","darkorchid1","plum1","lightsalmon1","indianred1","lightpink","chocolate1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","olivedrab1","palegreen","mediumseagreen","plum","darkorchid1","plum1","lightsalmon1","indianred1","lightpink","chocolate1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,14,15)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,16:18)]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()

################################################################################
################
edges1 <- 2
edges2 <- 1
##### Loop count
Anticor10 <- loops_n.count(topo.anticor10, edges1 = 2, edges2 = 1)
Corr10 <- loops_n.count(topo.corr10, edges1 = 2, edges2 = 1)
No_sel10 <- loops_n.count(topo.no_sel10, edges1 = 2, edges2 = 1) #
Sel10 <- loops_n.count(topo.sel10, edges1 = 2, edges2 = 1)

pdf(paste0("figures/",filename,"_N_loops",".pdf"), width=4, height=6)
par(mar = c(5,4, 1,1))
layout(matrix(c(1:2), 2, 1, byrow = TRUE))
hist(rbind(Anticor10,Corr10)$Loop_number, main = "Plastic genes", freq = FALSE, xlim = c(0,30), ylim=c(0,1), breaks = c(0:max(rbind(Anticor10,Corr10)$Loop_number)), xlab = "Number of FFL per gene")
hist(Sel10$Loop_number, main = "Stable genes", freq = FALSE, ylim=c(0,1), xlim = c(0,30), breaks = c(0:max(rbind(Anticor10,Corr10)$Loop_number)), xlab = "Number of FFL per gene")
dev.off()
################################################################################
#FBL
################
edges <- 2
#####
Anticor10 <- FBL.type(topo.anticor10, edges = edges)
Corr10 <- FBL.type(topo.corr10, edges = edges)
No_sel10 <- FBL.type(topo.no_sel10, edges = edges) #
Sel10 <- FBL.type(topo.sel10, edges = edges)

df <- as.data.frame(rbind(colSums(rbind(Anticor10,Corr10))/2*100/nrow(Anticor10), colSums(Anticor10)*100/nrow(Anticor10), colSums(Corr10)*100/nrow(Anticor10), colSums(No_sel10)*100/nrow(Anticor10), colSums(Sel10)*100/nrow(Sel10)))
rownames(df) <- c("All_plast","Anticor","Cor","No_sel","Sel")
write.csv(df, paste0("scripts/data/",filename,"_FBL",".csv"))

pdf(paste0("figures/", filename,"_FBL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:4]), col=c("grey","olivedrab1","palegreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","olivedrab1","palegreen"),
       legend= colnames(df[,2:4]) )
dev.off()