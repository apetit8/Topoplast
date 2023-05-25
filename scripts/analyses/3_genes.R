source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/3g", recursive = TRUE)
genes <- 3
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/3g"
#####################
#DATA
##################
df.3 <- df.simul(sims.dirs1, all.gen = TRUE)
df.3$envir <- str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,3]
df.3$anc_id <- str_split(str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]


################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
#############

topo.anticor3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Anticorrelated"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))

topo.corr3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Correlated"),
                           treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))

topo.no_sel3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Control_no_sel"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))

topo.sel3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Control_sel"),
                           treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))

pdf("figures/3g_sel.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(topo.sel3)) {
  mainT <- length(which(sapply(1:length(topo.sel3),function(x) length(which(paste0(topo.sel3[[x]],
                                                                                       collapse = "") == paste0(unique(topo.sel3)[[j]], collapse = "")))) == 1))/length(topo.sel3)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()

pdf("figures/3g_corr.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(topo.corr3)) {
  mainT <- length(which(sapply(1:length(topo.corr3),function(x) length(which(paste0(topo.corr3[[x]],
                                                                                   collapse = "") == paste0(unique(topo.corr3)[[j]], collapse = "")))) == 1))/length(topo.corr3)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()

################
cutoff.max <- 3
cutoff.min <- 1
##### Loop count
Anticor3_n <- loops_n.count(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
Corr3_n <- loops_n.count(topo.corr3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
No_sel3_n <- loops_n.count(topo.no_sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2) #
Sel3_n <- loops_n.count(topo.sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)

pdf(paste0(pdfname,"_N_FFL",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor3_n$Loop_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr3_n$Loop_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel3_n$Loop_number, main = "No Sel", ylim = c(0, 300))
hist(Sel3_n$Loop_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()

################
##### Coherence
Anticor3 <- FFL.coherence(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)
Corr3 <- FFL.coherence(topo.corr3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)
No_sel3 <- FFL.coherence(topo.no_sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2, randomFF=TRUE) #
Sel3 <- FFL.coherence(topo.sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c("Anticor3", "Corr3", "No_sel3", "Sel3")

pdf(paste0(pdfname,"_coherence_FFL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL","No Loop"))
dev.off()


################
##### Homogeneity
Anticor3 <- homog.count(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
Corr3 <- homog.count(topo.corr3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
No_sel3 <- homog.count(topo.no_sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2, randomFF=FALSE) #
Sel3 <- homog.count(topo.sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor3_n$Loop_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr3_n$Loop_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel3_n$Loop_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel3_n$Loop_number), 1)))

pdf(paste0(pdfname,"_homog_FFL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/300, col=c(2,4,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(2,4,"grey"),
       legend=c("Het FFL", "Hom FFL","No FFL"))
dev.off()

################
##### Type
Anticor3 <- FFL.type(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, frequencies=FALSE)
Corr3 <- FFL.type(topo.corr3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, frequencies=FALSE)
No_sel3 <- FFL.type(topo.no_sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, frequencies=FALSE) #
Sel3 <- FFL.type(topo.sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, frequencies=FALSE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))

pdf(paste0(pdfname,"_type_FFL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:5])*100/300, col=c("darkolivegreen","yellowgreen","darkblue","lightblue","grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkolivegreen","yellowgreen","darkblue","lightblue","grey"),
       legend=c( "Activating","Inhibiting","Z_act","Z_inh","No FFL"))
dev.off()


################
##### Type
Anticor3 <- FFL.type2(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
Corr3 <- FFL.type2(topo.corr3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
No_sel3 <- FFL.type2(topo.no_sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min) #
Sel3 <- FFL.type2(topo.sel3, cutoff.max = cutoff.max, cutoff.min = cutoff.min)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c("Anticor","Cor","No_sel","Sel")

pdfname <- "figures/3g"
pdf(paste0(pdfname,"_type_FFL2_sign",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,2:10])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","Input_Dep_Amplifying_neg","Input_Dep_Amplifying_pos","Input_Dep_Disruptive_neg","Input_Dep_Disruptive_pos","Input_Ind_Amplifying_neg","Input_Ind_Amplifying_pos","Input_Ind_Disruptive_neg","Input_Ind_Disruptive_pos") )
dev.off()


################################################################################
#FEEDBACKS

#Renforcing = same sign
#Defusing = opposite sign

################
cutoff <- 3
##### Loop count
Anticor3_n <- FBL_n.count(topo.anticor3, cutoff = cutoff, target = 2)
Corr3_n <- FBL_n.count(topo.corr3, cutoff = cutoff, target = 2)
No_sel3_n <- FBL_n.count(topo.no_sel3, cutoff = cutoff, target = 2) #
Sel3_n <- FBL_n.count(topo.sel3, cutoff = cutoff, target = 2)

pdf(paste0(pdfname,"_N_FBL",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor3_n$Loop_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr3_n$Loop_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel3_n$Loop_number, main = "No Sel", ylim = c(0, 300))
hist(Sel3_n$Loop_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()


################
##### Type
Anticor3 <- FBL.type(topo.anticor3, cutoff = cutoff, target = 2, randomFF=FALSE)
Corr3 <- FBL.type(topo.corr3, cutoff = cutoff, target = 2, randomFF=FALSE)
No_sel3 <- FBL.type(topo.no_sel3, cutoff = cutoff, target = 2, randomFF=FALSE) #
Sel3 <- FBL.type(topo.sel3, cutoff = cutoff, target = 2, randomFF=FALSE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor3_n$FBL_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr3_n$FBL_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel3_n$FBL_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel3_n$FBL_number), 1)))

pdf(paste0(pdfname,"_type_FBL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/300, col=c("sienna1","orchid3","grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("sienna1","orchid3","grey"),
       legend=c("Defusing", "Renforcing FFL","No FFL"))
dev.off()



