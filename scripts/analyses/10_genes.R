source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/10g", recursive = TRUE)
genes <- 10
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/10g.pdf"
#####################
#DATA
##################
df.10 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10$envir <- str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,3]
df.10$anc_id <- str_split(str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]


#Add Slope
for (i in 1:nrow(df.10)) {
  W <- t(matrix(as.numeric(df.10[i,7:(genes^2+6)]), ncol = genes))
  df.10[i,(genes^2+9)] <- getSlope.ALR(W=W, n.env=21, target.gene=target, min=min, max=max, a=0.5)
}


################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
#############

topo.anticor10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Anticorrelated"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Correlated"),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.no_sel10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Control_no_sel"),
                               treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Control_sel"),
                            treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

core.sel10 <- core_topo.alt(topo.sel10)
core.anticor10 <- core_topo.alt(topo.anticor10)
core.corr10 <- core_topo.alt(topo.corr10)
core.no_sel10 <- core_topo.alt(topo.no_sel10)

pdf(pdfname, width=6, height=6)
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
#
plot(core.anticor10, layout=layout_in_circle, edge.color=ifelse(E(core.anticor10)$weight > .8, "black",ifelse(E(core.anticor10)$weight < -.8, "red","grey")), vertex.size=30, main="Anticor", vertex.color=c("green","orange", "yellow"))
#
plot(core.corr10, layout=layout_in_circle, edge.color=ifelse(E(core.corr10)$weight > .8, "black",ifelse(E(core.corr10)$weight < -.8, "red","grey")), vertex.size=30, main="Corr", vertex.color=c("green","orange", "yellow"))
# 
plot(core.no_sel10, layout=layout_in_circle, edge.color=ifelse(E(core.no_sel10)$weight > .8, "black",ifelse(E(core.no_sel10)$weight < -.8, "red","grey")), vertex.size=30, main="No_sel", vertex.color=c("green","orange", "yellow"))
#
plot(core.sel10, layout=layout_in_circle, edge.color=ifelse(E(core.sel10)$weight > .8, "black",ifelse(E(core.sel10)$weight < -.8, "red","grey")), vertex.size=30, main="Sel", vertex.color=c("green","orange", "yellow"))
# 
dev.off()


pdf("figures/10g_anticorr.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(topo.anticor10)) {
  mainT <- length(which(sapply(1:length(topo.anticor10),function(x) length(which(paste0(topo.anticor10[[x]],
                                                                                       collapse = "") == paste0(unique(topo.anticor10)[[j]], collapse = "")))) == 1))/length(topo.anticor10)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()
pdf("figures/10g_cor.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(topo.corr10)) {
  mainT <- length(which(sapply(1:length(topo.corr10),function(x) length(which(paste0(topo.corr10[[x]],
                                                                                    collapse = "") == paste0(unique(topo.corr10)[[j]], collapse = "")))) == 1))/length(topo.corr10)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()


###

Anticor10 <- coherence_count(topo.anticor10, cutoff.max = 9, cutoff.min = 1)
if(is.na(Anticor10["No_FF"])) Anticor10["No_FF"] <- 0
if(is.na(Anticor10["FF_Coherent"])) Anticor10["FF_Coherent"] <- 0
if(is.na(Anticor10["FF_Incoherent"])) Anticor10["FF_Incoherent"] <- 0
Corr10 <- coherence_count(topo.corr10, cutoff.max = 9, cutoff.min = 1)
if(is.na(corr10["No_FF"])) corr10["No_FF"] <- 0
if(is.na(corr10["FF_Coherent"])) corr10["FF_Coherent"] <- 0
if(is.na(corr10["FF_Incoherent"])) corr10["FF_Incoherent"] <- 0
#Enrichment in incoherent FFL when selected for plasticity !
No_sel10 <- coherence_count(topo.no_sel10, cutoff.max = 9, cutoff.min = 1)
if(is.na(No_sel10["No_FF"])) No_sel10["No_FF"] <- 0
if(is.na(No_sel10["FF_Coherent"])) No_sel10["FF_Coherent"] <- 0
if(is.na(No_sel10["FF_Incoherent"])) No_sel10["FF_Incoherent"] <- 0
Sel10 <- coherence_count(topo.sel10, cutoff.max = 9, cutoff.min = 1)
if(is.na(sel10["No_FF"])) sel10["No_FF"] <- 0
if(is.na(sel10["FF_Coherent"])) sel10["FF_Coherent"] <- 0
if(is.na(sel10["FF_Incoherent"])) sel10["FF_Incoherent"] <- 0


l <- list("Anticor10"=Anticor10, "Corr10"=Corr10, "No_sel10"=No_sel10, "Sel10"=Sel10)
df.coherence <- do.call(rbind, lapply(l, function(row) row[order(names(row))]))



pdf("figures/FFL_distrib_10g_cut_off.pdf", width=4, height=4)
layout(matrix(c(1), 1, 1, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
barplot(t(df.coherence)*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL", "No FFL"))
dev.off()





Anticor10 <- coherence_count(topo.anticor10, cutoff.max = 3, cutoff.min = 2)
if(is.na(Anticor10["No_FF"])) Anticor10["No_FF"] <- 0
if(is.na(Anticor10["FF_Coherent"])) Anticor10["FF_Coherent"] <- 0
if(is.na(Anticor10["FF_Incoherent"])) Anticor10["FF_Incoherent"] <- 0
Corr10 <- coherence_count(topo.corr10, cutoff.max = 3, cutoff.min = 2)
if(is.na(corr10["No_FF"])) corr10["No_FF"] <- 0
if(is.na(corr10["FF_Coherent"])) corr10["FF_Coherent"] <- 0
if(is.na(corr10["FF_Incoherent"])) corr10["FF_Incoherent"] <- 0
#Enrichment in incoherent FFL when selected for plasticity !
No_sel10 <- coherence_count(topo.no_sel10, cutoff.max = 3, cutoff.min = 2)
if(is.na(No_sel10["No_FF"])) No_sel10["No_FF"] <- 0
if(is.na(No_sel10["FF_Coherent"])) No_sel10["FF_Coherent"] <- 0
if(is.na(No_sel10["FF_Incoherent"])) No_sel10["FF_Incoherent"] <- 0
Sel10 <- coherence_count(topo.sel10, cutoff.max = 3, cutoff.min = 2)
if(is.na(sel10["No_FF"])) sel10["No_FF"] <- 0
if(is.na(sel10["FF_Coherent"])) sel10["FF_Coherent"] <- 0
if(is.na(sel10["FF_Incoherent"])) sel10["FF_Incoherent"] <- 0


l <- list("Anticor10"=Anticor10, "Corr10"=Corr10, "No_sel10"=No_sel10, "Sel10"=Sel10)
df.coherence <- do.call(rbind, lapply(l, function(row) row[order(names(row))]))

barplot(t(df.coherence)*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL", "No FFL"))







topo.anticor10 <- essential.topo(df=subset(df.10, Gen==3000 & envir=="Anticorrelated"),
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10 <- essential.topo(df=subset(df.10, Gen==3000 & envir=="Correlated"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.no_sel10 <- essential.topo(df=subset(df.10, Gen==3000 & envir=="Control_no_sel"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10 <- essential.topo(df=subset(df.10, Gen==3000 & envir=="Control_sel"),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))



Anticor10 <- c.count(topo.anticor10, cutoff.max = 3, cutoff.min = 1)
Corr10 <- c.count(topo.corr10, cutoff.max = 3, cutoff.min = 1)
No_sel10 <- c.count(topo.no_sel10, cutoff.max = 3, cutoff.min = 1, target = 2) #
Sel10 <- c.count(topo.sel10, cutoff.max = 3, cutoff.min = 1)

df <- rbind(colSums(Anticor10), colSums(Corr10), colSums(No_sel10), colSums(Sel10))
rownames(df) <- c("Anticor10", "Corr18", "No_sel10", "Sel10")



layout(matrix(c(1:2), 1, 2, byrow = TRUE))

barplot(t(df[,3:4])*100/300, col=c("grey","yellow2"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","yellow2"),
       legend=c("No Loop", "Loop(s)"))

barplot(t(df[,1:2]), col=c(7,3))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3),
       legend=c("Coherent FFL", "Incoherent FFL"))


#####
Anticor10 <- c.count(topo.anticor10, cutoff.max = 3, cutoff.min = 1, randomFF=TRUE)
Corr10 <- c.count(topo.corr10, cutoff.max = 3, cutoff.min = 1, randomFF=TRUE)
No_sel10 <- c.count(topo.no_sel10, cutoff.max = 3, cutoff.min = 1, target = 2, randomFF=TRUE) #
Sel10 <- c.count(topo.sel10, cutoff.max = 3, cutoff.min = 1, randomFF=TRUE)

df <- rbind(colSums(Anticor10), colSums(Corr10), colSums(No_sel10), colSums(Sel10))
rownames(df) <- c("Anticor10", "Corr18", "No_sel10", "Sel10")

layout(matrix(c(1:1), 1, 1, byrow = TRUE))

barplot(t(df[,1:3])*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL","No Loop"))

