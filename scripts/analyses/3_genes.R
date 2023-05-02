source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/3g", recursive = TRUE)
genes <- 3
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/3g.pdf"
#####################
#DATA
##################
df.3 <- df.simul(sims.dirs1, all.gen = TRUE)
df.3$envir <- str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,3]
df.3$anc_id <- str_split(str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]


#Add Slope
for (i in 1:nrow(df.3)) {
  W <- t(matrix(as.numeric(df.3[i,7:(genes^2+6)]), ncol = genes))
  df.3[i,(genes^2+9)] <- getSlope.ALR(W=W, n.env=21, target.gene=target, min=min, max=max, a=0.5)
}


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

core.sel3 <- core_topo.alt(topo.sel3)
core.anticor3 <- core_topo.alt(topo.anticor3)
core.corr3 <- core_topo.alt(topo.corr3)
core.no_sel3 <- core_topo.alt(topo.no_sel3)

pdf(pdfname, width=6, height=6)
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
  par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
  #
  plot(core.anticor3, layout=layout_in_circle, edge.color=ifelse(E(core.anticor3)$weight > .8, "black",ifelse(E(core.anticor3)$weight < -.8, "red","grey")), vertex.size=30, main="Anticor", vertex.color=c("green","orange", "yellow"))
  #
  plot(core.corr3, layout=layout_in_circle, edge.color=ifelse(E(core.corr3)$weight > .8, "black",ifelse(E(core.corr3)$weight < -.8, "red","grey")), vertex.size=30, main="Corr", vertex.color=c("green","orange", "yellow"))
  # 
  plot(core.no_sel3, layout=layout_in_circle, edge.color=ifelse(E(core.no_sel3)$weight > .8, "black",ifelse(E(core.no_sel3)$weight < -.8, "red","grey")), vertex.size=30, main="No_sel", vertex.color=c("green","orange", "yellow"))
  #
  plot(core.sel3, layout=layout_in_circle, edge.color=ifelse(E(core.sel3)$weight > .8, "black",ifelse(E(core.sel3)$weight < -.8, "red","grey")), vertex.size=30, main="Sel", vertex.color=c("green","orange", "yellow"))
  # 
dev.off()


pdf("figures/3g_anticorr.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(topo.anticor3)) {
  mainT <- length(which(sapply(1:length(topo.anticor3),function(x) length(which(paste0(topo.anticor3[[x]],
            collapse = "") == paste0(unique(topo.anticor3)[[j]], collapse = "")))) == 1))/length(topo.anticor3)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()
pdf("figures/3g_cor.pdf", width=6, height=6)
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


#####

Anticor3 <- c.count(topo.anticor3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)
Corr3 <- c.count(topo.corr3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)
No_sel3 <- c.count(topo.no_sel3, cutoff.max = 3, cutoff.min = 1, target = 2, randomFF=FALSE) #
Sel3 <- c.count(topo.sel3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c("Anticor3", "Corr3", "No_sel3", "Sel3")

layout(matrix(c(1:1), 1, 1, byrow = TRUE))

barplot(t(df[,1:3])*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFl", "Incoherent FFL","No FFL"))

#####

#####

Anticor3 <- homog.count(topo.anticor3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)
Corr3 <- homog.count(topo.corr3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)
No_sel3 <- homog.count(topo.no_sel3, cutoff.max = 3, cutoff.min = 1, target = 2, randomFF=FALSE) #
Sel3 <- homog.count(topo.sel3, cutoff.max = 3, cutoff.min = 1, randomFF=FALSE)

df <- rbind(colSums(Anticor3), colSums(Corr3), colSums(No_sel3), colSums(Sel3))
rownames(df) <- c("Anticor3", "Corr3", "No_sel3", "Sel3")

layout(matrix(c(1:1), 1, 1, byrow = TRUE))

barplot(t(df[,1:3])*100/300, col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFl", "Incoherent FFL","No FFL"))

#####













