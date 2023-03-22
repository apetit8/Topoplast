source("scripts/functions/functions.R")
#####################
sims.dirs1 <- list.dirs("simul/4_random_RN", recursive = TRUE)
genes <- 4
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/4g_core_topo.random.pdf"
#####################
#DATA
##################
df.4 <- df.simul(sims.dirs1, all.gen = TRUE)
df.4$envir <- str_split(df.4$data.dir, "/", n=8, simplify = TRUE)[,3]
df.4$anc_id <- str_split(str_split(df.4$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]

#Add Slope
for (i in 1:nrow(df.4)) {
  W <- t(matrix(as.numeric(df.4[i,7:(genes^2+6)]), ncol = genes))
  #First version :
  #basal <- if(grepl("Down", df.4[i,23])) 0.8 else if(grepl("Up", df.4[i,23])) 0.2 else 0.5
  #Alt :
  # max <- if(grepl("Down", df.4[i,23])) 0.5 else if(grepl("Up", df.4[i,23])) 0.8 else 0.8
  # min <- if(grepl("Down", df.4[i,23])) 0.2 else if(grepl("Up", df.4[i,23])) 0.5 else 0.2
  #
  df.4[i,(genes^2+9)] <- getSlope.ALR(W=W, n.env=21, target.gene=target, min=min, max=max, a=0.5)
}

ggplot(df.3, aes(Gen, W_1_2))+ geom_point()

################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
#############

#Need : a way to sort genes in essential.topo, for genes inhibiting target to always be in position X and genes activating to be in Y
#Necessary for the core_topo.alt algo to function correctly !
#Can be put in core_topo.alt ! in a lapply

anticor_down4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Anticorrelated_Down" & Fitness >= 0.98),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_AD4 <- core_topo.alt(anticor_down4, sorting_4g = TRUE)

anticor_up4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Anticorrelated_Up" & Fitness >= 0.98),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_AU4 <- core_topo.alt(anticor_up4, sorting_4g = TRUE)

anticor_UD4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Anticorrelated_UD" & Fitness >= 0.98),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_AUD4 <- core_topo.alt(anticor_UD4, sorting_4g = TRUE)

#
noncor_down4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Noncorrelated_UD" & Fitness >= 0.98 & P_mean_2 <= 0.35),
                               treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_ND4 <- core_topo.alt(noncor_down4, sorting_4g = TRUE)

noncor_up4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Noncorrelated_UD" & Fitness >= 0.98 & P_mean_2 >= 0.65),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_NU4 <- core_topo.alt(noncor_up4, sorting_4g = TRUE)

noncor_UD4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Noncorrelated_UD" & Fitness >= 0.98),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_NUD4 <- core_topo.alt(noncor_UD4, sorting_4g = TRUE)

#
corr_down4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Correlated_Down" & Fitness >= 0.98),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_CD4 <- core_topo.alt(corr_down4, sorting_4g = TRUE )

corr_up4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Correlated_Up" & Fitness >= 0.98),
                           treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_CU4 <- core_topo.alt(corr_up4, sorting_4g = TRUE)

corr_UD4 <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Correlated_UD" & Fitness >= 0.98),
                           treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:4), min=min, max=max)
G_CUD4 <- core_topo.alt(corr_UD4, sorting_4g = TRUE)




pdf(pdfname, width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
grid <- matrix(c( 1,1,4,4,1,4,4,1),ncol=2,byrow=TRUE)
  #
  plot(G_AD4, layout=grid, edge.color=ifelse(E(G_AD4)$weight > .9, "black",ifelse(E(G_AD4)$weight < -.9, "red","grey")), vertex.size=30, main="AD", vertex.color=c("green","orange", "yellow","yellow"))
  plot(G_AU4, layout=grid, edge.color=ifelse(E(G_AU4)$weight > .9, "black",ifelse(E(G_AU4)$weight < -.9, "red","grey")), vertex.size=30, main="AU", vertex.color=c("green","orange", "yellow", "yellow"))
  plot(G_AUD4, layout=grid, edge.color=ifelse(E(G_AUD4)$weight > .9, "black",ifelse(E(G_AUD4)$weight < -.9, "red","grey")), vertex.size=30, main="AUD", vertex.color=c("green","orange", "yellow", "yellow"))
  #
  plot(G_ND4, layout=grid, edge.color=ifelse(E(G_ND4)$weight > .9, "black",ifelse(E(G_ND4)$weight < -.9, "red","grey")), vertex.size=30, main="ND", vertex.color=c("green","orange", "yellow", "yellow"))
  plot(G_NU4, layout=grid, edge.color=ifelse(E(G_NU4)$weight > .9, "black",ifelse(E(G_NU4)$weight < -.9, "red","grey")), vertex.size=30, main="NU", vertex.color=c("green","orange", "yellow", "yellow"))
  plot(G_NUD4, layout=grid, edge.color=ifelse(E(G_NUD4)$weight > .9, "black",ifelse(E(G_NUD4)$weight < -.9, "red","grey")), vertex.size=30, main="NUD", vertex.color=c("green","orange", "yellow", "yellow"))
  # 
  plot(G_CD4, layout=grid, edge.color=ifelse(E(G_CD4)$weight > .9, "black",ifelse(E(G_CD4)$weight < -.9, "red","grey")), vertex.size=30, main="CD", vertex.color=c("green","orange", "yellow", "yellow"))
  plot(G_CU4, layout=grid, edge.color=ifelse(E(G_CU4)$weight > .9, "black",ifelse(E(G_CU4)$weight < -.9, "red","grey")), vertex.size=30, main="CU", vertex.color=c("green","orange", "yellow", "yellow"))
  plot(G_CUD4, layout=grid, edge.color=ifelse(E(G_CUD4)$weight > .9, "black",ifelse(E(G_CUD4)$weight < -.9, "red","grey")), vertex.size=30, main="CUD", vertex.color=c("green","orange", "yellow", "yellow"))
  # 
dev.off()


pdf(pdfname, width=6, height=3)
layout(matrix(c(1:2), 1, 2, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
grid <- matrix(c( 1,1,4,4,1,4,4,1),ncol=2,byrow=TRUE)
  #
  plot(G_AUD4, layout=grid, edge.color=ifelse(E(G_AUD4)$weight > .8, "black",ifelse(E(G_AUD4)$weight < -.8, "red","grey")), vertex.size=30, main="AUD", vertex.color=c("green","orange", "yellow", "yellow"))
  #
  plot(G_CUD4, layout=grid, edge.color=ifelse(E(G_CUD4)$weight > .8, "black",ifelse(E(G_CUD4)$weight < -.8, "red","grey")), vertex.size=30, main="CUD", vertex.color=c("green","orange", "yellow", "yellow"))
  # 
dev.off()


# pdf("../figures/4g_noncorr.pdf", width=6, height=6)
# layout(matrix(c(1:9), 3, 3, byrow = TRUE))
# par(mar=c(1.5, 1.5, 1.5, 0), mgp = c(1.75, 0.75, 0), las=0)
# for (i in unique(noncor)) {
#   #W matrix as a graph : 
#   G <- as.directed(graph.adjacency(t(i), weighted = T))
#   V(G)$color <- c("green","orange", "yellow", "yellow")
#   plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30)
# }
# dev.off()
# 
pdf("figures/4g_CU.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 0.5, 0), mgp = c(1.75, 0.75, 0), las=0)
for (i in unique(corr_up4)) {
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=grid, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30) #, layout=layout_in_circle
}
dev.off()
# 
# pdf("../figures/4g_corr.pdf", width=6, height=6)
# layout(matrix(c(1:9), 3, 3, byrow = TRUE))
# par(mar=c(0.5, 0.5, 0.5, 0), mgp = c(1.75, 0.75, 0), las=0)
# for (i in unique(corr)) {
#   #W matrix as a graph : 
#   G <- as.directed(graph.adjacency(t(i), weighted = T))
#   V(G)$color <- c("green","orange", "yellow", "yellow")
#   plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30)
# }
# dev.off()




#Would it be possible to extract the common regulations between all networks ?


G <- as.directed(graph.adjacency(t(core_topo.alt(corr)), weighted = T))
V(G)$color <- c("green","orange", "yellow", "yellow")
plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red"), vertex.size=30)


G <- as.directed(graph.adjacency(t(core_topo.alt(anticor)), weighted = T))
V(G)$color <- c("green","orange", "yellow", "yellow")
plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red"), vertex.size=30)



#
noncor <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Noncorrelated" & P_mean_2 > 0.6),
                         treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)

G <- as.directed(graph.adjacency(t(core_topo.alt(noncor)), weighted = T))
V(G)$color <- c("green","orange", "yellow", "yellow")
plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red"), vertex.size=30)


noncor <- essential.topo(df=subset(df.4, Gen==max(df.4$Gen) & envir=="Noncorrelated" & P_mean_2 < 0.2),
                         treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)

mean_W(noncor)

G <- as.directed(graph.adjacency(t(core_topo.alt(noncor, treshold=0.95)), weighted = T))
V(G)$color <- c("green","orange", "yellow", "yellow")
plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red"), vertex.size=30)
