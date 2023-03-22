source("scripts/functions/functions.R")
#####################
sims.dirs1 <- list.dirs("simul/3_random_RN", recursive = TRUE)
genes <- 3
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/3g_core_topo.connect.pdf"
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

anticor_down3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Anticorrelated_Down"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_AD3 <- core_topo.alt(anticor_down3)

anticor_up3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Anticorrelated_Up"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_AU3 <- core_topo.alt(anticor_up3)

anticor_UD3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Anticorrelated_UD"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_AUD3 <- core_topo.alt(anticor_UD3)

#
noncor_down3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Noncorrelated_UD" & P_mean_2 <= 0.35),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_ND3 <- core_topo.alt(noncor_down3)

noncor_up3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Noncorrelated_UD"& P_mean_2 >= 0.6),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_NU3 <- core_topo.alt(noncor_up3)

noncor_UD3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Noncorrelated_UD"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_NUD3 <- core_topo.alt(noncor_UD3)

#
corr_down3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Correlated_Down"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_CD3 <- core_topo.alt(corr_down3 )

corr_up3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Correlated_Up"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_CU3 <- core_topo.alt(corr_up3)

corr_UD3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Correlated_UD"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))
G_CUD3 <- core_topo.alt(corr_UD3)


pdf(pdfname, width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
  # #
  # plot(G_AD3, layout=layout_in_circle, edge.color=ifelse(E(G_AD3)$weight > 0, "black","red"), vertex.size=30, main="AD", vertex.color=c("green","orange", "yellow"))
  # plot(G_AU3, layout=layout_in_circle, edge.color=ifelse(E(G_AU3)$weight > 0, "black","red"), vertex.size=30, main="AU", vertex.color=c("green","orange", "yellow"))
  # plot(G_AUD3, layout=layout_in_circle, edge.color=ifelse(E(G_AUD3)$weight > 0, "black","red"), vertex.size=30, main="AUD", vertex.color=c("green","orange", "yellow"))
  # #
  # plot(G_ND3, layout=layout_in_circle, edge.color=ifelse(E(G_ND3)$weight > 0, "black","red"), vertex.size=30, main="ND", vertex.color=c("green","orange", "yellow"))
  # plot(G_NU3, layout=layout_in_circle, edge.color=ifelse(E(G_NU3)$weight > 0, "black","red"), vertex.size=30, main="NU", vertex.color=c("green","orange", "yellow"))
  # plot(G_NUD3, layout=layout_in_circle, edge.color=ifelse(E(G_NUD3)$weight > 0, "black","red"), vertex.size=30, main="NUD", vertex.color=c("green","orange", "yellow"))
  # #
  # plot(G_CD3, layout=layout_in_circle, edge.color=ifelse(E(G_CD3)$weight > 0, "black","red"), vertex.size=30, main="CD", vertex.color=c("green","orange", "yellow"))
  # plot(G_CU3, layout=layout_in_circle, edge.color=ifelse(E(G_CU3)$weight > 0, "black","red"), vertex.size=30, main="CU", vertex.color=c("green","orange", "yellow"))
  # plot(G_CUD3, layout=layout_in_circle, edge.color=ifelse(E(G_CUD3)$weight > 0, "black","red"), vertex.size=30, main="CUD", vertex.color=c("green","orange", "yellow"))
  # #
  #
  plot(G_AD3, layout=layout_in_circle, edge.color=ifelse(E(G_AD3)$weight > 0.75, "black",ifelse(E(G_AD3)$weight < -0.75, "red","grey")), vertex.size=30, main="AD", vertex.color=c("green","orange", "yellow"))
  plot(G_AU3, layout=layout_in_circle, edge.color=ifelse(E(G_AU3)$weight > 0.75, "black",ifelse(E(G_AU3)$weight < -0.75, "red","grey")), vertex.size=30, main="AU", vertex.color=c("green","orange", "yellow"))
  plot(G_AUD3, layout=layout_in_circle, edge.color=ifelse(E(G_AUD3)$weight > 0.75, "black",ifelse(E(G_AUD3)$weight < -0.75, "red","grey")), vertex.size=30, main="AUD", vertex.color=c("green","orange", "yellow"))
  #
  plot(G_ND3, layout=layout_in_circle, edge.color=ifelse(E(G_ND3)$weight > 0.75, "black",ifelse(E(G_ND3)$weight < -0.75, "red","grey")), vertex.size=30, main="ND", vertex.color=c("green","orange", "yellow"))
  plot(G_NU3, layout=layout_in_circle, edge.color=ifelse(E(G_NU3)$weight > 0.75, "black",ifelse(E(G_NU3)$weight < -0.75, "red","grey")), vertex.size=30, main="NU", vertex.color=c("green","orange", "yellow"))
  plot(G_NUD3, layout=layout_in_circle, edge.color=ifelse(E(G_NUD3)$weight > 0.75, "black",ifelse(E(G_NUD3)$weight < -0.75, "red","grey")), vertex.size=30, main="NUD", vertex.color=c("green","orange", "yellow"))
  #
  plot(G_CD3, layout=layout_in_circle, edge.color=ifelse(E(G_CD3)$weight > 0.75, "black",ifelse(E(G_CD3)$weight < -0.75, "red","grey")), vertex.size=30, main="CD", vertex.color=c("green","orange", "yellow"))
  plot(G_CU3, layout=layout_in_circle, edge.color=ifelse(E(G_CU3)$weight > 0.75, "black",ifelse(E(G_CU3)$weight < -0.75, "red","grey")), vertex.size=30, main="CU", vertex.color=c("green","orange", "yellow"))
  plot(G_CUD3, layout=layout_in_circle, edge.color=ifelse(E(G_CUD3)$weight > 0.75, "black",ifelse(E(G_CUD3)$weight < -0.75, "red","grey")), vertex.size=30, main="CUD", vertex.color=c("green","orange", "yellow"))
  #
dev.off()



pdf(pdfname, width=6, height=3)
layout(matrix(c(1:2), 1, 2, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
  #
  plot(G_AUD3, layout=layout_in_circle, edge.color=ifelse(E(G_AUD3)$weight > .8, "black",ifelse(E(G_AUD3)$weight < -.8, "red","grey")), vertex.size=30, main="AUD", vertex.color=c("green","orange", "yellow"))
  #
  plot(G_CUD3, layout=layout_in_circle, edge.color=ifelse(E(G_CUD3)$weight > .8, "black",ifelse(E(G_CUD3)$weight < -.8, "red","grey")), vertex.size=30, main="CUD", vertex.color=c("green","orange", "yellow"))
  # 
dev.off()


pdf("figures/3g_AU_connect.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 0.5, 0), mgp = c(1.75, 0.75, 0), las=0)
for (i in anticor_down3) {
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30) #, layout=layout_in_circle
}
dev.off()




