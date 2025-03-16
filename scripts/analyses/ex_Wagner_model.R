source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(igraph)

pdf(paste0("Def_Wagner.pdf"), width=6, height=3.8)
par(mgp=c(2.5, 1.2, 0))

tt <- phenotype_plastic_loop_R(matrix(data=c(0,1,0,0), nrow=2), a=0.2, S0=c(0.2,0.2) , steps = 24, sensors = NULL, measure=4)$full
plot(c(0:24), tt[1,], ylim = c(0,1), pch=19, col="sandybrown")
points(c(0:24), tt[2,], pch=19, col="lightskyblue")
lines(c(0:24), tt[2,], pch=19, col="lightskyblue")
lines(c(0:24), tt[1,], pch=19, col="sandybrown")


tt <- phenotype_plastic_loop_R(matrix(data=c(0,-1,0,0), nrow=2), a=0.2, S0=c(0.2,0.2) , steps = 24, sensors = NULL, measure=4)$full
plot(c(0:24), tt[1,], ylim = c(0,1), pch=19, col="sandybrown")
points(c(0:24), tt[2,], pch=19, col="lightskyblue")
lines(c(0:24), tt[2,], pch=19, col="lightskyblue")
lines(c(0:24), tt[1,], pch=19, col="sandybrown")


tt <- pheno.from.W(matrix(data=c(0,-1,0,0,0,1,0,0,0), nrow=3), a=0.2, S0=c(0,0.2,0.2) , steps = 24, sensors = 0, measure=4, full=TRUE)$full
plot(c(0:24), tt[2,], ylim = c(0,1), pch=19, col="sandybrown")
points(c(0:24), tt[3,], pch=19, col="lightskyblue")
lines(c(0:24), tt[3,], pch=19, col="lightskyblue")
lines(c(0:24), tt[2,], pch=19, col="sandybrown")
points(c(0:24), tt[1,], pch=19, col="olivedrab1")
lines(c(0:24), tt[1,], pch=19, col="olivedrab1")


tt <- phenotype_plastic_loop_R(matrix(data=c(0,-1,0,0,0,1,0,0,0), nrow=3), a=0.2, S0=c(0.8,0.2,0.2) , steps = 24, sensors = 0.8, measure=4)$full
plot(c(0:24), tt[2,], ylim = c(0,1), pch=19, col="sandybrown")
points(c(0:24), tt[3,], pch=19, col="lightskyblue")
lines(c(0:24), tt[3,], pch=19, col="lightskyblue")
lines(c(0:24), tt[2,], pch=19, col="sandybrown")
points(c(0:24), tt[1,], pch=19, col="olivedrab1")
lines(c(0:24), tt[1,], pch=19, col="olivedrab1")

dev.off()

data <- rnorm(100)
mat <- matrix(data=data, nrow=10)
diag(mat) <- 0

pdf(paste0("Def_Wagner_complex.pdf"), width=6, height=3.8)
par(mgp=c(2.5, 1.2, 0))
tt <- phenotype_plastic_loop_R(mat, a=0.2, S0=rep(0.2, 10) , steps = 24, measure=4, sensors = NULL)$full
plot(c(0:24), tt[1,], ylim = c(0,1), pch=19, col="sandybrown")
lines(c(0:24), tt[1,], pch=19, col="sandybrown")
points(c(0:24), tt[4,], pch=19, col="blue")
lines(c(0:24), tt[4,], pch=19, col="blue")
points(c(0:24), tt[5,], pch=19, col="red")
lines(c(0:24), tt[5,], pch=19, col="red")
points(c(0:24), tt[6,], pch=19, col="purple")
lines(c(0:24), tt[6,], pch=19, col="purple")
points(c(0:24), tt[7,], pch=19, col="darkgreen")
lines(c(0:24), tt[7,], pch=19, col="darkgreen")
points(c(0:24), tt[8,], pch=19, col="lightpink")
lines(c(0:24), tt[8,], pch=19, col="lightpink")
points(c(0:24), tt[9,], pch=19, col="maroon")
lines(c(0:24), tt[9,], pch=19, col="maroon")
points(c(0:24), tt[10,], pch=19, col="yellow")
lines(c(0:24), tt[10,], pch=19, col="yellow")
points(c(0:24), tt[2,], pch=19, col="lightskyblue")
lines(c(0:24), tt[2,], pch=19, col="lightskyblue")
points(c(0:24), tt[3,], pch=19, col="olivedrab1")
lines(c(0:24), tt[3,], pch=19, col="olivedrab1")

g1 <- as.directed(graph.adjacency(t(mat), weighted = T))
deg <- degree(g1, mode = "all")

plot(g1, layout=layout.circle, edge.color=ifelse(E(g1)$weight > 0.5, "black",ifelse(E(g1)$weight < -0.5, "red","white")),
     vertex.size=15, edge.curved=TRUE,
     vertex.color="grey")

dev.off()

round(mat,2)
