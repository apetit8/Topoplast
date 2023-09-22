#Script to : trace stability plot to decide at which generation to look 
#And then be happy about it
library(matrixStats)
ifelse(!dir.exists("figures/supp"), dir.create("figures/supp"), FALSE)
################################################################################
df.nbr.reg <- data.frame()
for (gen in seq(500, 20000, 500)) {
  ntw <- c(readRDS(paste0("scripts/data/Stab/list_plastic_topo_full_drift_",gen,".Rds")), readRDS(paste0("scripts/data/Stab/list_nnplast_topo_full_drift_",gen,".Rds")))
  df <- data.frame( Nbr_reg = mean(unlist(lapply(ntw, function(i){rowSums2(abs(i))[2]} ))), Gen = gen  )
  df$Gen <- gen
  df.nbr.reg <- rbind(df.nbr.reg, df)  
}


df.sign.reg <- data.frame()
for (gen in seq(500, 20000, 500)) {
  ntw <- c(readRDS(paste0("scripts/data/Stab/list_plastic_topo_full_drift_",gen,".Rds")), readRDS(paste0("scripts/data/Stab/list_nnplast_topo_full_drift_",gen,".Rds")))
  reg_pos <- unlist(lapply(ntw, function(i){ i[(i == -1)] <- 0
  return(rowSums2(abs(i))[2])} )) # [2] : the target of the network is in position #2
  reg_neg <- unlist(lapply(ntw, function(i){ i[(i == 1)] <- 0
  return(rowSums2(abs(i))[2])} ))
  df <- data.frame( Nbr_reg =  mean(reg_pos)/mean(reg_neg), Gen = gen  )
  df$Gen <- gen
  df.sign.reg <- rbind(df.sign.reg, df)  
}


pdf("figures/supp/Stab_drift_reg10000.pdf", width=10, height=4)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar = c(5,2, 2,1))
plot(df.nbr.reg$Gen, df.nbr.reg$Nbr_reg, ylab = "Number of regulations", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
plot(df.sign.reg$Gen, df.sign.reg$Nbr_reg, ylab = "Pos reg / Neg reg", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
dev.off()

################################################################################
df.gen.FFL <- data.frame()
for (gen in seq(500, 20000, 500)) {
  df <- read.csv(paste0("scripts/data/Stab/full_drift_",gen,"_FFL.csv"))[1,]
  df$Gen <- gen
  df.gen.FFL <- rbind(df.gen.FFL, df)  
}
plot(df.gen.FFL$Gen, df.gen.FFL$FFL)

df.gen.DMD <- data.frame()
for (gen in seq(500, 20000, 500)) {
  df <- read.csv(paste0("scripts/data/Stab/full_drift_",gen,"_DMD.csv"))[1,]
  df$Gen <- gen
  df.gen.DMD <- rbind(df.gen.DMD, df)  
}
plot(df.gen.DMD$Gen, df.gen.DMD$FFL)

df.gen.FFL2 <- data.frame()
for (gen in seq(500, 20000, 500)) {
  df <- read.csv(paste0("scripts/data/Stab/full_drift_",gen,"_FFL.csv"))[2,]
  df$Gen <- gen
  df.gen.FFL2 <- rbind(df.gen.FFL2, df)  
}
plot(df.gen.FFL$Gen, df.gen.FFL$FFL)

df.gen.DMD2 <- data.frame()
for (gen in seq(500, 20000, 500)) {
  df <- read.csv(paste0("scripts/data/Stab/full_drift_",gen,"_DMD.csv"))[2,]
  df$Gen <- gen
  df.gen.DMD2 <- rbind(df.gen.DMD2, df)  
}
plot(df.gen.DMD2$Gen, df.gen.DMD$FFL)

pdf("figures/supp/Stab_drift_loop10000.pdf", width=8, height=5.5)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar = c(5,5, 1,1))
plot(df.gen.FFL$Gen, df.gen.FFL$FFL, ylab = "Plastic genes : Number of FFL", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
plot(df.gen.DMD$Gen, df.gen.DMD$FFL, ylab = "Plastic genes : Number of DMD", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
plot(df.gen.FFL2$Gen, df.gen.FFL2$FFL, ylab = "Stable genes : Number of FFL", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
plot(df.gen.DMD2$Gen, df.gen.DMD2$FFL, ylab = "Stable genes : Number of DMD", xlab = "Generation")
abline(v = 4000, col="red", lwd=3, lty=2)
dev.off()


