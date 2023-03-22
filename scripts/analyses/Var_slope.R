source("functions/functions.R")
library(ggpubr)
#####################
sims.dirs <- list.dirs("../simul/Canalization", recursive = TRUE)
#####################
# RN genetic SD
df.canal.sp2 <- df.sampling2(sims.dirs)
#OR
df.canal.sp2 <- read.csv("data/canal_sampling2.csv", header=TRUE, sep = ",")

df.canal.sp2$Var_Slope2 <- sqrt(df.canal.sp2$Var_Slope2)
df.canal.sp2$Var_Slope3 <- sqrt(df.canal.sp2$Var_Slope3)

### RN Mutational SD
df.canal.sp3 <- df.sampling3(c("../simul/Canalization/Canal/rep-0001"))
#OR
df.canal.sp3 <- read.csv("data/canal_sampling3.csv", header=TRUE, sep = ",")


g1 <- ggline(df.canal.sp2, x="Before_Gen", y="Var_Slope2",add = c("mean_sd"),
             ylab = "RN Genetic SD: correlated", xlab = "Generation", ylim=c(0,0.15), color = "darkblue")
g2 <- ggline(df.canal.sp2, x="Before_Gen", y="Var_Slope3",add = c("mean_sd"), 
             ylab = "RN Genetic SD: constant", xlab = "Generation", ylim=c(0,0.15), color = "lightblue2")

g1.1 <- ggline(df.canal.sp3, x="Before_Gen", y="Sd_ind_2", add = c("mean_sd"),
               ylab = "RN Mutational SD : correlated", xlab = "Generation", ylim=c(0,0.15), color = "darkblue")

g2.2 <- ggline(df.canal.sp3, x="Before_Gen", y="Sd_ind_3", add = c("mean_sd"),
               ylab = "RN Mutational SD : constant", xlab = "Generation", ylim=c(0,0.15), color = "lightblue2")


cairo_pdf("../figures/fig_sd_genet.pdf", width=8.5, height=6)
#Plot 1
grid.arrange(
  g1,g1.1,g2,g2.2,
  ncol = 2,
  nrow = 2,
  clip = FALSE
)
dev.off()



 
# #####################
# rep <- 1
# sd <- 0.5
# sims.dirs <- list.dirs("../simul/GenFluct", full.names = TRUE, recursive = TRUE)
# sims.dirs <- sims.dirs[grep('new_envir', sims.dirs, invert = TRUE)]
# #####################
# 
# df.var_ind <- df.sampling3(sims.dirs)
# 
# df.var_ind$Before_Gen <- as.numeric(df.var_ind$Before_Gen)
# 
# ggplot(df.var_ind, aes(x = as.numeric(Before_Gen), y = Sd_ind_2))+
#   geom_point()+geom_smooth(span=0.1, alpha=0)+ theme_bw()+ ylim(c(0,0.25))
# 
# ggplot(df.var_ind, aes(x = as.numeric(Before_Gen), y = Sd_ind_3))+
#   geom_point()+geom_smooth(span=0.1, alpha=0)+ theme_bw()+ ylim(c(0,0.25))
# 
# library(ggpubr)
# ggerrorplot(df.var_ind, x="Before_Gen", y="Sd_ind_2", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean")
# 
# g1 <- ggline(df.var_ind, x="Before_Gen", y="Sd_ind_2",add = c("mean_sd"), ylab = "RN mutational variance of plastic gene")
# g2 <- ggline(df.var_ind, x="Before_Gen", y="Sd_ind_3",add = c("mean_sd"), ylab = "RN mutational variance of stable gene")
# #New Envir
# ################################################################################
# sims.dirs <- list.dirs("../simul/GenFluct", full.names = TRUE, recursive = TRUE)
# sims.dirs <- c(grep("new_envir", sims.dirs, value = TRUE))
# 
# #####################
# df.gf_new_envir <- df.simul(sims.dirs, all.gen = FALSE)
# df.gf_new_envir$Before_Gen <- as.numeric( str_split(str_split(df.gf_new_envir$data.dir, "/", n=8, simplify = TRUE)[,4], "g", simplify = TRUE)[,2] )
# df.gf_new_envir$anc_id <- str_split(str_split(df.gf_new_envir$data.dir, "/", n=8, simplify = TRUE)[,5], "-", simplify = TRUE)[,2]
# df.slope.stab <- slope.stablefrom(df.gf_new_envir)
# 
# for (i in 1:nrow(df.gf_new_envir)) {
#   W <- t(matrix(as.numeric(df.gf_new_envir[i,7:106]), ncol = 10))
#   df.gf_new_envir[i,109] <- getSlope.ALR(W=W, n.env=21, target.gene=2)
# }
# 
# ggplot(df.gf_new_envir)+
#   # geom_point(aes(x = Gen, y=Fitness, col=b_envir))+
#   geom_point(aes(x =Before_Gen, y=V109), show.legend = FALSE)+
#   ylim(c(-1,1))+ylab("Reaction Norm")+xlab("Generations")+theme_bw()
# 
# g3 <- ggline(df.slope.stab, x="V1", y="V2",add = c("mean_sd"), ylab = "Gen before RN stability of plastic G", xlab = "Generation")
# 
# 
# cairo_pdf("../figures/fig_RN_sd.pdf", width=10, height=10.5)
# grid.arrange(
#   g2,g1,g3,
#   ncol = 1,
#   nrow = 3,
#   clip = FALSE
# )
# dev.off()