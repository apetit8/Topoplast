
##RESULTS DATA
control_np <- read.csv("scripts/data/nonplast_E_coli_random2_FFL.csv", sep = ",")
control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_FFL.csv", sep = ",")
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_005_FFL.csv", sep = ",")[,1:11]
drift <- read.csv("scripts/data/full_netw_005_drift_FFL.csv", sep = ",")[,1:11]
s_shuffled <- read.csv("scripts/data/full_netw_005_shuffled_FFL.csv", sep = ",")[,1:11]

df_ffl <- as.data.frame(rbind(colSums(control_np[,5:12])*100/sum(control_np[,3]), colSums(control_pl[,5:12])*100/sum(control_pl[,3]),
                              colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                              theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2],
                              s_shuffled[2,4:11]*100/s_shuffled[2,2],s_shuffled[1,4:11]*100/s_shuffled[1,2],
                              drift[2,4:11]*100/drift[2,2],drift[1,4:11]*100/drift[1,2]))
rownames(df_ffl) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nShuffled","Plastic\nShuffled", "Non-plastic\nDrift","Plastic\nDrift")
df_ffl <- df_ffl[,order(as.character(colnames(df_ffl)), method = c("radix"))]
 
control_np <- read.csv("scripts/data/nonplast_E_coli_random2_diamond.csv", sep = ",")
control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_diamond.csv", sep = ",")
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/full_netw_005_DMD.csv", sep = ",")[,1:13]
drift <- read.csv("scripts/data/full_netw_005_drift_DMD.csv", sep = ",")[,1:13]
s_shuffled <- read.csv("scripts/data/full_netw_005_shuffled_DMD.csv", sep = ",")[,1:13]

df1 <- as.data.frame(rbind(colSums(control_np[,5:14])*100/sum(control_np[,3]), colSums(control_pl[,5:14])*100/sum(control_pl[,3]),
                           colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                           theory[2,4:13]*100/theory[2,2],theory[1,4:13]*100/theory[1,2],
                           s_shuffled[2,4:13]*100/s_shuffled[2,2],s_shuffled[1,4:13]*100/s_shuffled[1,2],
                           drift[2,4:13]*100/drift[2,2],drift[1,4:13]*100/drift[1,2]))
rownames(df1) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory", "Non-plastic\nShuffled","Plastic\nShuffled", "Non-plastic\nDrift","Plastic\nDrift")
df1 <- df1[,order(as.character(colnames(df1)), method = c("radix"))]


# df.full <- df.simul(list.dirs(paste0("simul/Full_netw"), recursive = TRUE), all.gen = TRUE)
# df.full.drift <- df.simul(list.dirs(paste0("simul/Full_netw_drift"), recursive = TRUE), all.gen = TRUE)


pdf(paste0("figures/Figsupp_corr",".pdf"), width=4, height=4)
par(mar = c(3.,3., 1.2,0.1))

##DEFAULT
plot(unlist(c(df_ffl[3,], df1[3,], df_ffl[4,], df1[4,])), unlist(c(df_ffl[5,], df1[5,], df_ffl[6,], df1[6,])),
     col=c(rep("red", 18), rep("black", 18)), pch=c(rep(16,8),rep(18,10),rep(16,8),rep(18,8)), ylab="", xlab = "",main="Default", ylim=c(0,100))
abline(lm( unlist(c( df_ffl[6,], df1[6,])) ~  unlist(c(df_ffl[4,], df1[4,] ))), col="black")
abline(lm( unlist(c(df_ffl[5,], df1[5,])) ~  unlist(c(df_ffl[3,], df1[3,] ))), col="red")
legend(-1.3, 104, legend=c("Plastic genes", "Non plastic genes", "FFL","DMD"),
       col=c("black", "red", "black", "black"), bty="o", cex=0.8, bg="white", pch = c(15,15,1,5))
title(ylab="% in simulations", line=2, cex.lab=1)
title(xlab="% in E. coli", line=2, cex.lab=1)

#plot(df.full$Gen, df.full$Fitness, main="Default", ylim=c(0,1), col=alpha("black", 0.2), pch=19)


##DRIFT

plot(unlist(c(df_ffl[3,], df1[3,], df_ffl[4,], df1[4,])), unlist(c(df_ffl[9,], df1[9,], df_ffl[10,], df1[10,])),
     col=c(rep("red", 18), rep("black", 18)), pch=c(rep(16,8),rep(18,10),rep(16,8),rep(18,8)), ylab="", xlab = "", main="Drift", ylim=c(0,100))
abline(lm( unlist(c( df_ffl[10,], df1[10,])) ~  unlist(c(df_ffl[4,], df1[4,] ))), col="black")
abline(lm( unlist(c(df_ffl[9,], df1[9,])) ~  unlist(c(df_ffl[3,], df1[3,] ))), col="red")
title(ylab="% in simulations", line=2, cex.lab=1)
title(xlab="% in E. coli", line=2, cex.lab=1)
#plot(df.full.drift$Gen, df.full.drift$Fitness, main="Drift", ylim=c(0,1), col=alpha("black", 0.2), pch=19)



##SUPPLEMENTAL
for (simu in c("a_less","a_more","pop_200","pop_2000","genes_more","genes_less","TF_more","TF_less")) {
  
  ################################################################################
  ##FFL motifs Data###############################################################*
  if(simu == "a_less") main <- "Constitutive expression K = 0.05"
  if(simu == "a_more") main <- "Constitutive expression K = 0.5"
  if(simu == "pop_200") main <- "Population size = 200"
  if(simu == "pop_2000") main <- "Population size = 2000"
  if(simu == "genes_more") main <- "Number of genes = 66"
  if(simu == "genes_less") main <- "Number of genes = 21"
  if(simu == "TF_more") main <- "Number of TF = 20"
  if(simu == "TF_less") main <- "Number of TF = 5"
  
  control_np <- read.csv("scripts/data/nonplast_E_coli_random2_FFL.csv", sep = ",")
  control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_FFL.csv", sep = ",")
  non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
  all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
  theory <- read.csv(paste0("scripts/data/full_netw_",simu,"_FFL.csv"), sep = ",")[,1:11]
  
  df_ffl <- as.data.frame(rbind(colSums(control_np[,5:12])*100/sum(control_np[,3]), colSums(control_pl[,5:12])*100/sum(control_pl[,3]),
                                colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                                theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2]))
  rownames(df_ffl) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
  df_ffl <- df_ffl[,order(as.character(colnames(df_ffl)), method = c("radix"))]
  
  control_np <- read.csv("scripts/data/nonplast_E_coli_random2_diamond.csv", sep = ",")
  control_pl <- read.csv("scripts/data/plast_genes_E_coli_random2_diamond.csv", sep = ",")
  non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
  all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
  theory <- read.csv(paste0("scripts/data/full_netw_",simu,"_DMD.csv"), sep = ",")[,1:13]
  
  df1 <- as.data.frame(rbind(colSums(control_np[,5:14])*100/sum(control_np[,3]), colSums(control_pl[,5:14])*100/sum(control_pl[,3]),
                             colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                             theory[2,4:13]*100/theory[2,2],theory[1,4:13]*100/theory[1,2]))
  rownames(df1) <- c("\n\n\nNon-Plastic\nrandomized\nE. coli", "\nPlastic\nrandomized\nE. coli", "Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
  df1 <- df1[,order(as.character(colnames(df1)), method = c("radix"))]
  
  
  plot(unlist(c(df_ffl[3,], df1[3,], df_ffl[4,], df1[4,])), unlist(c(df_ffl[5,], df1[5,], df_ffl[6,], df1[6,])),
       col=c(rep("red", 18), rep("black", 18)), pch=c(rep(16,8),rep(18,10),rep(16,8),rep(18,8)), ylab="", xlab = "", main=main, ylim=c(0,100))
  abline(lm( unlist(c( df_ffl[6,], df1[6,])) ~  unlist(c(df_ffl[4,], df1[4,] ))), col="black")
  abline(lm( unlist(c(df_ffl[5,], df1[5,])) ~  unlist(c(df_ffl[3,], df1[3,] ))), col="red")
  
  title(ylab="% in simulations", line=2, cex.lab=1)
  title(xlab="% in E. coli", line=2, cex.lab=1)
  
  # #Add here plot with fitness
  # df.full.gen <- df.simul(list.dirs(paste0("simul/Full_netw_",simu), recursive = TRUE), all.gen = TRUE)
  # plot(df.full.gen$Gen, df.full.gen$Fitness, main=simu, ylim=c(0,1), col=alpha("black", 0.2), pch=19)
  # 
  # 
}

dev.off()
