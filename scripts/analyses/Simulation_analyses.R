source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/Full_netw", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_005"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
#DATA
#####################
# df.full.gen <- df.simul(sims.dirs1, all.gen = TRUE)
# .pdf(paste0("figures/",filename,"_fitness",".pdf"), width=5, height=4)
# plot(df.full.gen$Gen, df.full.gen$Fitness)
# plot(subset(df.full.gen, Gen >=1000)$Optimum, subset(df.full.gen, Gen >=1000)$Fitness)
# dev.off()
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#df.full <- subset(df.full, Gen ==4000)
# df.full2 <- subset(df.full, Gen ==4000)
#plot(df.full2$Optimum, df.full2$Fitness, xlim=c(0,1), ylim=c(0,1))

#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)
plot(df.last50$Pmean_1, df.last50$Pmean_16)
abline(a = 0.1, b = 0.5, col=2, lwd=5)
abline(a = 0.1, b = 0.505, col=2, lwd=5)
# abline(a = 0.5, b = 0.51, col=2, lwd=5)
# abline(a = 0.5, b = 0.52, col=2, lwd=5)
# abline(a = 0.5, b = 0.05, col=3, lwd=5)

# Rprof()
topo_all <- mclapply(unique(df.last50[1:(nrow(df.last50)),1]), function(netw){
  df <- subset(df.last50, data.dir==netw)
  topo <- lapply(2:genes, function(i){ #start at 12 to not include TFs
    reg50 <- .lm.fit(cbind(rep(1, length(df$Pmean_1)), df$Pmean_1), df[,i+2])$coefficients # x9 times faster than lm()
    output <- essential.topo(df=subset(df.full, data.dir==netw), target=i, basal=basal,
                                                      treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)
    #Reorder network so the target is in position #2 in rows and columns. Necessary in further analyses, simpler than keeping "i" in topo_all object
    output1 <- output
    output1[[1]][,2] <- output[[1]][,i]
    output1[[1]][,i] <- output[[1]][,2]
    output2 <- output1
    output2[[1]][2,] <- output1[[1]][i,]
    output2[[1]][i,] <- output1[[1]][2,]
    ifelse(abs(reg50[2])<=0.025, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=0.5, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
    return(output2)
  })
  gc(verbose = FALSE)
  return(topo)
}, mc.cores = 3)
# summaryRprof()


#Saving outputs for plastic and non plastic genes separatly
topo_plastic <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "plastic"], recursive = FALSE)
saveRDS(topo_plastic, paste0("scripts/data/list_plastic_topo_",filename,".Rds"))
#
topo_nnplast <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "non_plastic"], recursive = FALSE)
saveRDS(topo_nnplast, paste0("scripts/data/list_nnplast_topo_",filename,".Rds"))

################################################################################
############FFL
topo_plastic <- readRDS(paste0("scripts/data/list_plastic_topo_",filename,".Rds"))
topo_nnplast <- readRDS(paste0("scripts/data/list_nnplast_topo_",filename,".Rds"))
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FFL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FFL",".csv"))

############FFL count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFFL",".csv"))
################################################################################
############Diamond
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_DMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_DMD",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_DMD",".csv"))

############DMD count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrDMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrDMD",".csv"))
################################################################################
#############FBL
edges <- c(2:5)
Plastic <- FBL.type(topo_plastic, edges = edges, target= 2)
Nnplast <- FBL.type(topo_nnplast, edges = edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FBL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FBL",".csv"))

#############FBL count
Plastic <- FBL_n.count(topo_plastic, edges=edges, target= 2)
Nnplast <- FBL_n.count(topo_nnplast, edges=edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFBL",".csv"))

#####################
################################################################################
#CONTROL: SHUFFLED
################################################################################
#####################
sims.dirs1 <- list.dirs("simul/Full_netw", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
filename <- "full_netw_005_shuffled"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
#DATA
#####################
# df.full.gen <- df.simul(sims.dirs1, all.gen = TRUE)
# pdf(paste0("figures/",filename,"_fitness",".pdf"), width=5, height=4)
# plot(df.full.gen$Gen, df.full.gen$Fitness)
# plot(subset(df.full.gen, Gen >=1000)$Optimum, subset(df.full.gen, Gen >=1000)$Fitness)
# dev.off()
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#df.full <- subset(df.full, Gen ==4000)
#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)
# plot(df.last50$Pmean_1, df.last50$Pmean_26)
# abline(a = 0.5, b = 0.3, col=2, lwd=5)
# abline(a = 0.5, b = 0.05, col=3, lwd=5)

# Rprof()
topo_all <- mclapply(unique(df.last50[2:(nrow(df.last50)),1]), function(netw){ #CAREFULL only 1/5 simuls for test
  df <- subset(df.last50, data.dir==netw)
  topo <- lapply(2:genes, function(i){ #starts at 12 to not include TFs
    reg50 <- .lm.fit(cbind(rep(1, length(df$Pmean_1)), df$Pmean_1), df[,i+2])$coefficients # x9 times faster than lm()
    output <- essential.topo(df=subset(df.full, data.dir==netw), target=i, basal=basal,
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)
    #Reorder network so the target is in position #2 in rows and columns. Necessary in further analyses, simpler than keeping "i" in topo_all object
    output1 <- output
    output1[[1]][,2] <- output[[1]][,i]
    output1[[1]][,i] <- output[[1]][,2]
    output2 <- output1
    output2[[1]][2,] <- output1[[1]][i,]
    output2[[1]][i,] <- output1[[1]][2,]
    ifelse(abs(reg50[2])<=0.025, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=0.5, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
    return(output2)
  })
  gc(verbose = FALSE)
  return(topo)
}, mc.cores = 3)
# summaryRprof()


#Saving outputs for plastic and non plastic genes separatly
topo_plastic <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "plastic"], recursive = FALSE)

# topo_plastic <- readRDS(paste0("scripts/data/list_plastic_topo_",filename1,".Rds"))
# topo_nnplast <- readRDS(paste0("scripts/data/list_nnplast_topo_",filename1,".Rds"))

for (i in 1:length(topo_plastic)) {
  for (j in 1:length(topo_plastic[[1]])) {
    if(topo_plastic[[i]][j] !=0 ) topo_plastic[[i]][j] <- sample(topo_plastic[[i]][(topo_plastic[[i]] != 0)], 1)
  }
}  

saveRDS(topo_plastic, paste0("scripts/data/list_plastic_topo_",filename,".Rds"))
#
topo_nnplast <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "non_plastic"], recursive = FALSE)

for (i in 1:length(topo_nnplast)) {
  for (j in 1:length(topo_nnplast[[1]])) {
    if(topo_nnplast[[i]][j] !=0 ) topo_nnplast[[i]][j] <- sample(topo_nnplast[[i]][(topo_nnplast[[i]][] != 0)], 1)
  }
} 

saveRDS(topo_nnplast, paste0("scripts/data/list_nnplast_topo_",filename,".Rds"))

################################################################################
############FFL
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
Drift   <- FFL.type2(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FFL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_FFL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FFL",".csv"))

############FFL count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
Drift <- loops_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFFL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_Dr_nbrFFL",".csv"))
################################################################################
############Diamond
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
Drift <- FFL.type2(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_DMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_DMD",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_DMD",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_DMD",".csv"))

############DMD count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
Drift <- loops_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrDMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrDMD",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_Dr_nbrDMD",".csv"))
################################################################################
#############FBL
edges <- c(2:5)
Plastic <- FBL.type(topo_plastic, edges = edges, target= 2)
Nnplast <- FBL.type(topo_nnplast, edges = edges, target= 2)
Drift <- FBL.type(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges = edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FBL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_FBL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FBL",".csv"))

#############FBL count
Plastic <- FBL_n.count(topo_plastic, edges=edges, target= 2)
Nnplast <- FBL_n.count(topo_nnplast, edges=edges, target= 2)
Drift <- FBL_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges=edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFBL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_nbrFBL",".csv"))


################################################################################
#CONTROL: DRIFT
################################################################################
#####################
sims.dirs1 <- list.dirs("simul/Full_netw_drift", recursive = TRUE)
genes <- 36
basal <- 0.15
min <- 0.15
max <- 0.85
filename <- "full_netw_005_drift"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
#DATA
#####################
# df.full.gen <- df.simul(sims.dirs1, all.gen = TRUE)
# plot(df.full.gen$Gen, df.full.gen$Fitness)
# plot(subset(df.full.gen, Gen >=1000)$Optimum, subset(df.full.gen, Gen >=1000)$Fitness)
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)
plot(df.last50$Pmean_1, df.last50$Pmean_3)
#plot(subset(df.last50, Gen==500)$Pmean_1, subset(df.last50, Gen==500)$Pmean_8)
abline(a = 0.5, b = 0.3, col=2, lwd=5) 
abline(a = 0.5, b = 0.025, col=3, lwd=5) 

topo_all <- mclapply(unique(df.last50[1:(nrow(df.last50)),1]), function(netw){
  df <- subset(df.last50, data.dir==netw)
  topo <- lapply(2:genes, function(i){ #start at 12 to not include TFs
    reg50 <- .lm.fit(cbind(rep(1, length(df$Pmean_1)), df$Pmean_1), df[,i+2])$coefficients # x9 times faster than lm()
    output <- essential.topo(df=subset(df.full, data.dir==netw), target=i, basal=basal,
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)
    #Reorder network so the target is in position #2 in rows and columns. Necessary in further analyses, simpler than keeping "i" in topo_all object
    output1 <- output
    output1[[1]][,2] <- output[[1]][,i]
    output1[[1]][,i] <- output[[1]][,2]
    output2 <- output1
    output2[[1]][2,] <- output1[[1]][i,]
    output2[[1]][i,] <- output1[[1]][2,]
    ifelse(abs(reg50[2])<=0.025, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=0.5, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
    return(output2)
  })
  gc(verbose = FALSE)
  return(topo)
}, mc.cores = 3)


#Saving outputs for plastic and non plastic genes separatly
topo_plastic <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "plastic"], recursive = FALSE)
saveRDS(topo_plastic, paste0("scripts/data/list_plastic_topo_",filename,".Rds"))
#
topo_nnplast <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "non_plastic"], recursive = FALSE)
saveRDS(topo_nnplast, paste0("scripts/data/list_nnplast_topo_",filename,".Rds"))

################################################################################
############FFL
topo_plastic <- readRDS(paste0("scripts/data/list_plastic_topo_",filename,".Rds"))
topo_nnplast <- readRDS(paste0("scripts/data/list_nnplast_topo_",filename,".Rds"))
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
Drift   <- FFL.type2(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FFL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_FFL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FFL",".csv"))

############FFL count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 1, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 1, target= 2)
Drift <- loops_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 1, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFFL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_Dr_nbrFFL",".csv"))
################################################################################
############Diamond
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
Drift <- FFL.type2(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_DMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_DMD",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_DMD",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_DMD",".csv"))

############DMD count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 2, target= 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 2, target= 2)
Drift <- loops_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges1 = 2, edges2 = 2, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrDMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrDMD",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_Dr_nbrDMD",".csv"))
################################################################################
#############FBL
edges <- c(2:5)
Plastic <- FBL.type(topo_plastic, edges = edges, target= 2)
Nnplast <- FBL.type(topo_nnplast, edges = edges, target= 2)
Drift <- FBL.type(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges = edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FBL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_FBL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FBL",".csv"))

#############FBL count
Plastic <- FBL_n.count(topo_plastic, edges=edges, target= 2)
Nnplast <- FBL_n.count(topo_nnplast, edges=edges, target= 2)
Drift <- FBL_n.count(unlist(unlist(topo_all, recursive = FALSE), recursive = FALSE), edges=edges, target= 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFBL",".csv"))
write.csv(Drift, paste0("scripts/data/",filename,"_drift_nbrFBL",".csv"))


