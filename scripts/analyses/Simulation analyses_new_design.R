source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/Full_netw", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
filename <- "full"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.01 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.01 # difference accepted in the RN linear regression intercept
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = TRUE)
pdf(paste0("figures/",filename,"_fitness",".pdf"), width=5, height=4)
plot(df.full$Gen, df.full$Fitness)
plot(subset(df.full, Gen >=9000)$Optimum, subset(df.full, Gen >=9000)$Fitness)
dev.off()
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)
plot(df.last50$Pmean_1, df.last50$Pmean_26)
abline(a = 0.5, b = 0.3, col=2, lwd=5) 
abline(a = 0.5, b = 0.05, col=3, lwd=5) 

# Rprof()
topo_all <- mclapply(unique(df.last50[,1]), function(netw){
  df <- subset(df.last50, data.dir==netw)
  topo <- lapply(1:genes, function(i){
    reg50 <- .lm.fit(cbind(rep(1, length(df$Pmean_1)), df$Pmean_1), df[,i+2])$coefficients # x9 times faster than lm()
    output <- essential.topo(df=subset(df.full, data.dir==netw), target=i,
                                                      treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes)
    ifelse(abs(reg50[2])<=0.01, attr(output, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=0.5, attr(output, "is.plastic") <- "plastic", attr(output, "is.plastic") <- "weak_plastic"))
    return(output)
  })
  gc(verbose = FALSE)
  return(topo)
}, mc.cores = 3)
# summaryRprof()


#Saving outputs for plastic and non plastic genes separatly
topo_plastic <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "plastic"], recursive = FALSE)
saveRDS(topo_plastic, paste0("scripts/data/list_plastic_topo",filename,".Rds"))
#
topo_nnplast <-  unlist(unlist(topo_all, recursive = FALSE)[sapply(unlist(topo_all, recursive = FALSE), FUN=attr, "is.plastic") == "non_plastic"], recursive = FALSE)
saveRDS(topo_nnplast, paste0("scripts/data/list_nnplast_topo.Rds",filename,".Rds"))

################################################################################
############FFL
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 1)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 1)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_FFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_FFL",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
df
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FFL",".csv"))

############FFL count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 1)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 1)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFFL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFFL",".csv"))
################################################################################
############Diamond
Plastic <- FFL.type2(topo_plastic, edges1 = 2, edges2 = 2)
Nnplast <- FFL.type2(topo_nnplast, edges1 = 2, edges2 = 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_plast_DMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_control_DMD",".csv"))

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_DMD",".csv"))

############DMD count
Plastic <- loops_n.count(topo_plastic, edges1 = 2, edges2 = 2)
Nnplast <- loops_n.count(topo_nnplast, edges1 = 2, edges2 = 2)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrDMD",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrDMD",".csv"))
################################################################################
#############FBL
edges <- c(2:6)
Plastic <- FBL.type(topo_plastic, edges = edges)
Nnplast <- FBL.type(topo_nnplast, edges = edges)

df <- as.data.frame(rbind(colSums(Plastic)*100/nrow(Plastic), colSums(Nnplast)*100/nrow(Nnplast)))
rownames(df) <- c("Plastic\ngenes","Non-plastic\ngenes")
write.csv(df, paste0("scripts/data/",filename,"_FBL",".csv"))

#############FBL count
Plastic <- FBL_n.count(topo_plastic, edges=edges)
Nnplast <- FBL_n.count(topo_nnplast, edges=edges)
write.csv(Plastic, paste0("scripts/data/",filename,"_Pl_nbrFBL",".csv"))
write.csv(Nnplast, paste0("scripts/data/",filename,"_NP_nbrFBL",".csv"))





