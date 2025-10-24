source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/Full_netw_pop_200", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_pop_200"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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

source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")

################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_pop_2000", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_pop_2000"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
#Loop to identify plastic genes
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_a_less", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_a_less"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_a_more", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_a_more"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
    return(output2)
  })
  gc(verbose = FALSE)
  return(topo)
}, mc.cores = 3)
# summaryRprof()
#sum(unlist(unlist(topo_all, recursive = FALSE))==1)/(sum(unlist(unlist(topo_all, recursive = FALSE))==-1)+sum(unlist(unlist(topo_all, recursive = FALSE))==1))*100
#57% positive regulation

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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_TF_more", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_TF_more"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_TF_less", recursive = TRUE)
genes <- 36
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_TF_less"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_genes_more", recursive = TRUE)
genes <- 66
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_genes_more"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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


################################################################################
################################################################################
################################################################################
sims.dirs1 <- list.dirs("simul/Full_netw_genes_less", recursive = TRUE)
genes <- 21
min <- 0.15
max <- 0.85
basal <- 0.15 #Basal expression
filename <- "full_netw_genes_less"
################################################################################
#Parameters when keeping "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.005 # difference accepted in the Reaction Norm linear regression slope
treshold_og    <- 0.005 # difference accepted in the RN linear regression intercept
#####################
slopemin_plast <- 0.65
slopemax_nnplast <- 0.1
#####################
#DATA
#####################
df.full <- df.simul(sims.dirs1, all.gen = FALSE)
#############
#############
df.last50 <- df.last.gens(sims.dirs1)

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
    ifelse(abs(reg50[2])<=slopemax_nnplast, attr(output2, "is.plastic") <- "non_plastic", ifelse(abs(reg50[2])>=slopemin_plast, attr(output2, "is.plastic") <- "plastic", attr(output2, "is.plastic") <- "weak_plastic"))
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

