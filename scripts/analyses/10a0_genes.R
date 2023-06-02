source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/10g_a0.15", recursive = TRUE)
genes <- 10
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/10g_a0.15"
#####################
#DATA
##################
df.10_a0_15 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10_a0_15$envir <- str_split(df.10_a0_15$data.dir, "/", n=8, simplify = TRUE)[,3]
df.10_a0_15$anc_id <- str_split(str_split(df.10_a0_15$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]

################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
gen <- 3000 #max(df.10_a0_15$Gen)
#############

topo.anticor10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Anticorrelated"),
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Correlated"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.no_sel10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Control_no_sel"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Control_sel"),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

ifelse(!dir.exists(file.path("scripts/data/essential_motif")), dir.create("scripts/data/essential_motif"), FALSE)#burning_andreas folder
lapply(topo.anticor10_a0_15, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.anticor10_a0_15.csv'  , append= T, sep=',' ))
lapply(topo.anticor10_a0_15, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.corr10_a0_15.csv'  , append= T, sep=',' ))
lapply(topo.anticor10_a0_15, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.no_sel10_a0_15.csv'  , append= T, sep=',' ))
lapply(topo.anticor10_a0_15, function(x) write.table( data.frame(x), 'scripts/data/essential_motif/topo.sel10_a0_15.csv'  , append= T, sep=',' ))
################
##### Type
Anticor10_a0_15 <- FFL.type2(topo.anticor10_a0_15, edges1 = 2, edges2 = 1)
Corr10_a0_15 <- FFL.type2(topo.corr10_a0_15, edges1 = 2, edges2 = 1)
No_sel10_a0_15 <- FFL.type2(topo.no_sel10_a0_15, edges1 = 2, edges2 = 1) #
Sel10_a0_15 <- FFL.type2(topo.sel10_a0_15, edges1 = 2, edges2 = 1)

df <- as.data.frame(rbind(colSums(Anticor10_a0_15), colSums(Corr10_a0_15), colSums(No_sel10_a0_15), colSums(Sel10_a0_15)))
rownames(df) <- c("Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/10g_a0_15"
pdf(paste0(pdfname,"_type_FFL2_cutoff_3",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:10])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","Input_Dep_Amplifying_neg","Input_Dep_Amplifying_pos","Input_Dep_Disruptive_neg","Input_Dep_Disruptive_pos","Input_Ind_Amplifying_neg","Input_Ind_Amplifying_pos","Input_Ind_Disruptive_neg","Input_Ind_Disruptive_pos") )
#Coherence
barplot(t(df[,c(2,11,12)])*100/300, col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,13,14)])*100/300, col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

#
Anticor10_a0_15 <- FFL.type2(topo.anticor10_a0_15, edges1 = 2, edges2 = 2)
Corr10_a0_15 <- FFL.type2(topo.corr10_a0_15, edges1 = 2, edges2 = 2)
No_sel10_a0_15 <- FFL.type2(topo.no_sel10_a0_15, edges1 = 2, edges2 = 2) #
Sel10_a0_15 <- FFL.type2(topo.sel10_a0_15, edges1 = 2, edges2 = 2)

df <- as.data.frame(rbind(colSums(Anticor10_a0_15), colSums(Corr10_a0_15), colSums(No_sel10_a0_15), colSums(Sel10_a0_15)))
rownames(df) <- c("Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3,5,6,8,10,11)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(4,7,9,12)]))
df$homogenous_pos <- rowSums2(as.matrix(df[,c(3,6,9)]))
df$homogenous_neg <- rowSums2(as.matrix(df[,c(5,8,12)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(4,7,10,11)]))

pdfname <- "figures/10g_a0_15"
pdf(paste0(pdfname,"_type_FFL2_diamond",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,13,14)])*100/300, col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,15:17)])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()


################################################################################
################################################################################
#FEEDBACKS

#Renforcing = same sign
#Defusing = opposite sign

################
cutoff <- 3
##### Loop count
Anticor10_a0_15_n <- FBL_n.count(topo.anticor10_a0_15, cutoff = cutoff, target = 2)
Corr10_a0_15_n <- FBL_n.count(topo.corr10_a0_15, cutoff = cutoff, target = 2)
No_sel10_a0_15_n <- FBL_n.count(topo.no_sel10_a0_15, cutoff = cutoff, target = 2) #
Sel10_a0_15_n <- FBL_n.count(topo.sel10_a0_15, cutoff = cutoff, target = 2)

pdf(paste0(pdfname,"_N_FBL",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor10_a0_15_n$FBL_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr10_a0_15_n$FBL_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel10_a0_15_n$FBL_number, main = "No Sel", ylim = c(0, 300))
hist(Sel10_a0_15_n$FBL_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()

################
##### Homogeneity
Anticor10_a0_15 <- FBL.type(topo.anticor10_a0_15, cutoff = cutoff, target = 2, randomFF=FALSE)
Corr10_a0_15 <- FBL.type(topo.corr10_a0_15, cutoff = cutoff, target = 2, randomFF=FALSE)
No_sel10_a0_15 <- FBL.type(topo.no_sel10_a0_15, cutoff = cutoff, target = 2, randomFF=FALSE) #
Sel10_a0_15 <- FBL.type(topo.sel10_a0_15, cutoff = cutoff, target = 2, randomFF=FALSE)

df <- rbind(colSums(Anticor10_a0_15), colSums(Corr10_a0_15), colSums(No_sel10_a0_15), colSums(Sel10_a0_15))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor10_a0_15_n$FBL_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr10_a0_15_n$FBL_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel10_a0_15_n$FBL_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel10_a0_15_n$FBL_number), 1)))

pdf(paste0(pdfname,"_type_FBL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/300, col=c("sienna1","orchid3","grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("sienna1","orchid3","grey"),
       legend=c("Defusing", "Renforcing FFL","No FFL"))
dev.off()


