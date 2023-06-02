source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/30g", recursive = TRUE)
genes <- 30
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/30g"
#####################
#DATA
##################
df.30 <- df.simul(sims.dirs1, all.gen = TRUE)
df.30$envir <- str_split(df.30$data.dir, "/", n=8, simplify = TRUE)[,3]
df.30$anc_id <- str_split(str_split(df.30$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]

################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
gen <- 3000 #max(df.30$Gen)
#############

topo.anticor30 <- essential.topo(df=subset(df.30, Gen==gen & envir=="Anticorrelated"),
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:30))

topo.corr30 <- essential.topo(df=subset(df.30, Gen==gen & envir=="Correlated"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:30))

topo.no_sel30 <- essential.topo(df=subset(df.30, Gen==gen & envir=="Control_no_sel"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:30))

topo.sel30 <- essential.topo(df=subset(df.30, Gen==gen & envir=="Control_sel"),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:30))


################
# Cutoffs are number of NODES
edges1 <- 2
edges2 <- 1
##### Loop count
Anticor30_n <- loops_n.count(topo.anticor30, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
Corr30_n <- loops_n.count(topo.corr30, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
No_sel30_n <- loops_n.count(topo.no_sel30, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2) #
Sel30_n <- loops_n.count(topo.sel30, cutoff.max = cutoff.max, cutoff.min = cutoff.min)

pdf(paste0(pdfname,"_N_loops",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor30_n$Loop_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr30_n$Loop_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel30_n$Loop_number, main = "No Sel", ylim = c(0, 300))
hist(Sel30_n$Loop_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()


################
##### Type
Anticor30 <- FFL.type2(topo.anticor30, edges1 = 2, edges2 = 1)
Corr30 <- FFL.type2(topo.corr30, edges1 = 2, edges2 = 1)
No_sel30 <- FFL.type2(topo.no_sel30, edges1 = 2, edges2 = 1)
Sel30 <- FFL.type2(topo.sel30, edges1 = 2, edges2 = 1)

df <- as.data.frame(rbind(colSums(Anticor30), colSums(Corr30), colSums(No_sel30), colSums(Sel30)))
rownames(df) <- c("Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:30)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,30)]))


pdfname <- "figures/30g"
pdf(paste0(pdfname,"_type_FFL2_2_1_edges",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:30])*100/150, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","Input_Dep_Amplifying_neg","Input_Dep_Amplifying_pos","Input_Dep_Disruptive_neg","Input_Dep_Disruptive_pos","Input_Ind_Amplifying_neg","Input_Ind_Amplifying_pos","Input_Ind_Disruptive_neg","Input_Ind_Disruptive_pos") )
#Coherence
barplot(t(df[,c(2,11,12)])*100/150, col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,13,14)])*100/150, col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

#
Anticor30 <- FFL.type2(topo.anticor30, edges1 = 2, edges2 = 2)
Corr30 <- FFL.type2(topo.corr30, edges1 = 2, edges2 = 2)
No_sel30 <- FFL.type2(topo.no_sel30, edges1 = 2, edges2 = 2) #
Sel30 <- FFL.type2(topo.sel30, edges1 = 2, edges2 = 2)

df <- as.data.frame(rbind(colSums(Anticor30), colSums(Corr30), colSums(No_sel30), colSums(Sel30)))
rownames(df) <- c("Anticor","Cor","No_sel","Sel")
df$coherent <- rowSums2(as.matrix(df[,c(3,5,6,8,30,11)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(4,7,9,12)]))
df$homogenous_pos <- rowSums2(as.matrix(df[,c(3,6,9)]))
df$homogenous_neg <- rowSums2(as.matrix(df[,c(5,8,12)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(4,7,30,11)]))

pdfname <- "figures/30g"
pdf(paste0(pdfname,"_type_FFL2_diamond",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12])*100/150, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,13,14)])*100/150, col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,15:17)])*100/150, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()

# df <- as.data.frame(dplyr::bind_rows(table(Anticor30$Motif_sign), table(Corr30$Motif_sign), table(No_sel30$Motif_sign), table(Sel30$Motif_sign)))
# rownames(df) <- c("Anticor","Cor","No_sel","Sel")
# df[is.na(df)] <- 0
# df <- df[, sort(colnames(df))]
# 
# pdf(paste0(pdfname,"_type_FFL2_sign",".pdf"), width=6.5, height=6)
# layout(matrix(c(1:1), 1, 1, byrow = TRUE))
# barplot(t(df)*300/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
# legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
#        legend=c( "No_FFL","Input_Dep_Amplifying_neg","Input_Dep_Amplifying_pos","Input_Dep_Disruptive_neg","Input_Dep_Disruptive_pos","Input_Ind_Amplifying_neg","Input_Ind_Amplifying_pos","Input_Ind_Disruptive_neg","Input_Ind_Disruptive_pos") )
# dev.off()



################
##### Homogeneity
Anticor30 <- homog.count(topo.anticor30, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
Corr30 <- homog.count(topo.corr30, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
No_sel30 <- homog.count(topo.no_sel30, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2, randomFF=FALSE) #
Sel30 <- homog.count(topo.sel30, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)

df <- rbind(colSums(Anticor30), colSums(Corr30), colSums(No_sel30), colSums(Sel30))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor30_n$Loop_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr30_n$Loop_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel30_n$Loop_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel30_n$Loop_number), 1)))

pdf(paste0(pdfname,"_homog_loops",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/150, col=c(2,4,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(2,4,"grey"),
       legend=c("Het FFL", "Hom FFL","No FFL"))
dev.off()

################################################################################
################################################################################
#FEEDBACKS

#Renforcing = same sign
#Defusing = opposite sign

################
cutoff <- 3
##### Loop count
Anticor30_n <- FBL_n.count(topo.anticor30, cutoff = cutoff, target = 2)
Corr30_n <- FBL_n.count(topo.corr30, cutoff = cutoff, target = 2)
No_sel30_n <- FBL_n.count(topo.no_sel30, cutoff = cutoff, target = 2) #
Sel30_n <- FBL_n.count(topo.sel30, cutoff = cutoff, target = 2)

pdf(paste0(pdfname,"_N_FBL",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor30_n$FBL_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr30_n$FBL_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel30_n$FBL_number, main = "No Sel", ylim = c(0, 300))
hist(Sel30_n$FBL_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()

################
##### Homogeneity
Anticor30 <- FBL.type(topo.anticor30, cutoff = cutoff, target = 2, randomFF=FALSE)
Corr30 <- FBL.type(topo.corr30, cutoff = cutoff, target = 2, randomFF=FALSE)
No_sel30 <- FBL.type(topo.no_sel30, cutoff = cutoff, target = 2, randomFF=FALSE) #
Sel30 <- FBL.type(topo.sel30, cutoff = cutoff, target = 2, randomFF=FALSE)

df <- rbind(colSums(Anticor30), colSums(Corr30), colSums(No_sel30), colSums(Sel30))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor30_n$FBL_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr30_n$FBL_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel30_n$FBL_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel30_n$FBL_number), 1)))

pdf(paste0(pdfname,"_type_FBL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/150, col=c("sienna1","orchid3","grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("sienna1","orchid3","grey"),
       legend=c("Defusing", "Renforcing FFL","No FFL"))
dev.off()


