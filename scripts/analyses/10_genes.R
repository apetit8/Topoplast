source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
sims.dirs1 <- list.dirs("simul/10g", recursive = TRUE)
genes <- 10
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname <- "figures/10g"
#####################
#DATA
##################
df.10 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10$envir <- str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,3]
df.10$anc_id <- str_split(str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]

################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
#############

topo.anticor10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Anticorrelated"),
                                 treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Correlated"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.no_sel10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Control_no_sel"),
                                treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10 <- essential.topo(df=subset(df.10, Gen==max(df.10$Gen) & envir=="Control_sel"),
                             treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))


################
cutoff.max <- 3
cutoff.min <- 1
##### Loop count
Anticor10_n <- loops_n.count(topo.anticor10, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
Corr10_n <- loops_n.count(topo.corr10, cutoff.max = cutoff.max, cutoff.min = cutoff.min)
No_sel10_n <- loops_n.count(topo.no_sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2) #
Sel10_n <- loops_n.count(topo.sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min)

pdf(paste0(pdfname,"_N_loops",".pdf"))
  layout(matrix(c(1:4), 2, 2, byrow = TRUE))
  hist(Anticor10_n$Loop_number, main = "Anticor Plast", ylim = c(0, 300))
  hist(Corr10_n$Loop_number, main = "Cor Plast", ylim = c(0, 300))
  hist(No_sel10_n$Loop_number, main = "No Sel", ylim = c(0, 300))
  hist(Sel10_n$Loop_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()

################
##### Coherence
Anticor10 <- c.count(topo.anticor10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)
Corr10 <- c.count(topo.corr10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)
No_sel10 <- c.count(topo.no_sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2, randomFF=TRUE) #
Sel10 <- c.count(topo.sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=TRUE)

df <- rbind(colSums(Anticor10), colSums(Corr10), colSums(No_sel10), colSums(Sel10))
rownames(df) <- c("Anticor10", "Corr10", "No_sel10", "Sel10")

pdf(paste0(pdfname,"_coherence_loops",".pdf"), width=6.5, height=6)
  layout(matrix(c(1:1), 1, 1, byrow = TRUE))
  barplot(t(df[,1:3])*100/300, col=c(7,3,"grey"))
  legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
         legend=c("Coherent FFL", "Incoherent FFL","No Loop"))
dev.off()


################
##### Homogeneity
Anticor10 <- homog.count(topo.anticor10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
Corr10 <- homog.count(topo.corr10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)
No_sel10 <- homog.count(topo.no_sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = 2, randomFF=FALSE) #
Sel10 <- homog.count(topo.sel10, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)

df <- rbind(colSums(Anticor10), colSums(Corr10), colSums(No_sel10), colSums(Sel10))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor10_n$Loop_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr10_n$Loop_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel10_n$Loop_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel10_n$Loop_number), 1)))

pdf(paste0(pdfname,"_homog_loops",".pdf"), width=6.5, height=6)
  layout(matrix(c(1:1), 1, 1, byrow = TRUE))
  barplot(t(df[,1:3])*100/300, col=c(2,4,"grey"))
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
Anticor10_n <- FBL_n.count(topo.anticor10, cutoff = cutoff, target = 2)
Corr10_n <- FBL_n.count(topo.corr10, cutoff = cutoff, target = 2)
No_sel10_n <- FBL_n.count(topo.no_sel10, cutoff = cutoff, target = 2) #
Sel10_n <- FBL_n.count(topo.sel10, cutoff = cutoff, target = 2)

pdf(paste0(pdfname,"_N_FBL",".pdf"))
layout(matrix(c(1:4), 2, 2, byrow = TRUE))
hist(Anticor10_n$FBL_number, main = "Anticor Plast", ylim = c(0, 300))
hist(Corr10_n$FBL_number, main = "Cor Plast", ylim = c(0, 300))
hist(No_sel10_n$FBL_number, main = "No Sel", ylim = c(0, 300))
hist(Sel10_n$FBL_number, main = "Sel Stable", ylim = c(0, 300))
dev.off()

################
##### Homogeneity
Anticor10 <- FBL.type(topo.anticor10, cutoff = cutoff, target = 2, randomFF=FALSE)
Corr10 <- FBL.type(topo.corr10, cutoff = cutoff, target = 2, randomFF=FALSE)
No_sel10 <- FBL.type(topo.no_sel10, cutoff = cutoff, target = 2, randomFF=FALSE) #
Sel10 <- FBL.type(topo.sel10, cutoff = cutoff, target = 2, randomFF=FALSE)

df <- rbind(colSums(Anticor10), colSums(Corr10), colSums(No_sel10), colSums(Sel10))
rownames(df) <- c(paste0("Anticor\n m loop = ", round(mean(Anticor10_n$FBL_number), 1)),
                  paste0("Corr\n m loop =  ", round(mean(Corr10_n$FBL_number), 1)),
                  paste0("No_sel\n m loop =  ", round(mean(No_sel10_n$FBL_number), 1)),
                  paste0("Sel\n m loop =  ", round(mean(Sel10_n$FBL_number), 1)))

pdf(paste0(pdfname,"_type_FBL",".pdf"), width=6.5, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:3])*100/300, col=c("sienna1","orchid3","grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("sienna1","orchid3","grey"),
       legend=c("Defusing", "Renforcing FFL","No FFL"))
dev.off()


