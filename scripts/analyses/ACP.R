library(FactoMineR)
library(factoextra)
library(gridExtra)
library(ggfortify)
library(ggplot2)
pdfname <- "figures/fig_netw"
##ACP on BOTH DMD and FFL #################################
non_plast <- cbind(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,c(5:12)],  read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")[,c(5:14)])
all_plast <- cbind(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,c(5:12)], read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")[,c(5:14)])

nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"
df.ACPe <- rbind(all_plast, non_plast)
df.ACPe <- df.ACPe[rowSums(df.ACPe[,1:18])>0,]
#For the next line to be usefull: FFlcolumns * number of FFL AND same for DMD
# df.ACPe[,1:18] <- df.ACPe[,1:18]*100/rowSums(df.ACPe[,1:18])

#normalize total number of loop again ?

pca_rese <- prcomp(df.ACPe[,1:18], scale. = FALSE)


g1 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'grey',alpha=0.7,
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("E. coli")+theme_bw()+ theme(legend.position = "none")+scale_color_manual(values = c("orange","cyan3"))

g1b <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.7,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("E. coli")+theme_bw()+ theme(legend.position = "none")+scale_color_manual(values = c("orange","cyan3"))

autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 2, y = 3,
         loadings = TRUE, loadings.colour = 'grey',alpha=0.7,
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("E. coli")+theme_bw()+ theme(legend.position = "none")+scale_color_manual(values = c("orange","cyan3"))

##
non_plast <- cbind(read.csv("scripts/data/full_netw_control_FFL.csv", sep = ",")[,c(4:11)],  read.csv("scripts/data/full_netw_control_DMD.csv", sep = ",")[,c(4:13)])
all_plast <- cbind(read.csv("scripts/data/full_netw_plast_FFL.csv", sep = ",")[,c(4:11)], read.csv("scripts/data/full_netw_plast_DMD.csv", sep = ",")[,c(4:13)])

nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")[,2]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
nbrloop_empl <- read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")[,2]
nbrloop_emnp <- read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")[,2]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"
#all_plast <- all_plast[sample(nrow(all_plast), 500), ]
df.ACPt <- rbind(all_plast, non_plast)
df.ACPt <- df.ACPt[rowSums(df.ACPt[,1:18])>0,]
# df.ACPt[,1:18] <- df.ACPt[,1:18]*100/rowSums(df.ACPt[,1:18])

pca_rest <- prcomp(df.ACPt[,1:18], scale. = FALSE)


g2 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'grey',alpha=0.7, 
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))

g2b <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.7, 
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))

autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 2, y = 3,
         loadings = TRUE, loadings.colour = 'grey',alpha=0.7, 
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))

  # ggplot(df.ACPt,aes(x = pca_rest$x[,1], y = pca_rest$x[,2], color = Genes))+
  #   geom_point() +
  #   stat_ellipse()+
  #   ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))
  


pdf(paste0(pdfname,"_ACP_all",".pdf"), width=10, height=10)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
grid.arrange(
  g1,g2,g1b,g2b,
  ncol = 2,
  nrow = 2,
  widths=c(0.8,1),
  clip = FALSE
)
dev.off()

################################################################################
##ACP on BOTH DMD and FFL and FBL###############################################
non_plast <- cbind(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")[,c(5:12)],
                   read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")[,c(5:14)],
                   read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")[,c(7:10)])
all_plast <- cbind(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")[,c(5:12)],
                   read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")[,c(5:14)],
                   read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")[,c(7:10)])

nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,3]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,3]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,3]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,3]
non_plast[,19:22] <- non_plast[,19:22]*nbrloop_emnp
all_plast[,19:22] <- all_plast[,19:22]*nbrloop_empl

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"
df.ACPe <- rbind(all_plast, non_plast)
df.ACPe <- df.ACPe[rowSums(df.ACPe[,1:22])>0,]
#For the next line to be usefull: FFlcolumns * number of FFL AND same for DMD
# df.ACPe[,1:22] <- df.ACPe[,1:22]*100/rowSums(df.ACPe[,1:22])

#normalize total number of loop again ?

pca_rese <- prcomp(df.ACPe[,1:22], scale. = FALSE)


g1 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.7,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("E. coli")+theme_bw()+ theme(legend.position = "none")+scale_color_manual(values = c("orange","cyan3"))

##
non_plast <- cbind(read.csv("scripts/data/full_netw_control_FFL.csv", sep = ",")[,c(4:11)],
                   read.csv("scripts/data/full_netw_control_DMD.csv", sep = ",")[,c(4:13)],
                   read.csv("scripts/data/full_netw_control_FBL.csv", sep = ",")[,c(6:9)])
all_plast <- cbind(read.csv("scripts/data/full_netw_plast_FFL.csv", sep = ",")[,c(4:11)],
                   read.csv("scripts/data/full_netw_plast_DMD.csv", sep = ",")[,c(4:13)],
                   read.csv("scripts/data/full_netw_plast_FBL.csv", sep = ",")[,c(6:9)])

nbrloop_empl <- subset(read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ","), Loop_number >=1)[,2]
nbrloop_emnp <- subset(read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ","), Loop_number >=1)[,2]
non_plast[,1:8] <- non_plast[,1:8]*nbrloop_emnp
all_plast[,1:8] <- all_plast[,1:8]*nbrloop_empl
nbrloop_empl <- subset(read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ","), Loop_number >=1)[,2]
nbrloop_emnp <- subset(read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ","), Loop_number >=1)[,2]
non_plast[,9:18] <- non_plast[,9:18]*nbrloop_emnp
all_plast[,9:18] <- all_plast[,9:18]*nbrloop_empl
nbrloop_empl <- subset(read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ","), FBL_number >=1)[,2]
nbrloop_emnp <- subset(read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ","), FBL_number >=1)[,2]
non_plast[,19:22] <- non_plast[,19:22]*nbrloop_emnp
all_plast[,19:22] <- all_plast[,19:22]*nbrloop_empl

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"
#all_plast <- all_plast[sample(nrow(all_plast), 500), ]
df.ACPt <- rbind(all_plast, non_plast)
df.ACPt <- df.ACPt[rowSums(df.ACPt[,1:22])>0,]
df.ACPt[,1:22] <- df.ACPt[,1:22]*100/rowSums(df.ACPt[,1:22])

pca_rest <- prcomp(df.ACPt[,1:22], scale. = FALSE)


g2 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.7, 
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+
  ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))

# ggplot(df.ACPt,aes(x = pca_rest$x[,1], y = pca_rest$x[,2], color = Genes))+
#   geom_point() +
#   stat_ellipse()+
#   ggtitle("Simulations")+theme_bw() +scale_color_manual(values = c("orange","cyan3"))



pdf(paste0(pdfname,"_ACP_all",".pdf"), width=10, height=5)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
grid.arrange(
  g1,g2,
  ncol = 2,
  nrow = 1,
  widths=c(0.8,1),
  clip = FALSE
)
dev.off()




##FBL###########################################################################
non_plast <- subset(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ","), FFL==1)[,c(2,5:12)]
all_plast <- subset(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ","), FFL==1)[,c(2,5:12)]

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"

df.ACPe <- rbind(non_plast,all_plast)


# res.pca <- PCA(df.ACPe[,2:9], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_rese <- prcomp(df.ACPe[,2:9], scale. = FALSE)

############
non_plast <- subset(read.csv("scripts/data/full_netw_control_FFL.csv", sep = ","), FFL==1)[,c(4:11)]
all_plast <- subset(read.csv("scripts/data/full_netw_plast_FFL.csv", sep = ","), FFL==1)[,c(4:11)]


non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"

df.ACPt <- rbind(non_plast,all_plast)

# 
# res.pca <- PCA(df.ACPt[,1:10], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_rest <- prcomp(df.ACPt[,1:8], scale. = FALSE)

g1 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+ ggtitle("E. coli")+theme_bw()#+ theme(legend.position = "none")
g2 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+ ggtitle("E. coli")+theme_bw()
g3 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey', alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+theme_bw()+ ggtitle("Simulations")
g4 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey', alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black")+theme_bw()+ ggtitle("Simulations")

pdf(paste0(pdfname,"_ACP_FFL",".pdf"), width=10, height=8)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
grid.arrange(
  g1,g2,g3,g4,
  ncol = 2,
  nrow = 2,
  clip = FALSE
)
dev.off()

pdf(paste0(pdfname,"_ACP_FFL_ecoli",".pdf"), width=5, height=4)
g1
dev.off()

pdf(paste0(pdfname,"_ACP_FFL_simul",".pdf"), width=5, height=4)
g3
dev.off()

pdf(paste0(pdfname,"_ACP_FFL_PC1",".pdf"), width=8, height=4)
grid.arrange(
  g1,g3,
  ncol = 2,
  nrow = 1,
  widths = c(0.6,1),
  clip = FALSE
)
dev.off()

##DMD###########################################################################
non_plast <- subset(read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ","), FFL==1)[,c(2,5:14)]
all_plast <- subset(read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ","), FFL==1)[,c(2,5:14)]


non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"

df.ACPe <- rbind(non_plast,all_plast)

# res.pca <- PCA(df.ACP[,2:11], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_rese <- prcomp(df.ACPe[,2:11], scale. = FALSE)


############
non_plast <- subset(read.csv("scripts/data/full_netw_control_DMD.csv", sep = ","), FFL==1)[,c(4:13)]
all_plast <- subset(read.csv("scripts/data/full_netw_plast_DMD.csv", sep = ","), FFL==1)[,c(4:13)]

non_plast$Genes <- "Non plastic"
all_plast$Genes <- "Plastic"

df.ACPt <- rbind(non_plast,all_plast)

# res.pca <- PCA(df.ACP[,1:10], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_rest <- prcomp(df.ACPt[,1:10], scale. = FALSE)


g1 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black", frame = TRUE, frame.type = 'norm')+theme_bw()+ ggtitle("E. coli")
g2 <- autoplot(pca_rese, data = df.ACPe, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey',alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black", frame = TRUE, frame.type = 'norm')+theme_bw()+ ggtitle("E. coli")
g3 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 1, y = 2,
               loadings = TRUE, loadings.colour = 'grey', alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black", frame = TRUE, frame.type = 'norm')+theme_bw()+ ggtitle("Simulations")
g4 <- autoplot(pca_rest, data = df.ACPt, colour = 'Genes', label = FALSE, x = 3, y = 4,
               loadings = TRUE, loadings.colour = 'grey', alpha=0.5,
               loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black", frame = TRUE, frame.type = 'norm')+theme_bw()+ ggtitle("Simulations")

pdf(paste0(pdfname,"_ACP_DMD",".pdf"), width=10, height=8)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
grid.arrange(
  g1,g2,g3,g4,
  ncol = 2,
  nrow = 2,
  clip = FALSE
)
dev.off()


pdf(paste0(pdfname,"_ACP_DMD_ecoli",".pdf"), width=5, height=4)
g1
dev.off()

pdf(paste0(pdfname,"_ACP_DMD_simul",".pdf"), width=5, height=4)
g3
dev.off()






