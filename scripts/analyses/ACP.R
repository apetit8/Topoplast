library(FactoMineR)
library(factoextra)
library(gridExtra)
library(ggfortify)
library(ggplot2)
##FBL###########################################################################
non_plast <- subset(read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ","), FFL==1)[,c(2,5:12)]
all_plast <- subset(read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ","), FFL==1)[,c(2,5:12)]

non_plast$Plast <- "NP"
all_plast$Plast <- "Pl"

df.ACP <- rbind(non_plast,all_plast)


# res.pca <- PCA(df.ACP[,2:9], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_res <- prcomp(df.ACP[,2:9], scale. = FALSE)
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'black',alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 3, y = 4,
         loadings = TRUE, loadings.colour = 'black',alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()

################################################################################
non_plast <- subset(read.csv("scripts/data/full_20k_0-01_control_FFL.csv", sep = ","), FFL==1)[,c(4:11)]
all_plast <- subset(read.csv("scripts/data/full_20k_0-01_plast_FFL.csv", sep = ","), FFL==1)[,c(4:11)]


non_plast$Plast <- "NP"
all_plast$Plast <- "Pl"

df.ACP <- rbind(non_plast,all_plast)

# 
# res.pca <- PCA(df.ACP[,1:10], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_res <- prcomp(df.ACP[,1:8], scale. = FALSE)
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'black', alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 3, y = 4,
         loadings = TRUE, loadings.colour = 'black', alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()

##DMD###########################################################################
non_plast <- subset(read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ","), FFL==1)[,c(2,5:14)]
all_plast <- subset(read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ","), FFL==1)[,c(2,5:14)]


non_plast$Plast <- "NP"
all_plast$Plast <- "Pl"

df.ACP <- rbind(non_plast,all_plast)

# res.pca <- PCA(df.ACP[,2:11], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_res <- prcomp(df.ACP[,2:11], scale. = FALSE)
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'black',alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 3, y = 4,
         loadings = TRUE, loadings.colour = 'black',alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()

################################################################################
non_plast <- subset(read.csv("scripts/data/full_20k_0-01_control_DMD.csv", sep = ","), FFL==1)[,c(4:13)]
all_plast <- subset(read.csv("scripts/data/full_20k_0-01_plast_DMD.csv", sep = ","), FFL==1)[,c(4:13)]


non_plast$Plast <- "NP"
all_plast$Plast <- "Pl"

df.ACP <- rbind(non_plast,all_plast)


# res.pca <- PCA(df.ACP[,1:10], scale.unit=FALSE, graph = FALSE)
# eig.val <- res.pca$eig
# cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

pca_res <- prcomp(df.ACP[,1:10], scale. = FALSE)
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 1, y = 2,
         loadings = TRUE, loadings.colour = 'black', alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
autoplot(pca_res, data = df.ACP, colour = 'Plast', label = FALSE, x = 3, y = 4,
         loadings = TRUE, loadings.colour = 'black', alpha=0.5,
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
