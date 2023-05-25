source("scripts/functions/functions.R")
#####################
sims.dirs1 <- list.dirs("simul/4g", recursive = TRUE)
sims.dirs2 <- list.dirs("simul/2g", recursive = TRUE)
sims.dirs3 <- list.dirs("simul/2g_selfreg", recursive = TRUE)
pdfname <- "figures/4g"
gen <- 3000
#####################
#DATA
##################
df.4 <- df.simul(sims.dirs1, all.gen = TRUE)
df.4$envir <- str_split(df.4$data.dir, "/", n=8, simplify = TRUE)[,3]

df.2 <- df.simul(sims.dirs2, all.gen = TRUE)
df.2$envir <- str_split(df.2$data.dir, "/", n=8, simplify = TRUE)[,3]

df.2sr <- df.simul(sims.dirs3, all.gen = TRUE, size = 1000)
df.2sr$envir <- str_split(df.2sr$data.dir, "/", n=8, simplify = TRUE)[,3]

df.2lg <- subset(df.2, Gen==gen & envir %in% c("Correlated","Anticorrelated"))
df.2srlg<- subset(df.2sr, Gen==gen & envir %in% c("Correlated","Anticorrelated"))
df.3lg <- subset(df.3, Gen==gen & envir %in% c("Correlated","Anticorrelated"))
df.4lg <- subset(df.4, Gen==gen & envir %in% c("Correlated","Anticorrelated"))
df.10lg <- subset(df.10, Gen==gen & envir %in% c("Correlated","Anticorrelated"))

df.2lg <- cbind(Genes=rep("2g", nrow(df.2lg)), df.2lg )
df.2srlg <- cbind(Genes=rep("2sr", nrow(df.2srlg)), df.2srlg )
df.3lg <- cbind(Genes=rep("3g", nrow(df.3lg)), df.3lg )
df.4lg <- cbind(Genes=rep("4g", nrow(df.4lg)), df.4lg )
df.10lg <- cbind(Genes=rep("10g", nrow(df.10lg)), df.10lg )

df.g <- rbind(df.2lg[,1:7], df.2srlg[,1:7], df.3lg[,1:7], df.4lg[,1:7], df.10lg[,1:7])
df.g$Genes <- as.factor(df.g$Genes)

ggplot(df.g, aes(x = Genes, y=Fitness))+ ylim(0,1)+
  geom_boxplot()+theme_bw()
  
pdf("figures/fitness_gene_nbr.pdf", width=6.5, height=6)
ggplot(df.g, aes(x = Genes, y=Fitness))+ ylim(0,1)+
  geom_boxplot()+theme_bw()
dev.off()



