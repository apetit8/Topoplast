source("scripts/functions/functions.R")
#####################
sims.dirs1 <- list.dirs("simul/4g", recursive = TRUE)
genes <- 4
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
pdfname_suffixe <- ""
#####################
#DATA
##################
df.4 <- df.simul(sims.dirs1, all.gen = TRUE)
df.4$envir <- str_split(df.4$data.dir, "/", n=8, simplify = TRUE)[,3]
df.4$anc_id <- str_split(str_split(df.4$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]

#Add Slope
for (i in 1:nrow(df.4)) {
  W <- t(matrix(as.numeric(df.4[i,7:(genes^2+6)]), ncol = genes))
  # basal <- if(grepl("Down", df.4[i,23])) 0.8 else if(grepl("Up", df.4[i,23])) 0.2 else 0.5
  df.4[i,(genes^2+9)] <- getSlope.ALR(W=W, n.env=21, target.gene=target, min=min, max=max)
}
#####################
sims.dirs1 <- list.dirs("simul/3g", recursive = TRUE)
genes <- 3
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
#####################
#DATA
##################
df.3 <- df.simul(sims.dirs1, all.gen = TRUE)
df.3$envir <- str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,3]
df.3$anc_id <- str_split(str_split(df.3$data.dir, "/", n=8, simplify = TRUE)[,4], "-", simplify = TRUE)[,2]


#Add Slope
for (i in 1:nrow(df.3)) {
  W <- t(matrix(as.numeric(df.3[i,7:(genes^2+6)]), ncol = genes))
  df.3[i,(genes^2+9)] <- getSlope.ALR(W=W, n.env=21, target.gene=target, min=min, max=max)
}

#Figure#########################################################################


pdf(paste0("figures/fig_RN", pdfname_suffixe, ".pdf"), width=15, height=7)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
#Plot 1
par(mar=c(2.8, 3, 2.5, 0.5), mgp = c(1.75, 0.75, 0), las=0)
plot(df.3$Gen, NULL, type="n", xlim = c(min(df.3$Gen), 9500),ylim = c(-1.2,1.1),
     main="3 genes", ylab = "Reaction Norm", xlab = "Generation")
#Trace gene category by gene category, for both df
bymodel <- by(df.3$V18, list(df.3$Gen, df.3$envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_Down"], col="darkblue", lwd = 3, lty=3)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_Up"], col="darkblue", lwd = 3, lty=2)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_UD"], col="darkblue", lwd = 3, lty=1)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_Down"], col="maroon2", lwd = 3, lty=3)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_Up"], col="maroon2", lwd = 3, lty=2)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_UD"], col="maroon2", lwd = 3, lty=1)
lines(as.numeric(rownames(bymodel)), bymodel[,"Noncorrelated_UD"], col="yellowgreen", lwd = 3, lty=1)


legend("topleft", lty=c(1,1,1,1,2,3), box.lty=0,  bg="transparent", col=c("darkblue","maroon2","yellowgreen", "black", "black", "black"),
       legend=c("Anticorrelated", "Correlated", "Noncorrelated", "Down and Up regulated","Up regulated", "Down regulated"), lwd = 3)

plot(df.4$Gen, NULL, type="n", xlim = c(min(df.4$Gen), 9500),ylim = c(-1.2,1.1),
     main="4 genes", ylab = "Reaction Norm", xlab = "Generation")
#Trace gene category by gene category, for both df
bymodel <- by(df.4$V25, list(df.4$Gen, df.4$envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_Down"], col="darkblue", lwd = 3, lty=3)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_Up"], col="darkblue", lwd = 3, lty=2)
lines(as.numeric(rownames(bymodel)), bymodel[,"Anticorrelated_UD"], col="darkblue", lwd = 3, lty=1)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_Down"], col="maroon2", lwd = 3, lty=3)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_Up"], col="maroon2", lwd = 3, lty=2)
lines(as.numeric(rownames(bymodel)), bymodel[,"Correlated_UD"], col="maroon2", lwd = 3, lty=1)
lines(as.numeric(rownames(bymodel)), bymodel[,"Noncorrelated_UD"], col="yellowgreen", lwd = 3, lty=1)

dev.off()


pdf(paste0("figures/fig_RN_envir", pdfname_suffixe, ".pdf"), width=7, height=4)
ggplot(subset(df.3, envir!="Noncorrelated_UD" & Gen >=2500), aes(Optimum, P_mean_2, col=envir))+
  geom_point()+ggtitle("3 genes")+
  theme_bw()+ geom_smooth(se = FALSE, method = lm)
ggplot(subset(df.4, envir!="Noncorrelated_UD"& Gen >=2500), aes(Optimum, P_mean_2, col=envir))+
  geom_point()+theme_bw()+ggtitle("4 genes")+geom_smooth(se = FALSE, method = lm)
dev.off()


ggplot(subset(df.4, envir!="Noncorrelated_UD"& Gen >=2500), aes(Gen, Fitness, col=envir))+
  geom_point()+theme_bw()+ggtitle("4 genes")+geom_smooth(se = FALSE, method = lm)


