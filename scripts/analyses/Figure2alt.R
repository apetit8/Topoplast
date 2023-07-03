source("scripts/functions/functions.R")
pdfname <- "figures/fig_alt"
#Plot combining both empirical data and theoretical data
################################################################################
##FFL motifs####################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
theory <- read.csv("scripts/data/10g_FFL.csv", sep = ",")[c(1,5),1:11]

df_ffl <- as.data.frame(rbind(theory[1,4:11]*100/theory[1,2],theory[2,4:11]*100/theory[2,2], colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                              colSums(non_plast[,5:12])*100/sum(non_plast[,3]) ))
rownames(df_ffl) <- c("Plastic_Prediction", "Nonplastic_Prediction","Plastic_Empiric","NonPlastic_Empiric")

df <- rbind(df_ffl[1,]/df_ffl[2,], df_ffl[3,]/df_ffl[4,])
df$group <- c(1,2)

pdf(paste0(pdfname,"_FFL_motifs",".pdf"), width=5, height=5)
plot(df$group, df$C1, frame=FALSE, type = "b", pch = 19,
     col = "darkseagreen", xlab = "x", ylab = "y", 
     lty = 1, lwd = 1, ylim = c(0,max(df)))
abline(h=1)
lines(df$group, df$C2, pch = 19, col = "yellowgreen", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$C3, pch = 19, col = "dodgerblue", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$C4, pch = 19, col = "deepskyblue3", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$I1, pch = 19, col = "indianred1", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$I2, pch = 19, col = "lightpink", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$I3, pch = 19, col = "orange", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$I4, pch = 19, col = "lightgoldenrod1", type = "b", 
      lty = 1, lwd = 1)
legend("topright", box.lty=0,  bg="transparent", fill = c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend = c("C1","C2","C3","C4","I1","I2","I3","I4"))
dev.off()


plot(as.data.frame(t(df))$Plastic_Prediction-1, as.data.frame(t(df))$Plastic_Empiric-1, frame=FALSE, type = "p", pch = 19,
     col = "yellowgreen", xlab = "x", ylab = "y",
     lty = 1, lwd = 1)
abline(lm(Plastic_Prediction-1 ~ Plastic_Empiric-1, data = as.data.frame(t(df))), col = "black")


plot(t(df_ffl)[,c(1,3)], frame=FALSE, pch = 19, col = c("grey"), ylab = "Empiric", xlab="Prediction", ylim = c(0, max(df_ffl)), xlim = c(0, max(df_ffl)))
abline(lm(Plastic_Prediction ~ Plastic_Empiric, data = as.data.frame(t(df_ffl))), col = "grey")
lines(t(df_ffl)[,c(2,4)], pch = 19, type = "p", col = c("black"))
abline(lm(Nonplastic_Prediction ~ NonPlastic_Empiric, data = as.data.frame(t(df_ffl))), col = "black")
legend("topright", box.lty=0,  bg="transparent", fill = c("grey","black"), legend = c("Plastic","Non-Plastic"))


pdf(paste0(pdfname,"_FFL_lm",".pdf"), width=5, height=4.5)
par(mar = c(4,4, 1,1))
plot(NULL, NULL, type = "p", pch = 19,
     col = "darkseagreen", xlab = "Plastic/Non-plastic Prediction", ylab = "Plastic/Non-plastic Empiric", lty = 1, lwd = 1,
     ylim = c(min(as.data.frame(t(df))[,2]-1), max(as.data.frame(t(df))[,2]-1)),
     xlim = c(min(as.data.frame(t(df))[,1]-1), max(as.data.frame(t(df))[,1]-1)))
polygon(x=c(min(as.data.frame(t(df))[,1]-1)-0.15,min(as.data.frame(t(df))[,1]-1)-0.15,0,0),
        y=c(min(as.data.frame(t(df))[,2]-1)-0.075,0,0,min(as.data.frame(t(df))[,2]-1)-0.075),
        col="honeydew2", border=NA)
polygon(x=c(0,0,max(as.data.frame(t(df))[,1]-1)+0.15,max(as.data.frame(t(df))[,1]-1)+0.15),
        y=c(0,max(as.data.frame(t(df))[,2]-1)+0.075,max(as.data.frame(t(df))[,2]-1)+0.075,0),
        col="honeydew2", border=NA)
lines(as.data.frame(t(df))[1,1]-1, as.data.frame(t(df))[1,2]-1, pch = 22, cex=2, bg = "darkseagreen", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[2,1]-1, as.data.frame(t(df))[2,2]-1, pch = 22, cex=2, bg = "yellowgreen", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[3,1]-1, as.data.frame(t(df))[3,2]-1, pch = 22, cex=2, bg = "dodgerblue", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[4,1]-1, as.data.frame(t(df))[4,2]-1, pch = 22, cex=2, bg = "deepskyblue3", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[5,1]-1, as.data.frame(t(df))[5,2]-1, pch = 22, cex=2, bg = "indianred1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[6,1]-1, as.data.frame(t(df))[6,2]-1, pch = 22, cex=2, bg = "lightpink", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[7,1]-1, as.data.frame(t(df))[7,2]-1, pch = 22, cex=2, bg = "orange", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[8,1]-1, as.data.frame(t(df))[8,2]-1, pch = 22, cex=2, bg = "lightgoldenrod1", type = "p", lty = 1, lwd = 1)

abline(lm(Plastic_Prediction-1 ~ Plastic_Empiric-1, data = as.data.frame(t(df))), col = "black")
# legend("topright", box.lty=0,  bg="transparent", fill = c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
#        legend = colnames(df)[1:8])
dev.off()



################################################################################
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
theory <- read.csv("scripts/data/10g_diamond.csv", sep = ",")[c(1,5),1:13]

df1 <- as.data.frame(rbind(theory[1,4:13]*100/theory[1,2], theory[2,4:13]*100/theory[2,2], 
                           colSums(all_plast[,5:14])*100/sum(all_plast[,3]), colSums(non_plast[,5:14])*100/sum(non_plast[,3]) ))
rownames(df1) <- c("Plastic_Prediction", "Nonplastic_Prediction","Plastic_Empiric","NonPlastic_Empiric")

df <- rbind(df1[1,]/df1[2,], df1[3,]/df1[4,])
df$group <- c(1,2)

pdf(paste0(pdfname,"_DMD_motifs",".pdf"), width=5, height=5)
plot(df$group,df$PP, frame=FALSE, type = "b", pch = 19,
     col = "olivedrab1", xlab = "x", ylab = "y", 
     lty = 1, lwd = 1, ylim = c(0,max(df)))
abline(h=1)
lines(df$group, df$PM, pch = 19, col = "palegreen", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$PN, pch = 19, col = "mediumseagreen", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$NP, pch = 19, col = "plum", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$NM, pch = 19, col = "darkorchid1", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$NN, pch = 19, col = "plum1", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$MP, pch = 19, col = "lightsalmon1", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$MM1, pch = 19, col = "indianred1", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$MM2, pch = 19, col = "lightpink", type = "b", 
      lty = 1, lwd = 1)
lines(df$group, df$MN, pch = 19, col = "chocolate1", type = "b", 
      lty = 1, lwd = 1)
legend("topright", box.lty=0,  bg="transparent", fill = c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1",
               "plum1","lightsalmon1","indianred1","lightpink","chocolate1"), legend = c("PP","PM","PP","NP","NM","NN","MP","MM1","MM2","MN"))
dev.off()
#BUG : CALCULATING %

plot(t(df1)[,c(1,3)], frame=FALSE, pch = 19, col = c("grey"), ylab = "Empiric", xlab="Prediction")
abline(lm(Plastic_Prediction ~ Plastic_Empiric, data = as.data.frame(t(df1))), col = "grey", ylim = c(0, max(df1)), xlim = c(0, max(df1)))
lines(t(df1)[,c(2,4)], pch = 19, type = "p", col = c("black"))
abline(lm(Nonplastic_Prediction ~ NonPlastic_Empiric, data = as.data.frame(t(df1))), col = "black")
legend("topright", box.lty=0,  bg="transparent", fill = c("grey","black"), legend = c("Plastic","Non-Plastic"))



pdf(paste0(pdfname,"_DMD_lm",".pdf"), width=5, height=4.5)
par(mar = c(4,4, 1,1))
plot(NULL, NULL, type = "p",
     col = "olivedrab1", xlab = "Plastic/Non-plastic Prediction", ylab = "Plastic/Non-plastic Empiric", lty = 1, lwd = 1,
     ylim = c(min(as.data.frame(t(df))[,2]-1), max(as.data.frame(t(df))[,2]-1)),
     xlim = c(min(as.data.frame(t(df))[,1]-1), max(as.data.frame(t(df))[,1]-1)))
polygon(x=c(min(as.data.frame(t(df))[,1]-1)-0.2,min(as.data.frame(t(df))[,1]-1)-0.2,0,0),
        y=c(min(as.data.frame(t(df))[,2]-1)-0.045,0,0,min(as.data.frame(t(df))[,2]-1)-0.045),
        col="honeydew2", border=NA)
polygon(x=c(0,0,max(as.data.frame(t(df))[,1]-1)+0.2,max(as.data.frame(t(df))[,1]-1)+0.2),
        y=c(0,max(as.data.frame(t(df))[,2]-1)+0.05,max(as.data.frame(t(df))[,2]-1)+0.05,0),
        col="honeydew2", border=NA)
lines(as.data.frame(t(df))[1,1]-1, (as.data.frame(t(df))[1,2])-1, pch = 22, cex=2, bg = "olivedrab1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[2,1]-1, as.data.frame(t(df))[2,2]-1, pch = 22, cex=2, bg = "palegreen", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[3,1]-1, as.data.frame(t(df))[3,2]-1, pch = 22, cex=2, bg = "mediumseagreen", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[4,1]-1, as.data.frame(t(df))[4,2]-1, pch = 22, cex=2, bg = "plum", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[5,1]-1, as.data.frame(t(df))[5,2]-1, pch = 22, cex=2, bg = "darkorchid1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[6,1]-1, as.data.frame(t(df))[6,2]-1, pch = 22, cex=2, bg = "plum1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[7,1]-1, as.data.frame(t(df))[7,2]-1, pch = 22, cex=2, bg = "lightsalmon1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[8,1]-1, as.data.frame(t(df))[8,2]-1, pch = 22, cex=2, bg = "indianred1", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[9,1]-1, as.data.frame(t(df))[9,2]-1, pch = 22, cex=2, bg = "lightpink", type = "p", lty = 1, lwd = 1)
lines(as.data.frame(t(df))[10,1]-1, as.data.frame(t(df))[10,2]-1, pch = 22, cex=2, bg = "chocolate1", type = "p", lty = 1, lwd = 1)

abline(lm(Plastic_Prediction-1 ~ Plastic_Empiric-1, data = as.data.frame(t(df))), col = "black")
# legend("topright", box.lty=0,  bg="transparent", fill = c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1",
#                                                           "plum1","lightsalmon1","indianred1","lightpink","chocolate1"), legend = c("PP","PM","PN","NP","NM","NN","MP","MM2","MM1","MN"))
dev.off()









