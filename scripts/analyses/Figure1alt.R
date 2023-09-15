source("scripts/functions/functions.R")
library(viridis)
pdfname <- "figures/fig_1"
#Plot combining both empirical data and theoretical data
################################################################################
pval <- 0.05
################################################################################
#Function to rbind tables, modified from https://www.pogol.net/rbind-fill-for-1-dimensional-tables-in-r
rbind1dtable=function(tab1,tab2,tab3,tab4,fail=0){
  sapply(unique(c(names(tab1),names(tab2),names(tab3),names(tab4))),function(n){
    t1=tab1[which(names(tab1)==n)]
    t2=tab2[which(names(tab2)==n)]
    t3=tab3[which(names(tab3)==n)]
    t4=tab4[which(names(tab4)==n)]
    c(ifelse(length(t1)==0,fail,t1),ifelse(length(t2)==0,fail,t2),ifelse(length(t3)==0,fail,t3),ifelse(length(t4)==0,fail,t4))
  })
}
################################################################################
#FFL###############################
nbrloop_thpl <- read.csv("scripts/data/full_netw_Pl_nbrFFL.csv", sep = ",")
nbrloop_thnp <- read.csv("scripts/data/full_netw_NP_nbrFFL.csv", sep = ",")
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]

tt <- as.data.frame <- rbind1dtable(table(nbrloop_emnp[,2]), table(nbrloop_empl[,2]), 
                   table(nbrloop_thnp[,2]), table(nbrloop_thpl[,2]) )
rownames(tt) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
tt <- tt[,order(as.numeric(colnames(tt)), method="radix")]
for(i in 1:nrow(tt)) tt[i,] <- tt[i,]*100/sum(tt[i,])


pdf(paste0(pdfname,"_alt_FFL",".pdf"), width=5, height=4)
par(mar=c(3.3,3.3,2,2.8), xpd=TRUE)
par(lty = 0)
barplot(t(tt), main = "Number of FFL per gene", col=c("grey",rev(viridis(ncol(tt)-1))), space=c(0.0,0.1,0.35,0.1), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
title(xlab = "Genes", ylab = "Frequency", cex.lab=1.2, mgp=c(2.2,2.2,0))
text(1.05, 103, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
lgd_ = rep(NA, ncol(tt))
lgd_[c(1, round(ncol(tt)/4), round(ncol(tt)/2), round(ncol(tt)*3/4), ncol(tt))] = c(0, colnames(tt)[round(ncol(tt)/4)], colnames(tt)[round(ncol(tt)/2)], colnames(tt)[round(ncol(tt)*3/4)], colnames(tt)[ncol(tt)])
legend("topright", inset=c(-0.15,0.1), legend = lgd_, fill = c("grey",rev(viridis(ncol(tt)-1))), border = NA, x.intersp = 0.15,
       y.intersp = 0.3, cex = 1.1, pt.cex = 1)
dev.off()

print(paste0("E. coli non plastic mean number of FFL : ",mean(subset(nbrloop_emnp, Loop_number!=0)$Loop_number)," ; ",
             "E. coli plastic mean number of FFL : ",mean(subset(nbrloop_empl, Loop_number!=0)$Loop_number)," ; ",
             "Simulations non plastic mean number of FFL : ",mean(subset(nbrloop_thnp, Loop_number!=0)$Loop_number)," ; ",
             "simulations plastic mean number of FFL : ",mean(subset(nbrloop_thpl, Loop_number!=0)$Loop_number)," ; "
             ))

#DMD###############################
nbrloop_thpl <- read.csv("scripts/data/full_netw_Pl_nbrDMD.csv", sep = ",")
nbrloop_thnp <- read.csv("scripts/data/full_netw_NP_nbrDMD.csv", sep = ",")
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ",")[,c(1,3,4,5)]

tt <- as.data.frame <- rbind1dtable(table(nbrloop_emnp[,2]), table(nbrloop_empl[,2]), 
                                    table(nbrloop_thnp[,2]), table(nbrloop_thpl[,2]) )
rownames(tt) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
tt <- tt[,order(as.numeric(colnames(tt)), method="radix")]
for(i in 1:nrow(tt)) tt[i,] <- tt[i,]*100/sum(tt[i,])


pdf(paste0(pdfname,"_alt_DMD",".pdf"), width=5, height=4)
par(mar=c(3.3,3.3,2,2.8), xpd=TRUE)
par(lty = 0)
barplot(t(tt), main = "Number of DMD per gene", col=c("grey",rev(viridis(ncol(tt)-1))), space=c(0.0,0.1,0.35,0.1), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
title(xlab = "Genes", ylab = "Frequency", cex.lab=1.2, mgp=c(2.2,2.2,0))
text(1.05, 103, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
lgd_ = rep(NA, ncol(tt))
lgd_[c(1, round(ncol(tt)/4), round(ncol(tt)/2), round(ncol(tt)*3/4), ncol(tt))] = c(0, colnames(tt)[round(ncol(tt)/4)], colnames(tt)[round(ncol(tt)/2)], colnames(tt)[round(ncol(tt)*3/4)], colnames(tt)[ncol(tt)])
legend("topright", inset=c(-0.15,0.1), legend = lgd_, fill = c("grey",rev(viridis(ncol(tt)-1))), border = NA, x.intersp = 0.15,
       y.intersp = 0.13, cex = 1.1, pt.cex = 1)
dev.off()

print(paste0("E. coli non plastic mean number of FFL : ",mean(subset(nbrloop_emnp, Loop_number!=0)$Loop_number)," ; ",
             "E. coli plastic mean number of FFL : ",mean(subset(nbrloop_empl, Loop_number!=0)$Loop_number)," ; ",
             "Simulations non plastic mean number of FFL : ",mean(subset(nbrloop_thnp, Loop_number!=0)$Loop_number)," ; ",
             "simulations plastic mean number of FFL : ",mean(subset(nbrloop_thpl, Loop_number!=0)$Loop_number)," ; "
))

#FBL###############################
nbrloop_thpl <- read.csv("scripts/data/full_netw_Pl_nbrFBL.csv", sep = ",")
nbrloop_thnp <- read.csv("scripts/data/full_netw_NP_nbrFBL.csv", sep = ",")
nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]
nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ",")[,c(1,3,4,5)]

tt <- as.data.frame <- rbind1dtable(table(nbrloop_emnp[,2]), table(nbrloop_empl[,2]), 
                                    table(nbrloop_thnp[,2]), table(nbrloop_thpl[,2]) )
rownames(tt) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
tt <- tt[,order(as.numeric(colnames(tt)), method="radix")]
for(i in 1:nrow(tt)) tt[i,] <- tt[i,]*100/sum(tt[i,])


pdf(paste0(pdfname,"_alt_FBL",".pdf"), width=5, height=4)
par(mar=c(3.3,3.3,2,2.8), xpd=TRUE)
par(lty = 0)
barplot(t(tt), main = "Number of FBL per gene", col=c("grey",rev(viridis(ncol(tt)-1))), space=c(0.0,0.1,0.35,0.1), 
        args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
title(xlab = "Genes", ylab = "Frequency", cex.lab=1.2, mgp=c(2.2,2.2,0))
text(1.05, 103, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
lgd_ = rep(NA, ncol(tt))
lgd_[c(1, round(ncol(tt)/4), round(ncol(tt)/2), round(ncol(tt)*3/4), ncol(tt))] = c(0, colnames(tt)[round(ncol(tt)/4)], colnames(tt)[round(ncol(tt)/2)], colnames(tt)[round(ncol(tt)*3/4)], colnames(tt)[ncol(tt)])
legend("topright", inset=c(-0.15,0.1), legend = lgd_, fill = c("grey",rev(viridis(ncol(tt)-1))), border = NA, x.intersp = 0.15,
       y.intersp = 0.06, cex = 1.1, pt.cex = 1)
dev.off()

print(paste0("E. coli non plastic mean number of FFL : ",mean(subset(nbrloop_emnp, FBL_number!=0)$FBL_number)," ; ",
             "E. coli plastic mean number of FFL : ",mean(subset(nbrloop_empl, FBL_number!=0)$FBL_number)," ; ",
             "Simulations non plastic mean number of FFL : ",mean(subset(nbrloop_thnp, FBL_number!=0)$FBL_number)," ; ",
             "simulations plastic mean number of FFL : ",mean(subset(nbrloop_thpl, FBL_number!=0)$FBL_number)," ; "
))





# ##TEST##########################################################################
# nbrloop_thpl <- read.csv("scripts/data/full_a0-7_Pl_nbrDMD.csv", sep = ",")
# nbrloop_thnp <- read.csv("scripts/data/full_a0-7_NP_nbrDMD.csv", sep = ",")
# nbrloop_empl <- read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
# nbrloop_emnp <- read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ",")[,c(1,3,4,5)]
# 
# tt <- as.data.frame <- rbind1dtable(table(nbrloop_emnp[,2]), table(nbrloop_empl[,2]), 
#                                     table(nbrloop_thnp[,2]), table(nbrloop_thpl[,2]) )
# rownames(tt) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
# tt <- tt[,order(as.numeric(colnames(tt)), method="radix")]
# for(i in 1:nrow(tt)) tt[i,] <- tt[i,]*100/sum(tt[i,])
# 
# 
# par(mar=c(3.3,3.3,2,2.8), xpd=TRUE)
# par(lty = 0)
# barplot(t(tt), main = "Number of FFL per gene", col=c("grey",rev(viridis(ncol(tt)-1))), space=c(0.0,0.1,0.35,0.1), 
#         args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
# title(xlab = "Genes", ylab = "Frequency", cex.lab=1.2, mgp=c(2.2,2.2,0))
# text(1.05, 103, paste0( ifelse(t.test(nbrloop_empl[,2], nbrloop_emnp[,2], var.equal=F)$p.value <=0.01, "***", ""  )), cex=2)
# lgd_ = rep(NA, ncol(tt))
# lgd_[c(1, round(ncol(tt)/4), round(ncol(tt)/2), round(ncol(tt)*3/4), ncol(tt))] = c(0, colnames(tt)[round(ncol(tt)/4)], colnames(tt)[round(ncol(tt)/2)], colnames(tt)[round(ncol(tt)*3/4)], colnames(tt)[ncol(tt)])
# legend("topright", inset=c(-0.15,0.2), legend = lgd_, fill = c("grey",rev(viridis(ncol(tt)-1))), border = NA, x.intersp = 0.15,
#        y.intersp = 0.3, cex = 1.1, pt.cex = 1)
# 
# 

