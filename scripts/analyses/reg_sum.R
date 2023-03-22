#Analyses ?
source("functions/functions.R")
# source("functions/plotnet.R")
#####################
sims.dirs1 <- list.dirs("../simul/Evoplast/burning/Correlated", recursive = TRUE)
sims.dirs2 <- list.dirs("../simul/Evoplast/new_envir/B_Correlated", recursive = TRUE)
sims.dirs2 <- sims.dirs2[grep('N_Noncorrelated', sims.dirs2)]
sims.dirs3 <- list.dirs("../simul/Evoplast/burning/Noncorrelated", recursive = TRUE)
sims.dirs4 <- list.dirs("../simul/Evoplast/new_envir/B_Noncorrelated", recursive = TRUE)
sims.dirs4 <- sims.dirs4[grep('N_Correlated', sims.dirs4)]
#####################

df.part1 <- df.simul(c(sims.dirs1, sims.dirs3), all.gen = TRUE)
df.part1$n_envir <- str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,5]
df.part1$b_envir <- paste0("B_",str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,5])
df.part1$anc_id <- str_split(str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,6], "-", simplify = TRUE)[,2]

df.part2 <- df.simul(c(sims.dirs2, sims.dirs4), all.gen = TRUE)
df.part2$n_envir <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,7]
df.part2$b_envir <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,5]
df.part2$anc_id <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,6]
df.part2$Gen <- df.part2$Gen + max(df.part1$Gen)

df.gensum <- rbind(df.part1, df.part2)

#Sum of reg on Plastic gene (=2) and Stable gene (=3)
df.gensum$V110 <- rowSums(cbind(abs(df.gensum$W_4_2), abs(df.gensum$W_5_2), abs(df.gensum$W_6_2), 
                       abs(df.gensum$W_7_2), abs(df.gensum$W_8_2), abs(df.gensum$W_9_2),  abs(df.gensum$W_10_2)))
df.gensum$V111 <- rowSums(cbind(abs(df.gensum$W_4_3), abs(df.gensum$W_5_3), abs(df.gensum$W_6_3), 
                       abs(df.gensum$W_7_3), abs(df.gensum$W_8_3), abs(df.gensum$W_9_3),  abs(df.gensum$W_10_3)))



cairo_pdf("../figures/fig_sum_gen.pdf", width=7, height=5)
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
#Plot 1
# par(mar=c(3, 3, 3, 0.5), mgp = c(1.75, 0.75, 0), las=0)
plot(df.gensum$Gen, NULL, type="n", xlim = c(min(df.gensum$Gen), max(df.gensum$Gen)),ylim = c(0,3),
     ylab = "Sum of regulations toward selected genes ", xlab = "Generation")
abline(v = 10000, lwd = 2, lty = 3)
#Trace cat by cat, indicating color manually each time (3 lines)
bymodel <- by(df.gensum$V111, list(df.gensum$Gen, df.gensum$b_envir), FUN=mean)
# lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="yellowgreen", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="lightblue2", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="thistle", lwd = 3)
bymodel <- by(df.gensum$V110, list(df.gensum$Gen, df.gensum$b_envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="maroon2", lwd = 3)
# legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("yellowgreen","darkblue","maroon2"),
#        legend=c(paste0("Stable"), paste0("Plastic to Stable"), paste0("Stable to Plastic")))
legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("darkblue","lightblue2","maroon2","thistle"),
       legend=c("Changing 1: Plastic to Stable","Constant with Changing 1","Changing 2: Stable to Plastic","Constant with Changing 2"))
dev.off()

# plot(subset(df.gensum, Gen !=0)$V110, df.gensum2$Sd_ind_2)
# plot(subset(df.gensum, Gen !=0)$V111, df.gensum2$Sd_ind_3)

df.test <- data.frame()
test <- subset(subset(df.gensum, b_envir=="B_Correlated"), Gen !=0)
for (i in 1:nrow(test)) {
  W <- t(matrix(as.numeric(test[i,7:106]), ncol = 10))
  df.test[i,1] <- getSlope.ALR(W=W, n.env=21, target.gene=2)
}

plot(df.test[,1], subset(df.gensum2, b_envir=="B_Correlated")$Sd_ind_2)

df.test[,2] <- rowSums(cbind(abs(test$W_1_4), abs(test$W_1_5), abs(test$W_1_6), 
                                abs(test$W_1_7), abs(test$W_1_8), abs(test$W_1_9)))

plot(df.test[,1], df.test[,2])

plot(test$Gen, df.test[,2])

################################################################################

#Sum of reg on Plastic gene (=2) and Stable gene (=3)
df.gensum$V112 <- rowSums(cbind(abs(df.gensum$W_2_4), abs(df.gensum$W_2_5), abs(df.gensum$W_2_6), 
                                abs(df.gensum$W_2_7), abs(df.gensum$W_2_8), abs(df.gensum$W_2_9),  abs(df.gensum$W_2_10)))
df.gensum$V113 <- rowSums(cbind(abs(df.gensum$W_3_4), abs(df.gensum$W_3_5), abs(df.gensum$W_3_6), 
                                abs(df.gensum$W_3_7), abs(df.gensum$W_3_8), abs(df.gensum$W_3_9),  abs(df.gensum$W_3_10)))



cairo_pdf("../figures/fig_sum_gen_out.pdf", width=7, height=5)
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
#Plot 1
# par(mar=c(3, 3, 3, 0.5), mgp = c(1.75, 0.75, 0), las=0)
plot(df.gensum$Gen, NULL, type="n", xlim = c(min(df.gensum$Gen), max(df.gensum$Gen)),ylim = c(0,3),
     ylab = "Sum of regulations from selected genes", xlab = "Generation")
abline(v = 10000, lwd = 2, lty = 3)
#Trace cat by cat, indicating color manually each time (3 lines)
bymodel <- by(df.gensum$V113, list(df.gensum$Gen, df.gensum$b_envir), FUN=mean)
# lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="yellowgreen", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="lightblue2", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="thistle", lwd = 3)
bymodel <- by(df.gensum$V112, list(df.gensum$Gen, df.gensum$b_envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="maroon2", lwd = 3)
# legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("darkblue","lightblue2","maroon2","thistle"),
#        legend=c("Changing 1: Plastic to Stable","Constant with Changing 1","Changing 2: Stable to Plastic","Constant with Changing 2"))

dev.off()

################################################################################

#Sum of reg on Plastic gene (=2) and Stable gene (=3)
df.gensum$V115 <- rowSums(cbind(abs(df.gensum$W_1_4), abs(df.gensum$W_1_5), abs(df.gensum$W_1_6), 
                                abs(df.gensum$W_1_7), abs(df.gensum$W_1_8), abs(df.gensum$W_1_9),  abs(df.gensum$W_1_10)))


cairo_pdf("../figures/fig_sum_s_to_tf.pdf", width=7, height=5)
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
#Plot 1
# par(mar=c(3, 3, 3, 0.5), mgp = c(1.75, 0.75, 0), las=0)
plot(df.gensum$Gen, NULL, type="n", xlim = c(min(df.gensum$Gen), max(df.gensum$Gen)),ylim = c(0,5),
     ylab = "Sum of regulations from Sensor to RF", xlab = "Generation")
abline(v = 10000, lwd = 2, lty = 3)
#Trace cat by cat, indicating color manually each time (3 lines)
bymodel <- by(df.gensum$V115, list(df.gensum$Gen, df.gensum$b_envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="maroon2", lwd = 3)
# legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("yellowgreen","darkblue","maroon2"),
#        legend=c(paste0("Stable"), paste0("Plastic to Stable"), paste0("Stable to Plastic")))
legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("darkblue","maroon2"),
       legend=c("TF with Changing 1", "TF with Changing 2"))

dev.off()

################################################################################
#Var slope but for df = for pop mean W 
mut <- 100 #number of mutation to apply
###
sims.dirs1 <- list.dirs("../simul/Evoplast/burning/Correlated", recursive = TRUE)
sims.dirs2 <- list.dirs("../simul/Evoplast/new_envir/B_Correlated", recursive = TRUE)
sims.dirs2 <- sims.dirs2[grep('N_Noncorrelated', sims.dirs2)]
sims.dirs3 <- list.dirs("../simul/Evoplast/burning/Noncorrelated", recursive = TRUE)
sims.dirs4 <- list.dirs("../simul/Evoplast/new_envir/B_Noncorrelated", recursive = TRUE)
sims.dirs4 <- sims.dirs4[grep('N_Correlated', sims.dirs4)]
sims.dirs5 <- list.dirs("../simul/Evoplast/burning/Anticorrelated", recursive = TRUE)
sims.dirs6 <- list.dirs("../simul/Evoplast/new_envir/B_Correlated", recursive = TRUE)
sims.dirs6 <- sims.dirs6[grep('N_Anticorrelated', sims.dirs6)]
#####################


df.part1 <- df.simul(c(sims.dirs1, sims.dirs3, sims.dirs5), all.gen = TRUE)
df.part1$n_envir <- str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,5]
df.part1$b_envir <- paste0("B_",str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,5])
df.part1$anc_id <- str_split(str_split(df.part1$data.dir, "/", n=8, simplify = TRUE)[,6], "-", simplify = TRUE)[,2]

df.part2 <- df.simul(c(sims.dirs2, sims.dirs4, sims.dirs6), all.gen = TRUE)
df.part2$n_envir <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,7]
df.part2$b_envir <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,5]
df.part2$anc_id <- str_split(df.part2$data.dir, "/", n=8, simplify = TRUE)[,6]
df.part2$Gen <- df.part2$Gen + max(df.part1$Gen)


df.gensum <- rbind(df.part1, df.part2)
df.gensum2 <- subset(df.gensum, Gen !=0)

df.mut_rn_sd <- mclapply(1:nrow(df.gensum2), function(i) { 
  slopes2 <- c()
  slopes3 <- c()
  Ws <- t(matrix(as.numeric(df.gensum2[i,7:106]), ncol = 10))
  Wprob <- ifelse(Ws==0,0,1)
  for (m in 1:mut) {
    W_mut <- Ws
    #Get W cell that is different from 0
    Wij <- which(W_mut == sample(W_mut, 1, prob = Wprob))  
    #Change a random value in a W matrix 
    W_mut[Wij] <- rtruncnorm(n=1, mean=W_mut[Wij], sd=sd)
    slopes2 <- c(slopes2, getSlope.ALR(W=W_mut, n.env=10, target.gene=2))
    slopes3 <- c(slopes3, getSlope.ALR(W=W_mut, n.env=10, target.gene=3))
  }
  data.gen <- data.frame(
    row  = i, 
    data.dir <- df.gensum2[i,1],
    Gen = df.gensum2[i,2],
    b_envir = df.gensum2[i,108],
    anc_id = df.gensum2[i,109],
    Sd_ind_2 = sqrt(var(as.numeric(slopes2))),
    Sd_ind_3 = sqrt(var(as.numeric(slopes3))))
  return(data.gen)
}, mc.cores=4)
df.gensum2 <- do.call(rbind, df.mut_rn_sd)


cairo_pdf("../figures/fig_rob2.pdf", width=7, height=5.5)
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
#Plot 1
# par(mar=c(3, 3, 3, 0.5), mgp = c(1.75, 0.75, 0), las=0)
plot(c(0:max(df.gensum2$Gen)), NULL, type="n", xlim = c(min(df.gensum2$Gen), max(df.gensum2$Gen)),ylim = c(0,(0.06)),
     ylab = "RN mut var sd", xlab = "Generation")
abline(v = 10000, lwd = 2, lty = 3)
#Trace cat by cat, indicating color manually each time (3 lines)
bymodel <- by(df.gensum2$Sd_ind_3, list(df.gensum2$Gen, df.gensum2$b_envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="lightblue2", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="thistle", lwd = 3)
bymodel <- by(df.gensum2$Sd_ind_2, list(df.gensum2$Gen, df.gensum2$b_envir), FUN=mean)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Correlated"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"B_Noncorrelated"], col="maroon2", lwd = 3)
legend("topleft", lty=1, box.lty=0,  bg="transparent", col=c("darkblue","lightblue2","maroon2","thistle"),
       legend=c("Changing 1: Plastic to Stable","Constant with Changing 1","Changing 2: Stable to Plastic","Constant with Changing 2"))
dev.off()


