source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
#####################
genes <- 10
min <- 0.15
max <- 0.85
target <- 2 # Target gene in the network
################################################################################
#Keep "essential" connections. Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
treshold_coeff <- 0.001 # difference accepted in the Reaction Norm linear regression slope
treshold_og <- 0.001    # difference accepted in the RN linear regression intercept
gen <- 10000 #max(df.10_a0_15$Gen)
#####################
sims.dirs1 <- list.dirs("simul/10g_a0.15", recursive = TRUE)
#####################
#DATA
##################
df.10_a0_15 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10_a0_15$envir <- str_split(df.10_a0_15$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Anticorrelated"),
                                       treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Correlated"),
                                    treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10_a0_15 <- essential.topo(df=subset(df.10_a0_15, Gen==gen & envir=="Control_sel"),
                                   treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

#####################
sims.dirs1 <- list.dirs("simul/10g", recursive = TRUE)
#####################
#DATA
##################
df.10 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10$envir <- str_split(df.10$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Anticorrelated"),
                                       treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Correlated"),
                                    treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10 <- essential.topo(df=subset(df.10, Gen==gen & envir=="Control_sel"),
                                   treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))


#####################
sims.dirs1 <- list.dirs("simul/10g_a0.2", recursive = TRUE)
#####################
#DATA
##################
df.10_a0_2 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10_a0_2$envir <- str_split(df.10_a0_2$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10_a0_2 <- essential.topo(df=subset(df.10_a0_2, Gen==gen & envir=="Anticorrelated"),
                                       treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10_a0_2 <- essential.topo(df=subset(df.10_a0_2, Gen==gen & envir=="Correlated"),
                                    treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10_a0_2 <- essential.topo(df=subset(df.10_a0_2, Gen==gen & envir=="Control_sel"),
                                   treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

#####################
sims.dirs1 <- list.dirs("simul/10g_a0.3", recursive = TRUE)
#####################
#DATA
##################
df.10_a0_3 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10_a0_3$envir <- str_split(df.10_a0_3$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10_a0_3 <- essential.topo(df=subset(df.10_a0_3, Gen==gen & envir=="Anticorrelated"),
                                      treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10_a0_3 <- essential.topo(df=subset(df.10_a0_3, Gen==gen & envir=="Correlated"),
                                   treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10_a0_3 <- essential.topo(df=subset(df.10_a0_3, Gen==gen & envir=="Control_sel"),
                                  treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

#####################
sims.dirs1 <- list.dirs("simul/10g_a0.4", recursive = TRUE)
#####################
#DATA
##################
df.10_a0_4 <- df.simul(sims.dirs1, all.gen = TRUE)
df.10_a0_4$envir <- str_split(df.10_a0_4$data.dir, "/", n=8, simplify = TRUE)[,3]

topo.anticor10_a0_4 <- essential.topo(df=subset(df.10_a0_4, Gen==gen & envir=="Anticorrelated"),
                                      treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.corr10_a0_4 <- essential.topo(df=subset(df.10_a0_4, Gen==gen & envir=="Correlated"),
                                   treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))

topo.sel10_a0_4 <- essential.topo(df=subset(df.10_a0_4, Gen==gen & envir=="Control_sel"),
                                  treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3:10))


#Plot of ratio between positive regulations and negative regulations for different constitutive expression
################################################################################
################################################################################

ff <- table(unlist(E_coli_mat))
ff[1]*100/(ff[1]+ff[3])
ff[3]*100/(ff[1]+ff[3])



ff2 <- table(unlist(c(topo.anticor10_a0_15,topo.corr10_a0_15,topo.sel10_a0_15)))
ff2[1]*100/(ff2[1]+ff2[3])
ff2[3]*100/(ff2[1]+ff2[3])


ff3 <- table(unlist(c(topo.anticor10,topo.corr10,topo.sel10)))
ff3[1]*100/(ff3[1]+ff3[3])
ff3[3]*100/(ff3[1]+ff3[3])

ff4 <- table(unlist(c(topo.anticor10_a0_2,topo.corr10_a0_2,topo.sel10_a0_2)))
ff4[1]*100/(ff4[1]+ff4[3])
ff4[3]*100/(ff4[1]+ff4[3])

ff5 <- table(unlist(c(topo.anticor10_a0_3,topo.corr10_a0_3,topo.sel10_a0_3)))
ff5[1]*100/(ff5[1]+ff5[3])
ff5[3]*100/(ff5[1]+ff5[3])

ff6 <- table(unlist(c(topo.anticor10_a0_4,topo.corr10_a0_4,topo.sel10_a0_4)))
ff6[1]*100/(ff6[1]+ff6[3])
ff6[3]*100/(ff6[1]+ff6[3])

df <- data.frame(A=c(0.15,0.2,0.3,0.4,0.5),
           Ratio=c(ff2[1]*100/(ff2[1]+ff2[3]), ff4[1]*100/(ff4[1]+ff4[3]), ff5[1]*100/(ff5[1]+ff5[3]), ff6[1]*100/(ff6[1]+ff6[3]), ff3[1]*100/(ff3[1]+ff3[3])))


dd <- rbind(ff2, ff4, ff5, ff6, ff3)


plot(df)




