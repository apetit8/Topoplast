#FUNCTIONS

#Packages
suppressPackageStartupMessages(library(abind))
library(parallel)
library(plyr)
library(stringr)
library(matrixStats)
library(grid)
library(scales)
library(Rcpp)
library(tidyr)
library(rlist)
library(purrr)
library(data.table)
library(igraph)

##TOOLS#######################################################################################

#Matrix extraction from output files
extract.P.mean <- function(tt, what="MPhen",  gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  ex
}

extract.matrix <- function(tt, what="CovPhen", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  if (sqrt(length(ex)) %% 1 != 0) stop("No way to make a square matrix out of ", length(ex), " elements.")
  matrix(ex, ncol=sqrt(length(ex)))
}

extract.fitness <- function(tt, what="Mfit", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  return(ex)
}

extract.opt <- function(tt, what="FitOpt1", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  return(ex)
}

extract.P.matrix <- function(tt, gen=tt[nrow(tt),"Gen"]) {
  extract.matrix(tt, "CovPhen", gen=gen)
}

extract.M.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
  extract.matrix(tt, "MGenCov", gen=gen)
}

extract.W.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
  t(extract.matrix(tt, "MeanAll", gen=gen))
}

##
#function to read a parameter file
read.param <- function(parfile) {
  ff <- readLines(parfile)
  ss <- strsplit(ff, split="\\s+")
  ans <- sapply(ss, function(x) {cx <- suppressWarnings(as.numeric(x[-1])); if(any(is.na(cx))) x[-1] else cx})
  param.names <- sapply(ss, "[", 1)
  if (any(duplicated(param.names))) warning("Duplicated param names in ", parfile)
  names(ans) <- param.names
  ans
}

#MATRIX FEATURES############################################################################
matrix.features <- function(M, n.genes=ncol(M)) {
  stopifnot(ncol(M) > 1, n.genes <= ncol(M))
  if (!isSymmetric(M)) warning("Non-symmetric matrix: check your input, results are probably meaningless.")
  M <- M[1:n.genes, 1:n.genes]
  ee <- try(eigen(M), silent=TRUE)
  ee$vectors <- t (t(ee$vectors * sign(ee$vectors[2,]))) # by convention, the second eigenvector is always positive
  
  if (class(ee) == "try-error") {
    return(list(cor=NA, angle=NA, size=NA, eccentricity=NA))
  }
  ans <- list()
  ans$cor <- cov2cor(M)[upper.tri(M)]
  ans$angle <- modulopi(acos(ee$vectors[1,]))
  ans$size <- sum(diag(M))
  ans$eccentricity <- ee$values[2:ncol(M)]/ee$values[1]
  ans
}

#TRIGONOMETRY##########################################################################

# returns the angle of an ellipse (between -pi/2 and pi/2)
modulopi <- function(angle) {
  ans <- angle %% pi
  ifelse(ans > pi/2, ans-pi, ans)
}

#DATA EXTRACTION################################################################


df.simul <- function(sims.dir, all.gen=TRUE, size=50000, M=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$")
  files     <- files[grep('out-rep\\d+\\.txt', files)]
  files <- files[sapply(files, file.size) > size]
  mytopos <- lapply(files, function(ff) {
    print(ff)
    tt <- fread(ff, data.table=FALSE)
    mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
    for (gen in mygens) {
      #M data
      gen <-gen
      phen.mean <- extract.P.mean(tt, gen=gen)
      fitness <- extract.fitness(tt, "MFit", gen=gen)
      opt <- extract.fitness(tt, "FitOpt1", gen=gen)
      if(M==TRUE){M.mat <- extract.M.matrix(tt, gen=gen)
      msize <- matrix.features(M.mat,n.genes=c(3:3))[["size"]][1] #M matrix size, used as mutational robustness
      }else (msize <- 0)
      Wmat <- extract.W.matrix(tt, gen=gen) #W data
      data.gen <- c(ff, gen, opt[1], phen.mean[2], msize, fitness, t(Wmat) )
      filedata <- rbind(filedata, data.gen)
    }
    return(filedata)
  })
  simul.df <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
  nbrc <- ncol(simul.df)-6
  data <- list(let = c("W"), id = c(seq(1,sqrt(nbrc) ,1)), greeting = c(seq(1,sqrt(nbrc),1)), sep = c("_"))
  Wcol <- data %>% cross() %>% map(lift(paste)) 
  setnames(simul.df, 1:(6+nbrc), c("data.dir","Gen", "Optimum", "P_mean_2", "M_size", "Fitness", unlist(Wcol)))
  simul.df[,2:ncol(simul.df)] <- lapply( simul.df[,2:ncol(simul.df)], as.numeric)
  return(simul.df)
}

df.last.gens <- function(sims.dir, n.gen=50, size=50000){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$")
  files     <- files[grep('out-rep\\d+\\.txt', files)]
  files <- files[sapply(files, file.size) > size]
  mytopos <- lapply(files, function(ff) {
    print(ff)
    tt <- fread(ff, data.table=FALSE)
    mygens <- (tt[nrow(tt),"Gen"]-n.gen):tt[nrow(tt),"Gen"]
    for (gen in mygens) {
      #M data
      gen <-gen
      phen.mean <- extract.P.mean(tt, gen=gen)
      data.gen <- c(ff, gen, phen.mean)
      filedata <- rbind(filedata, data.gen)
    }
    return(filedata)
  })
  simul.df <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
  nbrc <- ncol(simul.df)-2
  data <- list(let = c("Pmean"), id = c(seq(1,nbrc ,1)), sep = c("_"))
  Wcol <- data %>% cross() %>% map(lift(paste)) 
  setnames(simul.df, 1:(2+nbrc), c("data.dir","Gen", unlist(Wcol)))
  simul.df[,2:ncol(simul.df)] <- lapply( simul.df[,2:ncol(simul.df)], as.numeric)
  return(simul.df)
}


#PHENOTYPES################################################################
# netW.R
#
# The version of the Wagner model used in 
# Rünneburger & Le Rouzic 2016 BMC Evol. Biol. 
#
# 
# Copyright Arnaud Le Rouzic / CNRS 2015-2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
#
library(Rcpp)
library(inline, quietly=TRUE)

sigma.M2 <- function(x, a) {
  1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
}

sigma.M2p <- function(x, lambda, mu) {
  1. / (1. + lambda * exp(-mu*x))
}

suppressMessages(library(compiler))
sigma.M2c <- cmpfun(sigma.M2p)

cppFunction('
	List internal_loop_cpp(const NumericMatrix &W, const NumericVector &S0, double a, unsigned int steps, unsigned int measure) {
		double lambda = (1-a)/a;
		double mu     = 1/(a*(1-a));
		NumericMatrix sto (S0.size(), steps+1);
		NumericVector sumx (S0.size());
		NumericVector sumx2 (S0.size()); 
		for (unsigned int i = 0; i < S0.size(); i++)
			sto(i,0) = S0(i);
		for (unsigned int t = 1; t <= steps; t++) {
			for (unsigned int i = 0; i < S0.size(); i++) {
				double tmp = 0.;
				for (unsigned int j = 0; j < S0.size(); j++) {
					tmp += sto(j,t-1) * W(i,j);
				}
				tmp =  1. / (1. + lambda * exp(-mu*tmp));
				sto(i,t) = tmp;
				if (t > steps-measure) {
					sumx(i) += tmp;
					sumx2(i) += tmp*tmp;
				}
			}
		}
		for (unsigned int i = 0; i < S0.size(); i++) {
			sumx(i) /= static_cast<double>(measure); // sumx(i) now contains the mean
			sumx2(i) /= static_cast<double>(measure);
			sumx2(i) += -sumx(i)*sumx(i); // sumx2(i) now contains the variance
		}
		return List::create(Named("full")=sto, Named("mean")=sumx, Named("var")=sumx2);
	}')

#OG
internal_loop_R <- function(W, S0, a, steps, measure, sensors=NULL) {
  lambda <- (1-a)/a
  mu <- 1/(a*(1-a))
  sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
  sto[,1] <- S0
  for (i in 1:steps) {
    S0 <- sigma.M2c((W %*% S0), lambda=lambda, mu=mu)
    if(is.null(sensors)==FALSE){S0[1:length(sensors)] <- sensors}
    sto[,i+1] <- S0
  }
  list(full=sto, sumx=rowSums(sto[,(steps-measure+1):steps]), sumx2=rowSums(sto[,(steps-measure+1):steps]^2))
}

#Modified
phenotype_plastic_loop_R <- function(W, S0, a, steps, measure, sensors=NULL) {
  # Simulates development with a plastic signal
  # sensors = value taken by the sensor gene (=first gene)
  # a = basal gene expression
  # Steps = nbr of cycle ; measure = mean phenotype obtained over the last x measure
  lambda <- (1-a)/a
  mu <- 1/(a*(1-a))
  sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
  sto[,1] <- S0
  for (i in 1:steps) {
    S0 <- sigma.M2p((W %*% S0), lambda=lambda, mu=mu) 
    if(is.null(sensors)==FALSE){S0[1:length(sensors)] <- sensors}
    sto[,i+1] <- S0
  }
  list(full=sto, mean=rowSums(sto[,(steps-measure+1):steps])/measure, var=(rowSums(sto[,(steps-measure+1):steps]^2)/measure - (rowSums(sto[,(steps-measure+1):steps])/measure)^2) )
}

#
pheno.from.W <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=20, measure=4, full=FALSE, loopFUN=phenotype_plastic_loop_R, sensors=NULL) {
  ans <- loopFUN(W, S0, a, steps, measure, sensors=sensors)
  if (!full) ans$full <- NULL
  return(ans)
}

# pheno.from.W <- cmpfun(pheno.from.W)

#Function that gives the number of plastic genes fom W matrix from a dataframe
plasticity <- function(dfplast, genes=51, treshold=0.01, envir1=c(0.15), envir2=c(0.85)){
  df <- data.frame()
  for (i in 1:nrow(dfplast)) {
    W <-  matrix(as.numeric(dfplast[i,25:(ncol(dfplast)-2)]), ncol = genes) 
    # sensors <- c(dfplast[i,4],dfplast[i,5])  #Take Sensors genes values dataframe
    plast <- length( which( abs(pheno.from.W(W, a=0.5, sensors = envir1)$mean -  pheno.from.W(W, a=0.5, sensors = envir2)$mean) > treshold) )
    #
    df[i, 1] <- dfplast[i,(ncol(dfplast))]
    df[i, 2] <- dfplast[i,(ncol(dfplast)-1)]
    df[i, 3] <- plast
    df[i, 4] <- dfplast[i,2]
  }
  setnames(df, 1:4, c("anc","envir","nbr_plg","Gen") )
  
  return(df)
}
#Add reaction norm slope by gene category ?
#Fitness ?

getSlope.ALR <- function(W, n.env=21, target.gene=2, min=0.15, max=0.85, giveback=2, a=0.5 ) {
  #giveback 2 = reg coeff ; 1 = reg origin
  envs  <- seq(min, max, length.out=n.env)
  phens <- sapply(envs, function(env) pheno.from.W(W, sensors = env, a=a)$mean[target.gene])
  # reg   <- lm(phens ~ envs)
  reg <- .lm.fit(cbind(rep(1, length(envs)), envs), phens)$coefficients  #much faster than lm()
  return(reg[giveback])  # the first coefficient is the regression intercept, the second is the slope
}


#TOPOLOGIES#####################################################################

# redefining the combinat::permn function, for exactly the same reason
permn <- function (x, fun = NULL, ...) 
{
  if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == 
      x) 
    x <- seq(x)
  n <- length(x)
  nofun <- is.null(fun)
  out <- vector("list", gamma(n + 1))
  p <- ip <- seqn <- 1:n
  d <- rep(-1, n)
  d[1] <- 0
  m <- n + 1
  p <- c(m, p, m)
  i <- 1
  use <- -c(1, n + 2)
  while (m != 1) {
    out[[i]] <- if (nofun) 
      x[p[use]]
    else fun(x[p[use]], ...)
    i <- i + 1
    m <- n
    chk <- (p[ip + d + 1] > seqn)
    m <- max(seqn[!chk])
    if (m < n) 
      d[(m + 1):n] <- -d[(m + 1):n]
    index1 <- ip[m] + 1
    index2 <- p[index1] <- p[index1 + d[m]]
    p[index1 + d[m]] <- m
    tmp <- ip[index2]
    ip[index2] <- ip[m]
    ip[m] <- tmp
  }
  out
}
my.permn <- function(x, fun=NULL, ...) {
  if (length(x) == 1 ) return(x)
  permn(x, fun, ...)
}

#Topo analyses 
library("combinat")


#Function that gives the essential topology
#Keep "essential" connections. Test every connection to see effect on target gene RN.
#Inspired by Burda et al., 2011
essential.topo <- function(df, min=0.15, max=0.85, target=2, treshold_coeff=0.05,
                           treshold_og=0.05, genes=4, basal=0.15, cores=2){
  #
  # target = Target gene in the network which RN will be tested
  # treshold_coeff = difference accepted in the Reaction Norm linear regression slope
  # treshold_og = difference accepted in the RN linear regression intercept
  # genes = number of genes in the network
  #
  essential_Ws <- mclapply(1:nrow(df), function(i){ #df : 
    W <-  t(matrix(as.numeric(df[i,7:(genes*genes+6)]), ncol = genes))
    W2 <- W
    RN_W_coeff <- getSlope.ALR(W=W, n.env=15, target.gene=target, min=min, max=max, a=basal)
    RN_W_og <- getSlope.ALR(W=W, n.env=15, target.gene=target, min=min, max=max, giveback=1, a=basal)
    for(Wij in 1:length(W)){  #(but diagonal)
      W_test <- W
      if(W_test[Wij]!=0){
        W_test[Wij] <- 0
        RN_Wij_coeff <- getSlope.ALR(W=W_test, n.env=15, target.gene=target, min=min, max=max, a=basal)
        RN_Wij_og <- getSlope.ALR(W=W_test, n.env=15, target.gene=target, min=min, max=max, giveback=1, a=basal)
        W2[Wij] <- ifelse(treshold_coeff > (abs(RN_W_coeff-RN_Wij_coeff)) ||
                            treshold_og > (abs(RN_W_og-RN_Wij_og)), 0, 1) * sign(W2[Wij])
        }}
    return(W2)}, mc.cores=cores)
  return(essential_Ws)
}


mean_W <- function(list_topo){
  #list_topo = result from essential.topo()
  W <- Reduce("+",list_topo)/length(list_topo)
  return(W)
}






