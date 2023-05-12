#FUNCTIONS

#Packages
suppressPackageStartupMessages(library(abind))
library(parallel)
library(plyr)
library(stringr)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
library(ggplot2)
library(ggExtra)
library(grid)
library(ggstance)
library(gtable)
library(scales)
library(Rcpp)
library(tidyr)
library(rlist)
library(tidyverse)
library(data.table)
library(truncnorm)
library(igraph)

##TOOLS#######################################################################################

#Variances with sliding window
slidingvar <- function(x, window=3) { # x= array of values
  ans <- rep(0, 1+length(x)-window)
  for (i in 1:length(ans))
    ans[i] <- var(x[i:(i+window-1)])
  ans
}

replicate.mean <- function(tt) {
  stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
  arr <- do.call(abind, c(tt, list(along=3)))
  rowMeans(arr, dims=2)
}

replicate.var <- function(tt) {
  # This rowVars function comes from https://stat.ethz.ch/pipermail/r-help/2006-April/103001.html
  # Author: David Brahm
  .rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
    if (SumSquares) return(rowSums(x^2, na.rm, dims))
    N <- rowSums(!is.na(x), FALSE, dims)
    Nm1 <- if (unbiased) N-1 else N
    if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
    (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
  }
  stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
  arr <- do.call(abind, c(tt, list(along=3)))
  .rowVars(arr, dims=2)
}

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

extract.matrix.all.pop <- function(tt, what="CovPhen", index=1) {
  # index can be a vector -> returns a list of matrices
  mycols <- grep(colnames(tt), pattern=what)
  if (length(mycols) == 0) stop("No match for ", what, " at generation ", index,".")
  if (sqrt(length(mycols)) %% 1 != 0) stop("No way to make a square matrix out of ", length(mycols), " elements.")
  
  if (length(index)== 1) {
    return(matrix(unlist(tt[index, mycols]), ncol=sqrt(length(mycols))))
  } else {
    return(lapply(index, function(i) matrix(unlist(tt[i, mycols]), ncol=sqrt(length(mycols)))))
  }
}

extract.W.matrix.all.pop <- function(tt, index=1) {
  t(extract.matrix.all.pop(tt, "MeanAll", index=index))
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

##
#function to write in a parameter file
write.param <- function(parlist, parfile) {
  unlink(parfile)
  for (nn in names(parlist)) 
    cat(nn, "\t", parlist[[nn]], "\n", file=parfile, append=TRUE)
}

#function to write multiple parameter files for multiple changing parameters
generate.param.list <- function(template, param.list, reps=1, launchfile=NA, sep=c("-", "-", "-R"), vec.indx=rep(1, length(param.list)), variable = 1, name="simu",
                                simevolv="$SIMEVOLV", param.dir=dirname(template), simu.dir=dirname(template), param.ext=".par", simu.ext=".txt", nextpar=FALSE) {
  
  launch <- NULL
  templ <- read.param(template)
  param.grid <- do.call(expand.grid, param.list)
  if (nextpar==1) param.grid[,1] <- as.character(param.grid[,1] ) #the number = number of the argument : nextpar
  if (nextpar==2) param.grid[,2] <- as.character(param.grid[,2] ) #the number = number of the argument : nextpar
  names.grid <- do.call(expand.grid,  list(lapply(param.list, function(pp) if (is.null(names(pp))) as.character(pp) else  as.character(pp))))
  colnames(names.grid) <- paste0(colnames(names.grid), ifelse(vec.indx == 1, "", vec.indx))
  rownames(param.grid) <- do.call(paste, c(lapply(colnames(names.grid), function(nn) paste(nn, names.grid[,nn], sep=sep[1])), list(sep=sep[2])))
  for (i in seq_len(nrow(param.grid))) {
    mypar <- templ
    for (j in seq_along(names(param.list))) {
      if (!names(param.list)[j] %in% names(mypar))
        stop("Error: parameter ", names(param.list)[j], " is not in file ", template, ".")
      # if (length(mypar[[names(param.list)[j]]]) < vec.indx[j])
      #   stop("Error: parameter ", names(param.list)[j], " does not have ", vec.indx[j], "elements.")
      mypar[[names(param.list)[j]]][vec.indx[j]] <- param.grid[i,j]
    }
    # parfile <- paste0(param.dir, "/", rownames(param.grid)[i], param.ext)
    parfile <- paste0(param.dir, "/", name, variable, param.ext)
    # browser()
    write.param(parfile=parfile, parlist=mypar)
    for (rr in seq_len(reps)) {
      launch <- c(launch, paste(simevolv, "-p", paste0(rownames(param.grid)[i], param.ext), "-o", paste0(rownames(param.grid)[i], sep[3], formatC(rr, width = 1+floor(log10(reps)), format = "d", flag = "0"), simu.ext), sep=" "))
    }
  }
  # browser()
  if (!is.na(launchfile)) {
    dir <- as.character(dirname(template))
    cat(launch, file=print(sprintf("%s/%s", dir , launchfile)), sep="\n")
  }
  launch
}

print.launcher <- function(param.list, dir="../simul/", reps=30, simevolv="$SIMEVOLV", simu.ext=".txt", launchfile="launcher.sh"){
  launch <- NULL
  for (i in param.list) {
        for (rr in seq_len(reps)) {
          launch <- c(launch, paste(simevolv, "-p", i, "-o", paste0(i,rr,simu.ext), sep=" "))
        }
        if (!is.na(launchfile)) {
          cat(launch, file=print(sprintf("%s/%s", dir , launchfile)), sep="\n")
        }
  }
}

# print.launcher.P <- function(param.list, dir="../simul/", reps=30, simevolv="$SIMEVOLV", simu.ext=".txt", launchfile="launcher.sh", prev_pop="pop.pop"){
#   launch <- NULL
#   jj <- 1
#   for (prev in prev_pop) {
#     for (i in param.list) {
#     
#       for (rr in seq_len(reps)) {
#         launch <- c(launch, paste(simevolv, "-p", i, "-o", paste0(i,"_anc",jj,"/",rr,simu.ext), sep=" "," -P ", prev))
#         # dir.create(paste0(i,"_anc",jj,"/",rr,simu.ext))
#         ifelse(!dir.exists(file.path(paste0(i,"_anc",jj))), dir.create(paste0(i,"_anc",jj,"/")), FALSE)
#         }
#       if (!is.na(launchfile)) {
#         cat(launch, file=print(sprintf("%s", launchfile)), sep="\n")
#       }
#     }
#     jj <- jj + 1
#   }
# }


print.launcher.P <- function(param.list, dir="../simul/", reps=30, simevolv="$SIMEVOLV",
                    simu.ext=".txt", launchfile="launcher.sh", prev_pop="pop.pop", directory="_amp/"){
  launch <- NULL
  jj <- 1
  for (prev in prev_pop) {
    for (i in param.list) {
      for (rr in seq_len(reps)) {
        parfiledir <- str_split(i, "param", n=2, simplify = TRUE)[1]
        parfile_id <- str_split(i, "param", n=2, simplify = TRUE)[2]
        launch <- c(launch, paste(simevolv, "-p", i, "-o", paste0(parfiledir, directory,parfile_id,"_anc",jj,"/",rr,simu.ext), sep=" "," -P ", prev))
        # dir.create(paste0(i,"_anc",jj,"/",rr,simu.ext))
        ifelse(!dir.exists(file.path(paste0(parfiledir, directory,parfile_id,"_anc",jj))), dir.create(paste0(parfiledir, directory,parfile_id,"_anc",jj)), FALSE)
      }
      if (!is.na(launchfile)) {
        cat(launch, file=print(sprintf("%s", launchfile)), sep="\n")
      }
    }
    jj <- jj + 1
  }
}

print.launcher.J <- function(param.list, reps=30, simevolv="$SIMEVOLV", simu.ext=".txt", launchfile="launcher.sh", prev_pop="pop.pop"){
  launch <- NULL
  for (prev in prev_pop) {
    for (i in param.list) {
      for (rr in seq_len(reps)) {
        parfiledir <- str_split(i, "optimum", n=2, simplify = TRUE)[1]
        parfile_id <- str_split(i, "optimum", n=2, simplify = TRUE)[2]

        ancpop <- str_split(prev, "../simul/Population_", n=2, simplify = TRUE)
        ancpop <- str_split(ancpop[1,2], "/rep", n=2, simplify = TRUE)
        prevpop_id <- str_split(ancpop[1,2], "/out", n=2, simplify = TRUE)[1]
        launch <- c(launch, paste(simevolv, "-p", i, "-o", paste0(parfiledir, ancpop[,1],"/",prevpop_id,"/",parfile_id,"_",rr,simu.ext), sep=" "," -P ", prev))
        # browser()
        }
    }
  }
  if (!is.na(launchfile)) {
    cat(launch, file=print(sprintf("%s", launchfile)), sep="\n", append = FALSE)
  }
}

#Function to print parameter files of the wanted S matrix
param.from.sel.features <- function(param.template, param.out="param.par", cor=NA, angle=NA, size=NA, eccentricity=NA, reps=40) {
  pp <- read.param(param.template)
  S.mat <- matrix2.from.features(cor=cor, angle=angle, size=size, eccentricity=eccentricity)
  n.genes <- if ("GENET_NBPHEN" %in% names(pp)) pp$GENET_NBPHEN else pp$GENET_NBLOC
  pS <- param.S.matrix(S.mat, n.genes=n.genes)
  pp[names(pS)] <- pS
  browser()
  write.param(pp, parfile=param.out)
}

extract.nbphen <- function(parfile) {
  rp <- read.param(parfile)
  if (rp$TYPE_ARCHI %in% c("additive", "multilinear")) {
    np <- if("GENET_NBPHEN" %in% names(rp)) rp$GENET_NBPHEN else 1
  } else if (rp$TYPE_ARCHI %in% c("wagner", "siegal", "m2")) {
    np <- rp$GENET_NBLOC
  }
  np
}

extract.theta <- function(parfile) {
  rp <- read.param(parfile)
  np <- extract.nbphen(parfile)
  th <- rp$FITNESS_OPTIMUM
  if (length(th) == 1) return(rep(th, np))
  if (length(th) == np) return(th)
  stop("Param file: ", parfile, ", The number of optima (", length(th), ") does not match the number of traits (", np, ").")
}


extract.S.matrix <- function(parfile) {
  .cor2cov <- function(cc, sd) t(cc*sd)*sd
  rp <- read.param(parfile)
  np <- extract.nbphen(parfile)
  fs <- rp$FITNESS_STRENGTH
  if (length(fs)==1)
    fs <- rep(fs, np)
  if (length(fs) != np)
    stop("Param file: ", parfile, ", The number of fitness strengths (", length(fs), ") does not match the number of traits (", np, ").")
  if (rp$FITNESS_TYPE == "gaussian") {
    ans <- diag(1/(2*fs))
  } else if (rp$FITNESS_TYPE == "multivar_gaussian") {
    vv  <- 1/(2*fs)
    rr <- rp$FITNESS_CORRELATION
    if (length(rr) != np*(np-1)/2)
      stop("Param file: ", parfile, ", The number of fitness correlations (", length(rr), ") does not match the number of traits (", np, ").")
    ansr <- diag(np)
    ansr[upper.tri(ansr)] <- rr
    ansr <- as.matrix(Matrix::forceSymmetric(ansr))
    ans <- .cor2cov(ansr, sqrt(vv))
  } else {
    stop("Asking the S matrix for non-gaussian selection is meaningless.")
  }
  ans
}

extract.S.matrix.fig3 <- function(parfile) {
  .cor2cov <- function(cc, sd) t(cc*sd)*sd
  rp <- read.param(parfile)
  fs <- rp$FITNESS_STRENGTH
  np <- 4
  vv  <- 1/(2*fs)
  rr <- rp$FITNESS_CORRELATION
  ansr <- diag(np)
  ansr[upper.tri(ansr)] <- rr
  ansr <- as.matrix(Matrix::forceSymmetric(ansr))
  ans <- .cor2cov(ansr, sqrt(vv))
  ans
}

##
#Take properties of S and convert it in S parameters for Simevolv 
param.S.matrix <- function(S.mat, n.genes=ncol(S.mat)) {
  stopifnot(is.matrix(S.mat), ncol(S.mat) == nrow(S.mat))
  stopifnot(n.genes > 0)
  if (n.genes < ncol(S.mat))
    S.mat <- S.mat[1:n.genes, 1:n.genes]
  if (n.genes > ncol(S.mat)) {
    rec <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=nrow(S.mat))
    dg <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=n.genes-nrow(S.mat))
    diag(dg) <- Inf
    S.mat <- rbind(cbind(S.mat, rec), cbind(t(rec), dg))
  }
  list(
    FITNESS_STRENGTH = 1/(2*diag(S.mat)),
    FITNESS_CORRELATION = cov2cor(S.mat)[upper.tri(S.mat)])
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

#Functions to find a matrix that gives an ellipse of given parameters (Numeric methods)
matrix2.targetdist <- function(pp, target) {
  if (pp["v1"] <= 0 || pp["v2"] <= 0) return(Inf)
  mf <- matrix.features(matrix(c(pp["v1"], pp["c"], pp["c"], pp["v2"]), ncol=2))
  ans <- sum(sapply(names(target), function(nn) mf[[nn]][1]-target[nn])^2)
  ans
}

matrix2.optim <- function(target) {
  stopifnot(all(names(target) %in% c("cor","angle","size","eccentricity")))
  stopifnot(length(target) == 3)
  
  st <- c(target["size"]/2, target["size"]/2, sign(target["angle"])*target["size"]/10)
  names(st) <- c("v1","v2","c")
  ans <- optimx::optimx(par=st, matrix2.targetdist, method="Nelder-Mead", target=target, itnmax=2000)
  if (ans$convcode != 0) 
    warning("Convergence failed for ", paste0(names(target), "=", target, collapse="  "))
  matrix(c(ans$v1, ans$c, ans$c, ans$v2), ncol=2)
}

# Common wrapper function
matrix2.from.features <- function(cor=NA, angle=NA, size=NA, eccentricity=NA, method=c("numeric", "analytic")[1]) {
  if (method == "numeric") {
    tt <- c(cor=cor, angle=modulopi(angle), size=size, eccentricity=eccentricity)
    tt <- tt[!is.na(tt)]
    return(matrix2.optim(tt))
  } else if (method == "analytic") {
    if (is.na(cor) || is.na(size) || is.na(eccentricity) ||!is.na(angle))
      stop("Analytic method not available for this combination.")
    return(matrix2.noangle(r=cor, S=size, e=eccentricity))
  } else {
    stop("Method ", method, " not available.")
  }
}


#TRIGONOMETRY##########################################################################

# returns the angle of an ellipse (between -pi/2 and pi/2)
modulopi <- function(angle) {
  ans <- angle %% pi
  ifelse(ans > pi/2, ans-pi, ans)
}

#Diff between angle, allows to pick the modulo
modulo.all <- function(angle, modulo=pi) {
  angle <- angle %% modulo
  ifelse(angle > modulo/2, angle - modulo, angle)
}

#Calculate an angular mean between pi and -pi
mean.angle.2pi <- function(data) {
  atan2(mean(sin(data)), mean(cos(data)))
}

#Calculate an angular mean between pi/2 and -pi/2
mean.angle.pi <- function(data) {
  modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
}

#Put mean M direction close to the S defined between -pi/2 and pi/2
mean.angle.pi.byS <- function(data, ang_S) {
  df <- modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
  if ((df - mean(ang_S, 7)/2) < -pi/2){ df = df + pi}
  if ((df - mean(ang_S, 7)/2) > pi/2){ df = df - pi}
  return(df)
} 

#DATA EXTRACTION################################################################

#For 10 genes

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


#For 10 genes
df.sampling <- function(sims.dir, gen=1, size=50000){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$") #Check if this works
  files <- files[sapply(files, file.size) > size]
  for (i in files) { 
    print(i)
    mytopos <- lapply(i, function(ff) {
      tt <- fread(ff, data.table=FALSE) 
      #Data
      b.gen <- str_split(ff, "/", n=8, simplify = TRUE)[,4]
      Wmat <- extract.W.matrix(tt, gen=gen) #W data
      W <- as.matrix(Wmat, ncol=10)
      slope <- getSlope.ALR(W=W, n.env=21, target.gene=2)
      prevpop <- str_split(ff, "/sampling/", n=2, simplify = TRUE)[,1]
      data.gen <- c(ff, prevpop, b.gen, slope )
      filedata <- rbind(filedata, data.gen)
      return(as.list(filedata))
    })
    newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
    setnames(newt, 1:4, c("data.dir","Before_Pop","Before_Gen", "Slope"))
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df$Before_Gen <- as.numeric( str_split(simul.df$Before_Gen, "g", simplify = TRUE)[,2] )
  simul.df[,3:4] <- lapply( simul.df[,3:4], as.numeric)
  return(simul.df)
}


#To use every gen sampling info
#Problem !! : There is no first column that can be used as an index (as "Gen" for regular output files).
#Solution : replicate "extract W" but for row index instead of Gen.
df.sampling2 <- function(sims.dir, size=50000){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\_G") #Check if this works
  files <- files[sapply(files, file.size) > size]
  file.results <- mclapply(files, function(ff) { 
    print(ff)
    tt <- fread(ff, data.table=FALSE)
    #Data
    gen <- str_split(str_split(str_split(ff, "/", n=8, simplify = TRUE)[,6], "_G", simplify = TRUE)[,2], ".txt", simplify = TRUE)
    
    Ws <- extract.W.matrix.all.pop(tt, index=1:nrow(tt)) #Is it in the right order ?
    slopes2 <- sapply(Ws, function(W) getSlope.ALR(W=t(W), n.env=10, target.gene=2))#WHY t() ??Why two times, one in the extract matrix function and one now ? It's not the case anywhere else
    slopes3 <- sapply(Ws, function(W) getSlope.ALR(W=t(W), n.env=10, target.gene=3))#WHY t() ??Why two times, one in the extract matrix function and one now ? It's not the case anywhere else
    
    var.slope2 <- var(slopes2)
    mean.slope2 <- mean(slopes2)
    var.slope3 <- var(slopes3)
    mean.slope3 <- mean(slopes3)
    prevpop <- str_split(ff, "/", simplify = TRUE)[,4]
    data.gen <- data.frame(
      data.dir   = ff, 
      Before_Pop = prevpop, 
      Before_Gen = gen[,1], 
      Var_Slope2  = sqrt(var.slope2), 
      Mean_Slope2 = mean.slope2,
      Var_Slope3  = sqrt(var.slope3), 
      Mean_Slope3 = mean.slope3)
  }, mc.cores=4)
  
  simul.df <- do.call(rbind, file.results)
  return(simul.df)
}
      

#Function that gives the mutational robustness of the reaction norm of a pop.
df.sampling3 <- function(sims.dir, size=50000, rep=6){
  #rep has to be => 2 to calculate a variance
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\_G") #Check if this works
  files <- files[grep('G00000', files, invert = TRUE)]
  files <- files[sapply(files, file.size) > size]
  file.results <- lapply(files, function(ff) { #mclapply(files, function(ff) { 
    print(ff)
    tt <- fread(ff, data.table=FALSE)
    #Data
    gen <- str_split(str_split(str_split(ff, "/", n=8, simplify = TRUE)[,6], "_G", simplify = TRUE)[,2], ".txt", simplify = TRUE)
    #
    slopes2 <- c()
    slopes3 <- c()
    varslp2 <- c()
    varslp3 <- c()
    W <- extract.W.matrix.all.pop(tt, index=1:nrow(tt)) #Is it in the right order ?
    for (j in 1:length(W)) {
      Ws <- as.matrix(W[[j]], nc)
      Wprob <- ifelse(Ws==0,0,1)
      #loop for each mut
      for (i in 1:rep) {
        W_mut <- Ws
        #Get W cell that is different from 0
        Wij <- which(W_mut == sample(W_mut, size=1, replace=TRUE, prob = Wprob))  
        #Change a random value in a W matrix 
        W_mut[Wij] <- rtruncnorm(n=1, mean=W_mut[Wij], sd=sd)
        slopes2 <- c(slopes2, getSlope.ALR(W=t(W_mut), n.env=5, target.gene=2))
        slopes3 <- c(slopes3, getSlope.ALR(W=t(W_mut), n.env=5, target.gene=3))
      }
      varslp2 <- c(var(as.numeric(slopes2)))
      varslp3 <- c(var(as.numeric(slopes3)))
    }
    sd_rep2 <- mean(varslp2)
    sd_rep3 <- mean(varslp3)
    prevpop <- str_split(ff, "/", simplify = TRUE)[,4]
    #
    data.gen <- data.frame(
      data.dir   = ff, 
      Before_Pop = prevpop, 
      Before_Gen = gen[,1], 
      Sd_ind_2 = sqrt(sd_rep2),
      Sd_ind_3 = sqrt(sd_rep3))
    
    return(data.gen)
  })#, mc.cores=4)
  
  simul.df <- do.call(rbind, file.results)
  return(simul.df)
}

#DF to compare mean centrality (every genes of every pop) for the 5 gene categories
centrality.df <- function(df, Ss=1, Pg=1, Sg=1, Tf=7, start=7){
  # Ss = number of sensor genes ; Pg = plastic genes ; Sg = stables genes ; Tf = transcription factors
  # start = column in df from which the network begun
  df.centrality <- data.frame() 
  j<-1
  for (i in 1:nrow(df)) {
    W <- t(matrix(as.numeric(df[i,start:((Ss+Tf+Pg+Sg)*(Ss+Tf+Pg+Sg)+start-1)]), ncol = (Ss+Tf+Pg+Sg)))
    #W matrix as a graph : 
    G <- as.directed(graph.adjacency(t(W), weighted = T))
    G <- delete.edges(G, E(G)[ abs(weight) < 0.1 ])
    ec <- eigen_centrality(G, directed=T, weights=NA, options=list(maxiter=1000000))$vector
    deg <- degree(G, mode = "all")
    degin <- degree(G, mode = "in")
    degout <- degree(G, mode = "out")
    #By gene category 
    df.centrality[j:(j+Ss-1), 1] <-"Sensor"
    df.centrality[(j+Ss):(j+Pg+Ss-1), 1] <- "Plastic"
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 1] <- "Stable"
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 1] <- "Transcriptor"
    #
    df.centrality[j:(j+Ss-1), 2] <- ec[1:Ss]
    df.centrality[(j+Ss):(j+Pg+Ss-1), 2] <- ec[(Ss+1):(Ss+Pg)]
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 2] <- ec[(Ss+Pg+1):(Ss+Pg+Sg)]
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 2] <- ec[(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf)]
    #
    df.centrality[j:(j+Ss-1), 3] <- deg[1:Ss]
    df.centrality[(j+Ss):(j+Pg+Ss-1), 3] <- deg[(Ss+1):(Ss+Pg)]
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 3] <- deg[(Ss+Pg+1):(Ss+Pg+Sg)]
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 3] <- deg[(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf)]
    #
    df.centrality[j:(j+Ss-1), 4] <- degin[1:Ss]
    df.centrality[(j+Ss):(j+Pg+Ss-1), 4] <- degin[(Ss+1):(Ss+Pg)]
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 4] <- degin[(Ss+Pg+1):(Ss+Pg+Sg)]
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 4] <- degin[(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf)]
    #
    df.centrality[j:(j+Ss-1), 5] <- degout[1:Ss]
    df.centrality[(j+Ss):(j+Pg+Ss-1), 5] <- degout[(Ss+1):(Ss+Pg)]
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 5] <- degout[(Ss+Pg+1):(Ss+Pg+Sg)]
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 5] <- degout[(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf)]
    #
    df.centrality[j:(j+Ss-1), 6] <- sum(abs(W[1:Ss,]))/(1:Ss)
    df.centrality[(j+Ss):(j+Pg+Ss-1), 6] <- sum(abs(W[(Ss+1):(Ss+Pg),]))/((Ss+1):(Ss+Pg)) 
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 6] <- sum(abs(W[(Ss+Pg+1):(Ss+Pg+Sg),]))/((Ss+Pg+1):(Ss+Pg+Sg))
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 6] <- sum(abs(W[(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf),]))/((Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf))
    #    
    df.centrality[j:(j+Ss-1), 7] <- sum(abs(W[,1:Ss]))/(1:Ss)
    df.centrality[(j+Ss):(j+Pg+Ss-1), 7] <- sum(abs(W[,(Ss+1):(Ss+Pg)]))/((Ss+1):(Ss+Pg)) 
    df.centrality[(j+Pg+Ss):(j+Pg+Ss+Sg-1), 7] <- sum(abs(W[,(Ss+Pg+1):(Ss+Pg+Sg)]))/((Ss+Pg+1):(Ss+Pg+Sg))
    df.centrality[(j+Pg+Ss+Sg):(j+Tf+Pg+Ss+Sg-1), 7] <- sum(abs(W[,(Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf)]))/((Ss+Pg+Sg+1):(Ss+Pg+Sg+Tf))
    
    #I can add as many indexes as i want
    df.centrality[j:(j+Tf+Pg+Ss+Sg-1), 8] <- df[i,2] #Gen
    df.centrality[j:(j+Tf+Pg+Ss+Sg-1), 9] <- df[i,(107)] #Environment
    
    j<- j+(Ss+Tf+Pg+Sg)
  }
  df.centrality$V10 <- paste0(df.centrality$V1,"_", df.centrality$V9)
  setnames(df.centrality, 1:10, c("Gene_cat","Centrality","Connect","ConnectIn","ConnectOut","Sum_in","Sum_out","Gen","Envir","Cat_envir" ) )
  
  return(df.centrality)
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
    S0 <- sigma.M2c((W %*% S0), lambda=lambda, mu=mu) 
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

pheno.from.W <- cmpfun(pheno.from.W)

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
  reg   <- lm(phens ~ envs)
  return(coef(reg)[giveback])  # the first coefficient is the regression intercept, the second is the slope
}

slope.stablefrom <- function(df, thresh=0.001, window=5, ...) {
  df.stab <- data.frame()
  for (file in unique(df$data.dir)) {
    slopes <- c()
    df2 <- subset(df, data.dir==file) #Get all gen for one output file
    #Get slope for all gen
    for (i in 1:nrow(df2)) {
      W <- t(matrix(as.numeric(df2[i,7:106]), ncol = 10))
      slope <- getSlope.ALR(W, n.env = 5)
      slopes <- c(slopes, slope)
    }
    varslopes <- slidingvar(slopes, window=window)
    #Find at which gen the slope is stabilized
    Genstab <- which(varslopes < thresh)
    Genstab <- max(Genstab)
    stabslope <- c(unique(df2$Before_Gen), df2[Genstab,2]) #Column 2 is the generation.
    
    df.stab[nrow(df.stab)+1,1:2] <- stabslope
  }
  return(df.stab)  #return : df with Before_Gen and Gen to reach stable slope
}



#TOPOLOGIES#####################################################################

# redefining the combinat::permn function, for exactly the same reason
my.permn <- function(x, fun=NULL, ...) {
  if (length(x) == 1 ) return(x)
  combinat::permn(x, fun, ...)
}

netw.group <- function(df, Ss=1, Pg=1, Sg=1, Tf=7, start=7, target=2){
  # Ss = number of sensor genes ; Pg = plastic genes ; Sg = stables genes ; Tf = transcription factors
  # start = column in df from which the network begun
  #Target : 2="Changing gene" ; 3="Constant" gene ; 1="Signal" gene
  file.results <- lapply(1:nrow(df), function(i) { #mclapply(files, function(ff) {
    W <- t(matrix(as.numeric(df[i,start:((Ss+Tf+Pg+Sg)*(Ss+Tf+Pg+Sg)+start-1)]), ncol = (Ss+Tf+Pg+Sg)))
    #W matrix as a graph : 
    G <- as.directed(graph.adjacency(t(W), weighted = T))
    G <- delete.edges(G, E(G)[ abs(weight) < mean(abs(weight)) ]) #delete.edges(G, E(G)[ abs(weight) < 0.1 ]) ; mean(abs(weight))
    
    V(G)$color <- c("darkred","orange", "green","yellow",  "yellow", "yellow", "yellow", "yellow", "yellow", "yellow")
    
    GG <- induced_subgraph(G, c(1,target, neighbors(G,target, mode="all")))
    GG <- t(as.matrix(as_adj(GG, attr="weight")))
    
    return(GG)
  })
  return(file.results)
}


#Topo analyses 
library("combinat")

unique.topo <- function(topo, groups=as.list(1:ncol(topo))) {
  # topo is a n x n square matrix (forced to be a signed matrix)
  # groups is a list of groups of genes, numbered from 1 to the n
  # Just in case the user is sloppy, smart guesses:
  topo <- sign(topo)
  if (!all(1:ncol(topo) %in% unlist(groups)))
    groups <- c(groups, as.list((1:ncol(topo))[! 1:ncol(topo) %in% unlist(groups)]))
  stopifnot(
    is.matrix(topo), ncol(topo) > 0, ncol(topo) == nrow(topo),
    all(topo %in% c(-1,0,1)),
    length(unlist(groups)) == ncol(topo), all(unlist(groups) %in% 1:ncol(topo))
  )
  equiv.topos(topo, groups)[[1]]
}


equiv.topos <- function(topo, groups=as.list(1:ncol(topo)), sorted=TRUE, unique=TRUE) {
  # topo is a n x n square matrix (forced to be a signed matrix)
  # groups is a list of groups of genes, numbered from 1 to the n
  # Just in case the user is sloppy, smart guesses:
  topo <- sign(topo)
  if (!all(1:ncol(topo) %in% unlist(groups)))
    groups <- c(groups, as.list((1:ncol(topo))[! 1:ncol(topo) %in% unlist(groups)]))
  #
  stopifnot(
    is.matrix(topo), ncol(topo) > 0, ncol(topo) == nrow(topo),
    all(topo %in% c(-1,0,1)),
    length(unlist(groups)) == ncol(topo), all(unlist(groups) %in% 1:ncol(topo))
  )
  #
  group.perms <- lapply(groups, function(gr) my.permn(gr))
  glob.perms   <- do.call(expand.grid, lapply(sapply(group.perms, length), seq_len))
  all.perms  <- apply(glob.perms, 1, function(x) unlist(lapply(seq_along(x), function(i) group.perms[[i]][x[i]])))
  topo.list <- lapply(as.data.frame(all.perms), function(p) topo[p,p])
  tokeep <- 1:length(topo.list)
  if (unique || sorted) {
    topo.df <- as.data.frame(do.call(rbind, lapply(topo.list, c)))
    if (sorted)
      tokeep <- do.call(order, c(as.list(topo.df[tokeep,]), list(decreasing=TRUE)))
    if (unique)
      tokeep <- tokeep[!duplicated(topo.df[tokeep,])]			
  }
  topo.list[tokeep] 
}


#Function that gives the essential topology
#Keep "essential" connections. Test every connection to see effect on target gene RN.
#Inspired by Burda et al., 2011
essential.topo <- function(df, min=0.15, max=0.85, target=2, treshold_coeff=0.05,
                           treshold_og=0.05, genes=4, groups=list(1,2,3:4), basal=0.5){
  #
  # target = Target gene in the network which RN will be tested
  # treshold_coeff = difference accepted in the Reaction Norm linear regression slope
  # treshold_og = difference accepted in the RN linear regression intercept
  # genes = number of genes in the network
  #
  essential_Ws <- lapply(1:nrow(df), function(i){ #df : 
    W <-  t(matrix(as.numeric(df[i,7:(genes*genes+6)]), ncol = genes))
    W2 <- W
    #basal <- if(grepl("Down", df[i,(ncol(df)-1)])) 0.8 else if(grepl("Up", df[i,(ncol(df)-1)])) 0.2 else 0.5
    #The line above : retrieve basal value from file direction. Only works for MY CURRENT design.
    RN_W_coeff <- getSlope.ALR(W=W, n.env=30, target.gene=target, min=min, max=max, a=basal)
    RN_W_og <- getSlope.ALR(W=W, n.env=30, target.gene=target, min=min, max=max, giveback=1, a=basal)
    for(Wij in 1:length(W)){  #(but diagonal)
      W_test <- W
      W_test[Wij] <- 0
      RN_Wij_coeff <- getSlope.ALR(W=W_test, n.env=30, target.gene=target, min=min, max=max, a=basal)
      RN_Wij_og <- getSlope.ALR(W=W_test, n.env=30, target.gene=target, min=min, max=max, giveback=1, a=basal)
      # W2[Wij] <- ifelse(RN_Wij_coeff >= (RN_W_coeff-treshold_coeff) &
      #                     RN_Wij_coeff <= (RN_W_coeff+treshold_coeff) &
      #                     RN_Wij_og >= (RN_W_og-treshold_og)&
      #                     RN_Wij_og <= (RN_W_og+treshold_og), 0, 1) * sign(W2[Wij])
      W2[Wij] <- ifelse(treshold_coeff > (abs(RN_W_coeff-RN_Wij_coeff)) ||
                          treshold_og > (abs(RN_W_og-RN_Wij_og)), 0, 1) * sign(W2[Wij])
    }
    return(W2)})
  
   simul.topo <- lapply(1:length(essential_Ws), function(i){
     untopo <- unique.topo(essential_Ws[[i]], groups=groups)
     # untopo <- equiv.topos(essential_Ws[[i]], groups)[[1]]
     return(untopo)
  })
  return(simul.topo)
}


mean_W <- function(list_topo){
  #list_topo = result from essential.topo()
  W <- Reduce("+",list_topo)/length(list_topo)
  return(W)
}


core_topo <- function(list_topo, treshold=0.90, asgraph=TRUE, sorting_4g=FALSE){
  #list_topo = result from essential.topo()
  if(sorting_4g==TRUE) list_topo <- sorting_topo_4g(list_topo)
  mean_topo <- mean_W(list_topo)
  for (cell in 1:length(mean_topo)) {
    mean_topo[cell] <- ifelse(abs(mean_topo[cell])>=treshold,1,0) * sign(mean_topo[cell])
  }
  if (asgraph==TRUE) return(as.directed(graph.adjacency(t(mean_topo), weighted = T))) else return(mean_topo)
}


core_topo.alt <- function(list_topo, treshold=0.9, asgraph=TRUE, sorting_4g=FALSE){
  #list_topo = result from essential.topo()
  if(sorting_4g==TRUE) list_topo <- sorting_topo_4g(list_topo)
  mean_topo <- mean_W(list_topo)
  freq_topo <- mean_W(lapply(list_topo, abs)) #Freq matrix
  print(mean_topo)
  print(freq_topo)
  # browser()
  #Keep only cells with freq > treshold
  for (cell in 1:length(freq_topo)) {
    freq_topo[cell] <- ifelse(freq_topo[cell]>=treshold,1,0)
  } 
  #Only keep the regulations that are always there ; if value != from 1 or -1 it means that the reg can be positive or negative
  for (cell in 1:length(mean_topo)) {
    mean_topo[cell] <- ifelse(freq_topo[cell]==1,1,0) * mean_topo[cell]
  }
  if (asgraph==TRUE) return(as.directed(graph.adjacency(t(mean_topo), weighted = T))) else return(mean_topo)
}


sorting_topo_4g <- function(list_topo, target=2, genes=c(3,4)){
  sort_topo <- lapply(list_topo, function(W){
    W2 <- W
    if(W[target,genes[1]]==1){
      #transfo matrix, put column 1 into position 2
      W1 <- W
      W1[,genes[1]] <- W[,genes[2]]
      W1[,genes[2]] <- W[,genes[1]]
      W2 <- W1
      W2[genes[1],] <- W1[genes[2],]
      W2[genes[2],] <- W1[genes[1],]
    }
    if(W[target,genes[2]]==-1){
      #transfo matrix
      W1 <- W
      W1[,genes[1]] <- W[,genes[2]]
      W1[,genes[2]] <- W[,genes[1]]
      W2 <- W1
      W2[genes[1],] <- W1[genes[2],]
      W2[genes[2],] <- W1[genes[1],]
    }
    return(W2)  })
}





