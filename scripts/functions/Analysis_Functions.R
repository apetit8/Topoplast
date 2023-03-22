#!/usr/bin/env Rscript

###########################################################
# 
# R functions for the analysis of simulations results
# produced by "simevolv"
# 
# Copyright Arnaud Le Rouzic / CNRS 2019
# <lerouzic@egce.cnrs-gif.fr>
# Copyright Andreas Odorico / CNRS 2019
# <andreas.odorico@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


# Functions are organized in the chronological order in which they are used in most Analysis Scripts.

##############################################
### Dependencies
##############################################

# The only package used in the analysis script is R "parallel"
library(parallel)

"mcsapply" <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
# mcsapply is a parallel-equivalent to sapply
    if (!require(parallel)) mclapply <- lapply
    FUN <- match.fun(FUN)
    answer <- mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!identical(simplify, FALSE) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

##############################################
### Miscellaneous functions
##############################################

InterruptFun <- function(msg) {
# Used in interactive mode to force the user to notice an error message... saving useless calculi that may last a couple hours...
	A <- readline(msg)
	while( ! (tolower(A) == "yes" | tolower(A) =="y") ) { A <- readline("\n\tPlease, type \"yes\" or \"y\" to Continue.\n") }
}


##############################################
### Function used to Obtain Data Files
##############################################

Getfiles <- function(path_to_sim, names, filePtrn="out_rep", files=TRUE, absolutePaths=TRUE) {
# arborescence should be "path_to_sim/PARAM=VALUE/names_of_subsimulation/outputs/out_repN
	if(!absolutePaths&&!files) { message("\nWhen outputting directories (files=FALSE) absolutePaths is irrelevant : full path is provided\n") }

	pathS <- list.dirs(path_to_sim, recursive=TRUE, full.names=FALSE)
		## full.names=FALSE to avoid grepping words present higher in the arborescence "outputs" from the results
			# this also allows to have a matrix 2 lines down further down instead of a list of vectors...
	pathS <- grep("output", unlist(strsplit(pathS, path_to_sim)), value=TRUE)
	pathS <- paste(path_to_sim, pathS, sep="/")
	# browser()
		for(ii in names) { if(! sum(grepl(ii, pathS))>0) { stop("Some \'names\' provided have not been found in the path to simulation. Stopping execution\n") } }
	rdr <- sapply(names, function(x) which(grepl(x, pathS)))
	pathS <- pathS[rdr]	# grouped by "VAR" then by "subcdt" -> grouped by "subcdt"
	if(!is.matrix(rdr)) {
		names(pathS) <- names
		if(! files) {
			return(pathS)
		} else {
			ans <- lapply(names, function(x) list.files(pattern=filePtrn, path=pathS[[x]], full.names=absolutePaths) )
			names(ans) <- names
			return(ans)
		}
	} else { #multiple "param=value"
		rdr <- sapply(names, function(x) which(grepl(x, pathS))) # After re-ordering pathS, rdr needs to be reobtained for further use...

		VarVal <- list.dirs(path_to_sim, recursive=FALSE)
		VarVal <- unname(sapply(VarVal, function(x) strsplit(basename(x), "=")[[1]] ))
		if(! all(VarVal[1,]==VarVal[1,1])) stop("multiple parameters on the first level of the arborescence") # if length(unique(VarVal[1,]))==1
			VAR <<- VarVal[1,1] ; VAL <<- VarVal[2,] #useless
		MMM <- (sapply(VAL, function(y) apply(rdr, 1, function(x) grepl(y, pathS[x]))))
			# matrix of "being in "VAR=VAL" folder"

# MMM is only used to get the right names (list elements) associated to the right paths
## the re-ordering is done after this list is produced

######## IF MMM contains errors :
if(length(pathS)!=nrow(MMM) | ncol(MMM)!=nrow(rdr) | !all( colSums(MMM) == (nrow(MMM)/ncol(MMM)) ) | !all(rowSums(MMM)==1) ) {
	## if some patterns of "VAL" overlap : ex "1000" is found in "10000"
		## all cells in which the rowSums is > 1 or colSums > (nrow(MMM)/ncol(MMM)) [[why not "length(subsim)" ?]] are due to this error

		## in these cells, set TRUE to false :
	ggg <-  expand.grid( which(rowSums(MMM)>1) , ccc <- which(colSums(MMM)>nrow(rdr)) )
		for(ii in 1:nrow(ggg)) { MMM[ggg[ii,1], ggg[ii,2]] <- FALSE }
# get "the right" diagonal (or else the names used later will be inversed)
	while( !all(diag(unique(MMM))) ) {
# identify which columns mapped multiple time, and thus need to be swapped to get "the right" diagonal
		toTest <- sapply(VAL, grepl, VAL)
		toChange <- which(colSums(toTest)>1)
		changeTo <- which(rowSums(toTest)>1)
		toTest <- expand.grid(toChange, changeTo)
	newMMM <- MMM
		for(sss in 1:nrow(toTest)) {
			# try to swap.
			newMMM[,unlist(toTest[sss,])] <- newMMM[,rev(unlist(toTest[sss,]))]
		# as soon as a new reordering has got "a better diagonal" : save it
			if(sum(diag(unique(newMMM)))>sum(diag(unique(MMM)))) { colnames(newMMM)[rev(unlist(toTest[sss,]))] <- colnames(newMMM)[unlist(toTest[sss,])] ; MMM <- newMMM ; break
#										^ MOST IMPORTANTLY : CHANGE THE NAMES 
		# as swappings are saved, if a swapping isn't better, it needs to be reversed
			} else { newMMM <- MMM }
		}
# restart until the diagonal is restored
	}
} ###### "REPAIR" MMM

	stopifnot(length(pathS)==nrow(MMM)) ; stopifnot(ncol(MMM)==nrow(rdr))
	# MMM's purpose is only to get the names in the right order.
	rownames(rdr) <- colnames(MMM)
	# Some of these tests might be useful in some conditions. As only the diagonal is important for name reordering, they have been deactivated
#~	stopifnot(all( colSums(MMM) == (nrow(MMM)/ncol(MMM)) )) ; stopifnot(all(rowSums(MMM)==1))

		if(! files) {
			ans <- lapply(1:length(rownames(rdr)), function(YY) {
				lapply(1:length(colnames(rdr)), function(XX) {
					pathS[rdr[YY,XX]] })})
		} else {
			ans <- lapply(1:length(rownames(rdr)), function(YY) {
				lapply(1:length(colnames(rdr)), function(XX) {
					list.files(pattern=filePtrn, path=pathS[rdr[YY,XX]], full.names=absolutePaths) })})
		}
		names(ans) <- rownames(rdr)
		for(i in (1:length(ans)) ) { names(ans[[i]]) <- colnames(rdr) }
		return(ans)
	}
}


##############################################
### Functions used to Extract Results
##############################################

#Give list of W for every or selected generations
read.W <- function(simfile, gen=NA, last.only=FALSE, first=TRUE, pattern="MeanAll", generReplace=FALSE) {
# Retrieve values as a matrix. Used to get the average interaction matrix of a population.
# The "generReplace" argument allows replacement of generation number with line number. It is used when multiple generations share the same identification number as a result of individual sampling.
	stopifnot(
		is.character(simfile),
		file.exists(simfile)
	)

	daata <- read.table(simfile, header=TRUE)
	if(generReplace) { daata$Gen <- 0:(length(daata$Gen)-1) }
	if (is.na(gen)) {
		gen <- daata$Gen
		if (!first) gen <- gen[-1]
		if (last.only) gen <- gen[length(gen)]
	}
	ans <- list()
	for (gg in gen) {
		mat <- as.numeric(daata[which(daata$Gen == gg), grep(pattern, colnames(daata))])
		if (length(mat) == 0) stop("No field called \"", pattern, "\" in file \"", simfile, "\".")
		mat <- matrix(mat, ncol=sqrt(length(mat)), byrow=TRUE) #byrow=FALSE to have t(W)
		ans[[as.character(gg)]] <- mat
	}
	return(ans)
}


values_txt <- function(fiile, pattern){
# Retrieve all the values of a simul file, pattern = name/pattern of the column wanted
	val <- read.table(fiile, header=TRUE)
	val <- as.matrix(val[, grep(pattern, colnames(val))])
	return(val)
}


##############################################
### Functions used to simulate (plastic) development
##############################################

sigma.M2 <- function(x, a) {
# Sigmoid scaling function
	1. / (1. + exp((-x/(a*(1-a)))+log(1/a-1)))
}


sigma.M2p <- function(x, lambda, mu) {
# Sigmoid scaling function optimized for calculus (use values pre-calculated from the fixed "a")
	# lamba = (1-a)/a	# mu = (1/(a*(1-a)))
	1. / (1. + lambda * exp(-mu*x))
}


model.M2 <- function(W, a=0.2, S0=rep(a, nrow(W)), steps=24, measure=4, full=FALSE, varFUN=function(x) mean((x-mean(x))^2)) {
# Simulates development from an interaction matrix and an initial state
    measure=measure-1	# so that "measure" and "steps" both 'start counting' at 1
    lambda <- (1-a)/a
    mu <- 1/(a*(1-a))
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0 <- sigma.M2p((W %*% S0), lambda=lambda, mu=mu) 			
		sto[,i+1] <- S0
	}
	ans <- list()
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, varFUN)
	if (full) ans$full <- sto
	return(ans)
}


model.plasti <- function(W, plastrate=c(0), plastsignal=c(0), a=0.2, S0=rep(a, nrow(W)), steps=24, measure=4, full=FALSE) {
# Simulates development with a plastic signal
	measure=measure-1	# so that "measure" and "steps" both 'start counting' at 1
		# note that S0 = sto[,1], hence the modification
	stopifnot(length(plastrate)==length(plastsignal))	# although annoying, it is good control.
	stopifnot(length(S0)==nrow(W))
	stopifnot( (length(which(plastrate>1)) | length(which(plastrate<0)) | length(which(plastsignal<0)) | length(which(plastsignal>1)) )==0) #rmq... on peut juste ne pas inputter de la merde
	plastic=which(plastrate!=0)

	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
		#Re-evaluate S0 from standard S0 with plasticity
		S0[plastic] <- S0[plastic] + ((plastsignal[plastic]-S0[plastic])*plastrate[plastic])
	St <- S0
	sto[,1] <- S0
	for (i in 1:steps) {
		St <- sigma.M2(W %*% St, a)
		if(length(plastic)!=0){
			St[plastic] <- St[plastic] + ((plastsignal[plastic]-St[plastic])*plastrate[plastic])
		}
		sto[,i+1] <- St
	}
	ans <- list()
	if (measure == 0) {
		ans$val <- sto[,(steps-measure)]
	} else {
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, function(x) mean((x-mean(x))^2))
	}
	if (full) ans$full <- sto
	return(ans)
}


plasticity_index <- function(W, opt, index, signalEnv=seq(0,1,0.05), Ssel=c(0, 10, 10, 10, 10, 10), Sprime=46000, outputPheno=TRUE, displayWarnings=TRUE, model=model.plasti, plastrate=c(1,0,0), ...) {
    # Wrapper for model.plasti.
    # Calculates development in multiple environments. May be informed with a selection function with a different optimum in each environment. May output fitnesses in every environment.
    			# multiple indexes :
    			#	"fit" uses strength of selection, optima and expression to assess a fitness score in each "signalEnv" value
    			#	"dist" outputs the distance to the optima
    			#	"none" outputs no specific index. An error if thrown if index="none" and "outputPheno=FALSE" ==> no output.
    
    if(index=="none" & !outputPheno) { stop("No output selected. Choose another \"index\" or set \"outputPheno\" to \"TRUE\"") }
    	### Multiple check-ups. No specific warnings.
    stopifnot(ncol(W)==nrow(W))
    N_loc <- ncol(W)
    N_env <- length(signalEnv)
    stopifnot( sum(index %in% c("dist", "fit", "none"))==1 )
    
    	### 1) make a matrix of [signal, optima]
    if(!is.character(opt)){ stop("opt should be an expression as a function of \'x\', given as a string. The optima of every genes should be separated by a \";\" ")  }
    OPTIMA <- unlist(strsplit(x=opt, split=";"))
    # OPTIMA will fill with 0 to the size of N_loc
    while(length(OPTIMA) < N_loc) { OPTIMA[length(OPTIMA)+1] <- 0 ; if(displayWarnings){message("length(opt) < nrow(W) : auto-filling optimum with 0")} }
    
    OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(signalEnv, function(x) eval(parse(text=OPTIMA[oo])) ) )
    	# outputs a matrix of [N_env, gene_optima]
    
    	### 2) Development with every environmental signal
    # plastrate autofills with 0.
    while(length(plastrate) < N_loc) { plastrate[length(plastrate)+1] <- 0 ; if(displayWarnings){message("length(plastrate) < nrow(W) : auto-filling plasticity rate with 0")} }
    
    Pheno <- sapply(1:N_env, function(ss) model(W=W, plastsignal=OPTIMA[ss,], plastrate=plastrate, ...) )
    
    Pheno <- list( t(sapply(1:N_env, function(x) unlist(Pheno[1,x]) )), t(sapply(1:N_env, function(x) unlist(Pheno[2,x]) )) )
    names(Pheno)=c("mean", "var")
    	# Pheno is a list of matrices which give : [[1]]=> mean, [[2]]=> var ; for every combination of [signal,gene]
    
    	# Strength of stabilizing selection autofills unspecified values to the size of N_loc with 0
    while(length(Ssel) < N_loc) { Ssel[length(Ssel)+1] <- 0 ; if(displayWarnings){message("length(Ssel) < nrow(W) : auto-filling strength of selection on expression with 0")} } 
    
    if(index=="fit") {
    	# Strength of selection on developmental stability autofills unspecified values to the size of N_loc with itself (the first value]
    	while(length(Sprime) < N_loc) { Sprime[length(Sprime)+1] <- Sprime[1] ; if(displayWarnings){message("length(Sprime) < nrow(W) : auto-filling strength of selection on stability with the first value")} } 
    
    	# Fitness calculation
    	ansVEC <- sapply(1:N_env, function(ee) { prod(exp(-Ssel*(Pheno$mean[ee,] - OPTIMA[ee,])^2)) * prod(exp(-Sprime*Pheno$var[ee,])) } )
    	names(ansVEC) <- signalEnv
    		# returns a vector of fitness in every "N_env" condition
    
    } else if (index=="dist") {
    
    	if(sum(Ssel!=0)==0) stop("No selected gene : can't calculate distance to optimum")
    	if(sum(OPTIMA[1,-1]!=0) != sum(Ssel!=0) && displayWarnings  ) {message("\nThere might be an inconsistency between opt[-1] and Ssel\n",
    							"This may be due to selection for null gene expression (gene expression extinction) ; PROCEEDING with distance calculation\n")}
    
    	#ansVEC <- colMeans(sapply(1:N_env, function(ee) { abs(Pheno$mean[ee,] - OPTIMA[ee,]) } ))
    	ansVEC <- sapply(1:N_env, function(ee) { abs(Pheno$mean[ee,(Ssel!=0)] - OPTIMA[ee,(Ssel!=0)]) } )
    
    		if(is.vector(ansVEC) && !is.matrix(ansVEC) ) { ansVEC <- t(as.matrix(ansVEC)) }
    	colnames(ansVEC) <- signalEnv
    	ansVEC <- colMeans(ansVEC)
    		# returns a vector of distance to optimum in every "N_env" condition
    }
    
    # Determination of Return
    if(outputPheno){
    	if(index != "none") {
    		ans <- list(ansVEC, Pheno$mean, Pheno$var)
    		names(ans)=c(index, "mean", "var")
    	} else {
    		ans <- list(Pheno$mean, Pheno$var)
    		names(ans)=c("mean", "var")
    	}
    }else{
    #~	ans <- ansVEC	# old output... added a list level for homogeneity of outputs.
    	ans <- list(ansVEC) ; names(ans) <- index
    }
    return(ans)
}


labile.plasti <- function(W, signalEnv, outputPheno=TRUE, plastrate=c(1, rep(0,19)), a=0.2, firstSzero=delayedAssign("firstSzero", rep(a, N_loc)), ...) {
# Used to develop a character as a labile character.
#   (i.e, starting the development in a next environment from the last value of the previous environment)

	# in the ellipse are the arguments of "model.plasti"
## This function can, in its state, only support one plastic gene which is the first gene.
	# signal should be a vector, of all used signals (and not "of the signal for every gene")

	stopifnot(exists("model.plasti"))
	stopifnot(ncol(W)==nrow(W))
	N_loc <- ncol(W)
stopifnot(length(plastrate)==N_loc)

	eval(firstSzero)
	Szero <- firstSzero

	ans <- list(c(), c(), c(), c())
	names(ans) <- c("envs", "means", "vars", "devs")

	randEnv <- signalEnv[sample(order(signalEnv))]
for(EE in 1:length(signalEnv)) {
	res <- model.plasti(W=W, plastrate=plastrate, plastsignal=c(randEnv[EE], rep(0, N_loc-1)), S0=Szero, a=a, full=outputPheno, ...)
	Szero <- res$mean

	if(outputPheno) {
		naames <- c()
		for(nnn in 0:(ncol(res$full)-1)) { naames <- append(naames, paste0("step", nnn, "_env", randEnv[EE])) }
			colnames(res$full) <- naames
			ans$devs <- cbind(ans$devs, res$full)
	}
		ans$means <- cbind(ans$means, res$mean)
			colnames(ans$means)[ncol(ans$means)] <- paste0("env",randEnv[EE])
		ans$vars <-  cbind(ans$vars, res$var)
			colnames(ans$vars)[ncol(ans$vars)] <- paste0("env",randEnv[EE])

		ans$envs <- append(ans$envs, randEnv[EE])

}
	if(!outputPheno) { ans[["devs"]] <- NULL }
#~	if(!outputPheno) { ans <- list(res$mean, res$var) ; names(ans) <- c("lastMeans", "lastVars") } #else { names(ans$envs) <- ans$envs }
return(ans)
}


labileVSdev <- function(W, signalEnv, OptFenv, outputPheno=FALSE, thresh="no", optSTABstr, devSTABstr, plastrate, ...) {
# Used to compare independent development in multiple environments ; with sequential developments in these environments

	# In the ellipse are the arguments of plasticity_index
	# thresh can be "no" (output the difference between the two geometric fitnesses)
		# or a numeric value (will output a boolean depending on whether the fitness difference is or not under the threshold)

	# OUTPUT IS INDEPENDENT FITNESS - DEPDENDENT FITNESS
	# 	=> So positive if independent fitness is higher than fitness for a labile character. And vice versa.

	geoMean <- function(x) { stopifnot(is.vector(x) && length(x)>1) ; prod(x)^(1/length(x)) }
	stopifnot(exists("plasticity_index")) ; stopifnot(exists("labile.plasti"))
	stopifnot(is.numeric(thresh) | thresh=="no")

# Make a development :
	# In all environments INDEPENDENTLY
indepPhen <- plasticity_index(W=W, index="fit", opt=OptFenv, signalEnv=signalEnv, outputPheno=outputPheno, Ssel=optSTABstr, Sprime=devSTABstr, plastrate=plastrate, ...)
	# In all environments ONE AFTER THE OTHER
SequPhen <- labile.plasti(W=W, signalEnv=signalEnv, outputPheno=outputPheno, plastrate=plastrate)

# Obtain fitness optima (for the sequential environments)
OPTIMA <- unlist(strsplit(x=OptFenv, split=";"))
	while(length(OPTIMA)<length(as.vector(tail(t(SequPhen$means),1)))) { OPTIMA <- append(OPTIMA, 0) }
OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(signalEnv, function(x) eval(parse(text=OPTIMA[oo])) ) )

	### reorder fitness optima in the (randomly drawn) order used in seqPhen
if(all( order(OPTIMA[,1]) == (1:length(OPTIMA[,1])) )) OPTIMA <- OPTIMA[order(order(SequPhen$envs)),]

# produce a value of fitness in every environment (for an individual who went through all environments sequentially)
SequFit <- sapply(1:ncol(SequPhen$means), function(ee) { prod(exp(-optSTABstr*( SequPhen$means[,ee] - OPTIMA[ee,])^2)) * prod(exp(-devSTABstr*SequPhen$vars[,ee])) } )
	# calculate geometric fitness with these fitness values
SequFit <- geoMean(SequFit)

# produce a value of fitness in every environment (for an individual who went through all environments independently)
	# calculate a geometric fitness with the independent fitness values
if(outputPheno==FALSE) { indepFit <- geoMean(indepPhen) } else { indepFit <- geoMean(indepPhen[["fit"]]) }

	if(thresh=="no") { ans <- indepFit - SequFit
	} else {  ans <- (indepFit - SequFit) < thresh }

return(ans)
}


##############################################
### Functions used for Data Treatment
##############################################

getSlope <- function(envGeneMatrix, envirVal, targetGene, plastic=1) {
# Get a slope from the output of the function simulating a development with plasticity (i.e.: an "envGeneMatrix")

	# envGeneMatrix is a matrix giving expression of : [environment, gene]
	# envirVal should be all the environments of interest
	# targetGene should be all the genes of interest
	stopifnot(is.matrix(envGeneMatrix))
	stopifnot(length(envirVal)>=2)
	stopifnot(length(targetGene)>=1)

	envLines <- envGeneMatrix[,plastic] %in% envirVal 

	ans <- lm(envGeneMatrix[envLines, targetGene] ~ envGeneMatrix[envLines, plastic])$coefficients
# if length(targetGene)>1 => outputs a matrix of [(intercept;slope),Gene]
	if(is.matrix(ans)) { ans <- ans[2,] } else { ans <- ans[2] }

	names(ans) <- paste0("slope", targetGene)
return(ans) }


geoMean <- function(x) {
# Calculate geometric mean
	stopifnot(is.vector(x) && length(x)>1)
	prod(x)^(1/length(x))
}


eucliDist <- function(A, B){
# Outputs euclidean distance between two whole matrices (works fine with vectors or scalars), A and B
	sqrt(sum((A-B)^2))
}


##############################################
### Miscellaneous graphic functions
##############################################

lighten <- function(color, factor=0.5){
# Lightens a color
    col <- col2rgb(color)
    col <- col+(255-col)*factor
    col <- rgb(t(col), maxColorValue=255)
    col
}


alphize <- function(color, alpha=120){
# Makes a color transparent
	col <- col2rgb(color)
	col <- rbind(col, alpha, paste0("alphized_", color), 255)
		# rgb arguments, in the order, are : r, g, b, alpha, NAMES, maxColorValue
	col <- do.call(rgb, as.list(col))
return(col)
}

# Vectorized function : (It will not use all combinations of arguments, but a 1:1 correspondance)
alphize <- Vectorize(alphize, USE.NAMES=FALSE)


##############################################
### Functions used to Plot Results
##############################################

linesArea <- function(X,Y, rowNumAvg=11,rowNumBot=2,rowNumTop=20, colz,alpha, ELty=1, ELwd=1.25, lines=TRUE, polygone=TRUE) {
	colTrans <- alphize(colz, alpha)
	if(polygone) polygon(x=c(X, rev(X)), y=c(Y[rowNumBot,], rev(Y[rowNumTop,])), col=colTrans, border=NA)
	if(lines) lines(x=X, Y[rowNumAvg,], type="l", lty=ELty, lwd=ELwd, col=colz)
}

linesBars <- function(X,Y, rowNumAvg=11,rowNumBot=2,rowNumTop=20, Subset=NULL, colz,alphaz=120, ELty=1,ELwd=1, lines=TRUE, bars=TRUE) {
	if(bars) {
		colTrans <- alphize(colz, alphaz)
		if(!is.null(Subset)) { arrows(x0=X[Subset], y0=Y[rowNumBot,Subset], y1=Y[rowNumTop,Subset], col=colTrans, length=0.0, code=3, lwd=ELwd)	# by default x1=x0
		} else if(is.null(Subset)) { arrows(x0=X, y0=Y[rowNumBot,], y1=Y[rowNumTop,], col=colTrans, length=0.0, code=3, lwd=ELwd) }
	}
	if(lines) lines(x=X, Y[rowNumAvg,], type="l", col=colz, lty=ELty, lwd=ELwd)
}


