# To be "Sourced" from a R script file

suppressPackageStartupMessages(require(R.utils)) # Path manipulation 

library(parallel)
library(digest)
library(ade4)
library(igraph)
library(Rcpp)
library(inline)
library(truncnorm)

tar.param <- function(param.files, compressed.file, compression="gzip") {
  #~ 	tar(compressed.file, files=param.files, compression=compression)
  # There can be many param files. just store the list in a file to avoid too long command lines
  pfile <- file.path(dirname(param.files[1]), "allfiles.txt")
  cat(basename(param.files), file=pfile, sep="\n")
  command <- paste0("cd ", dirname(param.files[1]), "&& tar cfz ", path.expand(compressed.file), " -T ", basename(pfile))
  system(command)
  unlink(param.files)
  unlink(pfile)
}

untar.param <- function(compressed.file, remove.tar=FALSE) {
  #~ 	untar(compressed.file)
  system(paste0("cd ", dirname(compressed.file), "; tar xfz ", basename(compressed.file)))
  if (remove.tar)
    unlink(compressed.file)
}

create.launchfile <- function(prog.path, param.files, output.files, compressed.files, launch.file="./launchfile.sh", relative.paths=TRUE) {
  # relative.paths sets all paths relative to the launch.file directory
  stopifnot(
    length(param.files) > 0, 
    length(param.files) == length(output.files),
    # file.exists(prog.path), 
    all(file.exists(param.files) | file.exists(compressed.files)))
  
  
  if (relative.paths) {
    launch.dir <- dirname(launch.file)
    param.files <- getRelativePath(param.files, launch.dir)
    output.files<- getRelativePath(output.files, launch.dir) 
    compressed.files<- getRelativePath(compressed.files, launch.dir)
  }
  
  command <- paste0( prog.path, " -z ", compressed.files, " -p ", param.files, " -o ", output.files, " -Q ", output.files,".pop"  )
  writeLines(command, launch.file)
}

create.launchfile.alt <- function(prog.path, param.files, output.files, oldpop="pevpopfile", compressed.files, launch.file="./launchfile.sh", relative.paths=TRUE, prevpop=TRUE, sampling=" -i ") {
  # relative.paths sets all paths relative to the launch.file directory
  #sampling : " -I " for every generations ; " -i " for last gen ; "" for no sampling
  
  stopifnot(
    length(param.files) > 0, 
    length(param.files) == length(output.files),
    all(file.exists(param.files) | file.exists(compressed.files)))

  if (relative.paths) {
    launch.dir <- dirname(launch.file)
    param.files <- getRelativePath(param.files, launch.dir)
    output.files<- getRelativePath(output.files, launch.dir) 
    compressed.files<- getRelativePath(compressed.files, launch.dir)
  }
  if(prevpop==TRUE){ command <- paste0( prog.path, " -z ", compressed.files, " -p ", param.files, " -o ", output.files, sampling, " -P ",oldpop, " -Q ", output.files,".pop\n"  )} 
  else command <- paste0( prog.path, " -z ", compressed.files, sampling, " -p ", param.files, " -o ", output.files, "\n"  ) #" -Q ", output.files,".pop",
 
  cat(command, file=launch.file, append=TRUE)
}

create.launchfile.old <- function(prog.path, param.files, output.files, compressed.files, launch.file="./launchfile.sh", relative.paths=TRUE) {
  # relative.paths sets all paths relative to the launch.file directory
  stopifnot(
    length(param.files) > 0, 
    length(param.files) == length(output.files),
    file.exists(prog.path), 
    all(file.exists(param.files) | file.exists(compressed.files)))
  
  
  if (relative.paths) {
    launch.dir <- dirname(launch.file)
    prog.path   <- getRelativePath(prog.path, launch.dir)
    param.files <- getRelativePath(param.files, launch.dir)
    output.files<- getRelativePath(output.files, launch.dir) 
    compressed.files<- getRelativePath(compressed.files, launch.dir)
  }
  
  uncompress.command <- ifelse(is.na(compressed.files), "", paste0("tar xfz ", compressed.files, " -C ", dirname(compressed.files), " && "))
  #~ 	compress.command <- ifelse(is.na(compressed.files), "", paste0("tar cfz ", compressed.files, " ", dirname(compressed.files), "/param*.txt;"))
  compress.command <- ifelse(is.na(compressed.files), "", paste0("rm -f ", dirname(compressed.files), "/param*.txt;"))
  
  commands <- paste0(
    uncompress.command,
    prog.path, " -p ", param.files, " -o ", output.files, " && " ,
    compress.command)
  writeLines(commands, launch.file)
}

read.param <- function(param.file) {
  # The structure of the parameter file is assumed to be, for each line:
  # PARNAME x y
  # The function returns a list named by the parameter names
  
  lines <- readLines(param.file)
  lines <- lines[!grepl(lines, pattern="^#")]
  splines <- strsplit(lines, split="\\s+")
  nn <- sapply(splines, "[", 1)
  ans <- lapply(splines, function(ll) {
    if (length(ll) <= 1) return(NA)
    ll <- ll[-1]
    numll <- suppressWarnings(as.numeric(ll))
    if (any(is.na(numll))) 
      return(ll)
    else 
      return(numll)
  })
  names(ans) <- nn
  ans
}

write.param <- function(param.file, param) {
  content <- paste(names(param), sapply(param, paste, collapse="\t"), sep="\t")
  writeLines(content, param.file)
}


#make.randopt but to make random or custom optimum RN
make.randopt <- function(oldopt, pattern, begin.bottleneck=FALSE, random=TRUE, oldenv=0.5, sd=0.05, min=0.2, max=0.8, RN=list(c(0,0))) {
  # Returns a vector of random optima according to the provided pattern
  # (constant vs fluctuating optima)
  stopifnot(length(oldopt) == length(pattern),
            all(pattern %in% 0:5))
  newopt <- oldopt
  if(random==TRUE) env <- runif(1, min, max)
  # env <- rnorm(1, mean=oldenv, sd=0.05, max=0.2, min=0.8)
  if(random==FALSE) env <-  rtruncnorm(n=1, a=min, b=max, mean=oldenv, sd=sd)
  newopt[pattern == 0] <- 0
  newopt[pattern == 1] <- oldopt[pattern == 1] # does nothing, but clearer this way
  newopt[pattern == 2] <- env
  if (begin.bottleneck) newopt[pattern == 3] <- runif(sum(pattern == 3))
  newopt[pattern == 4] <- (1-env)
  #Loop for random plastic RN
  for (i in 1:length(newopt[pattern == 5])) {
    newopt[pattern == 5][i] <- (env*RN[[i]][1])+RN[[i]][2]
  }
  newopt
}


make.selstr <- function(pattern) {
  # Mirrors the make.random function, but for the selection strength. 
  # The algorithm is trivial, just sets zero for non-selected genes
  # (including the first one = signal). 
  
  ans <- pattern != 0
  ans[1] <- 0
  ans
}


#Remove every manipulation of the fitness strength as we only use the pre-burn part of this function.
create.paramseries <- function(param.template.file, extparam.file, simul.dir, overwrite=FALSE,
                               verbose=FALSE, allow.extrapar=c("GENET_MUTRATES"), tar.param=TRUE, mc.cores=detectCores()-2,
                               sd=0.05, min=0.15, max=0.85, random=FALSE, manual=NULL) 
  # This is the main algorithm that create the simulation structure. 
  # The function retruns the necessary information to make a launchfile
  # (it does not write the launchfile, because all information is not available here)
  # Two parameter files are necessary : 
  #  - the parameter template, which is a regular parameter file that will be used for the default settings
  #  - the "extended" parameter file, which provides "meta" information to run the simulation program
  #    this is specific for the domestication project, and involves three simulation stages:
  #     . before the bottleneck (burn-in)
  #     . during the bottleneck
  #     . after the bottleneck
  # It is expected that the information from the two parameter files match (e.g. number of genes), 
# consequences of a non-matching pattern are undefined. 
# The algorithm reads the two parameter files and create the necessary files for simulations
# For each simulation replicate (REPLICATES variable in the ext.par file), a simulation directory
# is created, and as many parameter files as generations are created within each directory
#  Note: this represents a very large number of (small) files, which can represent a burder for 
#  the system, especially when creating the files or compressing the directory
# The first parameter file of a simulation is complete, the next ones only change what is necessary:
#   . every generation: the environmental variable and the fluctuating optima
#   . before and after the bottleneck: the population size
#   . at the beginning of the bottleneck: the "domestication" selection pattern
# The directory and parameter naming pattern is specified in a few internal functions. 
{
  ndigits.rep <- 4
  ndigits.gen <- 5
  sf.format.rep <- paste0("%0", ndigits.rep, "d")
  sf.format.gen <- paste0("%0", ndigits.gen, "d")
  .repID <- function(rep) #, ndigits=4)
    sprintf(sf.format.rep, rep)
  .genID <- function(gen) #, ndigits=5)
    sprintf(sf.format.gen, gen)
  .repDir  <- function(rep, rep.template = "rep") 
    paste0(rep.template, "-", .repID(rep))
  .repFile <- function(rep, gen, rep.template = "rep", gen.template="gen")
    paste0("param-", rep.template, .repID(rep), "-", gen.template, .genID(gen), ".txt")
  .repTarFile <- function(rep, rep.template = "rep")
    paste0("param-", rep.template, .repID(rep), ".tar.gz")
  .outFile <- function(rep, rep.template = "rep")
    paste0("out-", rep.template, .repID(rep), ".txt")
  stopifnot(file.exists(param.template.file), file.exists(extparam.file))
  param.template <- read.param(param.template.file)
  extparam       <- read.param(extparam.file)
  myparam <- param.template
  extrapar <- extparam[names(extparam) %in% allow.extrapar]
  # A few shortcuts to make the code more readable
  totgen <- extparam$GENERATION_MAX
  # Checks for the consistency of both parameter files
  if (!"SIMUL_MAXGEN" %in% names(myparam)) 
    myparam$SIMUL_MAXGEN <- totgen

  stopifnot(
    length(extparam$SCENARIO_PART1) == myparam$GENET_NBLOC , 
    length(extparam$SCENARIO_PART1) == length(extparam$SCENARIO_PART2))
  stopifnot(extparam$REPLICATES > 0)
  # Now the parameter files:
  dir.create(file.path(simul.dir), showWarnings = FALSE)
  # It looks difficult to make it simpler: parameter files need to be updated
  # on a case-by-case basis, as the parameters to change depend on the 
  # simulation stage (before, during, of after the bottleneck)
  if (verbose)
    pb <- txtProgressBar(min = 1, max = extparam$REPLICATES, style = 3)
  # browser()
  mclapply(1:extparam$REPLICATES, function(rep) { #mclapply
    repdir <- file.path(simul.dir, .repDir(rep))
    compressed.file.name <- file.path(repdir, .repTarFile(rep))
    if (exists(compressed.file.name))
      untar.param(compressed.file.name)
    dir.create(repdir, showWarnings = FALSE)
    if (overwrite) { # This is quite powerful, use with caution (simulation results are deleted prior to launching a new sim)
      unlink(list.files(path=repdir, full.names=TRUE))
    }
    #Draw random RNs for "manual" plasticity
    RN <- lapply(1:length(which(extparam$SCENARIO_PART1==5)), function(i){
      RNslope <- sample(c(runif(1,0.5718,1.15), runif(1,-1.15,-0.5718)), 1) #(75% chance to have negative correlation)
      if(RNslope > 1) RNintercept <- runif(1, 1-RNslope, 0) else if(RNslope > 0) RNintercept <- runif(1, 0, 1-RNslope ) else if(RNslope < -1) RNintercept <- runif(1, 1, -RNslope) else RNintercept <- runif(1, -RNslope, 1)
      return(c(RNslope, RNintercept))
    })
    # First generation: use the full template file
    optim <- make.randopt(runif(myparam$GENET_NBLOC, min, max), extparam$SCENARIO_PART1, sd=sd, min=min, max=max, random=random, RN=RN)
    #Replace with chosen optimum if wanted
    if(is.null(manual)==FALSE){optim[(length(optim)-length(manual)+1):length(optim)] <- manual}
    myparam$FITNESS_OPTIMUM <- optim
    myparam$FILE_NEXTPAR <- suppressWarnings(normalizePath(file.path(repdir, .repFile(rep, 1))))
    par.file.name <- file.path(repdir, .repFile(rep, 0))
    if (!file.exists(par.file.name[1]) || overwrite)
      write.param(par.file.name[1], myparam)
    # From here, just update what is necessary
      for (gen in 1:(totgen-1)) {
        par.file.name <- c(par.file.name, file.path(repdir, .repFile(rep, gen)))
        if (!file.exists(par.file.name[gen+1]) || overwrite) {
          optim <- make.randopt(optim, extparam$SCENARIO_PART1, oldenv=optim[1], sd=sd, min=min, max=max, random=random, RN=RN) ## Here : give back the prev env.
          if(gen == totgen-50){ #More outputs before the end of the simulation (to infer reaction norms during results analyses)
            write.param(par.file.name[gen+1], list(FITNESS_OPTIMUM = optim, 
                             FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(repdir, .repFile(rep, gen+1)))),
                             SIMUL_OUTPUT = 1))    }
          else write.param(par.file.name[gen+1], list(FITNESS_OPTIMUM = optim, 
                                FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(repdir, .repFile(rep, gen+1))))))
        }}
    extparlast <- extparam$SCENARIO_PART1
    gen <- gen+1
    
    # Very last generation (no NEXTPAR)
    par.file.name <- c(par.file.name, file.path(repdir, .repFile(rep, gen)))
    if (!file.exists(par.file.name[gen+1]) || overwrite) {
      optim <- make.randopt(optim, extparlast, oldenv=optim[1], sd=sd, min=min, max=max, random=random, RN=RN)
      write.param(par.file.name[gen+1],
                  list(FITNESS_OPTIMUM = optim))
    }
    gc(verbose = FALSE)
    stopifnot(all(!duplicated(par.file.name))) # just to check if everything is all right... 
    if (tar.param)
      tar.param(par.file.name, compressed.file.name)
    if (verbose) setTxtProgressBar(pb, mc.cores * (rep %/% mc.cores))
  }, mc.cores=mc.cores)
  if (verbose) cat("\n")
  # Returns the parameter (and output) file names that will be necessary to make the launchfile
  return(list(param=file.path(simul.dir, .repDir(1:extparam$REPLICATES), .repFile(1:extparam$REPLICATES, 0)), 
              out =file.path(simul.dir, .repDir(1:extparam$REPLICATES), .outFile(1:extparam$REPLICATES)), 
              compressed = if (tar.param) file.path(simul.dir,  .repDir(1:extparam$REPLICATES), .repTarFile(1:extparam$REPLICATES)) else NA))
}
