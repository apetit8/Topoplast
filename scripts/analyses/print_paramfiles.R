#!/usr/bin/env Rscript
source(here::here("scripts/analyses/paths.R"))
source(here::here("scripts/functions/functions.R"))
source(here::here("scripts/functions/makeparam_funct.R"))

# This script generates the parameter files for the simulations. 

# It consists in two parts that are almost identical; 
# part 1 for the full network simulations
# part 2 for the control (drift) simulations
# note that path variables are overwritten between part 1 and part 2

overwrite <- TRUE
# options(warn=1)
prog.path <- path$prog

sd <- 0.25  #0.25
min <- 0.15 #min environment
max <- 0.85 #max environment

####### FULL NETWORK SIMULATIONS ################

param.dir      <- path$param.dir$fullnet
directory      <- path$sim.dir$fullnet
launchfilename <- path$launch.file$fullnet

ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))


####
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file (used to generate different environmental file each generations).
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)

all.sims <- rbind(
  pops = c("param1.txt", "extparam1.txt")
)
#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=FALSE, sd=sd, min=min, max=max)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=FALSE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

######### DRIFT SIMULATIONS #####################

param.dir      <- path$param.dir$drift
directory      <- path$sim.dir$drift
launchfilename <- path$launch.file$drift

ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file (used to generate different environmental file each generations).
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist

ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)

all.sims <- rbind(
  pops = c("param1.txt", "extparam1.txt")
)
#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=FALSE, sd=sd, min=min, max=max)
    
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=FALSE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

