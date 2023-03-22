#!/usr/bin/env Rscript
source("scripts/functions/functions.R")
source("scripts/functions/makeparam_funct.R")
##########################
script.dir <- getwd()
user.dir  <- getwd()
param.dir <- file.path(script.dir, "templates/4g")
ifelse(!dir.exists(file.path("simul/4g")), dir.create("simul/4g"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path(script.dir, "simul/4g"))
src.dir   <- file.path(script.dir, "functions")

overwrite <- FALSE
# options(warn=1)
prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
random=TRUE #New environment randomly drawn or not
min <- 0.15 #min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/4g-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)


all.sims <- rbind(
  Correlated_Down = c("Down/param1.txt", "Down/Corr_extparam1.txt"),
  Correlated_Up = c("Up/param1.txt", "Up/Corr_extparam1.txt"),
  Correlated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Corr_extparam1.txt"),
  Anticorrelated_Down = c("Down/param1.txt", "Down/Anticorr_extparam1.txt"),
  Anticorrelated_Up = c("Up/param1.txt", "Up/Anticorr_extparam1.txt"),
  Anticorrelated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Anticorr_extparam1.txt"),
  Noncorrelated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Noncorr_extparam1.txt")
)

for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


################################################################################
param.dir <- file.path(script.dir, "templates/3g")
ifelse(!dir.exists(file.path("simul/3g")), dir.create("simul/3g"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path(script.dir, "simul/3g"))
src.dir   <- file.path(script.dir, "functions")

overwrite <- FALSE
# options(warn=1)
prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
random=TRUE #New environment randomly drawn or not
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/3g-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)


all.sims <- rbind(
  Correlated_Down = c("Down/param1.txt", "Down/Corr_extparam1.txt"),
  Correlated_Up = c("Up/param1.txt", "Up/Corr_extparam1.txt"),
  Correlated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Corr_extparam1.txt"),
  Anticorrelated_Down = c("Down/param1.txt", "Down/Anticorr_extparam1.txt"),
  Anticorrelated_Up = c("Up/param1.txt", "Up/Anticorr_extparam1.txt"),
  Anticorrelated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Anticorr_extparam1.txt"),
  Noncorrelated_UD = c("Up_and_Down/param1.txt", "Up_and_Down/Noncorr_extparam1.txt")
)

for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}




