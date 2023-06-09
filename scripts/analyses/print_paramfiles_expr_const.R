#!/usr/bin/env Rscript
source("scripts/functions/functions.R")
source("scripts/functions/makeparam_funct.R")
################################################################################
param.dir <- file.path("templates/10g_a0.2")
ifelse(!dir.exists(file.path("simul/10g_a0.2")), dir.create("simul/10g_a0.2"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path("simul/10g_a0.2"))

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
launchfilename <- "launchers/10g_a0-2-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)


#Correlated
all.sims <- rbind(
  Correlated = c("Correlated/param1.txt", "Correlated/extparam1.txt")
)
#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,0.5,1)
  # RNintercept <- runif(1, 0, 1-RNslope )
  RNslope <- runif(1,0.5718,1.4282)
  if(RNslope > 1) RNintercept <- runif(1, 1-RNslope, 0) else RNintercept <- runif(1, 0, 1-RNslope )
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Correlated
all.sims <- rbind(
  Anticorrelated = c("Anticorrelated/param1.txt", "Anticorrelated/extparam1.txt")
)
#Negative RN = Anticorrelated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,-1,-0.5)
  # RNintercept <- runif(1, -RNslope, 1) 
  RNslope <- runif(1,-1.4282,-0.5718)
  if(RNslope < -1) RNintercept <- runif(1, 1, -RNslope) else RNintercept <- runif(1, -RNslope, 1)
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Controls
all.sims <- rbind(
  Control_sel = c("Control_sel/param1.txt", "Control_sel/extparam1.txt")
)
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)}

################################################################################
param.dir <- file.path("templates/10g_a0.3")
ifelse(!dir.exists(file.path("simul/10g_a0.3")), dir.create("simul/10g_a0.3"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path("simul/10g_a0.3"))

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
launchfilename <- "launchers/10g_a0-3-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)


#Correlated
all.sims <- rbind(
  Correlated = c("Correlated/param1.txt", "Correlated/extparam1.txt")
)
#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,0.5,1)
  # RNintercept <- runif(1, 0, 1-RNslope )
  RNslope <- runif(1,0.5718,1.4282)
  if(RNslope > 1) RNintercept <- runif(1, 1-RNslope, 0) else RNintercept <- runif(1, 0, 1-RNslope )
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Correlated
all.sims <- rbind(
  Anticorrelated = c("Anticorrelated/param1.txt", "Anticorrelated/extparam1.txt")
)
#Negative RN = Anticorrelated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,-1,-0.5)
  # RNintercept <- runif(1, -RNslope, 1) 
  RNslope <- runif(1,-1.4282,-0.5718)
  if(RNslope < -1) RNintercept <- runif(1, 1, -RNslope) else RNintercept <- runif(1, -RNslope, 1)
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Controls
all.sims <- rbind(
  Control_sel = c("Control_sel/param1.txt", "Control_sel/extparam1.txt")
)
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)}


################################################################################
param.dir <- file.path("templates/10g_a0.4")
ifelse(!dir.exists(file.path("simul/10g_a0.4")), dir.create("simul/10g_a0.4"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path("simul/10g_a0.4"))

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
launchfilename <- "launchers/10g_a0-4-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)


#Correlated
all.sims <- rbind(
  Correlated = c("Correlated/param1.txt", "Correlated/extparam1.txt")
)
#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,0.5,1)
  # RNintercept <- runif(1, 0, 1-RNslope )
  RNslope <- runif(1,0.5718,1.4282)
  if(RNslope > 1) RNintercept <- runif(1, 1-RNslope, 0) else RNintercept <- runif(1, 0, 1-RNslope )
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Correlated
all.sims <- rbind(
  Anticorrelated = c("Anticorrelated/param1.txt", "Anticorrelated/extparam1.txt")
)
#Negative RN = Anticorrelated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,-1,-0.5)
  # RNintercept <- runif(1, -RNslope, 1) 
  RNslope <- runif(1,-1.4282,-0.5718)
  if(RNslope < -1) RNintercept <- runif(1, 1, -RNslope) else RNintercept <- runif(1, -RNslope, 1)
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE, RNslope=RNslope, RNintercept=RNintercept)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}


#Controls
all.sims <- rbind(
  Control_sel = c("Control_sel/param1.txt", "Control_sel/extparam1.txt")
)
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)}

