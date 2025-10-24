#!/usr/bin/env Rscript
source("scripts/functions/functions.R")
source("scripts/functions/makeparam_funct.R")
################################################################################
param.dir <- file.path("templates/Full_netw_TF_more")
directory <- "simul/Full_netw_TF_more"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
# prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_TF_more-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_TF_less")
directory <- "simul/Full_netw_TF_less"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_TF_less-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_pop_2000")
directory <- "simul/Full_netw_pop_2000"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_pop_2000-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_pop_200")
directory <- "simul/Full_netw_pop_200"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
# prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_pop_200-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}
################################################################################
param.dir <- file.path("templates/Full_netw_genes_more")
directory <- "simul/Full_netw_genes_more"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
# prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_genes_more-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_genes_less")
directory <- "simul/Full_netw_genes_less"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_genes_less-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_a_more")
directory <- "simul/Full_netw_a_more"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
# prog.path <- "/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_a_more-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

################################################################################
param.dir <- file.path("templates/Full_netw_a_less")
directory <- "simul/Full_netw_a_less"
ifelse(!dir.exists(directory), dir.create(directory), FALSE)
cache.dir <- normalizePath(file.path(directory))

overwrite <- FALSE
# options(warn=1)
prog.path <- "../../simevolv/bin/Release/Simul_Prog"
# prog.path <- "../../simevolv/bin/Release/Simul_Prog"
#####
sd <- 0.25 #0.25
min <- 0.15#min environment
max <- 0.85 #max environment
#
launchfilename <- "launchers/Full_netw_a_less-launch.sh"
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
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE, sampling="")
}

