#Print_paramfile for Control
#!/usr/bin/env Rscript
source("scripts/functions/functions.R")
source("scripts/functions/makeparam_funct.R")
################################################################################
param.dir <- file.path("templates/3g_control")
ifelse(!dir.exists(file.path("simul/3g_control")), dir.create("simul/3g_control"), FALSE)#burning_andreas folder
cache.dir <- normalizePath(file.path("simul/3g_control"))

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
launchfilename <- "launchers/3g_control-launch.sh"
##########################
# Generate all simulations
# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

#Delete previous launch file if exist
ifelse(file.exists(launchfilename), unlink(launchfilename), FALSE)

all.sims <- rbind(
  No_sel = c("param1.txt", "extparam_no_sel.txt"),
  Sel = c("param1.txt", "extparam_sel.txt")
)

#Positive RN = Correlated
for (sim.name in rownames(all.sims)) {
  cat("Setting up simulation", sim.name, "...\n")
  # RNslope <- runif(1,0.5,1)
  # RNintercept <- runif(1, 0, 1-RNslope )
  pars <- create.paramseries(
    file.path(param.dir, all.sims[sim.name, 1]), 
    file.path(param.dir, all.sims[sim.name, 2]), 
    file.path(cache.dir, sim.name), 
    overwrite=overwrite, verbose=TRUE, sd=sd, min=min, max=max, bottleneck=FALSE)
  
  create.launchfile.alt(prog.path, pars$param, pars$out, pars$compressed, relative.paths=TRUE,
                        oldpop="none", file.path(launchfilename), prevpop=FALSE)
}

