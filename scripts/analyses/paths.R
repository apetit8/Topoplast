path <- list(
    prog      = here::here("simevolv/bin/Release/Simul_Prog"),
    param.dir = list(
        fullnet = here::here("templates/Full_netw"),
        drift   = here::here("templates/Full_netw_drift")),
    sim.dir   = list(
        fullnet = here::here("simul/Full_netw"), 
        drift   = here::here("simul/Full_netw_drift")),
    launch.file= list(
        fullnet = here::here("launchers/Full_netw-launch.sh"),
        drift   = here::here("launchers/Full_drift-launch.sh"))
)

