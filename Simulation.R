# Packages and functions
library(snowfall)
source("Function.R")

# Simulation
sfInit(parallel = T, cpus = 40)
sfExport("n", "rho", "time_num", "p_X", "p_Y", "comp", "p_ture")
sfSource("Function.R")
sfSource("LCCA_function.R")
Result <- sfLapply(seq(1, 10000, length.out = 100), simu_func,
                   n = n, rho = rho, time_num = time_num,
                   p_X = p_X, p_Y = p_Y, comp = comp, p_ture = p_ture)


sfStop()

save(Result, file = paste0("Result/Result_", n, "_", p_X, ".rda"), version = 2)

