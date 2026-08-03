# Read in data from simulation and generate the plots for presentation

library(ggplot2)
source("sim_utils.R")

# Load in data as from sim_arrays_job.R
sim_arrays_1 <- readRDS("outputs/sim_arrays_1.RData")

selex_arrays_1 <- omnibus_test(sim_arrays_1, 1, 0.05, 3)
selex_arrays_2 <- omnibus_test(sim_arrays_1, 2, 0.05, 3)
selex_arrays_3 <- omnibus_test(sim_arrays_1, 3, 0.05, 3)

viz_matrix <- viz_fitter(selex_arrays_1)
