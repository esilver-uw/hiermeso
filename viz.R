# Read in data from simulation and generate the plots for presentation

library(ggplot2)
source("sim_utils.R")

# Load in data as from sim_arrays_job.R
sim_arrays_1 <- readRDS("outputs/sim_arrays_1.RData")

# Globals
TEST_GROUP <- "res_2_group_2"

selex_arrays <- omnibus_test(sim_arrays_1, TEST_GROUP, 0.05, 3)

viz_matrix <- viz_fitter(selex_arrays)