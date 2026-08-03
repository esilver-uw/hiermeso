# Read in data from simulation and generate the plots for presentation

library(ggplot2)
source("sim_utils.R")

# Load in data as from sim_arrays_job.R
sim_arrays_1 <- readRDS("outputs/sim_arrays_1.RData")

selex_arrays_1 <- omnibus_test(sim_arrays_1, 1, 0.05, 3)
selex_arrays_2 <- omnibus_test(sim_arrays_1, 2, 0.05, 3)
selex_arrays_3 <- omnibus_test(sim_arrays_1, 3, 0.05, 3)

viz_matrix <- viz_fitter(selex_arrays_1)

# Plot n procedures
detex_plot <- function(viz_mat) {
  ggplot() +
    geom_path(mapping = aes(x = as.integer(viz_mat$size), y = viz_mat$detex, colour = viz_mat$method, group = viz_mat$method)) +
    xlab(label = "Signal Size") + ylab(label = "Proportion of Rejections") +
    scale_colour_discrete(name = "Procedure")
}

# Plot n resolutions
res_plot <- function(sizes, means) {
  ggplot() + 
    geom_path(mapping = aes(x = rep(sizes, 3), y = means[,1], colour = as.factor(means$res1))) +
    xlab(label = "Signal Size") +
    ylab(label = "Avg. Rejection Resolution") +
    scale_colour_discrete(name = "Perturbation Size")
}