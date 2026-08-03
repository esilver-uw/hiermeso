# Read in data from simulation and generate the plots for presentation

library(ggplot2)
source("sim_utils.R")

# Load in data as from sim_arrays_job.R
sim_arrays_1 <- readRDS("outputs/sim_arrays_1.RData")

selex_arrays_1 <- omnibus_test(sim_arrays_1, "res_1_group_2", 0.05, 3)
selex_arrays_2 <- omnibus_test(sim_arrays_1, "res_2_group_2", 0.05, 3)
selex_arrays_3 <- omnibus_test(sim_arrays_1, "res_3_group_2", 0.05, 3)

viz_matrix <- viz_fitter(selex_arrays_1)

# Plot n procedures
detex_plot <- function(viz_mat, methods) {
  viz_mat <- viz_mat[viz_mat$method %in% methods]
  ggplot() +
    geom_path(mapping = aes(x = viz_mat$Size, y = viz_mat$Detex, colour = viz_mat$Comp, group = viz_mat$Comp)) +
    xlab(label = "Signal Size") + ylab(label = "Proportion of Rejections") +
    scale_colour_discrete(name = "Procedure") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size)))
}

# Plot n resolutions
res_plot <- function(viz_mat) {
  ggplot() + 
    geom_col(mapping = aes(x = viz_mat$Size, y = viz_mat$Detex, fill = viz_mat$Comp), position = position_identity()) +
    xlab(label = "Signal Size") +
    ylab(label = "Avg. Rejection Resolution") +
    scale_colour_discrete(name = "Perturbation Size") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size)))
}
