# Read in data from simulation and generate the plots for presentation

library(ggplot2)

# TODO: Replace with local
res_mean <- read.csv2("projects/def-pwmacdon/e4silver/hiermeso/resolution_mean.csv")[,2:4]
detex_mean <- read.csv2("projects/def-pwmacdon/e4silver/hiermeso/detections_mean.csv")[,2:4]
comp_mean <- read.csv2("projects/def-pwmacdon/e4silver/hiermeso/detections_comparator_mean.csv")[,2:4]
benj_hoch_mean <- read.csv2("projects/def-pwmacdon/e4silver/hiermeso/detections_p_bh_comp_mean.csv")[,2:4]

# We want the following plots: for each resolution, across signal size, three line plots using the detections. 
# Then, a resolution plot as a power curve. 
# TODO: this will require recomputation.

detex_mean$Type <- rep("Hierarchical eLP", length(SIZES))
comp_mean$Type <- rep("Flat eLP", length(SIZES))
benj_hoch_mean$Type <- rep("BH Comparator", length(SIZES))
detex_means <- rbind(detex_mean, comp_mean, benj_hoch_mean)

# Plot three procedures
detex_plot <- function(sizes, means, resolution) {
  ggplot() +
    geom_path(mapping = aes(x = rep(sizes, 3), y = means[,resolution], colour = means$Type)) +
    xlab(label = "Signal Size") + ylab(label = "Proportion of Rejections") +
    scale_colour_discrete(name = "Procedure")
}

res_mean$res1 <- rep(1, length(SIZES))
res_mean$res2 <- rep(2, length(SIZES))
res_mean$res3 <- rep(3, length(SIZES))
res_means <- data.frame(rbind(as.matrix(res_mean[,c(1,4)]), as.matrix(res_mean[,c(2,5)]), as.matrix(res_mean[,c(3,6)])))

# Plot three resolutions
res_plot <- function(sizes, means) {
  ggplot() + 
    geom_path(mapping = aes(x = rep(sizes, 3), y = means[,1], colour = as.factor(means$res1))) +
    xlab(label = "Signal Size") +
    ylab(label = "Avg. Rejection Resolution") +
    scale_colour_discrete(name = "Perturbation Size")
}

detex_plot(SIZES, detex_means, 1)
detex_plot(SIZES, detex_means, 2)
detex_plot(SIZES, detex_means, 3)
res_plot(SIZES, res_means)
