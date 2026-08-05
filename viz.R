# Read in data from simulation and generate the plots for presentation

library(ggplot2)
library(gt)
source("sim_utils.R")

# THE LENSES: Two extremely basic apply wrappers.
# INPUT: 
# selex_array: an array of test group selections by iteration and size.
# OUTPUT:
# detex_array: whether or not any detection is made at a given method/iteration/size.
detex_selex <- function(selex_array) {
  norm <- dim(selex_array)[1]
  if (length(dim(selex_array)) > 2) {
    # Start by collapsing the test group dimension into 1 on whether or not a detection is made.
    selex_array <- apply(selex_array, c(1,2), sum)
    selex_array[selex_array >= 1] <- 1
    detex_array <- apply(selex_array, 2, sum)/norm
  } else {
    # Set margins based on signal sizes
    selex_array <- apply(selex_array, 1, sum)
    selex_array[selex_array >= 1] <- 1
    detex_array <- selex_array/norm
  }
  
  return(detex_array)
}

# OUTPUT:
# mean_array: mean detection resolution. NOTE: uses inverted resolutions because I was not at the time of coding this library thinking about resolutions in the most 
# intuitive way, 
mean_selex <- function(selex_array) {
  selex_array_wgt <- selex_array[]
  if (length(dim(selex_array)) > 2) {
    for (grp in dimnames(selex_array)[[3]]) {
      grp_num <- which(GROUPS[[4]]$res_Group == grp)
      wgt <- 4 - GROUPS[[4]]$Resolution[grp_num]
      selex_array_wgt[,,grp_num] <- selex_array[,,grp_num] * wgt
    }
    
    # Weight detections by resolution. Note the inversion.
    selex_array_wgt <- apply(selex_array_wgt, c(1,2), sum)
    detex_array <- apply(selex_array, c(1,2), sum)
    mean_array <- selex_array_wgt/detex_array
    mean_array[is.nan(mean_array)] <- 0
    
    mean_array <- apply(mean_array, 2, mean)
  } else {
    for (grp in dimnames(selex_array)[[2]]) {
      grp_num <- which(GROUPS[[4]]$res_Group == grp)
      wgt <- 4 - GROUPS[[4]]$Resolution[grp_num]
      selex_array_wgt[,grp_num] <- selex_array[,grp_num] * wgt
    }
    # Weight detections by resolution. Note the inversion.
    selex_array_wgt <- apply(selex_array_wgt, 1, sum)
    detex_array <- apply(selex_array, 1, sum)
    mean_array <- selex_array_wgt/detex_array
    mean_array[is.nan(mean_array)] <- 0
    
    mean_array <- mean(mean_array)
  }
  
  return(mean_array)
}

# Fitter function for visualization.
# INPUT: 
# selex_arrays: a list of arrays of test group selections by iteration and size, indexed by method or group.
# lens: sum_selex, mean_selex, or detex_selex. Default: detex_selex.
# OUTPUT: 
# viz_mat: a data.frame in the correct format for ggplot to group by method or group, with size as x and detections as y.
viz_fitter <- function(selex_arrays, lens = detex_selex) {
  viz_mat <- data.frame()
  for (i in dim(selex_arrays)[1]:1) {
    selex_mat <- lens(selex_arrays[i,,,])
    selex_mat <- cbind(selex_mat, dimnames(selex_arrays)[[3]], rep(dimnames(selex_arrays)[[1]][i], length(selex_mat)))
    viz_mat <- rbind(viz_mat, selex_mat)
  }
  colnames(viz_mat) <- c("Detex", "Size", "Comp")
  rownames(viz_mat) <- 1:dim(viz_mat)[1]
  viz_mat$Detex <- as.numeric(viz_mat$Detex)
  viz_mat$Size <- as.numeric(viz_mat$Size)
  return(viz_mat)
}

# Fitter function for FDR table.
# INPUT: 
# selex_arrays: a list of arrays of test group selections by iteration and size, indexed by method.
# OUTPUT: 
# fdr_mat: a data.frame of detections by method. 
fdr_fitter <- function(selex_arrays) {
  fdr_mat <- data.frame(t(apply(selex_arrays, 1, mean)))
  return(fdr_mat)
}

# THE PLOTS
# detex_plot plots curves of detections, res_plot plots bars of detections over resolutions.
# INPUT: 
# viz_mat: visualization matrix obtained by passing an omnibus or omnires test through viz_fitter.
# methods (detex plot): which methods to consider
# OUTPUT:
# none- plots a graph.
# Plot n procedures
detex_plot <- function(viz_mat, title, methods = c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_mean_1", "lr_mean_2", "lr_mean_3", "lr_prior_1", "lr_prior_2", "lr_prior_3")) {
  viz_mat <- viz_mat[viz_mat$Comp %in% methods,]
  ggplot() +
    geom_path(mapping = aes(x = viz_mat$Size, y = viz_mat$Detex, colour = viz_mat$Comp, group = viz_mat$Comp)) +
    xlab(label = "Signal Size") + ylab(label = "Proportion of Rejections") +
    scale_colour_discrete(name = "Procedure") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size))) +
    labs(title = title)
}

# Plot n resolutions
res_plot <- function(viz_mat, title) {
  ggplot() + 
    geom_col(mapping = aes(x = viz_mat$Size, y = viz_mat$Detex, fill = viz_mat$Comp), position = position_identity()) +
    xlab(label = "Signal Size") +
    ylab(label = "Avg. Rejection Resolution") +
    scale_fill_discrete(name = "Perturbation Size") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size))) +
    labs(title = title)
}

# Load in data as from sim_arrays_job.R

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
selex_array_1 <- readRDS("outputs/selex_array_1_vm.RData")
selex_array_2 <- readRDS("outputs/selex_array_2_vm.RData")
selex_array_3 <- readRDS("outputs/selex_array_3_vm.RData")

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
method_array_2 <- readRDS("outputs/method_array_2_vm.RData")
method_array_5 <- readRDS("outputs/method_array_5.RData")
method_array_6 <- readRDS("outputs/method_array_6_vm.RData")
method_array_7 <- readRDS("outputs/method_array_7.RData")
method_array_8 <- readRDS("outputs/method_array_8.RData")
method_array_9 <- readRDS("outputs/method_array_9_vm.RData")
method_array_10 <- readRDS("outputs/method_array_10.RData")

# Plots

viz_sa_1 <- viz_fitter(selex_array_1)
viz_sa_2 <- viz_fitter(selex_array_2)
viz_sa_3 <- viz_fitter(selex_array_3)

detex_plot(viz_sa_1, "Resolution 1 Detections")
detex_plot(viz_sa_1, "Resolution 1 Detections", c("p_value", "lr_mean_3", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_2, "Resolution 2 Detections")
detex_plot(viz_sa_2, "Resolution 2 Detections", c("p_value", "lr_mean_2", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_3, "Resolution 3 Detections")
detex_plot(viz_sa_3, "Resolution 3 Detections", c("p_value", "lr_mean_1", "cal_kappa_1", "lr_prior_2", "cal_mix"))

viz_ma_2 <- viz_fitter(method_array_2)
viz_ma_5 <- viz_fitter(method_array_5)
viz_ma_6 <- viz_fitter(method_array_6)
viz_ma_7 <- viz_fitter(method_array_7)
viz_ma_8 <- viz_fitter(method_array_8)
viz_ma_9 <- viz_fitter(method_array_9)
viz_ma_10 <- viz_fitter(method_array_10)

res_plot(viz_ma_2, "Calibrator, k = 0.25")
res_plot(viz_ma_5, "Calibrator, mixture")
res_plot(viz_ma_6, "Likelihood Ratio, mean = 2.5")
res_plot(viz_ma_7, "Likelihood Ratio, mean = 5")
res_plot(viz_ma_8, "Likelihood Ratio, mean = 7.5")
res_plot(viz_ma_9, "LR Mixture, prior sigma = 5")
res_plot(viz_ma_10, "LR Mixture, prior sigma = 20")

# False detex
fdr_mat_1 <- fdr_fitter(fdr_selex_array_1)
fdr_mat_2 <- fdr_fitter(fdr_selex_array_2)
fdr_mat_3 <- fdr_fitter(fdr_selex_array_3)

fdr_tbl <- rbind(fdr_mat_1, fdr_mat_2, fdr_mat_3)
rownames <- c("P. Res. 1", "P. Res. 2", "P. Res. 3")
fdr_tbl <- data.frame(cbind(rownames, fdr_tbl))
gt(fdr_tbl) |> fmt_number(decimals = 6) |> gt_split(col_slice_at = 6)
