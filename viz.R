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

# Over signal sizes, for fixed p_group and method, get proportion of highest-res-rejections per resolution.
# INPUT: 
# selex_arrays: a list of arrays of test group selections by iteration and size, indexed by group.
# OUTPUT: 
# res_mat: a data.frame in the correct format for ggplot to group by resolution, with size as x and detections as y.
res_fitter <- function(selex_arrays, comp_idx = 3) {
  selex_array <- selex_arrays[comp_idx,,,,drop=T]
  res_mat <- data.frame()
  for (j in 1:dim(selex_array)[2]) {
    res_mat <- rbind(res_mat, data.frame("Detex" = 0, "Size" = rep(dimnames(selex_array)[[2]][j], length(GROUP.SIZES)), "Comp" = 1:length(GROUP.SIZES)))
  }
  for (i in 1:dim(selex_array)[1]) {
    for (j in 1:dim(selex_array)[2]) {
      max_res <- max(selex_array[i,j,] * GROUPS[[4]]$Resolution)
      row <- which(res_mat$Size == dimnames(selex_array)[[2]][j] & res_mat$Comp == max_res)
      res_mat[row,1] <- res_mat[row,1] + 1
    }
  }
  res_mat$Detex <- res_mat$Detex/dim(selex_array)[1]
  return(res_mat)
}

# Fitter function for FDR table.
# INPUT: 
# selex_arrays: a list of arrays of test group selections by iteration and size, indexed by method.
# OUTPUT: 
# fdr_mat: a data.frame of detections by method. 
fdr_fitter <- function(selex_arrays, fd_arrays) {
  # Per iteration, if it makes a rejection, make it 1, if it doesn't, make it 0
  fdr_mat <- apply(fd_arrays, c(1,2), sum)/apply(selex_arrays, c(1,2), sum)
  fdr_mat[is.nan(fdr_mat)] <- 0
  fdr_vec <- apply(fdr_mat, 1, sum)/dim(selex_arrays)[2]
  return(fdr_vec)
}

# Fitter function for coverage rate plot.

# Get vector of edges in each high-resolution subgroup. Get vector of edges for all rejected hypotheses. 
# For each high resolution subgroup, 1 if it's completely contained in that vector, 0 if there's any point left out.
cvg_fitter <- function(selex_arrays, p_group) {
  res_grp <- GROUPS[[4]][GROUPS[[4]]$res_Group == p_group, 1:2]
  p_edges <- which(GROUPS[[1]][,res_grp[2]] == res_grp[1])
  cvg_mat <- data.frame(nrow = dim(selex_arrays)[2], ncol = dim(selex_arrays)[3])
  for (i in 1:dim(selex_arrays)[2]) {
    for (j in 1:dim(selex_arrays)[3]) {
      # The architecture on this one is complex. Consider a helper function.
      GROUPS[[1]][selex_arrays[p_group, i, j,] * dimnames(selex_arrays)[4]]
    }
  }
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
res_plot <- function(viz_mat, title, pos = position_stack(reverse = T)) {
  ggplot() + 
    geom_col(mapping = aes(x = viz_mat$Size, y = viz_mat$Detex, fill = as.factor(viz_mat$Comp)), position = pos) +
    xlab(label = "Signal Size") +
    ylab(label = "Avg. Rejection Resolution") +
    scale_fill_discrete(name = "Detection Size") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size))) +
    labs(title = title)
}

# Consider a plot of rejection coverage.

# Load in data as from sim_arrays_job.R

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
selex_array_1 <- readRDS("outputs/selex_array_1_vm_2.RData")
selex_array_2 <- readRDS("outputs/selex_array_2_vm_2.RData")
selex_array_3 <- readRDS("outputs/selex_array_3_vm_2.RData")

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
method_array_2 <- readRDS("outputs/method_array_2_vm_2.RData")
method_array_5 <- readRDS("outputs/method_array_5_vm_2.RData")
method_array_6 <- readRDS("outputs/method_array_6_vm_2.RData")
method_array_7 <- readRDS("outputs/method_array_7_vm_2.RData")
method_array_8 <- readRDS("outputs/method_array_8_vm_2.RData")
method_array_9 <- readRDS("outputs/method_array_9_vm_2.RData")
method_array_10 <- readRDS("outputs/method_array_10_vm_2.RData")

fdr_selex_array_1 <- readRDS("outputs/fdr_selex_array_1_vm_2.RData")
fdr_selex_array_2 <- readRDS("outputs/fdr_selex_array_2_vm_2.RData")
fdr_selex_array_3 <- readRDS("outputs/fdr_selex_array_3_vm_2.RData")

# Get true/false detections.

td_array_1 <- filter_true_selex(selex_array_1, "res_1_group_2", T)
fd_array_1 <- filter_true_selex(selex_array_1, "res_1_group_2", T, T)
td_array_2 <- filter_true_selex(selex_array_2, "res_2_group_2", T)
fd_array_2 <- filter_true_selex(selex_array_2, "res_2_group_2", T, T)
td_array_3 <- filter_true_selex(selex_array_3, "res_3_group_2", T)
fd_array_3 <- filter_true_selex(selex_array_3, "res_3_group_2", T, T)

viz_sa_1 <- viz_fitter(td_array_1)
viz_sa_2 <- viz_fitter(td_array_2)
viz_sa_3 <- viz_fitter(td_array_3)

detex_plot(viz_sa_1, "Resolution 1 Detections")
detex_plot(viz_sa_1, "Resolution 1 Detections", c("p_value", "lr_mean_3", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_2, "Resolution 2 Detections")
detex_plot(viz_sa_2, "Resolution 2 Detections", c("p_value", "lr_mean_2", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_3, "Resolution 3 Detections")
detex_plot(viz_sa_3, "Resolution 3 Detections", c("p_value", "lr_mean_1", "cal_kappa_1", "lr_prior_2", "cal_mix"))

# False detex
fdr_mat_1 <- fdr_fitter(selex_array_1, fd_array_1)
fdr_mat_2 <- fdr_fitter(selex_array_2, fd_array_2)
fdr_mat_3 <- fdr_fitter(selex_array_3, fd_array_3)

fdr_tbl <- rbind(fdr_mat_1, fdr_mat_2, fdr_mat_3)
rownames <- c("P. Res. 1", "P. Res. 2", "P. Res. 3")
fdr_tbl <- data.frame(cbind(rownames, fdr_tbl))
gt(fdr_tbl) |> fmt_number(decimals = 6) |> gt_split(col_slice_at = 6)

# Get true detections for methods.

td_mthd_2 <- filter_true_selex(method_array_2, semi_true = T)
td_mthd_5 <- filter_true_selex(method_array_5, semi_true = T)
td_mthd_6 <- filter_true_selex(method_array_6, semi_true = T)
td_mthd_7 <- filter_true_selex(method_array_7, semi_true = T)
td_mthd_8 <- filter_true_selex(method_array_8, semi_true = T)
td_mthd_9 <- filter_true_selex(method_array_9, semi_true = T)
td_mthd_10 <- filter_true_selex(method_array_10, semi_true = T)

res_mat_2_3 <- res_fitter(td_mthd_2)
res_mat_5_3 <- res_fitter(td_mthd_5)
res_mat_6_3 <- res_fitter(td_mthd_6)
res_mat_7_3 <- res_fitter(td_mthd_7)
res_mat_8_3 <- res_fitter(td_mthd_8)
res_mat_9_3 <- res_fitter(td_mthd_9)
res_mat_10_3 <- res_fitter(td_mthd_10)
res_mat_2_2 <- res_fitter(td_mthd_2, 2)
res_mat_5_2 <- res_fitter(td_mthd_5, 2)
res_mat_6_2 <- res_fitter(td_mthd_6, 2)
res_mat_7_2 <- res_fitter(td_mthd_7, 2)
res_mat_8_2 <- res_fitter(td_mthd_8, 2)
res_mat_9_2 <- res_fitter(td_mthd_9, 2)
res_mat_10_2 <- res_fitter(td_mthd_10, 2)
res_mat_2_1 <- res_fitter(td_mthd_2, 1)
res_mat_5_1 <- res_fitter(td_mthd_5, 1)
res_mat_6_1 <- res_fitter(td_mthd_6, 1)
res_mat_7_1 <- res_fitter(td_mthd_7, 1)
res_mat_8_1 <- res_fitter(td_mthd_8, 1)
res_mat_9_1 <- res_fitter(td_mthd_9, 1)
res_mat_10_1 <- res_fitter(td_mthd_10, 1)

res_plot(res_mat_2_1, "Res 2 at Res 1")
res_plot(res_mat_2_2, "Res 2 at Res 2")
res_plot(res_mat_2_3, "Res 2 at Res 3")

viz_ma_2 <- viz_fitter(method_array_2)
viz_ma_5 <- viz_fitter(method_array_5)
viz_ma_6 <- viz_fitter(method_array_6)
viz_ma_7 <- viz_fitter(method_array_7)
viz_ma_8 <- viz_fitter(method_array_8)
viz_ma_9 <- viz_fitter(method_array_9)
viz_ma_10 <- viz_fitter(method_array_10)

res_plot(viz_ma_2, "Calibrator, k = 0.25", position_identity())
res_plot(viz_ma_5, "Calibrator, mixture", position_identity())
res_plot(viz_ma_6, "Likelihood Ratio, mean = 2.5", position_identity())
res_plot(viz_ma_7, "Likelihood Ratio, mean = 5", position_identity())
res_plot(viz_ma_8, "Likelihood Ratio, mean = 7.5", position_identity())
res_plot(viz_ma_9, "LR Mixture, prior sigma = 5", position_identity())
res_plot(viz_ma_10, "LR Mixture, prior sigma = 20", position_identity())
