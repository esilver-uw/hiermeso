# Library for visualization functions

# Imports
library(ggplot2)
library(gt)

# Filters to include only true selections.
# INPUT:
# selex_arrays: an array of test group selections by method or group, iteration and size.
# p_group: the true perturbation group.
# semi-true: whether to also include 'semi-true' selections (selections which contain or are contained by the true perturbation group).
# false: reverse polarity-- filter out true selections.
# OUTPUT:
# filtered_selex_array: a filtered array of test group selections by iteration and size.
filter_true_selex <- function(selex_arrays, p_group = -1, semi_true = F, false = F) {
  if (p_group != -1) {
    accept <- GROUPS[[2]][p_group][[1]]
    if (semi_true) {
      accept <- union(accept, GROUPS[[4]]$res_Group[sapply(GROUPS[[2]], function(x) p_group %in% x)])
    }
    
    true_idxs <- which(GROUPS[[4]]$res_Group %in% accept)
    
    for (i in 1:dim(selex_arrays)[1]) {
      if (false) {
        selex_arrays[i,,,true_idxs] <- 0
      } else {
        selex_arrays[i,,,-true_idxs] <- 0
      }
    }
  } else {
    for (i in 1:dim(selex_arrays)[1]) {
      p_group <- dimnames(selex_arrays)[[1]][i]
      
      accept <- GROUPS[[2]][p_group][[1]]
      if (semi_true) {
        accept <- union(accept, GROUPS[[4]]$res_Group[sapply(GROUPS[[2]], function(x) p_group %in% x)])
      }
      
      true_idxs <- which(GROUPS[[4]]$res_Group %in% accept)
      
      if (false) {
        selex_arrays[i,,,true_idxs] <- 0
      } else {
        selex_arrays[i,,,-true_idxs] <- 0
      }
    }
  }
  
  return(selex_arrays)
}

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
  fdr_vec <- apply(fdr_mat, 1, sum)/200
  fdr_vec[is.nan(fdr_vec)] <- 0
  fdr_mat[is.nan(fdr_mat)] <- 0
  fdr_vec <- apply(fdr_mat, 1, sum)/dim(selex_arrays)[2]
  return(fdr_vec)
}

# Fitter function for coverage rate plot. 
# INPUT: 
# l_group: the group that may be inside the others. Intended to be highest-resolution.
# r_groups: vector of groups that may contain the others.
# OUTPUT: 
# contains: 1 if contained or 0
group_containment <- function(l_group, r_groups) {
  contains <- 0
  
  for (group in r_groups) {
    if (l_group %in% GROUPS[[2]][[group]]) {
      contains <- 1
    }
  }
  
  return(contains)
}

# Get vector of edges in each high-resolution subgroup. Get vector of edges for all rejected hypotheses. 
# For each high resolution subgroup, 1 if it's completely contained in that vector, 0 if there's any point left out.
cvg_fitter <- function(selex_arrays, p_group) {
  # Get resolution-1 subgroups of p_group
  p_subgroups <- GROUPS[[2]][[p_group]]
  res_1_groups <- GROUPS[[4]][GROUPS[[4]]$Resolution == 1,4]
  p_res_1_subgroups <- intersect(p_subgroups, res_1_groups)
  
  # Method and signal size.
  cvg_mat <- matrix(nrow = dim(selex_arrays)[1], ncol = dim(selex_arrays)[3])
  dimnames(cvg_mat) <- list(dimnames(selex_arrays)[[1]], dimnames(selex_arrays)[[3]])
  
  for (i in 1:dim(selex_arrays)[1]) {
    for (k in 1:dim(selex_arrays)[3]) {
      # We produce the running count over these indices
      running_rate <- NA
      for (j in 1:dim(selex_arrays)[2]) {
        ct <- 0
        selex_vec <- selex_arrays[i,j,k,]
        detex_groups <- names(selex_vec[selex_vec == 1])
        
        for (g in p_res_1_subgroups) {
          ct <- ct + group_containment(g, detex_groups)
        }
        ct <- ct/length(p_res_1_subgroups)
        if (ct > 0) {
          running_rate <- mean(c(running_rate, ct), na.rm = T)
        }
      }
      if (is.na(running_rate)) {
        running_rate <- 0
      }
      cvg_mat[i,k] <- running_rate
    }
  }
  
  cvg_fit <- data.frame()
  
  for (i in 1:dim(cvg_mat)[2]) {
    cvg_fit <- rbind(cvg_fit, cbind(as.numeric(cvg_mat[,i]), rep(dimnames(cvg_mat)[[2]][i], dim(cvg_mat)[1]), dimnames(cvg_mat)[[1]]))
  }
  
  colnames(cvg_fit) <- c("Cvg", "Size", "Comp")
  
  return(cvg_fit)
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

# Plot n procedures
cvg_plot <- function(viz_mat, title, methods = c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_mean_1", "lr_mean_2", "lr_mean_3", "lr_prior_1", "lr_prior_2", "lr_prior_3")) {
  viz_mat <- viz_mat[viz_mat$Comp %in% methods,]
  ggplot() +
    geom_path(mapping = aes(x = as.factor(viz_mat$Size), y = as.numeric(viz_mat$Cvg), colour = viz_mat$Comp, group = viz_mat$Comp)) +
    xlab(label = "Signal Size") + ylab(label = "Coverage Proportion") +
    scale_colour_discrete(name = "Procedure") +
    scale_x_discrete(limits = as.factor(unique(viz_mat$Size))) +
    scale_y_continuous(limits = c(0,1)) + 
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