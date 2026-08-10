# Library for simulation functions

# Imports
# library(network)
# library(igraph)
library(abind)
library(CVXR)

# Helper from Stack Overflow:
# Source - https://stackoverflow.com/a/8189441
# Posted by Ken Williams, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-18, License - CC BY-SA 4.0
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Create groups for upper triangle, inspired by KeLP architecture.
# INPUT:
# group_sizes: vector of group sizes.
# n: number of nodes.
# OUTPUT:
# for L length of GROUP_SIZES
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# group_subgroups: list of subgroups by group. (res_group notation)
# group_memberships: list of edges by group. (res_group notation)
# group_info: data frame containing group, level, and group-level for each group. (res_group notation)
generate_groups <- function(group_sizes, n) {
  # number of levels
  L <- length(group_sizes)
  
  # nodes_edge
  nodes_edge <- NULL
  for (i in 1:n) {
    for (j in i:n) {
      nodes_edge <- rbind(nodes_edge, c(i,j))
    }
  }
  
  # groups will be all combinations of node_groups. 
  upper_mask <- matrix(1:n^2, nrow = n)[upper.tri(matrix(1:n^2, nrow = n), T)]
  groups <- cbind(matrix(0, n^2, length(group_sizes)), expand.grid(1:n,1:n))
  
  # group memberships will be lists of nodes for each group and resolution, with same name as in group[[4]]
  group_memberships <- list()
  for (l in 1:L) {
    g <- 0
    group_size <- group_sizes[l]
    group_nodes <- split(1:n, ceiling(1:n/group_size))
    for (i in 1:ceiling(n/group_size)) {
      for (j in 1:ceiling(n/group_size)) {
        pairs <- as.matrix(expand.grid(group_nodes[[i]], group_nodes[[j]]))
        group <- intersect(matrix(1:n^2, nrow = n)[pairs], upper_mask)
        if (length(group)) {
          g <- g + 1
          name <- paste0("res_",l,"_group_",g)
          group_memberships[[name]] <- group
          
          for (k in group) {
            groups[k,l] <- g
          }
        }
      }
    }
  }
  
  # Remove empty rows
  groups <- groups[which(groups[,1] != 0),]
  
  group_info <- c()
  
  for (l in 1:L) {
    df_temp <- data.frame("Group_Number" = unique(groups[,l]))
    df_temp$Resolution <- l
    df_temp$group <- paste0("group_", df_temp$Group_Number)
    df_temp$res_Group <- paste0("res_",df_temp$Resolution, "_", df_temp$group)
    
    group_info <- rbind(group_info, df_temp)
  }
  
  group_subgroups <- NULL
  
  for (g in group_info$res_Group) {
    subgroups <- NULL
    group <- group_info$Group_Number[group_info$res_Group == g]
    res <- group_info$Resolution[group_info$res_Group == g]
    for (i in 1:res) {
      subgroups_i <- unique(paste0("group_",groups[groups[,res] == group, i]))
      subgroups <- append(subgroups, paste0("res_",i,"_",subgroups_i))
    }
    group_subgroups[[g]] <- subgroups
  }
  
  return(list(groups, group_subgroups, group_memberships, group_info))
}

# Returns observations from a given expected adjacency matrix.
# INPUT:
# theta: expected (mean) adjacency matrix
# n: sample size
# OUTPUT: 
# A: array of n observed adjacency matrices
sample_network <- function(theta, n) {
  A <- array(NA, dim = c(dim(theta)[1], dim(theta)[2], n))
  for (i in 1:n) {
    A[,,i] <-  matrix(rnorm(length(theta), theta, SIGMA), nrow = dim(theta)[1], ncol = dim(theta)[2])
  }
  return(A)
}

# TODO: Be able to perturb in a more targeted manner.
# Apply a perturbation to the parameter adjacency matrix.
# INPUT:
# THETA: parameter adjacency matrix.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# size: size and direction of perturbation to apply.
# OUTPUT:
# theta_prime: perturbed parameter adjacency matrix.
perturb_expected_matrix <- function(theta, g, size) {
  edges <- GROUPS[[3]][[g]]
  
  # Apply perturbation
  theta_prime <- theta[]
  theta_prime[edges] <- theta_prime[edges] + size
  return(theta_prime)
}

# Create Location Constraint Matrix.
# INPUT:
# groups: list of group attributes.
# n_base_level: number of base resolution groups (hypotheses)
# n_groups: total number of groups
# OUTPUT:
# location_constraint_matrix: a matrix with entries i,j = 1 if a base resolution group j is a subgroup of group i, else 0.
create_lcm <- function(groups, n_base_level, n_groups) {
  location_constraint_matrix <- matrix(0, n_groups, n_base_level)
  for (i in 1:dim(location_constraint_matrix)[1]) {
    # Get rows of the group at j
    res_Group <- groups[[4]][i,4]
    indices <- which(groups[[4]][groups[[4]]$Resolution == 1,4] %in% groups[[2]][[res_Group]])
    location_constraint_matrix[i,indices] <- 1
  }
  return(location_constraint_matrix)
}

# TODO: Implement weighting per Gablenz & Sabatti.
# Run eLP: Largely adapted from Gablenz & Sabatti.
# INPUT:
# e_vals: vector of e_values by aligned with groups
# t_groups: vector of group names used.
# alpha: alpha level of the test.
# OUTPUT: 
# detections: hypotheses rejected by the algorithm.
elp <- function(e_vals, t_groups, alpha) {
  # Get number of base level hypotheses and number of total hypotheses
  n_base_level <- length(intersect(GROUPS[[4]][GROUPS[[4]]$Resolution == 1,4], t_groups))
  n_groups <- length(t_groups)
  
  x <- CVXR::Variable(n_groups, integer = TRUE)
  objective <- CVXR::Maximize(sum(x))
  
  location_constraint_matrix <- create_lcm(GROUPS, n_base_level, n_groups)
  
  b <- rep(1, n_base_level)
  constraints <- list(x >= 0,
                      x <= 1,
                      t(location_constraint_matrix) %*% x <= b)
  
  constraints <- c(constraints, list(n_groups - e_vals * alpha * sum(x) <= n_groups * (1 - x)))
  
  problem <- CVXR::Problem(objective = objective, constraints = constraints)
  
  result <- psolve(problem)
  selections <- value(x)
  
  return(selections)
}

# Truncation function per Wang & Ramdas 2022
# INPUT: 
# K: K for the truncation function (number of hypotheses)
# x: x for the truncation.
trunc <- function(K, x) {
  if (x >= 1) {
    return(K/ceiling(K/x))
  }
  return(0)
}

# Create 4D simulation array, indexed by perturbation group, tested group, signal size, and method.
# INPUT: 
# p_groups: vector of groups to perturb. Default: all.
# OUTPUT: 
# sim_array: a single iteration simulation array.
array_step <- function(p_groups = GROUPS[[4]]$res_Group, t_groups = GROUPS[[4]]$res_Group) {
  sim_array <- array(dim = c(length(p_groups), length(t_groups), length(SIGNAL.SIZES), 11))
  dimnames(sim_array) <- list(p_groups, t_groups, SIGNAL.SIZES,
                              c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_mean_1", "lr_mean_2", "lr_mean_3", "lr_prior_1", "lr_prior_2", "lr_prior_3")) 
  A1 <- sample_network(THETA, N)
  # A1 Sample-wise Mean
  A1_sm <- apply(A1, c(1,2), mean)
  
  i = 0
  for (size in SIGNAL.SIZES) {
    i <- i + 1
    j = 0
    for (p_group in p_groups) {
      j <- j + 1
      k = 0
      
      theta_prime <- perturb_expected_matrix(THETA, p_group, size)
      A2 <- sample_network(theta_prime, N)
      # A2 Sample-wise Mean
      A2_sm <- apply(A2, c(1,2), mean)
      
      for (t_group in t_groups) {
        k <- k + 1
        
        # Get vector of edges in the group then yield node pairs.
        edges <- GROUPS[[3]][[t_group]]
        
        # This can happen once per group.
        m <- length(edges) * N
        
        A2_bar <- mean(A1_sm[edges])
        A1_bar <- mean(A2_sm[edges])
        
        d_bar <- A1_bar - A2_bar
        
        p_val <- p_value(d_bar, m)
        
        kappas <- c(0.25, 0.5, 0.75)
        cal_kappa_1 <- cal_kappa(p_val, kappas[1])
        cal_kappa_2 <- cal_kappa(p_val, kappas[2])
        cal_kappa_3 <- cal_kappa(p_val, kappas[3])
        
        cal_mix <- cal_mixture(p_val)
        
        pts <- c(2.5, 5, 7.5)
        lr_mean_1 <- lr_delta(d_bar, m, pts[1])
        lr_mean_2 <- lr_delta(d_bar, m, pts[2])
        lr_mean_3 <- lr_delta(d_bar, m, pts[3])
        
        priors <- c(5, 20, 35)
        lr_prior_1 <- lr_prior(d_bar, m, priors[1])
        lr_prior_2 <- lr_prior(d_bar, m, priors[2])
        lr_prior_3 <- lr_prior(d_bar, m, priors[3])
        
        sim_array[j,k,i,] <- c(p_val, cal_kappa_1, cal_kappa_2, cal_kappa_3, cal_mix, lr_mean_1, lr_mean_2, lr_mean_3, lr_prior_1, lr_prior_2, lr_prior_3)
      }
    }
  }
  return(sim_array)
}

# We want to parallellize the construction of a list of sim_arrays. Depending on ct. of cores c, divide total iteration ct. by c, assign it/c to each core, and join list.
# Creates multiple simulation arrays for parallellization purposes.
# INPUT: 
# ct: number of arrays to create.
# p_groups: vector of groups to perturb. Default: all.
# OUTPUT:
# sim_arrays: 5D array with first index iteration, then perturbation group, tested group, signal size, and method.
batchable_array <- function(ct, p_groups = GROUPS[[4]]$res_Group) {
  sim_arrays <- array(dim = c(ct, length(p_groups), dim(GROUPS[[4]])[1], length(SIGNAL.SIZES), 11))
  for (i in 1:ct) {
    sim_arrays[i,,,,] <- array_step(p_groups) 
  }
  
  dimnames(sim_arrays) <- list(1:ct, p_groups, GROUPS[[4]]$res_Group, SIGNAL.SIZES,
                               c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_mean_1", "lr_mean_2", "lr_mean_3", "lr_prior_1", "lr_prior_2", "lr_prior_3"))
  return(sim_arrays)
}

# Performs e-value testing on multiple iterations on a simulation array at a perturbation group and method combination.
# INPUT: 
# sim_arrays: a multiple-iteration simulation array.
# p_group: a single perturbation group index to consider.
# method_idx: a method to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_array: an array of test group selections by iteration and size.
test_step <- function(sim_arrays, p_group_num, method_idx, alpha, t_groups = GROUPS[[4]]$res_Group) {
  selex_array <- array(dim = c(dim(sim_arrays)[1], dim(sim_arrays)[4], length(t_groups)))
  dimnames(selex_array) <- list(dimnames(sim_arrays)[[1]], dimnames(sim_arrays)[[4]], t_groups)
  
  for (i in 1:dim(sim_arrays)[1]) {
    for (j in 1:dim(sim_arrays)[4]) {
      e_vals <- sim_arrays[i,p_group_num,1:length(t_groups),j,method_idx]
      selex_array[i,j,] <- elp(e_vals, t_groups, alpha)
    }
  }
  
  return(selex_array)
}

# Performs p_value testing on multiple iterations on a simulation array at a perturbation group.
# INPUT: 
# sim_arrays: a multiple-iteration simulation array.
# p_group_num: a single perturbation group index to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_array: an array of test group selections by iteration and size.
p_value_test <- function(sim_arrays, p_group, alpha, t_groups = GROUPS[[4]]$res_Group) {
  # Split t_groups into res 1 and res >1 
  r1_groups <- intersect(GROUPS[[4]]$res_Group[GROUPS[[4]]$Resolution == 1], t_groups)
  nr1_groups <- intersect(setdiff(GROUPS[[4]]$res_Group, r1_groups), t_groups)
  
  selex_array <- array(dim = c(dim(sim_arrays)[1], dim(sim_arrays)[4], length(r1_groups)))
  dimnames(selex_array) <- list(dimnames(sim_arrays)[[1]], dimnames(sim_arrays)[[4]], r1_groups)
  
  for (i in 1:dim(sim_arrays)[1]) {
    for (j in 1:dim(sim_arrays)[4]) {
      p_vals <- sim_arrays[i,p_group,1:length(r1_groups),j,1]
      selex_array[i,j,] <- as.integer(p.adjust(p_vals, method = "BH") <= alpha)
    }
  }
  
  # Array to conform p_value array to the rest.
  dummy_array <- array(0, dim = c(dim(sim_arrays)[1], dim(sim_arrays)[4], length(nr1_groups)))
  dimnames(dummy_array) <- list(dimnames(sim_arrays)[[1]], dimnames(sim_arrays)[[4]], nr1_groups)
  
  return(abind(selex_array, dummy_array))
}

# Performs e-value testing on multiple iterations on a simulation array on a single perturbation group on all methods.
# INPUT: 
# sim_arrays: a multiple-iteration simulation array.
# p_group: a single perturbation group to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_arrays: 4D array of test group selections by method, iteration, and size.
omnibus_test <- function(sim_arrays, p_group, alpha, t_groups = GROUPS[[4]]$res_Group) {
  
  selex_arrays <- array(dim = c(11, dim(sim_arrays)[1], dim(sim_arrays)[4], length(t_groups)))
  dimnames(selex_arrays) <- list(dimnames(sim_arrays)[[5]], dimnames(sim_arrays)[[1]], dimnames(sim_arrays)[[4]], t_groups)
  
  selex_array <- p_value_test(sim_arrays, p_group, alpha, t_groups)
  
  selex_arrays[1,,,] <- selex_array
  
  for (i in 2:(dim(sim_arrays)[5])) {
    selex_array <- test_step(sim_arrays, p_group, i, alpha, t_groups)
    
    selex_arrays[i,,,] <- selex_array
  }
  
  return(selex_arrays)
}

# Performs e-value testing on multiple iterations on a simulation array across perturbation groups on a single method.
# INPUT: 
# sim_arrays: a multiple-iteration simulation array.
# method_idx: method to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_arrays: 4D array of test group selections by group, iteration, and size.
omnires_test <- function(sim_arrays, method_idx, alpha, t_groups = GROUPS[[4]]$res_Group) {
  selex_arrays <- array(dim = c(dim(sim_arrays)[2], dim(sim_arrays)[1], dim(sim_arrays)[4], dim(sim_arrays)[3]))
  dimnames(selex_arrays) <- list(dimnames(sim_arrays)[[2]], dimnames(sim_arrays)[[1]], dimnames(sim_arrays)[[4]], dimnames(sim_arrays)[[3]])
  
  for (i in 1:(dim(sim_arrays)[2])) {
    selex_array <- test_step(sim_arrays, i, method_idx, alpha, t_groups)
    p_group <- dimnames(sim_arrays)[[2]][i]
    
    selex_arrays[i,,,] <- selex_array
  }
  return(selex_arrays)
}