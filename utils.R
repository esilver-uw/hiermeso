# Library for functions

# Imports
# library(network)
# library(igraph)
library(abind)
library(CVXR)
library(Rglpk)
library(parallel)

# Devise globals, including parameter adjacency matrix
SIGMA <- 50
N.SIZE <- 64
N <- 20
GROUP.SIZES <- c(8,16,32)

# Helper from Stack Overflow:
# Source - https://stackoverflow.com/a/8189441
# Posted by Ken Williams, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-18, License - CC BY-SA 4.0
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Create parameter adjacency matrix
set.seed(1970)
THETA <- matrix(0, nrow = N.SIZE, ncol = N.SIZE)
for (i in 1:N.SIZE) {
  for (j in i:N.SIZE) {
    THETA[i,j] <- runif(1, -15, 15)
  }
}
THETA[lower.tri(THETA)] = t(THETA)[lower.tri(THETA)]

# Create groups, inspired by KeLP architecture.
# INPUT:
# GROUP_SIZES: vector of group sizes.
# n: number of nodes.
# OUTPUT:
# for L length of GROUP_SIZES
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# node_groups: data frame of dimension n x L specifying node-group membership for each node.
# nodes_edge: data frame of dimension n^2 x 3 linking node pairs to edges.
# group_info: data frame containing group, level, and group-level for each group.
generate_groups <- function(GROUP_SIZES, n) {
  # number of levels
  L <- length(GROUP_SIZES)
  
  # nodes_edge
  nodes_edge <- cbind(expand.grid(1:n, 1:n), 1:n^2)
  
  # groups will be all combinations of node_groups. 
  groups <- matrix(0, n^2, L)
  for (l in 1:L) {
    index <- 1
    for (i in 1:(n/GROUP_SIZES[l])) {
      i_valid <- (i*GROUP_SIZES[l]) + 1 - 1:GROUP_SIZES[l]
      for (j in i:(n/GROUP_SIZES[l])) {
        # Group size is tolerance frame. If 
        j_valid <- (j*GROUP_SIZES[l]) + 1 - 1:GROUP_SIZES[l]
        groups[(nodes_edge[,1] %in% i_valid & nodes_edge[,2] %in% j_valid) | (nodes_edge[,1] %in% j_valid & nodes_edge[,2] %in% i_valid),l] <- index
        index <- index + 1
      }
    }
  }
  
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
  
  return(list(groups, group_subgroups, nodes_edge, group_info))
}

# Create groups for upper triangle, inspired by KeLP architecture.
# INPUT:
# GROUP_SIZES: vector of group sizes.
# n: number of nodes.
# OUTPUT:
# for L length of GROUP_SIZES
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# node_groups: data frame of dimension n x L specifying node-group membership for each node.
# nodes_edge: data frame of dimension n^2 x 3 linking node pairs to edges.
# group_info: data frame containing group, level, and group-level for each group.
generate_groups <- function(GROUP_SIZES, n) {
  # number of levels
  L <- length(GROUP_SIZES)
  
  # nodes_edge
  nodes_edge <- NULL
  for (i in 1:n) {
    for (j in i:n) {
      nodes_edge <- rbind(nodes_edge, c(i,j))
    }
  }
  
  # groups will be all combinations of node_groups. 
  groups <- NULL
  groups <- cbind(matrix(0, n^2, L), expand.grid(1:n,1:n))
  for (l in 1:L) {
    index <- 1
    for (i in 1:(n/GROUP_SIZES[l])) {
      i_valid <- (i*GROUP_SIZES[l]) + 1 - 1:GROUP_SIZES[l]
      for (j in i:(n/GROUP_SIZES[l])) {
        # Group size is tolerance frame. If
        j_valid <- (j*GROUP_SIZES[l]) + 1 - 1:GROUP_SIZES[l]
        # groups[(nodes_edge[,1] %in% i_valid & nodes_edge[,2] %in% j_valid) | (nodes_edge[,1] %in% j_valid & nodes_edge[,2] %in% i_valid),l] <- index
        groups[(nodes_edge[,1] %in% i_valid & nodes_edge[,2] %in% j_valid),l] <- index
        index <- index + 1
      }
    }
  }
  
  # me caveman. me de-duplicate up to ordering.
  groups <- groups[groups[,L+1] %in% nodes_edge[,1] & groups[,L+2] %in% nodes_edge[,2],1:3]
  
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
  
  return(list(groups, group_subgroups, nodes_edge, group_info))
}

# Returns an observation from a given parameter adjacency matrix.
sample_network <- function(theta, n) {
  A <- array(NA, c(nrow(theta), ncol(theta), n))
  for (k in 1:n) {
    for (i in 1:nrow(theta)) {
      for (j in 1:ncol(theta)) {
        A[i,j,k] = rnorm(1, theta[i,j], SIGMA)
      }
    }
  }
  return(A)
}

# Given a group and the two populations, return a p-value.
# INPUT:
# A.1: array of Adjacency Matrices from the first population.
# A.2: array of Adjacency Matrices from the second population.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# OUTPUT: 
# p_value: the p_value associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
p_value <- function(A1, A2, groups, g) {
  # Get vector of edges in the group
  res_grp <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_grp$Resolution] == res_grp$Group_Number)
  # Yield node pairs of the edges in question
  m <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  n <- m * dim(A1)[3]
  A1_bar <- mean(apply(A1, c(1,2), mean)[edges])
  A2_bar <- mean(apply(A2, c(1,2), mean)[edges]) # For some reason the old formulation didn't work. Probably some silly indexing.
  
  z.stat <- (A1_bar - A2_bar) / (SIGMA / sqrt(n))
  p_value <- pnorm(z.stat)
  
  return(p_value)
}

# Calibrate p_value to e_value.
# INPUT:
# p_value: a p_value, perhaps associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
# OUTPUT: 
# e_value: an e_value calibrated from that p_value.
e_value_cal <- function(p_value) {
  if (p_value > 0 & p_value < 1) {
    num <- 1 - p_value + p_value * log(p_value)
    denom <- p_value * (-log(p_value))^2
    e_value <- num/denom
  } else if (p_value == 0) {
    e_value <- 1e+15
  } else if (p_value == 1) {
    e_value <- 0.5
  }
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value = 1e+15
  }
  
  return(e_value)
}

e_value_dir <- function(A1, A2, groups, g, d) {
  # Get vector of edges in the group
  res_grp <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_grp$Resolution] == res_grp$Group_Number)
  # Yield node pairs of the edges in question
  m <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  n <- m * dim(A1)[3]
  A1_bar <- mean(apply(A1, c(1,2), mean)[edges])
  A2_bar <- mean(apply(A2, c(1,2), mean)[edges]) # For some reason the old formulation didn't work. Probably some silly indexing.
  
  z.stat <- (A1_bar - A2_bar) / (SIGMA / sqrt(n))
  
  print(d)
  print(z.stat)
  
  e_value <- exp(d * z.stat - d^2/2) + exp(-d * z.stat - d^2/2)
    
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value = 1e+15
  }
  
  return(e_value)
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
perturb_parameter_matrix <- function(theta, groups, g, size) {
  # Get vector of edges in the group
  res_Group <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_Group$Resolution] == res_Group$Group_Number)
  
  # Apply perturbation
  theta_prime <- theta
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
# groups: list of group attributes.
# alpha: alpha level of the test.
# OUTPUT: 
# detections: hypotheses rejected by the algorithm.
elp <- function(e_vals, groups, alpha) {
  # Get number of base level hypotheses and number of total hypotheses
  n_base_level <- length(groups[[4]][groups[[4]]$Resolution == 1,2])
  n_groups <- dim(groups[[4]])[1]
  
  x <- CVXR::Variable(n_groups, integer = TRUE)
  objective <- CVXR::Maximize(sum(x))
  
  location_constraint_matrix <- create_lcm(groups, n_base_level, n_groups)
  
  b <- rep(1, n_base_level)
  constraints <- list(x >= 0,
                      x <= 1,
                      t(location_constraint_matrix) %*% x <= b)
  
  constraints <- c(constraints, list(n_groups - e_vals * alpha * sum(x) <= n_groups * (1 - x)))
  
  problem <- CVXR::Problem(objective = objective, constraints = constraints)
  
  result <- psolve(problem, solver = 'GLPK_MI')
  selections <- value(x)
  
  return(selections)
}

# Helper function to get oracle delta
# INPUT: 
# size: expected size of the perturbation
# groups: list of group attributes
# treatment_g: group in question
get_delta <- function(size, groups, treatment_g) {
  group_num <- groups[[4]][groups[[4]]$res_Group == treatment_g,1]
  res <- groups[[4]][groups[[4]]$res_Group == treatment_g,2]
  group.size <- sum(groups[[1]][res] == group_num)
  
  
  delta <- (sqrt(group.size * N.SIZE) * size) / SIGMA
  return(delta)
}
