# Specification: We want a core 'simulation' function that takes in strength of signal and size of signal, that is to say, what level of scale the signal should be on. 
# The intent is to demonstrate the eLP algorithm working on different mesoSCALES. The only reason not to reject at the smallest possible scale is if there is
# insufficient power to do so.
# We will use a preset family of random networks derived from adjacency matrices A and A', where A' is an adjacency matrix derived from A by applying a perturbation with
# signal strength x and signal level l. 
# We will use hypotheses H_i^l where l=1 is the lowest level, defined recursively. For simplicity we will use binary splitting: that is, H_i^2 = [H_i^1, H_j^1], and 
# so forth. A single hypothesis H_i^1 = [N_1, N_2] where N_i is a node set. In other words, a base level hypothesis is the set of edges between those node sets.
# In this way, a single hypothesis is a region of the adjacency matrix. For now, we will be dealing with hierarchical (tree-like) hypothesis chains.
# The perturbation will be applied to the edges in a hypothesis. 
# The e-value function [will be unspecified for now.]

# Imports
library(network)
library(igraph)
library(abind)
library(CVXR)
library(Rglpk)

# Devise globals, including parameter adjacency matrix
SIGMA <- 10
N.SIZE <- 64

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
# L: number of levels.
# OUTPUT:
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# node_groups: data frame of dimension n x L specifying node-group membership for each node.
# nodes_edge: data frame of dimension n^2 x 3 linking node pairs to edges.
# group_info: data frame containing group, level, and group-level for each group.
generate_groups <- function(GROUP_SIZES, n, L) {
  node_groups <- matrix(0, n, L)
  for (l in 1:L) {
    node_groups[,l] <- rep(1:(n/GROUP_SIZES[l]), each = GROUP_SIZES[l])
  }
  
  # nodes_edge
  node_groups <- matrix(0, n, L)
  for (l in 1:L) {
    node_groups[,l] <- rep(1:(n/GROUP_SIZES[l]), each = GROUP_SIZES[l])
  }
  
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
  
  return(list(groups, node_groups, nodes_edge, group_info))
}

GROUPS <- generate_groups(c(1, 2, 4), 8, 3)

sample_network <- function(theta) {
  A <- matrix(NA, nrow = nrow(theta), ncol = ncol(theta))
  for (i in 1:nrow(theta)) {
    for (j in 1:ncol(theta)) {
      A[i,j] = rnorm(1, theta[i,j], SIGMA)
    }
  }
  return(A)
}

# Given a group and the two populations, return an e-value.
# INPUT:
# A.1: array of Adjacency Matrices from the first population.
# A.2: array of Adjacency Matrices from the second population.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# OUTPUT: 
# e_value: the e_value associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
e_value <- function(A1, A2, groups, g) {
  # Get vector of edges in the group
  res_Group <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_Group$Resolution] == res_Group$Group_Number)
  # Yield node pairs of the edges in question
  node_pairs <- groups[[3]][edges,1:2]
  m <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  n <- m * dim(A1)[3]
  A1_bar <- apply(A1[node_pairs], c(1,2), mean)
  A2_bar <- apply(A2[node_pairs], c(1,2), mean)
  s_A1 <- apply(A1[node_pairs], c(1,2), sd)
  s_A2 <- apply(A2[node_pairs], c(1,2), sd)
  s_pooled <- sqrt((s_A1^2 + s_A2^2)/n)
  t.stat <- (A1_bar - A2_bar) / s_pooled
  p_value <- pt(t.stat, 2*n - 2)
  
  # Calibrate p_value to e_value.
  num <- 1 - p_value + p_value * log(p_value)
  denom <- p_value * (-log(p_value))^2
  return(num/denom)
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
    res_Group <- data.frame(groups[[4]][i,1:2])
    indices <- which(groups[[1]][,res_Group$Resolution] == res_Group$Group_Number)
    print(groups[[1]][indices,])
    for (j in 1:dim(location_constraint_matrix)[2]) {
      if (j %in% groups[[1]][indices,1]) {
        location_constraint_matrix[i,j] <- 1
      }
    }
  }
  return(location_constraint_matrix)
}

# TODO: Implement weighting per Gablenz & Sabatti. 
# Run eLP
# INPUT:
# e_vals: vector of e_values by aligned with groups
# groups: list of group attributes.
# alpha: alpha level of the test.
# OUTPUT: 
# detections: hypotheses rejected by the algorithm.
elp <- function(e_vals, groups, alpha) {
  # Get number of base level hypotheses and number of total hypotheses
  n_base_level <- length(groups[[4]][groups[[4]]$Resolution == 1,2])
  n_groups <- dim(groups[[4]][,2:3])[1]
  
  x <- CVXR::Variable(n_groups, integer = TRUE)
  objective <- CVXR::Maximize(sum(x))
  
  location_constraint_matrix <- create_lcm(groups, n_base_level, n_groups)
  
  b <- rep(1, n_base_level)
  constraints <- list(x >= 0,
                      x <= 1,
                      t(location_constraint_matrix) %*% x <= b)
  
  constraints <- c(constraints, list(n_groups - e_vals * alpha * sum(x) <= n_groups * (1 - x)))
  
  problem <- CVXR::Problem(objective = objective, constraints = constraints)
  
  result <- CVXR::solve(problem, solver = 'GLPK')
  
  selections <- as.numeric(result$getValue(x))
  detections <- groups[[4]][which(selections == 1)]
  
  return(detections)
}

# Main simulation: iterate over a vector of perturbation SIZES, generating (with set seed) the two populations of networks and applying eLP as implemented above.
# Output rejection set for each, then display number of rejections and resolution of rejections. Permit selection of group. Try sprinkling equal little signals
# everywhere and see if increasing them yields different resolution rejections. We want the output to lead to a density plot per resolution. Experiment with weighting.
# INPUT: 
# groups: list of group attributes.
# alpha: alpha level of the test.
# THETA: unperturbed parameter adjacency matrix.
# perturb_g: res_Group identifier matching a group in groups.
# sizes: vector of sizes (including sign for direction) of perturbations.
# OUTPUT: 
# detections: list (or something) of sizes & resultant detections.
simulation <- function(groups, alpha, theta, perturb_g, sizes) {
  detections <- NA
  # Iterate over the sizes (magnitude + direction) of perturbations to apply
  for (size in sizes) {
    theta_prime <- perturb_parameter_matrix(theta = theta, groups = groups, g = perturb_g, size = size)
    A1 <- sample_network(theta = theta)
    A2 <- sample_network(theta = theta_prime)
    
    e_vals <- NULL
    for (g in groups) {
      e_vals <- append(e_vals, e_value(theta, theta_prime, groups, g))
    }
    
    detections <- append(detections, elp(e_vals, groups, alpha))
  }
  
  return(detections)
}

detections <- simulation(GROUPS, 0.05, THETA, "res_1_group_1", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
sample_network(THETA)
