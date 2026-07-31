# Library for functions

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
# node_groups: data frame of dimension n x L specifying node-group membership for each node.
# nodes_edge: data frame of dimension n^2 x 3 linking node pairs to edges.
# group_info: data frame containing group, level, and group-level for each group.
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
  
  return(list(groups, group_subgroups, nodes_edge, group_info))
}

# Returns observations from a given expected adjacency matrix.
# INPUT:
# theta: expected (mean) adjacency matrix
# n: sample size
# OUTPUT: 
# A: array of n observed adjacency matrices
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

# TODO: Be able to perturb in a more targeted manner.
# Apply a perturbation to the parameter adjacency matrix.
# INPUT:
# THETA: parameter adjacency matrix.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# size: size and direction of perturbation to apply.
# OUTPUT:
# theta_prime: perturbed parameter adjacency matrix.
perturb_expected_matrix <- function(theta, groups, g, size) {
  # Get vector of edges in the group
  res_grp <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_grp$Resolution] == res_grp$Group_Number)
  
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
