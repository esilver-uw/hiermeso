# Simulation to determine optimal choice of parameter d for the e-value to determine how the e-value behaves. Approach: approximate log-e power via sampling.
# Spec: same as ordinary simulation.

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

e_value <- function(A1, A2, groups, g, d) {
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
  
  return(exp(d * z.stat - d^2/2) + exp(-d * z.stat - d^2/2))
}

# Main simulation: iterate over sizes and 3 treatment resolution groups. Output mean log e_value.
simulation <- function(groups, theta, treatment_g, sizes, iterations, ds) {
  e_powers <- matrix(nrow = length(sizes), ncol = length(ds))
  # Get treated parameter matrix and network observations
  for (i in 1:length(sizes)) {
    for (j in 1:length(ds)) {
      theta_prime <- perturb_parameter_matrix(theta = theta, groups = groups, g = treatment_g, size = sizes[i])
      e_vals <- rep(0, iterations)
      
      for (k in 1:iterations) {
        A1 <- sample_network(theta = theta, N)
        A2 <- sample_network(theta = theta_prime, N)
        
        e_vals[k] <- e_value(A1, A2, groups, treatment_g, ds[j])
        
      }
      e_powers[i,j] <- mean(log(e_vals))
    }
  }
  
  return(e_powers)
}

simple_sim <- function(n, diff, sigma, iterations, d,) {
  e_val_obs <- rep(NA, length(iterations))
  for (j in 1:length(iterations)) {
    x1s <- rnorm(n, sd = sigma)
    x2s <- rnorm(n, mean = diff, sd = sigma)
    
    bar.x1 <- mean(x1s)
    bar.x2 <- mean(x2s)
    
    z.stat <- (bar.x2 - bar.x1)*sqrt(n)/sigma
    
    e_val_obs[j] <- (exp(d * z.stat - d^2/2) + exp(-d * z.stat - d^2/2))/2
  }
  e_val <- mean(e_val_obs)
  return(log(e_val))
}



sizes <- 20
hyp_ds <- 5:15/10
simulation(GROUPS, THETA, "res_1_group_2", sizes, 20, hyp_ds)

simmat3 <- simulation(GROUPS, THETA, "res_3_group_2", sizes, 20, hyp_ds)
rownames(simmat3) <- sizes
colnames(simmat3) <- hyp_ds

simmat3.ftr <- simmat3
simmat3.ftr[simmat3 <= log(20)] <- NA

# Consider a pilot study: Use some data to learn delta, then compute with delta. Do comparisons for different deltas. Try a simpler simulation. 
