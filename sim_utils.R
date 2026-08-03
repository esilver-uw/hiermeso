# Create the data object for the simulation. 
# Goal: e-values for each perturbation group, tested group, signal size, method.
# Methods: Three values of kappa, three values of delta, three values of prior variance.

# Imports
# library(network)
# library(igraph)
library(abind)
library(CVXR)
library(Rglpk)
library(parallel)

# Devise globals, including expected adjacency matrix
SIGMA <- 25
N.SIZE <- 64
N <- 20
GROUP.SIZES <- c(4,8,16)
SIGNAL.SIZES <- 1:10
source("utils.R")
source("e_procedures.R")

GROUPS <- generate_groups(GROUP.SIZES, N.SIZE)
# Create expected adjacency matrix
set.seed(1970)
THETA <- matrix(0, nrow = N.SIZE, ncol = N.SIZE)
for (i in 1:N.SIZE) {
  for (j in i:N.SIZE) {
    THETA[i,j] <- runif(1, -15, 15)
  }
}
THETA[lower.tri(THETA)] = t(THETA)[lower.tri(THETA)]

# Create 4D simulation array, indexed by perturbation group, tested group, signal size, and method.
# INPUT: 
# p_groups: vector of groups to perturb. Default: all.
# OUTPUT: 
# sim_array: a single iteration simulation array.
array_step <- function(p_groups = GROUPS[[4]]$res_Group) {
  sim_array <- array(dim = c(length(p_groups), dim(GROUPS[[4]])[1], length(SIGNAL.SIZES), 11))
  dimnames(sim_array) <- list(p_groups, GROUPS[[4]]$res_Group, SIGNAL.SIZES,
                              c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_delta_1", "lr_delta_2", "lr_delta_3", "lr_prior_1", "lr_prior_2", "lr_prior_3")) 
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
      
      theta_prime <- perturb_expected_matrix(THETA, GROUPS, p_group, size)
      A2 <- sample_network(theta_prime, N)
      # A2 Sample-wise Mean
      A2_sm <- apply(A2, c(1,2), mean)
      
      for (t_group in GROUPS[[4]]$res_Group) {
        k <- k + 1
        
        # Get vector of edges in the group then yield node pairs.
        res_grp <- data.frame(GROUPS[[4]][GROUPS[[4]]$res_Group == t_group,1:2])
        edges <- which(GROUPS[[1]][,res_grp$Resolution] == res_grp$Group_Number)
        
        # This can happen once per group.
        m <- length(edges) * dim(A1)[3]
        
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
        lr_delta_1 <- lr_delta(d_bar, m, pts[1])
        lr_delta_2 <- lr_delta(d_bar, m, pts[2])
        lr_delta_3 <- lr_delta(d_bar, m, pts[3])
        
        priors <- c(10, 25, 50)
        lr_prior_1 <- lr_prior(d_bar, m, priors[1])
        lr_prior_2 <- lr_prior(d_bar, m, priors[2])
        lr_prior_3 <- lr_prior(d_bar, m, priors[3])
        
        sim_array[j,k,i,] <- c(p_val, cal_kappa_1, cal_kappa_2, cal_kappa_3, cal_mix, lr_delta_1, lr_delta_2, lr_delta_3, lr_prior_1, lr_prior_2, lr_prior_3)
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
# sims_array: 5D array with first index iteration, then perturbation group, tested group, signal size, and method.
batchable_array <- function(ct, p_groups = GROUPS[[4]]$res_Group) {
  sims_array <- array(dim = c(ct, length(p_groups), dim(GROUPS[[4]])[1], length(SIGNAL.SIZES), 11))
  for (i in 1:ct) {
    sims_array[i,,,,] <- array_step(p_groups) 
  }
  
  dimnames(sims_array) <- list(1:ct, p_groups, GROUPS[[4]]$res_Group, SIGNAL.SIZES,
                               c("p_value", "cal_kappa_1", "cal_kappa_2", "cal_kappa_3", "cal_mix", "lr_delta_1", "lr_delta_2", "lr_delta_3", "lr_prior_1", "lr_prior_2", "lr_prior_3"))
  return(sims_array)
}

# Performs e-value testing over a vector of method indices for one iteration simulation array at a perturbation group and size combination.
# INPUT: 
# sim_array: a single simulation array
# p_group: a single perturbation group to consider.
# size_idx: a single size index to consider.
# method_idxs: a vector of methods to consider (kappa calibrators (1-3), mixture calibrator (4), point likelihood ratios (1-3), prior likelihood ratios (1-3))
# alpha: alpha level
# OUTPUT: 
# methods_selex: array of selections indexed by method.
test_methods_step <- function(sim_array, p_group, size_idx, method_idxs = 2:11, alpha) {
  # Remove 1 if present.
  method_idxs <- method_idxs[!method_idxs == 1]
  p_group_num <- which(GROUPS[[4]]$res_Group == p_group)
  
  methods_selex <- matrix(nrow = dim(GROUPS[[4]])[1], ncol = length(method_idxs))
  
  # Apply boosting treatment
  
  i <- 0
  for (method_idx in method_idxs) {
    i <- i + 1
    e_vals <- sim_array[p_group,,size_idx,method_idx]
    methods_selex[,i] <- elp(e_vals, GROUPS, alpha)
  }
  rownames(methods_selex) <- GROUPS[[4]]$res_Group
  colnames(methods_selex) <- method_idxs
  
  return(methods_selex)
}

# Performs e-value testing on multiple iterations on a simulation array at a perturbation group and method combination.
# INPUT: 
# sims_array: a multiple-iteration simulation array.
# p_group: a single perturbation group index to consider.
# method_idx: a method to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_array: an array of test group selections by iteration and size.
test_step <- function(sims_array, p_group_num, method_idx, alpha, t_groups = GROUPS) {
  
  selex_array <- array(dim = c(dim(sims_array)[1], dim(sims_array)[4], dim(sims_array)[3]))
  
  for (i in 1:dim(sims_array)[1]) {
    for (j in 1:dim(sims_array)[4]) {
      e_vals <- sims_array[i,p_group_num,,j,method_idx]
      selex_array[i,j,] <- elp(e_vals, t_groups, alpha)
    }
  }
  dimnames(selex_array) <- list(1:dim(sims_array)[1], SIGNAL.SIZES, GROUPS[[4]]$res_Group)
  
  return(selex_array)
}

# Performs p_value testing on multiple iterations on a simulation array at a perturbation group.
# INPUT: 
# sims_array: a multiple-iteration simulation array.
# p_group_num: a single perturbation group index to consider.
# alpha: alpha level.
# t_groups: groups to consider.
# OUTPUT: 
# selex_array: an array of test group selections by iteration and size.
p_value_test <- function(sims_array, p_group_num, alpha, mode = 1, t_groups = GROUPS[[4]]$res_Group[GROUPS[[4]]$Resolution == 1]) {
  selex_array <- array(dim = c(dim(sims_array)[1], dim(sims_array)[4], length(t_groups)))
  dimnames(selex_array) <- list(dimnames(sims_array)[[1]], dimnames(sims_array)[[4]], t_groups)
  
  for (i in 1:dim(sims_array)[1]) {
    for (j in 1:dim(sims_array)[4]) {
      p_vals <- sims_array[i,p_group_num,1:length(t_groups),j,1]
      selex_array[i,j,] <- as.integer(p.adjust(p_vals, method = "BH") <= alpha)
    }
  }
  
  return(selex_array)
}

# Filters to include only true selections.
# INPUT:
# selex_array: an array of test group selections by iteration and size.
# p_group: the true perturbation group.
# semi-true: whether to also include 'semi-true' selections (selections which contain or are contained by the true perturbation group).
# OUTPUT:
# filtered_selex_array: a filtered array of test group selections by iteration and size.
filter_true_selex <- function(selex_array, p_group, semi_true = F) {
  accept <- GROUPS[[2]][p_group][[1]]
  if (semi_true) {
    accept <- union(accept, GROUPS[[4]]$res_Group[sapply(GROUPS[[2]], function(x) p_group %in% x)])
  }
  
  true_idxs <- which(GROUPS[[4]]$res_Group %in% accept)
  
  selex_array[,,-true_idxs] <- 0
  return(selex_array)
}

# Performs e-value testing on multiple iterations on a simulation array on a single perturbation group on all methods.
# INPUT: 
# sims_array: a multiple-iteration simulation array.
# p_group: a single perturbation group to consider.
# alpha: alpha level.
# mode: 1 for no filtering, 2 for true filtering, 3 for semi-true filtering. Default: 1.
# t_groups: groups to consider.
# OUTPUT: 
# selex_arrays: 4D array of test group selections by method, iteration, and size.
omnibus_test <- function(sims_array, p_group, alpha, mode = 1, t_groups = GROUPS) {
  p_group_num <- which(dimnames(sims_array)[[2]] == p_group)
  
  selex_arrays <- array(dim = c(11, dim(sims_array)[1], dim(sims_array)[4], dim(sims_array)[3]))
  dimnames(selex_arrays) <- list(dimnames(sims_array)[[5]], dimnames(sims_array)[[1]], dimnames(sims_array)[[4]], dimnames(sims_array)[[3]])
  
  # Array to conform p_value array to the rest.
  dummy_array <- array(0, dim = c(200, 10, 46))
  dimnames(dummy_array) <- list(dimnames(sims_array)[[1]], dimnames(sims_array)[[4]], GROUPS[[4]]$res_Group[GROUPS[[4]]$Resolution != 1])
    
  selex_array <- abind(p_value_test(sims_array, p_group_num, alpha, mode), dummy_array)
  
  if (mode == 2) {
    selex_array <- filter_true_selex(selex_array, p_group, F)
  } else if (mode == 3) {
    selex_array <- filter_true_selex(selex_array, p_group, T)
  }
  
  selex_arrays[1,,,] <- selex_array
  
  for (i in 2:(dim(sims_array)[5])) {
    selex_array <- test_step(sims_array, p_group_num, i, alpha, t_groups)
    
    if (mode == 2) {
      selex_array <- filter_true_selex(selex_array, p_group, F)
    } else if (mode == 3) {
      selex_array <- filter_true_selex(selex_array, p_group, T)
    }
    
    selex_arrays[i,,,] <- selex_array
  }
  return(selex_arrays)
}

# Performs e-value testing on multiple iterations on a simulation array across perturbation groups on a single method.
# INPUT: 
# sims_array: a multiple-iteration simulation array.
# method_idx: method to consider.
# alpha: alpha level.
# mode: 1 for no filtering, 2 for true filtering, 3 for semi-true filtering. Default: 1.
# t_groups: groups to consider.
# OUTPUT: 
# selex_arrays: 4D array of test group selections by group, iteration, and size.
omnires_test <- function(sims_array, method_idx, alpha, mode = 1, t_groups = GROUPS) {
  selex_arrays <- array(dim = c(dim(sims_array)[2], dim(sims_array)[1], dim(sims_array)[4], dim(sims_array)[3]))
  dimnames(selex_arrays) <- list(dimnames(sims_array)[[2]], dimnames(sims_array)[[1]], dimnames(sims_array)[[4]], dimnames(sims_array)[[3]])
  
  for (i in 1:(dim(sims_array)[2])) {
    selex_array <- test_step(sims_array, i, method_idx, alpha, t_groups)
    p_group <- dimnames(sims_array)[[2]][i]
    
    if (mode == 2) {
      selex_array <- filter_true_selex(selex_array, p_group, F)
    } else if (mode == 3) {
      selex_array <- filter_true_selex(selex_array, p_group, T)
    }
    
    selex_arrays[i,,,] <- selex_array
  }
  return(selex_arrays)
}