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
SIGMA <- 50
N.SIZE <- 64
N <- 5
GROUP.SIZES <- c(8,16,32)
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
  A1 <- sample_network(THETA, N)
  
  i = 0
  for (size in SIGNAL.SIZES) {
    i <- i + 1
    j = 0
    for (p_group in p_groups) {
      j <- j + 1
      k = 0
      
      theta_prime <- perturb_expected_matrix(THETA, GROUPS, p_group, size)
      A2 <- sample_network(theta_prime, N)
      
      for (t_group in GROUPS[[4]]$res_Group) {
        k <- k + 1
        
        p_val <- p_value(A1, A2, GROUPS, t_group)
        
        kappas <- c(0.25, 0.5, 0.75)
        cal_kappa_1 <- cal_kappa(p_val, kappas[1])
        cal_kappa_2 <- cal_kappa(p_val, kappas[2])
        cal_kappa_3 <- cal_kappa(p_val, kappas[3])
        
        cal_mix <- cal_mixture(p_val)
        
        pts <- c(2.5, 5, 7.5)
        lr_delta_1 <- lr_delta(A1, A2, GROUPS, t_group, pts[1])
        lr_delta_2 <- lr_delta(A1, A2, GROUPS, t_group, pts[2])
        lr_delta_3 <- lr_delta(A1, A2, GROUPS, t_group, pts[3])
        
        priors <- c(10, 25, 50)
        lr_prior_1 <- lr_prior(A1, A2, GROUPS, t_group, priors[1])
        lr_prior_2 <- lr_prior(A1, A2, GROUPS, t_group, priors[2])
        lr_prior_3 <- lr_prior(A1, A2, GROUPS, t_group, priors[3])
        
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
    e_vals <- sim_array[p_group_num,,size_idx,method_idx]
    methods_selex[,i] <- elp(e_vals, GROUPS, alpha)
  }
  rownames(methods_selex) <- GROUPS[[4]]$res_Group
  colnames(methods_selex) <- method_idxs
  
  return(methods_selex)
}

# Performs e-value testing on multiple iterations on a simulation array at a perturbation group and method combination.
# INPUT: 
# sims_array: a multiple-iteration simulation array.
# p_group: a single perturbation group to consider.
# method_idx: a method to consider.
# alpha: alpha level.
# OUTPUT: 
# selex_array: an array of test group selections by iteration and size.
test_step <- function(sims_array, p_group, method_idx, alpha) {
  p_group_num <- which(GROUPS[[4]]$res_Group == p_group)
  
  selex_array <- array(dim = c(dim(sims_array)[1], dim(sims_array)[4], dim(sims_array)[3]))
  
  for (i in 1:dim(sims_array)[1]) {
    for (j in 1:dim(sims_array)[4]) {
      e_vals <- sims_array[i,p_group,,j,method_idx]
      selex_array[i,j,] <- elp(e_vals, GROUPS, alpha)
    }
  }
  dimnames(selex_array) <- list(1:dim(sims_array)[1], SIGNAL.SIZES, GROUPS[[4]]$res_Group)
  
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
  return(filter_array)
}

# Three extremely basic apply wrappers.
# INPUT: 
# selex_array: an array of test group selections by iteration and size.
# OUTPUT:
# sum or mean or detection matrix, group rows, signal size columns.
sum_selex <- function(selex_array) {
  sum_array <- t(apply(selex_array, c(2,3), sum))
  return(sum_array)
}

# Ditto above.
mean_selex <- function(selex_array) {
  mean_array <- t(apply(selex_array, c(2,3), mean))
  return(mean_array)
}

detex_selex <- function(selex_array) {
  sum_array <- sum_array(selex_array)
  detex_array <- sum_array[sum_array >= 1] <- 1
  return(detex_array)
}

# need a function to stack methods vertically using rowbind. wrap filter and mean/sum to do that.