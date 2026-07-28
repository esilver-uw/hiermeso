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
SIGNAL.SIZES <- 1:5*5
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

# perturbation group, tested group, signal size, method
sim_array <- array(dim = c(dim(GROUPS[[4]])[1], dim(GROUPS[[4]])[1], length(SIGNAL.SIZES), 11))

i = 0
for (size in SIGNAL.SIZES) {
  i <- i + 1
  j = 0
  for (p_group in GROUPS[[4]]$res_Group) {
    j <- j + 1
    k = 0
    
    theta_prime <- perturb_expected_matrix(THETA, GROUPS, p_group, size)
    A1 <- sample_network(THETA, N)
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
      
      priors <- c(25, 50, 75)
      lr_prior_1 <- lr_prior(A1, A2, GROUPS, t_group, priors[1])
      lr_prior_2 <- lr_prior(A1, A2, GROUPS, t_group, priors[2])
      lr_prior_3 <- lr_prior(A1, A2, GROUPS, t_group, priors[3])
      
      sim_array[j,k,i,] <- c(p_val, cal_kappa_1, cal_kappa_2, cal_kappa_3, cal_mix, lr_delta_1, lr_delta_2, lr_delta_3, lr_prior_1, lr_prior_2, lr_prior_3)
    }
  }
}

sim_test <- function(sim_array, p_group, size_idx, method_idx, alpha) {
  e_vals <- sim_array[p_group,,size_idx,method_idx]
  names(e_vals) <- GROUPS[[4]]$res_Group
  
  print(e_vals)
  
  # Apply boosting treatment
  
  selections <- elp(e_vals, GROUPS, alpha)
}
