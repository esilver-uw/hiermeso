# Procedures for simulations

source("utils.R")

# Performs eLP procedure on a subset of the groups.
# INPUT:
# A1/A2: Observation adjacency matrix stacks
# groups: Attribute list of all groups
# groups_subset: List of groups to perform eLP on
# size: an estimate/oracle value of perturbation size. if not provided, calibrate a p_value, 
# if provided, use direct e-value method.
# OUTPUT: 
# result: a data frame with two columns for groups and their associated e-values
proc_elp <- function(A1, A2, groups, groups_subset, size = F) {
  # Get e_values
  n.groups <- length(groups_subset)
  
  e_vals <- rep(NA, n.groups)
  if (size) {
    for (i in 1:n.groups) {
      delta <- get_delta(size, groups, groups_subset[i])
      e_vals[i] <- e_value_dir(A1, A2, groups, groups_subset[i], delta)
    }
  } else {
    for (i in 1:n.groups) {
      p_val <- p_value(A1, A2, groups, groups_subset[i])
      e_vals[i] <- e_Value_cal(p_val)
    }
  }
  
  selections <- elp(e_vals, groups_subset, 0.05)
  result <- data.frame("Selections" = selections == 1, "e-values" = e_vals)
  
  return(result)
}

# Wrapper for groups = groups_subset. Inherits definition from above, sans groups_subset.
proc_omnibus_elp <- function(A1, A2, groups, size = F) {
  return(proc_elp(A1, A2, groups, groups[[4]]$res_Group, size))
}

# Performs eLP procedure starting with the largest resolution.
proc_seq_elp <- function(A1, A2, groups, size = F) {
  resolutions <- unique(groups[[4]]$Resolution)
  significant <- T
  i = 1
  selections_matrix <- data.frame()
  
  while(significant) {
    groups_subset <- groups[[4]][groups[[4]]$Resolution == resolutions[i]]
    
    result <- proc_elp(A1, A2, groups, groups_subset, size)
    selections <- result$Selections
    if (sum(selections) > 0) {
      selections_matrix[,resolutions[i]] <- selections
    }
  }
}

theta_prime <- perturb_parameter_matrix(theta = THETA, groups = GROUPS, g = "res_3_group_2", size = 8)

A1 <- sample_network(theta = THETA, N)
A2 <- sample_network(theta = theta_prime, N)

E_VALS <- proc_elp(A1, A2, GROUPS, GROUPS[[4]]$res_Group, 10)

selections <- elp(E_VALS, GROUPS, 0.05)
GROUPS[[4]][which(selections == 1),]
E_VALS[selections == 1]
