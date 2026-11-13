# Runs a simulation intended to replicate Rosenblatt et al.'s spatial specificity paradox in complement.

source("sim_utils.R")
source("e_procedures.R")

# Get test statistics for a given edge set
# INPUT:
# M: number of subjects per matrix
# control_mats: array of matrices for the controls
# treatment_mats: array of matrices for the treatments
# edges: vector of edges in the group
# OUTPUT: 
# group_vec: vector of test statistics (p_val, cal_kappa(s), cal_mix, lr_mean(s), lr_prior(s))

# Some thinking will be required: verify the expression of sample standard deviation in terms of V_n and S_n for the difference of means case.

e_vals <- function(M, A1, A2, groups, s, ps) {
  e_vals <- rep(NA, length(groups[[4]][,1]))
  for (i in 1:length(groups[[4]][,1])) {
    lab <- groups[[4]][i,4]
    
    edges <- GROUPS[[3]][[lab]]
    
    control_sm <- apply(A1, c(1,2), sum)
    treatment_sm <- apply(A2, c(1,2), sum)
    
    n.edges <- M * length(edges)
    
    control_mean <- sum(control_sm[edges])/n.edges
    treatment_mean <- sum(treatment_sm[edges])/n.edges
    
    d_bar <- treatment_mean - control_mean
    
    e_vals[i] <- lr_prior(d_bar, n.edges, s, ps)
  }
  
  return(e_vals)
}

# Try bigger networks
# Try weights
# Try building the nested powerset of 16x16 

# Nested powerset

# Source - https://stackoverflow.com/a/18718066
# Posted by Vincent Zoonekynd
# Retrieved 2026-11-11, License - CC BY-SA 3.0

powerset <- function(set) { 
  n <- length(set)
  masks <- 2^(1:n-1)
  lapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] )
}

get_powerset_groups <- function(mat, size = dim(mat)[1]) {
  tri <- matrix(1:16^2, nrow = 16)[upper.tri(matrix(1:16^2, nrow = 16), T)]
  pset_groups <- list()
  splits <- split(1:size, rep(1:2, each = size/2))
  g <- 1
  group_nuclei <- list()
  for (i in 1:2) {
    for (j in 1:2) {
      mat_sub <- mat[splits[[i]], splits[[j]]]
      if (length(intersect(mat_sub, tri)) > 0) {
        group_nuclei[[g]] <- intersect(mat_sub, tri)
        g <- g + 1
      }
    }
  }
  
  pset <- powerset(1:length(group_nuclei))
  pset <- pset[-c(1, length(pset))]
  
  gc <- 1
  for (set in pset) {
    if (length(set) > 0) {
      temp <- NULL
      for (i in set) {
        temp <- append(temp, group_nuclei[[i]])
      }
      pset_groups[[gc]] <- temp
      gc <- gc + 1
    }
  }
  return(pset_groups)
}

Mat <- matrix(1:16^2, nrow = 16)
tri <- matrix(1:16^2, nrow = 16)[upper.tri(matrix(1:16^2, nrow = 16), T)]
pset_groups <- list()
pset_groups[[1]] <- tri

for (size in c(8, 4, 2)) {
  splits <- split(1:16, rep(1:(16/size), each = size))
  for (i in 1:(16/size)) {
    for (j in 1:(16/size)) {
      Mat_split <- matrix(1:16^2, nrow = 16)[splits[[i]], splits[[j]]]
      if (length(intersect(Mat_split, tri)) > 0) {
        pset_groups <- c(pset_groups, get_powerset_groups(Mat_split))
      }
    }
  }
}

names(pset_groups) <- paste0("group_", 1:length(pset_groups))

lcm <- matrix(0, length(pset_groups), length(tri))
for (i in 1:dim(lcm)[1]) {
  indices <- which(tri %in% pset_groups[[i]])
  lcm[i, indices] <- 1
}

pset_e_vals <- function(M, A1, A2, pset_groups, s, ps) {
  e_vals <- rep(NA, length(pset_groups))
  for (i in 1:length(pset_groups)) {
    indices <- pset_groups[[i]]
    
    control_sm <- apply(A1, c(1,2), sum)
    treatment_sm <- apply(A2, c(1,2), sum)
    
    n.edges <- M * length(indices)
    
    control_mean <- sum(control_sm[indices])/n.edges
    treatment_mean <- sum(treatment_sm[indices])/n.edges
    
    d_bar <- treatment_mean - control_mean
    
    e_vals[i] <- lr_prior(d_bar, n.edges, s, ps)
  }
  
  return(e_vals)
}

pset_elp <- function(e_vals, pset_groups, lcm, alpha, lambda = 0) {
  n_groups <- dim(lcm)[1]
  n_base_level <- dim(lcm)[2]
  x <- CVXR::Variable(n_groups, integer = TRUE)
  
  weights <- get_weight(lambda, pset_groups)
  
  objective <- CVXR::Maximize(sum(weights * x))
  
  b <- rep(1, n_base_level)
  constraints <- list(x >= 0,
                      x <= 1,
                      t(lcm) %*% x <= b)
  
  constraints <- c(constraints, list(n_groups - e_vals * alpha * sum(x) <= n_groups * (1 - x)))
  
  problem <- CVXR::Problem(objective = objective, constraints = constraints)
  
  result <- psolve(problem)
  selections <- value(x)
  
  return(selections)
}

# Try log a as well.
# Weight parameter is an exponent from -1 to 1
# Write down some of these cases. "How many small groups do I want to need to beat a big group"
# Also add a half coverage case (2/4 present) and a hidden tiny signal case.
# Can we make the weights depend on the data?
get_weight <- function(lambda, groups) {
  weights <- rep(NA, length(groups))
  for (i in 1:length(groups)) {
    weights[i] <- length(groups[[i]])^lambda
  }
  return(weights)
}

SIZES <- c(2, 4, 8)
GROUPS <- generate_groups(SIZES, 16)
SIGMA <- 20
PRIOR <- 5
M <- 25

SIGNAL.1 <- 20
SIGNAL.2 <- 20
SIGNAL.3 <- 20

T1 <- matrix(0, nrow = 16, ncol = 16)

# Add the signals up
T2 <- T1
T2[1:8,-(1:8)] <- matrix(SIGNAL.1, nrow = 8, ncol = 8)
T2[1:4,-(1:12)] <- matrix(SIGNAL.2, nrow = 4, ncol = 4)
T2[1:2,-(1:14)] <- matrix(SIGNAL.3, nrow = 2, ncol = 2)

A1 <- sample_network(T1, M, SIGMA)

A2 <- sample_network(T2, M, SIGMA)
e_vals(M, A1, A2, GROUPS, SIGMA, PRIOR)
e_vec <- e_vals(M, A1, A2, GROUPS, SIGMA, PRIOR)
GROUPS[[4]]$res_Group[which(elp(e_vec, GROUPS, GROUPS[[4]]$res_Group, 0.05) == 1)]

demo_mat <- T1
for (g in GROUPS[[4]][GROUPS[[4]][,2] == 3,4]) {
  demo_mat[GROUPS[[3]][[g]]] <- as.numeric(strsplit(g, "_")[[1]][4])
}
demo_mat

props <- NULL
for (sigma in c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20)) {
  A1 <- sample_network(T1, M, sigma)
  A2 <- sample_network(T2, M, sigma)
  e_vec <- e_vals(M, A1, A2, GROUPS, sigma, PRIOR)
  selex <- GROUPS[[4]]$res_Group[which(elp(e_vec, GROUPS, GROUPS[[4]]$res_Group, 0.05) == 1)]
  edges <- NULL
  for (s in selex) {
    edges <- append(edges, GROUPS[[3]][[s]])
  }
  
  props <- append(props, length(edges)/length(GROUPS[[3]]$res_3_group_2))
}
names(props) <- c(80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20)
print(props)

# Get the sawtooth graph by continuously varying sigma with replications.

sigmas <- seq(100, 1, by = -1)

props_list <- list()
avg_res_list <- list()

iterations <- 2

for (i in 1:iterations) {
  props <- NULL
  avg_res <- NULL
  indices_mat <- matrix(0, nrow = 256, ncol = length(sigmas))
  j <- 1
  
  for (sigma in sigmas) {
    A1 <- sample_network(T1, M, sigma)
    A2 <- sample_network(T2, M, sigma)
    e_vec <- pset_e_vals(M, A1, A2, pset_groups, sigma, PRIOR)
    selex <- pset_groups[which(pset_elp(e_vec, pset_groups, lcm, 0.05, 0.25) == 1)]
    indices <- NULL
    ress <- NULL
    
    for (s in selex) {
      indices <- append(indices, s)
      ress <- append(ress, length(s))
    }
    
    indices_mat[indices, j] <- 1
    j <- j + 1
    
    props <- append(props, length(intersect(indices, GROUPS[[3]]$res_3_group_2))/length(GROUPS[[3]]$res_3_group_2))
    if (!is.null(ress)) {
      avg_res <- append(avg_res, 1/mean(ress))
    } else {
      avg_res <- append(avg_res, NA)
    }
    
  }
  props_list[[i]] <- props
  avg_res_list[[i]] <- avg_res
}

cond_mean <- \(x) mean(x[x != 0], na.rm = TRUE)

props_means <- apply(do.call(rbind, props_list), 2, cond_mean)
avg_res_means <- apply(do.call(rbind, avg_res_list), 2, cond_mean)

for (i in 1:length(props_means)) {
  if (is.nan(props_means[i])) props_means[i] <- NA
  if (is.nan(avg_res_means[i])) avg_res_means[i] <- NA
}

plot(avg_res_means, type = "l", col = "red")
lines(props_means)


# Get a 1:100 grid and go over it, and use a sample size. Also implement weights. 
# For writing the simulation, set in stone the e_value method. 
# Roadmap for Friday.