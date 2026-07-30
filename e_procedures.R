# Procedures for simulations

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
  n <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  m <- n * dim(A1)[3]
  A1_bar <- mean(apply(A1, c(1,2), mean)[edges])
  A2_bar <- mean(apply(A2, c(1,2), mean)[edges]) # For some reason the old formulation didn't work. Probably some silly indexing.
  
  z_stat <- (A1_bar - A2_bar) * sqrt(m) / (sqrt(2) * SIGMA)
  p_value <- pnorm(z_stat)
  
  return(p_value)
}

# Calibrate p_value to e_value using a kappa parameter.
# INPUT:
# p_value: a p_value, perhaps associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
# k: power form parameter
# OUTPUT: 
# e_value: an e_value calibrated from that p_value.
cal_kappa <- function(p_value, k) {
  if (p_value > 0 & p_value < 1) {
    e_value <- k*p_value^(k-1)
  } else if (p_value == 0) {
    e_value <- 1e+15
  } else if (p_value == 1) {
    e_value <- 0.5
  }
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value <- 1e+15
  }
  
  return(e_value)
}

# Calibrate p_value to e_value.
# INPUT:
# p_value: a p_value, perhaps associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
# OUTPUT: 
# e_value: an e_value calibrated from that p_value.
cal_mixture <- function(p_value) {
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
    e_value <- 1e+15
  }
  
  return(e_value)
}

# Elicit a direct e-value using a delta parameter derived from the mean.
# INPUT: 
# A.1: array of Adjacency Matrices from the first population.
# A.2: array of Adjacency Matrices from the second population.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# pt: prior point estimate for the mean.
# OUTPUT: 
# e_value: an e_value.
lr_delta <- function(A1, A2, groups, g, pt) {
  # Get vector of edges in the group
  res_grp <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_grp$Resolution] == res_grp$Group_Number)
  # Yield node pairs of the edges in question
  n <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  # Use to create d.
  m <- n * dim(A1)[3]
  d <- pt * sqrt(m)
  A1_bar <- mean(apply(A1, c(1,2), mean)[edges])
  A2_bar <- mean(apply(A2, c(1,2), mean)[edges]) # For some reason the old formulation didn't work. Probably some silly indexing.
  
  x_bar <- (A1_bar - A2_bar)
  obs_term <- (m * pt * x_bar)/SIGMA^2
  delta_term <- (m * pt^2)/(2 * SIGMA^2)
  
  e_value <- (exp(obs_term - delta_term) + exp(-obs_term - delta_term))/2
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value <- 1e+15
  }
  
  return(e_value)
}

# Elicit a direct e-value using a mixture over a prior on the mean.
# INPUT: 
# A.1: array of Adjacency Matrices from the first population.
# A.2: array of Adjacency Matrices from the second population.
# groups: list of group attributes.
# g: res_Group identifier matching a group in groups.
# ps: prior sigma; prior variance of the mean.
# OUTPUT: 
# e_value: an e_value.
lr_prior <- function(A1, A2, groups, g, ps) {
  # Get vector of edges in the group
  res_grp <- data.frame(groups[[4]][groups[[4]]$res_Group == g,1:2])
  edges <- which(groups[[1]][,res_grp$Resolution] == res_grp$Group_Number)
  # Yield node pairs of the edges in question
  n <- length(edges)
  
  # Because groups are homogeneous across adjacency matrices, we can simply pool both sample sizes (I think?).
  # Use to create d.
  m <- n * dim(A1)[3]
  A1_bar <- mean(apply(A1, c(1,2), mean)[edges])
  A2_bar <- mean(apply(A2, c(1,2), mean)[edges]) # For some reason the old formulation didn't work. Probably some silly indexing.
  
  x_bar <- (A1_bar - A2_bar)
  norm_term <- 1/(ps*sqrt(m/SIGMA^2 + 1/ps^2))
  expo_term <- (m^2*x_bar^2)/(2*SIGMA^4*sqrt(m/SIGMA^2 + 1/ps^2))
  
  # combine into e-value.
  e_value <- norm_term * exp(expo_term)
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value <- 1e+15
  }
  
  return(e_value)
}
