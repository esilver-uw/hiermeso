# Procedures for simulations

# Given a difference stat and the two populations, return a p-value.
# INPUT:
# d_bar: difference between A1 and A2.
# m: total number of edges averaged.
# OUTPUT: 
# p_value: the p_value associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
p_value <- function(d_bar, m) {
  z_stat <- abs(d_bar * sqrt(m) / (sqrt(2) * SIGMA))
  p_value <- pnorm(z_stat, lower.tail = F)
  
  return(p_value)
}

# Calibrate p_value to e_value using a kappa parameter.
# INPUT:
# p_value: a p_value, perhaps associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
# k: power form parameter.
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
# d_bar: difference between A1 and A2.
# m: total number of edges averaged.
# pt: prior point estimate for the mean.
# OUTPUT: 
# e_value: an e_value.
lr_delta <- function(d_bar, m, pt) {
  obs_term <- (m * pt * d_bar)/(2 * SIGMA^2)
  delta_term <- (m * pt^2)/(4 * SIGMA^2)
  
  e_value <- (exp(obs_term - delta_term) + exp(-obs_term - delta_term))/2
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value <- 1e+15
  }
  
  return(e_value)
}

# Elicit a direct e-value using a mixture over a prior on the mean.
# INPUT: 
# d_bar: difference between A1 and A2.
# m: total number of edges averaged.
# ps: prior sigma; prior variance of the mean.
# OUTPUT: 
# e_value: an e_value.
lr_prior <- function(d_bar, m, ps) {
  # Simply much more tractable in tau notation.
  t <- SIGMA^(-2)
  t_0 <- ps^(-2)
  
  norm_term <- sqrt(t_0)/sqrt(m/2 * t + t_0)
  exp_term_num <- m^2/4*t^2*d_bar^2
  exp_term_denom <- 2*(m/2*t + t_0)
  
  # combine into e-value.
  e_value <- norm_term * exp(exp_term_num/exp_term_denom)
  
  # Threshold cutoff
  if (e_value >= 1e+15) {
    e_value <- 1e+15
  }
  
  return(e_value)
}
