# Procedures for simulations

# Given a difference stat and the two populations, return a p-value.
# INPUT:
# d_bar: difference between A1 and A2.
# m: total number of edges averaged.
# OUTPUT: 
# p_value: the p_value associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
p_value <- function(d_bar, m) {
  d_bar <- abs(d_bar)
  p_value <- 2*pnorm(d_bar, sd = (sqrt(2) * SIGMA / sqrt(m)), lower.tail = F)
  
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
    e_value <- 1e+10
  } else if (p_value == 1) {
    e_value <- 0.5
  }
  
  # Threshold cutoff
  if (e_value >= 1e+10) {
    e_value <- 1e+10
  }
  if (e_value <= 1e-5) {
    e_value <- 1e-5
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
    e_value <- 1e+10
  } else if (p_value == 1) {
    e_value <- 0.5
  }
  
  # Threshold cutoff
  if (e_value >= 1e+10) {
    e_value <- 1e+10
  }
  if (e_value <= 1e-5) {
    e_value <- 1e-5
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
  if (e_value >= 1e+10) {
    e_value <- 1e+10
  }
  if (e_value <= 1e-5) {
    e_value <- 1e-5
  }
  
  return(e_value)
}

# Elicit a direct e-value using a mixture over a prior on the mean.
# INPUT: 
# d_bar: difference between A1 and A2.
# s: sigma, or in our bad behaviour case the estimate.
# m: total number of edges averaged.
# ps: prior sigma; prior variance of the mean.
# OUTPUT: 
# e_value: an e_value.
lr_prior <- function(d_bar, m, s, ps) {
  # Simply much more tractable in tau notation.
  t <- s^(-2)
  t_0 <- ps^(-2)
  
  norm_term <- sqrt(t_0)/sqrt(m/2 * t + t_0)
  exp_term_num <- m^2/4*t^2*d_bar^2
  exp_term_denom <- 2*(m/2*t + t_0)
  
  # combine into e-value.
  e_value <- norm_term * exp(exp_term_num/exp_term_denom)
  
  # Threshold cutoff
  if (e_value >= 1e+10) {
    e_value <- 1e+10
  }
  if (e_value <= 1e-5) {
    e_value <- 1e-5
  }
  
  return(e_value)
}

# Given a difference stat and the two populations, return a p-value, assuming t-distribution.
# INPUT:
# d_bar: difference between A1 and A2.
# s: standard deviation
# m: total number of edges averaged.
# OUTPUT: 
# p_value: the p_value associated with the group-wise null hypothesis that A.1 and A.2 are generated from the same parameter matrix.
p_t <- function(d_bar, s, m) {
  t <- abs(d_bar)/(2 * s / sqrt(m))
  p_value <- 2*pt(t, df = 2*m - 2, lower.tail = F)
  
  return(p_value)
}

# For now, patch using abs.

# Elicit a direct e-value using a mixture over a prior on the mean, assuming t-distribution.
# INPUT: 
# d: difference between A1 and A2 means.
# s: pooled standard deviation of both A1 and A2
# m: total number of edges averaged.
# c: prior sigma; prior variance of the mean.
# OUTPUT: 
# e_value: an e_value.
lr_t <- function(d, s, m, c) {
  # Simply much more tractable in tau notation.
  
  k <- 2*m-2
  
  norm <- (1 + c^2)
  num <- 2 * s^2 * (1 + c^2) * k + m * d^2
  denom <- 2 * s^2 * (1 + c^2) * k + m * (1 + c^2) * d^2
  
  # combine into e-value.
  e_value <- norm^(-1/2) * (num/denom)^(-(k+1)/2)
  
  # Threshold cutoff
  if (e_value >= 1e+10) {
    e_value <- 1e+10
  }
  if (e_value <= 1e-5) {
    e_value <- 1e-5
  }
  
  return(e_value)
}
