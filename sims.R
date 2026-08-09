# Run simulations and save data.

# Devise globals, including expected adjacency matrix
SIGMA <- 25
N.SIZE <- 64
N <- 20
GROUP.SIZES <- c(4,8,16)
# SIGNAL.SIZES <- 0:15
SIGNAL.SIZES <- 0:15
source("sim_utils.R")
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

saveRDS(THETA, file = "outputs/THETA.RData")

# Simulation Arrays
sim_arrays <- batchable_array(200, c("res_1_group_2", "res_2_group_2", "res_3_group_2"))

saveRDS(sim_arrays, file = "outputs/sim_arrays.RData")

# Selection Arrays
selex_array_1 <- omnibus_test(sims_array, "res_1_group_2", 0.05, 2)
selex_array_2 <- omnibus_test(sims_array, "res_2_group_2", 0.05, 2)
selex_array_3 <- omnibus_test(sims_array, "res_3_group_2", 0.05, 2)

saveRDS(selex_array_1, "outputs/selex_array_1.RData")
saveRDS(selex_array_2, "outputs/selex_array_2.RData")
saveRDS(selex_array_3, "outputs/selex_array_3.RData")

# Method Arrays (Pruned methods; not all are interesting.)
method_array_2 <- omnires_test(sims_array, 2, 0.05, 2)
# method_array_3 <- omnires_test(sims_array, 3, 0.05, 2)
# method_array_4 <- omnires_test(sims_array, 4, 0.05, 2)
method_array_5 <- omnires_test(sims_array, 5, 0.05, 2)
method_array_6 <- omnires_test(sims_array, 6, 0.05, 2)
method_array_7 <- omnires_test(sims_array, 7, 0.05, 2)
method_array_8 <- omnires_test(sims_array, 8, 0.05, 2)
method_array_9 <- omnires_test(sims_array, 9, 0.05, 2)
method_array_10 <- omnires_test(sims_array, 10, 0.05, 2)
# method_array_11 <- omnires_test(sims_array, 11, 0.05, 2)

saveRDS(method_array_2, "outputs/method_array_2.RData")
saveRDS(method_array_5, "outputs/method_array_5.RData")
saveRDS(method_array_6, "outputs/method_array_6.RData")
saveRDS(method_array_7, "outputs/method_array_7.RData")
saveRDS(method_array_8, "outputs/method_array_8.RData")
saveRDS(method_array_9, "outputs/method_array_9.RData")
saveRDS(method_array_10, "outputs/method_array_10.RData")
