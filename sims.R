# Run simulations and save data.

source("sim_utils.R")

JOB.NUM <- 2

sims_array <- readRDS(paste0("outputs/sim_arrays_", JOB.NUM, ".RData"))

# Selection Arrays
selex_array_1 <- omnibus_test(sims_array, "res_1_group_2", 0.05, 2)
selex_array_2 <- omnibus_test(sims_array, "res_2_group_2", 0.05, 2)
selex_array_3 <- omnibus_test(sims_array, "res_3_group_3", 0.05, 2)

saveRDS(selex_array_1, "outputs/selex_array_1.RData")
saveRDS(selex_array_2, "outputs/selex_array_2.RData")
saveRDS(selex_array_3, "outputs/selex_array_3.RData")

# Method Arrays (Pruned methods; not all are interesting.)
method_array_2 <- omnires_test(sims_array, 2, 0.05, 2)
# method_array_3 <- omnires_test(sims_array, 3, 0.05, 2)
# method_array_4 <- omnires_test(sims_array, 4, 0.05, 2)
method_array_5 <- omnires_test(sims_array, 5, 0.05, 2)
method_array_6 <- omnires_test(sims_array, 6, 0.05, 2)
# method_array_7 <- omnires_test(sims_array, 7, 0.05, 2)
method_array_8 <- omnires_test(sims_array, 8, 0.05, 2)
method_array_9 <- omnires_test(sims_array, 9, 0.05, 2)
method_array_10 <- omnires_test(sims_array, 10, 0.05, 2)
# method_array_11 <- omnires_test(sims_array, 11, 0.05, 2)

saveRDS(method_array_2, "outputs/method_array_2.RData")
saveRDS(method_array_5, "outputs/method_array_5.RData")
saveRDS(method_array_6, "outputs/method_array_6.RData")
saveRDS(method_array_8, "outputs/method_array_8.RData")
saveRDS(method_array_9, "outputs/method_array_9.RData")
saveRDS(method_array_10, "outputs/method_array_10.RData")