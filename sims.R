# Run simulations and save data.

source("sim_utils.R")

JOB.NUM <- 3

sim_arrays <- batchable_array(200, c("res_1_group_2", "res_2_group_2", "res_3_group_2"))

saveRDS(sim_arrays, file = paste("outputs/sim_arrays_", JOB.NUM, "_vm_2.RData", sep=""))

sims_array <- readRDS(paste0("outputs/sim_arrays_", JOB.NUM, "_vm_2.RData"))

# Selection Arrays
selex_array_1 <- omnibus_test(sims_array, "res_1_group_2", 0.05, 2)
selex_array_2 <- omnibus_test(sims_array, "res_2_group_2", 0.05, 2)
selex_array_3 <- omnibus_test(sims_array, "res_3_group_2", 0.05, 2)

saveRDS(selex_array_1, "outputs/selex_array_1_vm_2.RData")
saveRDS(selex_array_2, "outputs/selex_array_2_vm_2.RData")
saveRDS(selex_array_3, "outputs/selex_array_3_vm_2.RData")

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

saveRDS(method_array_2, "outputs/method_array_2_vm_2.RData")
saveRDS(method_array_5, "outputs/method_array_5_vm_2.RData")
saveRDS(method_array_6, "outputs/method_array_6_vm_2.RData")
saveRDS(method_array_7, "outputs/method_array_7_vm_2.RData")
saveRDS(method_array_8, "outputs/method_array_8_vm_2.RData")
saveRDS(method_array_9, "outputs/method_array_9_vm_2.RData")
saveRDS(method_array_10, "outputs/method_array_10_vm_2.RData")

# False detex rate: set signal sizes to 0 only
sims_array_fdr <- sims_array[,,,1,,drop=F]
# Choice of perturbation group is arbitrary since all perturbations are size 0. Do not filter any selections.
# Thus: all detections are false detections.
res_1_grps <- generate_groups(c(4), N.SIZE)
res_12_grps <- generate_groups(c(4,8), N.SIZE)

fdr_selex_array_1 <- omnibus_test(sims_array_fdr, "res_1_group_2", 0.05, 1, res_1_grps)
fdr_selex_array_2 <- omnibus_test(sims_array_fdr, "res_1_group_2", 0.05, 1, res_12_grps)
fdr_selex_array_3 <- omnibus_test(sims_array_fdr, "res_1_group_2", 0.05, 1, GROUPS)

saveRDS(fdr_selex_array_1, "outputs/fdr_selex_array_1_vm_2.RData")
saveRDS(fdr_selex_array_2, "outputs/fdr_selex_array_2_vm_2.RData")
saveRDS(fdr_selex_array_3, "outputs/fdr_selex_array_3_vm_2.RData")
