# We provide the following args: iterations, job #.
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
job_num <- args[2]

# ct <- 200
# job_num <- 1

source("sim_utils.R")

sim_arrays <- batchable_array(ct, c("res_1_group_2", "res_2_group_2", "res_3_group_3"))

saveRDS(sim_arrays, file = paste("outputs/sim_arrays_", job_num, ".RData", sep=""))