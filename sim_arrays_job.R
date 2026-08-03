# We provide the following args: iterations, job #.
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
job_num <- args[2]

source("sim_utils.R")

sim_arrays <- batchable_array(ct)

saveRDS(sim_arrays, file = paste("outputs/sim_arrays_", job_num, ".RData", sep=""))