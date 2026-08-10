# Run simulations and save data.

source("viz_utils.R")

# Devise globals, including expected adjacency matrix
SIGMA <- 25
N.SIZE <- 64
N <- 20
GROUP.SIZES <- c(4,8,16)
# SIGNAL.SIZES <- 0:15
SIGNAL.SIZES <- 0:15
source("utils.R")
source("e_procedures.R")

GROUPS <- generate_groups(GROUP.SIZES, N.SIZE)

# Read in data from simulation and generate the plots for presentation

sims_array <- readRDS("outputs/sim_arrays_3.RData")

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
selex_array_1 <- readRDS("outputs/selex_array_1.RData")
selex_array_2 <- readRDS("outputs/selex_array_2.RData")
selex_array_3 <- readRDS("outputs/selex_array_3.RData")

# You may note '_vm' appended to certain files. I carried out simulations using both
# a computer cluster and my laptop.
method_array_2 <- readRDS("outputs/method_array_2.RData")
method_array_5 <- readRDS("outputs/method_array_5.RData")
method_array_6 <- readRDS("outputs/method_array_6.RData")
method_array_7 <- readRDS("outputs/method_array_7.RData")
method_array_8 <- readRDS("outputs/method_array_8.RData")
method_array_9 <- readRDS("outputs/method_array_9.RData")
method_array_10 <- readRDS("outputs/method_array_10.RData")

fdr_selex_array_1 <- readRDS("outputs/fdr_selex_array_1.RData")
fdr_selex_array_2 <- readRDS("outputs/fdr_selex_array_2.RData")
fdr_selex_array_3 <- readRDS("outputs/fdr_selex_array_3.RData")

# Get true/false detections.

td_array_1 <- filter_true_selex(selex_array_1, "res_1_group_2", T)
fd_array_1 <- filter_true_selex(selex_array_1, "res_1_group_2", T, T)
td_array_2 <- filter_true_selex(selex_array_2, "res_2_group_2", T)
fd_array_2 <- filter_true_selex(selex_array_2, "res_2_group_2", T, T)
td_array_3 <- filter_true_selex(selex_array_3, "res_3_group_2", T)
fd_array_3 <- filter_true_selex(selex_array_3, "res_3_group_2", T, T)

viz_sa_1 <- viz_fitter(td_array_1)
viz_sa_2 <- viz_fitter(td_array_2)
viz_sa_3 <- viz_fitter(td_array_3)

detex_plot(viz_sa_1, "Resolution 1 Detections")
detex_plot(viz_sa_1, "Resolution 1 Detections", c("p_value", "lr_mean_3", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_2, "Resolution 2 Detections")
detex_plot(viz_sa_2, "Resolution 2 Detections", c("p_value", "lr_mean_2", "cal_kappa_1", "lr_prior_2", "cal_mix"))
detex_plot(viz_sa_3, "Resolution 3 Detections")
detex_plot(viz_sa_3, "Resolution 3 Detections", c("p_value", "lr_mean_1", "cal_kappa_1", "lr_prior_2", "cal_mix"))

# False detex
fdr_mat_1 <- fdr_fitter(selex_array_1, fd_array_1)
fdr_mat_2 <- fdr_fitter(selex_array_2, fd_array_2)
fdr_mat_3 <- fdr_fitter(selex_array_3, fd_array_3)

fdr_tbl <- rbind(as.numeric(fdr_mat_1), as.numeric(fdr_mat_2), as.numeric(fdr_mat_3))
rownames <- c("P. Res. 1", "P. Res. 2", "P. Res. 3")
fdr_tbl <- data.frame(cbind(rownames, fdr_tbl))
gt(fdr_tbl) |> fmt_number(decimals = 6)

# Get true detections for methods.

td_mthd_2 <- filter_true_selex(method_array_2, semi_true = T)
td_mthd_5 <- filter_true_selex(method_array_5, semi_true = T)
td_mthd_6 <- filter_true_selex(method_array_6, semi_true = T)
td_mthd_7 <- filter_true_selex(method_array_7, semi_true = T)
td_mthd_8 <- filter_true_selex(method_array_8, semi_true = T)
td_mthd_9 <- filter_true_selex(method_array_9, semi_true = T)
td_mthd_10 <- filter_true_selex(method_array_10, semi_true = T)

res_mat_2_3 <- res_fitter(td_mthd_2)
res_mat_5_3 <- res_fitter(td_mthd_5)
res_mat_6_3 <- res_fitter(td_mthd_6)
res_mat_7_3 <- res_fitter(td_mthd_7)
res_mat_8_3 <- res_fitter(td_mthd_8)
res_mat_9_3 <- res_fitter(td_mthd_9)
res_mat_10_3 <- res_fitter(td_mthd_10)
res_mat_2_2 <- res_fitter(td_mthd_2, 2)
res_mat_5_2 <- res_fitter(td_mthd_5, 2)
res_mat_6_2 <- res_fitter(td_mthd_6, 2)
res_mat_7_2 <- res_fitter(td_mthd_7, 2)
res_mat_8_2 <- res_fitter(td_mthd_8, 2)
res_mat_9_2 <- res_fitter(td_mthd_9, 2)
res_mat_10_2 <- res_fitter(td_mthd_10, 2)
res_mat_2_1 <- res_fitter(td_mthd_2, 1)
res_mat_5_1 <- res_fitter(td_mthd_5, 1)
res_mat_6_1 <- res_fitter(td_mthd_6, 1)
res_mat_7_1 <- res_fitter(td_mthd_7, 1)
res_mat_8_1 <- res_fitter(td_mthd_8, 1)
res_mat_9_1 <- res_fitter(td_mthd_9, 1)
res_mat_10_1 <- res_fitter(td_mthd_10, 1)

res_plot(res_mat_2_1, "Calibrator, k = 0.25, at Res 1")
res_plot(res_mat_2_2, "Calibrator, k = 0.25, at Res 2")
res_plot(res_mat_2_3, "Calibrator, k = 0.25, at Res 3")

res_plot(res_mat_5_1, "Calibrator, mixture, at Res 1")
res_plot(res_mat_5_2, "Calibrator, mixture, at Res 2")
res_plot(res_mat_5_3, "Calibrator, mixture, at Res 3")

res_plot(res_mat_6_1, "Likelihood Ratio, mean = 2.5, at Res 1")
res_plot(res_mat_6_2, "Likelihood Ratio, mean = 2.5, at Res 2")
res_plot(res_mat_6_3, "Likelihood Ratio, mean = 2.5, at Res 3")

res_plot(res_mat_7_1, "Likelihood Ratio, mean = 5, at Res 1")
res_plot(res_mat_7_2, "Likelihood Ratio, mean = 5, at Res 2")
res_plot(res_mat_7_3, "Likelihood Ratio, mean = 5, at Res 3")

res_plot(res_mat_8_1, "Likelihood Ratio, mean = 7.5, at Res 1")
res_plot(res_mat_8_2, "Likelihood Ratio, mean = 7.5, at Res 2")
res_plot(res_mat_8_3, "Likelihood Ratio, mean = 7.5, at Res 3")

res_plot(res_mat_9_1, "LR Mixture, prior sigma = 5, at Res 1")
res_plot(res_mat_9_2, "LR Mixture, prior sigma = 5, at Res 2")
res_plot(res_mat_9_3, "LR Mixture, prior sigma = 5, at Res 3")

res_plot(res_mat_10_1, "LR Mixture, prior sigma = 20, at Res 1")
res_plot(res_mat_10_2, "LR Mixture, prior sigma = 20, at Res 2")
res_plot(res_mat_10_3, "LR Mixture, prior sigma = 20, at Res 3")

viz_ma_2 <- viz_fitter(method_array_2)
viz_ma_5 <- viz_fitter(method_array_5)
viz_ma_6 <- viz_fitter(method_array_6)
viz_ma_7 <- viz_fitter(method_array_7)
viz_ma_8 <- viz_fitter(method_array_8)
viz_ma_9 <- viz_fitter(method_array_9)
viz_ma_10 <- viz_fitter(method_array_10)

res_plot(viz_ma_2, "Calibrator, k = 0.25", position_identity())
res_plot(viz_ma_5, "Calibrator, mixture", position_identity())
res_plot(viz_ma_6, "Likelihood Ratio, mean = 2.5", position_identity())
res_plot(viz_ma_7, "Likelihood Ratio, mean = 5", position_identity())
res_plot(viz_ma_8, "Likelihood Ratio, mean = 7.5", position_identity())
res_plot(viz_ma_9, "LR Mixture, prior sigma = 5", position_identity())
res_plot(viz_ma_10, "LR Mixture, prior sigma = 20", position_identity())

cvg_ma_3 <- cvg_fitter(td_array_3, "res_3_group_2")
cvg_plot(cvg_ma_3, "Coverage Proportion by Signal Size at Resolution 3", method = c("p_value", "cal_mix", "lr_prior_1", "lr_mean_1", "lr_mean_2", "lr_mean_3", "lr_prior_2"))
