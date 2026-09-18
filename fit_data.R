# Read in real data in MatLab format.

# Get ordinary libraries
library(R.matlab)
library(abind)
library(stringr)
source("e_procedures.R")

read_subjects <- function(block, type, lbls) {
  for (sub in list.files("./real_data/neurocon")) {
    filename <- paste0("real_data/neurocon/", sub, "/", sub, "_", type, "_correlation_matrix.mat")
    mat <- readMat(filename)$data
    colnames(mat) <- lbls
    rownames(mat) <- lbls
    
    if (substr(sub, 5, 11) == block) {
      if (exists("subject_mats")) {
        # Ensure it becomes an array
        subject_mats <- abind(subject_mats, mat, along = 3)
      } else {
        subject_mats <- mat
      }
    }
  }
  return(subject_mats)
}

library(brainGraph)

# Node groups-- LR: Hemispheres; MLR: Lobes; MHR: Hemisphere-Lobes; HR: Individual regions 
# Edge groups-- LRxLR, MLRxMLR, MHRxMHR, HRxHR
# We require (for fitting).
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# group_info: data frame containing group, level, and group-level for each group.
# group_info_mat: array of group memberships of edges by level
# group_subgroups: list linking groups to their subgroups
labs <- brainGraph::aal116

# Get node groups

node_groups <- matrix(NA, nrow = 116, ncol = 4)
for (i in 1:116) {
  node_groups[i,1:2] <- c(which(unique(labs[,6])[[1]] %in% labs[i,6][[1]]),
                          which(unique(labs[,5])[[1]] %in% labs[i,5][[1]]))
}
col_3_raw <- paste(node_groups[,1], node_groups[,2], sep = "")
for (i in 1:116) {
  node_groups[i,3] <- which(unique(col_3_raw) == col_3_raw[i])
  node_groups[i,4] <- i
}

# Get node groups (set topology)

node_groups <- list()

# Level 1
node_groups[[1]] <- list()
for (lab in unique(labs[,6])[[1]]) {
  node_groups[[1]][[lab]] <- which(labs[,6] == lab)
}

# Level 2
node_groups[[2]] <- list()
for (lab in unique(labs[,5])[[1]]) {
  node_groups[[2]][[lab]] <- which(labs[,5] == lab)
}

# Level 3
node_groups[[3]] <- list()
for (lab1 in unique(labs[,6])[[1]]) {
  for (lab2 in unique(labs[,5])[[1]]) {
    node_groups[[3]][[paste0(lab2, ".", lab1, sep = "")]] <- which(labs[,6] == lab1 & labs[,5] == lab2)
  }
}

node_groups[[4]] <- list()
for (i in 1:116) {
  node_groups[[4]][[labs[i,1][[1]]]] <- i
}

# Get edge groups (cartesian product topology) --> plan a function.
# yield the following.


# Generate (edge) groups from node groups, with group names: res_(l)_group_(g1-by-g2)
# INPUT:
# L: number of resolutions
# N: number of nodes
# node_groups: 2D list mapping, for each resolution, nodes to each node group
# OUTPUT: 
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# group_overlaps: list of overlaps by group. (res_group notation)
# group_memberships: list of edges by group. (res_group notation)
# group_info: data frame containing group, level, and group-level for each group. (res_group notation)
# eg_to_ng: list mapping each (edge) group to its composite node groups.
groups_from_ng <- function(L, N, node_groups) {
  # Remove lower triangle; these will be empty edges.
  ref_mat <- matrix(1:N^2, nrow = N, ncol = N)
  ref_mat[!upper.tri(ref_mat)] <- NA
  
  groups <- matrix(NA, nrow = N^2, ncol = L)
  group_memberships <- list()
  eg_to_ng <- list()
  for (idx in 1:N^2) {
    cds <- which(ref_mat == idx, arr.ind = TRUE)
    if (length(cds) > 0) {
      i <- cds[1]
      j <- cds[2]
      for (l in 1:L) {
        il <- names(which(sapply(node_groups[[l]], function(y) i %in% y)))
        jl <- names(which(sapply(node_groups[[l]], function(y) j %in% y)))
        
        group <- sort(c(il, jl))
        group_lab <- paste0(group[1], "-", group[2], sep = "")
        
        groups[idx,l] <- group_lab
        group_memberships[[group_lab]] <- append(group_memberships[[group_lab]], idx)
        
        if (is.null(eg_to_ng[[group_lab]])) {
          eg_to_ng[[group_lab]] <- group
        }
      }
    }
  }
  
  group_info <- c()
  
  for (l in 1:L) {
    df_temp <- data.frame("Group_Lab" = unique(na.omit(groups[,l])))
    df_temp$Resolution <- l
    df_temp$group <- paste0("group_", df_temp$Group_Lab)
    df_temp$res_Group <- paste0("res_",df_temp$Resolution, "_", df_temp$group)
    
    group_info <- rbind(group_info, df_temp)
  }
  
  # This one takes a while. But you only have to run it once.
  
  group_overlaps <- list()
  
  # Iterate over groups s_group for selected and look at other groups at lower resolutions
  for (i in 1:nrow(group_info)) {
    s_inf <- group_info[i,]
    
    s_group <- s_inf[1][[1]]
    
    s_res <- s_inf[2][[1]]
    
    # For each group of lower resolution c_group for comparator, s_group_1 s_group_2 p_group_1 p_group_2, get node groups
    p_infs <- group_info[group_info$Resolution > s_res,]
    s_overlaps <- s_group
    
    if (dim(p_infs)[1] > 0) {
      for (j in 1:nrow(p_infs)) {
        p_inf <- p_infs[j,]
        
        p_group <- p_inf[1][[1]]
        
        p_res <- p_inf[2][[1]]
        
        # Add p_group to group_overlaps$s_group
        # If node group on s_group_1 and p_group_1 have nonzero overlap AND same on s_group_2 and p_group_2 have nonzero overlap
        if (length(intersect(group_memberships[[s_group]], group_memberships[[p_group]])) > 0) {
          s_overlaps <- append(s_overlaps, p_group)
        }
      }
    }
    
    group_overlaps[[s_group]] <- s_overlaps
  }
  
  return(list(groups, group_overlaps, group_memberships, group_info, eg_to_ng))
}

# Get test statistics for a given edge set
# INPUT:
# M: number of subjects per matrix
# control_mats: array of matrices for the controls
# treatment_mats: array of matrices for the treatments
# edges: vector of edges in the group
# OUTPUT: 
# group_vec: vector of test statistics (p_val, cal_kappa(s), cal_mix, lr_mean(s), lr_prior(s))

stats_vec <- function(M, control_mats, treatment_mats, edges) {
  control_sm <- apply(control_mats, c(1,2), mean)
  treatment_sm <- apply(treatment_mats, c(1,2), mean)
  
  n.edges <- M * length(edges)
  
  control_bar <- mean(control_sm[edges])
  treatment_bar <- mean(treatment_sm[edges])
  
  d_bar <- control_bar - treatment_bar
  
  p_val <- p_value(d_bar, n.edges)
  
  kappas <- c(0.25, 0.5, 0.75)
  cal_kappa_1 <- cal_kappa(p_val, kappas[1])
  cal_kappa_2 <- cal_kappa(p_val, kappas[2])
  cal_kappa_3 <- cal_kappa(p_val, kappas[3])
  
  cal_mix <- cal_mixture(p_val)
  
  pts <- c(2.5, 5, 7.5)
  lr_mean_1 <- lr_delta(d_bar, n.edges, pts[1])
  lr_mean_2 <- lr_delta(d_bar, n.edges, pts[2])
  lr_mean_3 <- lr_delta(d_bar, n.edges, pts[3])
  
  priors <- c(5, 20, 35)
  lr_prior_1 <- lr_prior(d_bar, n.edges, priors[1])
  lr_prior_2 <- lr_prior(d_bar, n.edges, priors[2])
  lr_prior_3 <- lr_prior(d_bar, n.edges, priors[3])
  
  return(c(p_val, cal_kappa_1, cal_kappa_2, cal_kappa_3, cal_mix, lr_mean_1, 
           lr_mean_2, lr_mean_3, lr_prior_1, lr_prior_2, lr_prior_3))
}

fit_groups <- groups_from_ng(4, 116, node_groups)
control_mats <- read_subjects("control", "AAL116", labs[,1][[1]])
treatment_mats <- read_subjects("patient", "AAL116", labs[,1][[1]])

LAB <- "Frontal.R-Parietal.L"
RES <- 3
M <- dim(control_mats)[3]

edges <- which(fit_groups[[1]][,3] == LAB)
stats_vec(M,control_mats,treatment_mats,edges)

# Obvious problem! We need a t-test version.

# From here, we join the stats vecs into a matrix and use that for the e_procedure. All code should be relatively reusable.