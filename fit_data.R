# Read in real data in MatLab format.

# Get ordinary libraries
library(R.matlab)
library(abind)
library(stringr)

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

# Get edge groups

nodes_to_edges <- matrix(1:116^2, nrow = 116, ncol = 116)
groups_to_edges <- array(NA, dim = c(116, 116, 4))
groups <- matrix(NA, nrow = 116^2, 4)
for (idx in 1:116^2) {
  i <- which(nodes_to_edges == idx, arr.ind = TRUE)[1]
  j <- which(nodes_to_edges == idx, arr.ind = TRUE)[2]
  if (i <= j) {
    groups_ij <- data.frame(i = c(node_groups[i,1], node_groups[i,2], node_groups[i,3], node_groups[i,4]), 
                            j = c(node_groups[j,1], node_groups[j,2], node_groups[j,3], node_groups[j,4]))
    groups_to_edges[i,j,] <- paste(groups_ij[,1], groups_ij[,2], sep = "")
  }
}

# Construct group info

group_info_mat <- array(NA, dim = c(116, 116, 4))
groups_ct <- length(unique(na.omit(as.vector(groups_to_edges[,,4]))))
group_info <- matrix(NA, nrow = groups_ct, ncol = 4)
for (l in 1:4) {
  group_ns <- unique(na.omit(as.vector(groups_to_edges[,,l])))
  group_info[,l] <- append(1:length(group_ns), rep(NA, groups_ct - length(group_ns)))
  for (i in 1:116) {
    for (j in i:116) {
      group_info_mat[i,j,l] <- which(group_ns == groups_to_edges[i,j,l])
    }
  }
}

