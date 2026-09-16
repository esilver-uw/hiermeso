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
# Group names: res_(l)_group_(g1-by-g2)
# groups: data frame of dimension n^2 x L specifying group membership for each possible edge.
# group_overlaps: list of overlaps by group. (res_group notation)
# group_memberships: list of edges by group. (res_group notation)
# group_info: data frame containing group, level, and group-level for each group. (res_group notation)

# Remove lower triangle; these will be empty edges.
ref_mat <- matrix(1:116^2, nrow = 116, ncol = 116)
ref_mat[!upper.tri(ref_mat)] <- NA

groups <- matrix(NA, nrow = 116^2, ncol = 4)
for (idx in 1:116^2) {
  cds <- which(ref_mat == idx, arr.ind = TRUE)
  if (length(cds) > 0) {
    i <- cds[1]
    j <- cds[2]
    for (l in 1:4) {
      il <- names(which(sapply(node_groups[[l]], function(y) i %in% y)))
      jl <- names(which(sapply(node_groups[[l]], function(y) j %in% y)))
      
      group <- sort(c(il, jl))
      
      groups[idx,l] <- paste0(group[1], "-", group[2], sep = "")
    }
  }
}

group_info <- c()

for (l in 1:4) {
  df_temp <- data.frame("Group_Number" = unique(na.omit(groups[,l])))
  df_temp$Resolution <- l
  df_temp$group <- paste0("group_", df_temp$Group_Number)
  df_temp$res_Group <- paste0("res_",df_temp$Resolution, "_", df_temp$group)
  
  group_info <- rbind(group_info, df_temp)
}

groups_trimmed <- groups[complete.cases(groups),]

# Somewhat extensive in terms of search space, but simplified by a few rules:
# 1. Each resolution forms a partition
# 2. Each combination of (L,R) x (L,R) overlaps with each combination of lobe x lobe (verify this fact)
# 3. In terms of (L,R) x lobe groups, an overlap occurs with 2 same (L,R) or 2 same lobe labels.
# 4. Each region accords to a particular (L,R) and lobe, and thus will overlap with 2 same (L,R) or 2 same lobe labels.
# Again, verify this, but I think this schema should work.

