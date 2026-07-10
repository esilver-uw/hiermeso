library(DiagrammeR)
library(latex2exp)

graph <- create_graph(directed = FALSE) |>
  add_gnm_graph(8, 10, set_seed = 95, 
                node_aes = node_aes(fillcolor = "#EC008C"),
                edge_aes = edge_aes(penwidth = 2))
render_graph(graph)

adj_mat <- matrix(0, nrow = 8, ncol = 8)
for (i in 1:dim(graph$edges_df)[1]) {
  adj_mat[graph$edges_df[i,2], graph$edges_df[i,3]] = adj_mat[graph$edges_df[i,3], graph$edges_df[i,2]] = 1
}
