# causal model of spatial and cultural phylogenetic non-independence
plotDAG <- function() {
  out <-
    dagify(Y ~ X + U,
           X ~ U,
           U ~ L + G,
           coords = tibble(name = c("X", "Y", "U", "G", "L"),
                           x = c(1, 2, 1.5, 1, 2),
                           y = c(1, 1, 0, 0, 0))) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(data = . %>% filter(name %in% c("U")),
                   shape = 1, stroke = 2, color = "black") +
    geom_dag_text(color = "black", size = 10) +
    geom_dag_edges(edge_color = "black", edge_width = 2,
                   arrow_directed = grid::arrow(type = "closed")) +
    theme_dag_blank()
  # save
  ggsave(out, filename = "figures/dag.pdf", height = 3.5, width = 6)
  return(out)
}
