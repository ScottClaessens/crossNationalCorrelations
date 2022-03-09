# plot showing spatial autocorrelation example on map
plotSpatialExample <- function(data) {
  # sample countries
  d <-
    data %>%
    filter(country %in% c("GB","ES","RU","IT","FR",
                          "DE","CA","US","MX","BR",
                          "AR","CO","PE","CN","IN",
                          "UA","JP","DZ","NE","LY",
                          "AU","NZ","TH","PH","VN")) %>%
    mutate(Residual = residuals(lm(y ~ x, .)),
           line = y - Residual)
  # top plot
  pA <-
    ggplot(d) +
    geom_segment(aes(x = x, xend = x, y = y, yend = line, colour = Residual)) +
    geom_point(aes(x = x, y = y)) +
    geom_text_repel(aes(x = x, y = y, label = country), 
                    segment.colour = NA, size = 3.5, force = 1.5,
                    nudge_x = ifelse(d$country == "AU", -0.2, ifelse(d$country %in% c("GB","RU"), 0.1, 0)),
                    nudge_y = ifelse(d$country == "ES", 0.2, ifelse(d$country == "JP", -0.1, 0))) +
    geom_smooth(aes(x = x, y = y), method = "lm", 
                formula = y ~ x, se = FALSE, colour = "black") +
    scale_colour_gradient2(high = "red", low = "blue") +
    labs(y = "National-level outcome variable",
         x = "National-level predictor variable") +
    theme_classic() +
    theme(legend.position = "none")
  # bottom plot
  d$region <- countrycode(d$country, "iso2c", "country.name")
  d$region <- ifelse(d$region == "United States", "USA", d$region)
  d$region <- ifelse(d$region == "United Kingdom", "UK", d$region)
  world <- 
    map_data("world") %>% 
    filter(region != "Antarctica") %>%
    left_join(select(d, region, Residual), by = "region")
  pB <-
    ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      fill = "lightgrey"
    ) +
    geom_map(
      data = world[!is.na(world$Residual),], 
      map = world[!is.na(world$Residual),],
      aes(long, lat, map_id = region, fill = Residual)
    ) +
    scale_fill_gradient2(high = "red", low = "blue") +
    theme_void()
  # put together
  out <- plot_grid(pA, plot_grid(pB, NULL, rel_widths = c(1, 0.02)), ncol = 1)
  ggsave(out, filename = "figures/example.pdf", height = 6, width = 6.5)
  return(out)
}

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