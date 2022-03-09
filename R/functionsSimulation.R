# custom functions - simulation

# load isocode data
loadISO <- function(fileISO) {
  read.csv(fileISO) %>%
    mutate(iso3 = str_remove(Alpha.3.code, " "),
           iso2 = str_remove(Alpha.2.code, " ")) %>%
    distinct(iso2, iso3, Latitude..average., Longitude..average.)
}

# load geographic distance matrix
loadDistanceMatrix <- function(file, continent, iso, legal) {
  # load distance matrix
  out <- 
    read_excel(file, na = "") %>%
    select(-ISO) %>%
    as.matrix()
  # row and column names identical
  rownames(out) <- colnames(out)
  # log
  out <- log(out)
  # diagonal = 0
  diag(out) <- 0
  # keep only countries with longitude, latitude, continent, and legal origin data
  legalOriginISO <- iso$iso2[iso$iso3 %in% as.character(legal$countrycode)]
  out <- out[rownames(out)[rownames(out) %in% iso$iso2],
             rownames(out)[rownames(out) %in% iso$iso2]]
  out <- out[rownames(out)[rownames(out) %in% continent$country],
             rownames(out)[rownames(out) %in% continent$country]]
  out <- out[rownames(out)[rownames(out) %in% legalOriginISO],
             rownames(out)[rownames(out) %in% legalOriginISO]]
  # scale matrix
  out <- out / max(out)
  return(out)
}

# simulate data
simulateData <- function(distGeo, continent, iso, legal, seed, lambda, rho) {
  # set seed
  set.seed(seed)
  # parameters
  ### lambda = amount of autocorrelation for dependent variable
  ### rho = amount of autocorrelation for independent variable
  # simulate x and y, modified from https://www.jstor.org/stable/644404 (p. 761)
  n <- nrow(distGeo)
  x <- ((1 - lambda*distGeo)^-1) %*% rnorm(n)
  y <- ((1 -    rho*distGeo)^-1) %*% rnorm(n)
  # add some random noise
  x <- x + rnorm(n)
  y <- y + rnorm(n)
  # standardise x and y
  x <- as.numeric(scale(x))
  y <- as.numeric(scale(y))
  # produce data frame and match longitude, latitude, continent, and legal origin data
  out <- 
    tibble(x = x, y = y, country = rownames(distGeo), seed = seed) %>%
    left_join(iso, by = c("country" = "iso2")) %>%
    left_join(distinct(continent, country, .keep_all = TRUE), by = "country") %>%
    left_join(select(legal, countrycode, LO), by = c("iso3" = "countrycode")) %>%
    mutate(continent = factor(continent),
           LO = factor(get_labels(LO)[LO])) %>%
    rename(latitude = Latitude..average.,
           longitude = Longitude..average.)
  return(out)
}

# ols model
fitOLSModel <- function(formula, data) {
  # fit model
  m <- lm(formula, data = data)
  # get coefficient on x and 95% CIs
  tibble(
    bX = as.numeric(coef(m)["x"]),
    Q2.5 = confint(m)["x","2.5 %"],
    Q97.5 = confint(m)["x","97.5 %"]
  )
}

# initialise brms model
setupBrms <- function(data) {
  brm(y ~ 0 + Intercept + x + gp(latitude, longitude), data = data,
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(5), class = sdgp),
                prior(exponential(5), class = sigma)),
      chains = 0)
}

# fit brms model
fitBrmsModel <- function(brmsInitial, data) {
  # fit model
  m <- update(brmsInitial, newdata = data, chains = 4,
              control = list(adapt_delta = 0.90, max_treedepth = 15),
              cores = 4, seed = 2113, iter = 2000)
  # get coefficient on X and 95% CIs
  tibble(
    bX = fixef(m)["x","Estimate"],
    Q2.5 = fixef(m)["x","Q2.5"],
    Q97.5 = fixef(m)["x","Q97.5"],
    rhat = rhat(m)["b_x"],
    divergences = sum(rstan::get_divergent_iterations(m$fit))
  )
}

# fit conley se model
fitConleyModel <- function(data, dist_cutoff) {
  # fit model
  m <- conleyreg(y ~ x, data = data, dist_cutoff = dist_cutoff, 
                 lat = "latitude", lon = "longitude")
  # get coefficient on X and 95% CIs
  tibble(
    bX = m["x","Estimate"],
    Q2.5 = confint(m)["x","2.5 %"],
    Q97.5 = confint(m)["x","97.5 %"]
  )
}

# plot simulation results for medium autocorrelation levels
plotSimulation1 <- function(olsModel1, olsModel2, olsModel3, olsModel4, olsModel5, 
                            conleyModel1, conleyModel2, conleyModel3, brmsModel) {
  # plotting function
  plotFun <- function(results, title, ylab) {
    results %>%
      mutate(sig = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0),
             id = factor(1:nrow(.))) %>%
      ggplot(aes(x = fct_reorder(id, bX, .desc = TRUE), 
                 y = bX, ymin = Q2.5, ymax = Q97.5, colour = sig)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_pointrange(size = 0.03) +
      scale_y_continuous(name = ylab, limits = c(-1, 1)) +
      scale_colour_manual(values = c("#000000", "#FF0000")) +
      ggtitle(title) +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none",
            title = element_text(size = 8),
            axis.title.y = element_text(size = 9.5))
  }
  # plots
  pA <- plotFun(olsModel1, "Standard linear regression", "")
  pB <- plotFun(olsModel2, "Latitude control", "")
  pC <- plotFun(olsModel3, "Longitude control", "")
  pD <- plotFun(olsModel4, "Continent fixed effects", "Correlation")
  pE <- plotFun(olsModel5, "Legal origin fixed effects", "")
  pF <- plotFun(conleyModel1, "Conley SEs (100km)", "")
  pG <- plotFun(conleyModel2, "Conley SEs (1000km)", "")
  pH <- plotFun(conleyModel3, "Conley SEs (10000km)", "")
  pI <- plotFun(brmsModel, "Bayesian w/ Gaussian Process", "")
  # put together
  out <-
    plot_grid(
      plot_grid(pA, pB, pC, nrow = 1),
      plot_grid(pD, pE, pF, nrow = 1),
      plot_grid(pG, pH, pI, nrow = 1),
      nrow = 3
    )
  # save
  ggsave(out, filename = "figures/simulation1.pdf", width = 7.5, height = 5)
  return(out)
}

# combine simulation results
combineResults <- function(results_0.1_0.1, results_0.1_0.5, results_0.1_0.9,
                           results_0.5_0.1, results_0.5_0.5, results_0.5_0.9,
                           results_0.9_0.1, results_0.9_0.5, results_0.9_0.9) {
  # how many simulated datasets?
  n <- nrow(results_0.1_0.1)
  # combine results from different levels of autocorrelation
  rbind(
    results_0.1_0.1, results_0.1_0.5, results_0.1_0.9,
    results_0.5_0.1, results_0.5_0.5, results_0.5_0.9,
    results_0.9_0.1, results_0.9_0.5, results_0.9_0.9
  ) %>%
    # add variables
    mutate(id = rep(1:n, times = 9),
           lambda = rep(c(0.1, 0.5, 0.9), each = n*3),
           rho = rep(rep(c(0.1, 0.5, 0.9), each = n), times = 3),
           sig = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0))
}

# plot simulation results across all autocorrelation levels
plotSimulation2 <- function(olsModel1, olsModel2, olsModel3, olsModel4, olsModel5, 
                            conleyModel1, conleyModel2, conleyModel3, brmsModel) {
  # plotting function
  plotFun <- function(data, title, ylab, xlab) {
    data %>%
      group_by(lambda, rho) %>%
      summarise(prop = sum(sig) / n()) %>%
      mutate(lower = map(prop, function(x) getBootCIProp(x, 100, 1000)[1]), 
             upper = map(prop, function(x) getBootCIProp(x, 100, 1000)[2])) %>%
      unnest(c(lower, upper)) %>%
      ggplot(aes(x = lambda, y = prop, group = factor(rho), colour = factor(rho))) +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      geom_pointrange(aes(x = lambda, y = prop, ymin = lower, ymax = upper), 
                      position = position_dodge(0.05), size = 0.15) +
      geom_line(position = position_dodge(0.05)) +
      scale_x_continuous(name = xlab, limits = c(0, 1), breaks = c(0.1, 0.5, 0.9)) +
      scale_y_continuous(name = ylab, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      guides(colour = guide_legend(title = "Strength of spatial\nautocorrelation for\npredictor variable")) +
      ggtitle(title) +
      theme_classic()
  }
  # individual plots
  pA <- plotFun(olsModel1, "Standard linear regression", "", "")
  pB <- plotFun(olsModel2, "Latitude control", "", "")
  pC <- plotFun(olsModel3, "Longitude control", "", "")
  pD <- plotFun(olsModel4, "Continent fixed effects", "False positive rate", "")
  pE <- plotFun(olsModel5, "Legal origin fixed effects", "", "")
  pF <- plotFun(conleyModel1, "Conley SEs (100km)", "", "")
  pG <- plotFun(conleyModel2, "Conley SEs (1000km)", "", " \n ")
  pH <- plotFun(conleyModel3, "Conley SEs (10000km)", "", "Strength of spatial autocorrelation\nfor outcome variable")
  pI <- plotFun(brmsModel, "Bayesian w/ Gaussian Process", "", " \n ")
  # put together
  out <-
    plot_grid(
      pA + theme(legend.position = "none"),
      pB + theme(legend.position = "none"),
      pC + theme(legend.position = "none"),
      pD + theme(legend.position = "none"),
      pE + theme(legend.position = "none"),
      pF + theme(legend.position = "none"),
      pG + theme(legend.position = "none"),
      pH + theme(legend.position = "none"),
      pI + theme(legend.position = "none"),
      nrow = 3, rel_heights = c(1, 1, 1.1)
    )
  # add legend
  out <- plot_grid(out, get_legend(pA), nrow = 1, rel_widths = c(1, 0.2))
  # save
  ggsave(out, filename = "figures/simulation2.pdf", width = 10, height = 7)
  return(out)
}