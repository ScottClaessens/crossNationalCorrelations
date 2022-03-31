# custom functions - simulation

# load isocode data
loadISO <- function(fileISO) {
  read.csv(fileISO) %>%
    mutate(iso3 = str_remove(Alpha.3.code, " "),
           iso2 = str_remove(Alpha.2.code, " ")) %>%
    distinct(iso2, iso3, Latitude..average., Longitude..average.)
}

# load geographic and linguistic covariance matrices
loadSimCovMat <- function(file, continent, iso, langFam) {
  # load distance matrix
  out <- 
    read_excel(file, na = "") %>%
    dplyr::select(-ISO) %>%
    as.matrix()
  rownames(out) <- colnames(out)
  # keep only countries with longitude, latitude, continent, and language family data
  out <- out[rownames(out)[rownames(out) %in% iso$iso2],
             rownames(out)[rownames(out) %in% iso$iso2]]
  out <- out[rownames(out)[rownames(out) %in% continent$country],
             rownames(out)[rownames(out) %in% continent$country]]
  out <- out[rownames(out)[rownames(out) %in% langFam$iso],
             rownames(out)[rownames(out) %in% langFam$iso]]
  # distances between 0 and 1
  out <- out / max(out)
  diag(out) <- 0
  # 1 - distance = proximity (covariance)
  out <- 1 - out
  # ensure covariance matrix is positive definite by adding small constant to diagonal
  diag(out) <- diag(out) + 1e-06
  return(out)
}

# brms model to simulate data
fitSimulationModel <- function(covMat, lambda, rho, iter) {
  # signal for dependent variable = lambda
  # signal for independent variable = rho
  signalY <- lambda
  signalX <- rho
  # calculate sd and sigma values, assuming total variance sd^2 + sigma^2 = 1
  # signal = sd^2 / (sd^2 + sigma^2)
  # signal = sd^2 / 1
  # sd = sqrt(signal)
  sdY <- sqrt(signalY)
  sdX <- sqrt(signalX)
  # sd^2 + sigma^2 = 1
  # sigma = sqrt(1 - sd^2)
  sigmaY <- sqrt(1 - sdY^2)
  sigmaX <- sqrt(1 - sdX^2)
  # mock data for simulation (just to feed to stan so model can run)
  d <- data.frame(y = rnorm(nrow(covMat)), x = rnorm(nrow(covMat)), id = rownames(covMat))
  # simulation
  mSim <- brm(bf(y ~ 0 + (1 | gr(id, cov = covMat))) +
              bf(x ~ 0 + (1 | gr(id, cov = covMat))) + set_rescor(FALSE),
              data = d, data2 = list(covMat = covMat),
              prior = c(prior_string(paste0("constant(", sdY, ")"), class = "sd", resp = "y"),
                        prior_string(paste0("constant(", sigmaY, ")"), class = "sigma", resp = "y"),
                        prior_string(paste0("constant(", sdX, ")"), class = "sd", resp = "x"),
                        prior_string(paste0("constant(", sigmaX, ")"), class = "sigma", resp = "x")),
              sample_prior = "only", warmup = 0, iter = 1, chains = 1)
  return(mSim)
}

# simulate data
simulateData <- function(simModel, covMat, continent, iso, langFam, iter, rho, lambda) {
  # fit model again with seed
  mSim <- update(simModel, sample_prior = "only", warmup = 0, iter = 1, chains = 1, seed = iter)
  # get simulated data
  pred <- posterior_predict(mSim)
  # standardise simulated x and y
  y <- as.numeric(scale(pred[1,,"y"]))
  x <- as.numeric(scale(pred[1,,"x"]))
  # produce data frame and match longitude, latitude, continent, and language family data
  out <- 
    tibble(x = x, y = y, country = rownames(covMat), iter = iter) %>%
    left_join(iso, by = c("country" = "iso2")) %>%
    left_join(distinct(continent, country, .keep_all = TRUE), by = "country") %>%
    left_join(dplyr::select(langFam, iso, langFamily), by = c("country" = "iso")) %>%
    mutate(continent = factor(continent),
           langFamily = factor(langFamily),
           rho = rho, lambda = lambda) %>%
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
    Q97.5 = confint(m)["x","97.5 %"],
    sig = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0),
    seed = unique(data$seed),
    rho = unique(data$rho),
    lambda = unique(data$lambda)
  )
}

# initialise brms model
setupBrms <- function(data, covMat, type = "") {
  # formula
  if (type == "spatial")    f <- bf(y ~ 0 + Intercept + x + gp(latitude, longitude, k = 5, c = 5/4))
  if (type == "linguistic") f <- bf(y ~ 0 + Intercept + x + (1 | gr(country, cov = covMat)))
  if (type == "both")       f <- bf(y ~ 0 + Intercept + x + gp(latitude, longitude, k = 5, c = 5/4) +
                                      (1 | gr(country, cov = covMat)))
  # priors
  priors <- c(prior(normal(0, 0.5), class = b),
              prior(exponential(5), class = sigma))
  if (type == "spatial")    priors <- c(priors, prior(exponential(5), class = sdgp))
  if (type == "linguistic") priors <- c(priors, prior(exponential(5), class = sd))
  if (type == "both")       priors <- c(priors, prior(exponential(5), class = sdgp),
                                        prior(exponential(5), class = sd))
  # fit model
  brm(f, data = data, data2 = list(covMat = covMat), prior = priors, chains = 0)
}

# fit brms model
fitBrmsModel <- function(brmsInitial, data) {
  # fit model
  m <- update(brmsInitial, newdata = data, chains = 4,
              control = list(adapt_delta = 0.999, max_treedepth = 15),
              cores = 4, seed = 2113, iter = 2000)
  # get coefficient on X and 95% CIs
  tibble(
    bX = fixef(m)["x","Estimate"],
    Q2.5 = fixef(m)["x","Q2.5"],
    Q97.5 = fixef(m)["x","Q97.5"],
    sig = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0),
    rhat = rhat(m)["b_x"],
    divergences = sum(rstan::get_divergent_iterations(m$fit)),
    seed = unique(data$seed),
    rho = unique(data$rho),
    lambda = unique(data$lambda)
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
    Q97.5 = confint(m)["x","97.5 %"],
    sig = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0),
    seed = unique(data$seed),
    rho = unique(data$rho),
    lambda = unique(data$lambda)
  )
}

# plot individual simulation results at strong autocorrelation levels
plotSimInd <- function(olsModel1, olsModel2, olsModel3, olsModel4, olsModel5,
                       conleyModel1, conleyModel2, conleyModel3, brmsModel1,
                       brmsModel2, brmsModel3, type, file) {
  # plotting function
  plotFun <- function(results, title, ylab) {
    results %>%
      # strong autocorrelation levels
      filter(rho == 0.8 & lambda == 0.8) %>%
      # order by slope
      mutate(id = factor(1:nrow(.))) %>%
      ggplot(aes(x = fct_reorder(id, bX, .desc = TRUE), 
                 y = bX, ymin = Q2.5, ymax = Q97.5, colour = sig)) +
      # above or below is false positive
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_pointrange(size = 0.03) +
      scale_y_continuous(name = ylab, limits = c(-1, 1)) +
      # red = sig, black = not sig
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
  pA <- plotFun(olsModel1, "No control", "")
  pB <- plotFun(olsModel2, "Latitude", "")
  pC <- plotFun(olsModel3, "Longitude", "")
  pD <- plotFun(olsModel4, "Continent", "")
  pE <- plotFun(olsModel5, "Language family", "Correlation")
  pF <- plotFun(conleyModel1, "Conley SEs 100km", "")
  pG <- plotFun(conleyModel2, "Conley SEs 1000km", "")
  pH <- plotFun(conleyModel3, "Conley SEs 10000km", "")
  pI <- plotFun(brmsModel1, "Bayesian spatial", "")
  pJ <- plotFun(brmsModel2, "Bayesian linguistic", "")
  pK <- plotFun(brmsModel3, "Bayesian spatial and linguistic", "")
  # put together
  out <-
    plot_grid(
      plot_grid(ggdraw() + draw_label(paste0("Simulation with strong ", type, " non-independence"), 
                                      x = ifelse(type == "spatial", 0.295, 0.355), fontface = "bold")),
      plot_grid(pA, pB, pC, pD, nrow = 1),
      plot_grid(pE, pF, pG, pH, nrow = 1),
      plot_grid(NULL, pI, pJ, pK, NULL, nrow = 1, 
                rel_widths = c(0.5, 1, 1, 1, 0.5)),
      nrow = 4, rel_heights = c(0.3, 1, 1, 1)
    )
  # save
  ggsave(out, filename = file, width = 9.5, height = 5.5)
  return(out)
}

# plot simulation results across all autocorrelation levels
plotSimAll <- function(olsModel1, olsModel2, olsModel3, olsModel4, olsModel5,
                       conleyModel1, conleyModel2, conleyModel3, brmsModel1,
                       brmsModel2, brmsModel3, type, file) {
  # plotting function
  plotFun <- function(data, title, ylab, xlab) {
    data %>%
      group_by(lambda, rho) %>%
      summarise(prop = sum(sig) / n(), n = n()) %>%
      mutate(lower = map2(prop, n, function(x, y) getBootCIProp(x, y, 1000)[1]), 
             upper = map2(prop, n, function(x, y) getBootCIProp(x, y, 1000)[2])) %>%
      unnest(c(lower, upper)) %>%
      ggplot(aes(x = lambda, y = prop, group = factor(rho), colour = factor(rho))) +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      geom_pointrange(aes(x = lambda, y = prop, ymin = lower, ymax = upper), 
                      position = position_dodge(0.05), size = 0.15) +
      geom_line(position = position_dodge(0.05)) +
      scale_x_continuous(name = xlab, limits = c(0.1, 0.9), breaks = c(0.2, 0.5, 0.8)) +
      scale_y_continuous(name = ylab, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      guides(colour = guide_legend(title = "Strength of\nautocorrelation\nfor predictor\nvariable")) +
      ggtitle(title) +
      theme_classic()
  }
  # individual plots
  pA <- plotFun(olsModel1, "No control", "", "")
  pB <- plotFun(olsModel2, "Latitude", "", "")
  pC <- plotFun(olsModel3, "Longitude", "", "")
  pD <- plotFun(olsModel4, "Continent", "", "")
  pE <- plotFun(olsModel5, "Language family", "False positive rate", "")
  pF <- plotFun(conleyModel1, "Conley SEs 100km", "", "")
  pG <- plotFun(conleyModel2, "Conley SEs 1000km", "", "")
  pH <- plotFun(conleyModel3, "Conley SEs 10000km", "", "")
  pI <- plotFun(brmsModel1, "Bayesian spatial", "", " \n ")
  pJ <- plotFun(brmsModel2, "Bayesian linguistic", "", "Strength of autocorrelation\nfor outcome variable")
  pK <- plotFun(brmsModel3, "Bayesian spatial and linguistic", "", " \n ")
  # put together
  out <-
    plot_grid(
      plot_grid(ggdraw() + draw_label(paste0("Simulation with ", type, " non-independence"), 
                                      x = ifelse(type == "spatial", 0.27, 0.35), fontface = "bold")),
      plot_grid(
        pA + theme(legend.position = "none"),
        pB + theme(legend.position = "none"),
        pC + theme(legend.position = "none"),
        pD + theme(legend.position = "none"),
        nrow = 1
      ),
      plot_grid(
        pE + theme(legend.position = "none"),
        pF + theme(legend.position = "none"),
        pG + theme(legend.position = "none"),
        pH + theme(legend.position = "none"),
        nrow = 1
      ),
      plot_grid(
        NULL, pI + theme(legend.position = "none"),
        pJ + theme(legend.position = "none"),
        pK + theme(legend.position = "none"), NULL,
        nrow = 1, rel_widths = c(0.5, 1, 1, 1, 0.5)
      ),
      nrow = 4, rel_heights = c(0.3, 1, 1, 1.08)
    )
  # add legend
  out <- plot_grid(out, NULL, get_legend(pA), nrow = 1, rel_widths = c(1, 0.02, 0.15))
  # save
  ggsave(out, filename = file, width = 11, height = 7)
  return(out)
}