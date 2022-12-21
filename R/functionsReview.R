# custom functions - review

# load main review
loadReview <- function(fileReview, fileMethods) {
  read_excel(fileReview, sheet = "Review") %>%
    left_join(read_csv(fileMethods), by = "MethodNI")
}

# fit random effects model to retrieve corrected proportions
fitReviewMLM <- function(review, outcome) {
  # fill NAs (analyses that did not control for NI at all) with zeroes
  review <- 
    review %>%
    mutate_at(vars(starts_with("method_")), function(x) ifelse(is.na(x), 0, x))
  # fit model
  brm(
    formula = bf(paste0(outcome, " ~ 0 + Review + (1 | Review:Citation)")),
    data = review,
    family = bernoulli,
    prior = c(prior(normal(0, 1), class = b),
              prior(exponential(1), class = sd)),
    cores = 4, seed = 2113
  )
}

# analyis-level year spline model
fitYearAnalysis <- function(review) {
  # prepare data
  review <-
    review %>%
    mutate(
      # standardise year
      Year.std = as.numeric(scale(Year)),
      # as factor
      Review = factor(Review)
    )
  # fit spline model
  brm(ControlNI ~ 0 + Intercept + s(Year.std, by = Review) + (1 | Review:Citation),
      data = review, family = bernoulli,
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(2), class = sd),
                prior(exponential(2), class = sds)),
      cores = 4, seed = 2113, iter = 3000,
      control = list(adapt_delta = 0.99))
}

# analysis-level impact factor model
fitIFAnalysis <- function(review) {
  # prepare data
  review <- mutate(review, Review = factor(Review))
  # fit model
  brm(ControlNI ~ 0 + Review + Review:log(IF) + (1 | Review:Citation),
      data = review, family = bernoulli,
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(2), class = sd)),
      cores = 4, seed = 2113, iter = 3000,
      control = list(adapt_delta = 0.99))
}

# article-level year spline model
fitYearArticle <- function(review) {
  # prepare data
  review <-
    review %>%
    group_by(Review, Citation) %>%
    summarise(
      ControlNI = ifelse(sum(ControlNI == "Yes") > 0, 1, 0),
      Year = mean(Year)
      ) %>%
    mutate(
      # standardise year
      Year.std = as.numeric(scale(Year)),
      # as factor
      Review = factor(Review)
    )
  # fit spline model
  brm(ControlNI ~ 0 + Intercept + s(Year.std, by = Review),
      data = review, family = bernoulli,
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(2), class = sds)),
      cores = 4, seed = 2113, iter = 3000,
      control = list(adapt_delta = 0.99))
}

# article-level impact factor model
fitIFArticle <- function(review) {
  # prepare data
  review <-
    review %>%
    group_by(Review, Citation) %>%
    summarise(
      ControlNI = ifelse(sum(ControlNI == "Yes") > 0, 1, 0),
      IF = mean(IF)
    ) %>%
    mutate(
      # as factor
      Review = factor(Review)
    )
  # fit model
  brm(ControlNI ~ 0 + Review + Review:log(IF),
      data = review, family = bernoulli,
      prior = prior(normal(0, 0.5), class = b),
      cores = 4, seed = 2113, iter = 3000,
      control = list(adapt_delta = 0.99))
}

# function to get bootstrapped confidence intervals for proportion
getBootCIProp <- function(prop, n, nboot) {
  set.seed(1)
  p <- c()
  for (i in 1:nboot) {
    p[i] <- mean(rbinom(n, 1, prop))
  }
  return(quantile(p, c(0.025, 0.975)))
}

# summary of analysis-level review results
plotReviewAnalysis <- function(review, yearAnalysis, ifAnalysis, 
                               postRM1, postRM2, postRM3, postRM4, postRM5) {
  # prepare data
  review <-
    review %>%
    mutate_at(vars(starts_with("method_")), function(x) ifelse(is.na(x), 0, x)) %>%
    mutate(Review = ifelse(Review == "Economic development", "Economic\ndevelopment", "Cultural\nvalues"),
           Review = factor(Review, levels = c("Economic\ndevelopment", "Cultural\nvalues")))
  # proportion of analyses controlling for non-independence?
  pA1 <-
    tibble(prop = c(postRM1[,,1], postRM1[,,2]),
           Review = rep(c("Economic\ndevelopment", "Cultural\nvalues"), each = length(postRM1[,,1]))) %>%
    mutate(prop = inv_logit_scaled(prop),
           Review = factor(Review, levels = c("Economic\ndevelopment", "Cultural\nvalues"))) %>%
    group_by(Review) %>%
    summarise(avg   = median(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975)) %>%
    ggplot(aes(x = "Yes", y = avg, ymin = lower, ymax = upper, colour = Review)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    labs(x = " ",
         y = "Proportion of analyses\ncontrolling for non-independence") +
    ylim(c(0, 1)) +
    theme_classic()
  # proportion of analyses using methods of ni control?
  pA2 <-
    tibble(prop = c(postRM2[,,1], postRM3[,,1], postRM4[,,1], postRM5[,,1],
                    postRM2[,,2], postRM3[,,2], postRM4[,,2], postRM5[,,2]),
           Review = rep(c("Economic\ndevelopment", "Cultural\nvalues"), each = length(postRM2[,,1])*4),
           method = rep(rep(c("Region\nFEs", "Distance", "Cultural\nhistory", "Other"), 
                            each = length(postRM2[,,1])), times = 2)) %>%
    mutate(Review = factor(Review, levels = c("Economic\ndevelopment", "Cultural\nvalues")),
           prop = inv_logit_scaled(prop)) %>%
    group_by(Review, method) %>%
    summarise(avg   = median(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975)) %>%
    mutate(method = factor(method, levels = c("Region\nFEs", "Distance", "Cultural\nhistory", "Other"))) %>%
    ggplot(aes(x = method, y = avg, ymin = lower, ymax = upper, colour = Review)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    labs(x = "Control for non-independence?                   ", 
         y = "Proportion of analyses\ncontrolling for non-independence") +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  # impact factor results
  pBdata <- plot(conditional_effects(ifAnalysis, resolution = 500, prob = 0.50), plot = FALSE)[[2]]$data
  pB <-
    ggplot() +
    stat_dots(
      data = mutate(
        review, 
        ControlNI = ifelse(ControlNI == "Yes", 1, 0),
        side = ifelse(ControlNI == 0, "top", "bottom"),
        Review = ifelse(Review == "Cultural\nvalues", "Values", "Economic development")
        ),
      aes(y = ControlNI, x = IF, side = side, fill = Review, colour = Review),
      scale = 0.5
    ) +
    geom_ribbon(
      data = pBdata, 
      aes(x = IF, ymin = lower__, ymax = upper__, fill = Review),
      alpha = 0.4
      ) +
    geom_line(
      data = pBdata,
      aes(x = IF, y = estimate__, colour = Review)
      ) +
    scale_y_continuous(name = "Probability of analysis\ncontrolling for non-independence", 
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_log10(name = "Journal impact factor") +
    theme_classic()
  # years published
  years <- c(1995, 2000, 2005, 2010, 2015, 2020)
  # spline results
  pCdata <- plot(conditional_effects(yearAnalysis, prob = 0.50), plot = FALSE)[[3]]$data
  pC <-
    ggplot() +
    stat_dots(
      data = mutate(
        review, 
        ControlNI = ifelse(ControlNI == "Yes", 1, 0),
        side = ifelse(ControlNI == 0, "top", "bottom"),
        Review = ifelse(Review == "Cultural\nvalues", "Values", "Economic development"),
        Year.std = as.numeric(scale(Year))
      ),
      aes(y = ControlNI, x = Year.std, side = side, fill = Review, colour = Review),
      scale = 0.5
    ) +
    geom_ribbon(
      data = pCdata, 
      aes(x = Year.std, ymin = lower__, ymax = upper__, fill = Review),
      alpha = 0.4
    ) +
    geom_line(
      data = pCdata,
      aes(x = Year.std, y = estimate__, colour = Review)
    ) +
    scale_y_continuous(name = "Probability of analysis\ncontrolling for non-independence", 
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(name = "Year", labels = function(x) round((x*sd(review$Year)) + mean(review$Year), 0),
                       breaks = (years - mean(review$Year)) / sd(review$Year)) +
    theme_classic()
  # put together
  top <- plot_grid(
    pA1 + theme(legend.position = "none"), 
    pA2 + theme(legend.position = "none"),
    rel_widths = c(0.5, 1),
    nrow = 1,
    labels = c("a", ""),
    align = "h"
  )
  bottom <- plot_grid(
    pB + theme(legend.position = "none"), 
    pC + theme(legend.position = "none"), 
    nrow = 1,
    labels = c("b", "c"),
    align = "h"
  )
  out <- plot_grid(top, bottom, nrow = 2)
  legend <-
    get_legend(pA1 +
               theme(legend.spacing.y = unit(0.2, 'cm'))  +
               guides(colour = guide_legend(byrow = TRUE)))
  out <- plot_grid(out, legend, nrow = 1, rel_widths = c(1, 0.22))
  # save
  ggsave(out, filename = "figures/reviewAnalysis.pdf", height = 6, width = 7)
  return(out)
}

# summary of article-level review results
plotReviewArticle <- function(review, yearArticle, ifArticle) {
  # prepare data
  review <-
    review %>%
    mutate_at(vars(starts_with("method_")), function(x) ifelse(is.na(x), 0, x)) %>%
    mutate(Review = ifelse(Review == "Economic development", "Economic\ndevelopment", "Cultural\nvalues"),
           Review = factor(Review, levels = c("Economic\ndevelopment", "Cultural\nvalues")))
  # proportion of papers controlling for non-independence at least once?
  pA1 <-
    review %>%
    group_by(Review, Citation) %>%
    summarise(Yes = ifelse(sum(ControlNI == "Yes") > 0, 1, 0)) %>%
    group_by(Review) %>%
    summarise(prop = mean(Yes)) %>%
    # bootstrap cis
    mutate(lower = map(prop, function(x) getBootCIProp(x, 50, 1000)[1]), 
           upper = map(prop, function(x) getBootCIProp(x, 50, 1000)[2])) %>%
    unnest(c(lower, upper)) %>%
    ggplot(aes(x = "Yes", y = prop, ymin = lower, ymax = upper, colour = Review)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    labs(x = " ",
         y = "Proportion of articles\ncontrolling for non-independence") +
    ylim(c(0, 1)) +
    theme_classic()
  # proportion of papers using methods of ni control at least once?
  pA2 <-
    review %>%
    pivot_longer(cols = starts_with("method_"),
                 names_to = "method") %>%
    group_by(Review, Citation, method) %>%
    summarise(Yes = ifelse(sum(value == 1) > 0, 1, 0)) %>%
    group_by(Review, method) %>%
    summarise(prop = mean(Yes)) %>%
    # bootstrap cis
    mutate(lower = map(prop, function(x) getBootCIProp(x, 50, 1000)[1]), 
           upper = map(prop, function(x) getBootCIProp(x, 50, 1000)[2])) %>%
    unnest(c(lower, upper)) %>%
    mutate(method = str_remove(method, "method_"),
           method = ifelse(method == "SharedCulturalHistory", "Cultural\nhistory",
                           ifelse(method == "RegionalFixedEffects", "Region\nFEs",
                                  method)),
           method = factor(method, levels = c("Region\nFEs", "Distance", "Cultural\nhistory", "Other"))) %>%
    ggplot(aes(x = method, y = prop, ymin = lower, ymax = upper, colour = Review)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    labs(x = "Control for non-independence?                   ", 
         y = "Proportion of articles\ncontrolling for non-independence") +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  # prepare review raw data by article
  review100 <-
    review %>%
    group_by(Review, Citation) %>%
    summarise(
      ControlNI = ifelse(sum(ControlNI == "Yes") > 0, 1, 0),
      Year = mean(Year),
      IF = mean(IF)
    ) %>%
    mutate(
      Year.std = as.numeric(scale(Year)),
      Review = factor(Review)
    )
  # impact factor results
  pBdata <- plot(conditional_effects(ifArticle, resolution = 500, prob = 0.50), plot = FALSE)[[2]]$data
  pB <-
    ggplot() +
    stat_dots(
      data = mutate(
        review100, 
        side = ifelse(ControlNI == 0, "top", "bottom"),
        Review = ifelse(Review == "Cultural\nvalues", "Values", "Economic development")
      ),
      aes(y = ControlNI, x = IF, side = side, fill = Review, colour = Review),
      scale = 0.5, dotsize = 0.35
    ) +
    geom_ribbon(
      data = pBdata, 
      aes(x = IF, ymin = lower__, ymax = upper__, fill = Review),
      alpha = 0.4
    ) +
    geom_line(
      data = pBdata,
      aes(x = IF, y = estimate__, colour = Review)
    ) +
    scale_y_continuous(name = "Probability of article\ncontrolling for non-independence", 
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_log10(name = "Journal impact factor") +
    theme_classic()
  # years published
  years <- c(1995, 2000, 2005, 2010, 2015, 2020)
  meanYear <- review %>% group_by(Review, Citation) %>% summarise(Year = mean(Year)) %>% pull(Year) %>% mean()
  sdYear   <- review %>% group_by(Review, Citation) %>% summarise(Year = mean(Year)) %>% pull(Year) %>% sd()
  # spline results
  pCdata <- plot(conditional_effects(yearArticle, prob = 0.50), plot = FALSE)[[3]]$data
  pC <-
    ggplot() +
    stat_dots(
      data = mutate(
        review100, 
        side = ifelse(ControlNI == 0, "top", "bottom"),
        Review = ifelse(Review == "Cultural\nvalues", "Values", "Economic development")
      ),
      aes(y = ControlNI, x = Year.std, side = side, fill = Review, colour = Review),
      scale = 0.5, dotsize = 0.3
    ) +
    geom_ribbon(
      data = pCdata, 
      aes(x = Year.std, ymin = lower__, ymax = upper__, fill = Review),
      alpha = 0.4
    ) +
    geom_line(
      data = pCdata,
      aes(x = Year.std, y = estimate__, colour = Review)
    ) +
    scale_y_continuous(name = "Probability of article\ncontrolling for non-independence", 
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(name = "Year", labels = function(x) round((x*sdYear) + meanYear, 0),
                       breaks = (years - meanYear) / sdYear) +
    theme_classic()
  # put together
  top <- plot_grid(
    pA1 + theme(legend.position = "none"), 
    pA2 + theme(legend.position = "none"),
    rel_widths = c(0.5, 1),
    nrow = 1,
    labels = c("a", ""),
    align = "h"
  )
  bottom <- plot_grid(
    pB + theme(legend.position = "none"), 
    pC + theme(legend.position = "none"), 
    nrow = 1,
    labels = c("b", "c"),
    align = "h"
  )
  out <- plot_grid(top, bottom, nrow = 2)
  legend <-
    get_legend(pA1 +
                 theme(legend.spacing.y = unit(0.2, 'cm'))  +
                 guides(colour = guide_legend(byrow = TRUE)))
  out <- plot_grid(out, legend, nrow = 1, rel_widths = c(1, 0.22))
  # save
  ggsave(out, filename = "figures/reviewArticle.pdf", height = 6, width = 7)
  return(out)
}

# make table with all articles in review
makeTableArticlesReview <- function(review) {
  review %>%
    select(Review, Citation, CitationsPerYear) %>%
    distinct() %>%
    rename(Reference = Citation, `Citations per year` = CitationsPerYear) %>%
    mutate(Review = c("Cultural values", rep("", times = 49),
                      "Economic development", rep("", times = 49)))
}