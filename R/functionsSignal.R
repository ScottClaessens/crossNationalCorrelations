# custom functions - signal

# load covariance matrices
loadCovMat <- function(file, log) {
  # load distance matrix
  out <- 
    read_excel(file, na = "") %>%
    select(-ISO) %>%
    as.matrix()
  rownames(out) <- colnames(out)
  # log distances?
  if (log) out <- log(out)
  # distances between 0 and 1
  out <- out / max(out)
  diag(out) <- 0
  # 1 - distance = proximity (covariance)
  out <- 1 - out
  return(out)
}

# load HDI data
loadHDIData <- function(fileHDI, fileISOHDI) {
  # tidy up messy excel sheet
  read_xlsx(
    path = fileHDI,
    sheet = "Table 2",
    skip = 4
  ) %>%
    drop_na("HDI rank") %>%
    select(c("Country", "1990", "2000", "2010", 
             "2014", "2015", "2017", "2018", "2019")) %>%
    # get data in long form
    pivot_longer(cols = !Country, names_to = "year", values_to = "hdi") %>%
    filter(hdi != "..") %>%
    mutate(year = as.numeric(year),
           hdi = round(as.numeric(hdi), 3)) %>%
    # add iso-codes
    left_join(read.csv(fileISOHDI), by = "Country") %>%
    # complete missings
    complete(year, nesting(iso, Country)) %>%
    # fix cote d'ivoire error
    mutate(iso = ifelse(is.na(iso), "CI", iso))
}

# signal for hdi
fitHDISignal <- function(hdi, geoCov, linCov) {
  # add iso variables for modelling
  hdi$isoGeo <- hdi$iso
  hdi$isoLin <- hdi$iso
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% hdi$iso],
                   colnames(geoCov)[colnames(geoCov) %in% hdi$iso]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% hdi$iso],
                   colnames(linCov)[colnames(linCov) %in% hdi$iso]]
  # fit model
  brm(hdi ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)) + (1 | iso),
      data = hdi, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0.5, 0.1), class = b),
                prior(exponential(8), class = sd),
                prior(exponential(8), class = sigma)),
      seed = 2113, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.9))
}

# signal for wvs
fitWVSSignal <- function(wvs, outcome = "", geoCov, linCov) {
  wvs <-
    wvs %>%
    # select vars
    select(
      # country code
      S009,
      # traditional vs. secular-rational values
      F063, Y003, F120, G006, E018,
      # survival vs. self-expression values
      Y002, A008, F118, E025, A165
    ) %>%
      # remove missings
      drop_na() %>%
      filter(Y003 != -5)
  # get cultural map
  # https://www.worldvaluessurvey.org/WVSContents.jsp?CMSID=Findings
  m <- pca(select(wvs, -S009), nfactors = 2, rotate = "varimax")
  # summarise and prepare for modelling
  wvs <-
    wvs %>%
    mutate(trad = as.numeric(scale(m$scores[,2])),
           surv = as.numeric(scale(m$scores[,1])),
           isoGeo = S009, 
           isoLin = S009)
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% wvs$S009],
                   colnames(geoCov)[colnames(geoCov) %in% wvs$S009]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% wvs$S009],
                   colnames(linCov)[colnames(linCov) %in% wvs$S009]]
  # fit model
  bf1 <- bf(paste0(outcome, " ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)) + (1 | S009)"))
  brm(formula = bf1, data = wvs, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(8), class = sd),
                prior(exponential(8), class = sigma)),
      seed = 2113, iter = 2000, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.95, max_treedepth = 15))
}

# plot signal
plotGeoLinSignal <- function(postHDI, postTrad, postSurv) {
  # plotting function
  plotFun <- function(data, xlab, title) {
    data %>%
      pivot_longer(
        cols = everything(),
        names_to = "par",
        values_to = "signal"
      ) %>%
      mutate(par = factor(par, levels = c("Residual", "Cultural", "Geographic"))) %>%
      ggplot(aes(x = signal, y = par)) +
      stat_halfeye(normalize = "xy") +
      labs(x = xlab, y = NULL) +
      xlim(c(0, 1)) +
      ggtitle(title) +
      theme_classic() +
      theme(plot.title = element_text(size = 12))
  }
  # get data for plots
  pA <-
    tibble(
      Geographic = as.numeric(postHDI[,,2]^2 / (postHDI[,,1]^2 + postHDI[,,2]^2 + postHDI[,,3]^2)),
      Cultural   = as.numeric(postHDI[,,3]^2 / (postHDI[,,1]^2 + postHDI[,,2]^2 + postHDI[,,3]^2)),
      Residual   = as.numeric(postHDI[,,1]^2 / (postHDI[,,1]^2 + postHDI[,,2]^2 + postHDI[,,3]^2))
    )
  pB <-
    tibble(
      Geographic = as.numeric(postTrad[,,1]^2 / (postTrad[,,1]^2 + postTrad[,,2]^2 + postTrad[,,3]^2)),
      Cultural   = as.numeric(postTrad[,,2]^2 / (postTrad[,,1]^2 + postTrad[,,2]^2 + postTrad[,,3]^2)),
      Residual   = as.numeric(postTrad[,,3]^2 / (postTrad[,,1]^2 + postTrad[,,2]^2 + postTrad[,,3]^2))
    )
  pC <-
    tibble(
      Geographic = as.numeric(postSurv[,,1]^2 / (postSurv[,,1]^2 + postSurv[,,2]^2 + postSurv[,,3]^2)),
      Cultural   = as.numeric(postSurv[,,2]^2 / (postSurv[,,1]^2 + postSurv[,,2]^2 + postSurv[,,3]^2)),
      Residual   = as.numeric(postSurv[,,3]^2 / (postSurv[,,1]^2 + postSurv[,,2]^2 + postSurv[,,3]^2))
    )
  # individual plots
  pA <- plotFun(pA, xlab = " \n ", title = "Human Development\nIndex (HDI)")
  pB <- plotFun(pB, xlab = "Proportion of national-level\nvariance explained", 
                title = "Traditional vs.\nsecular values")
  pC <- plotFun(pC, xlab = " \n ", title = "Survival vs.\nself-expression values")
  # put together
  out <- plot_grid(pA, 
                   pB + theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank()), 
                   pC + theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank()), 
                   nrow = 1, rel_widths = c(1.35, 1, 1))
  # save
  ggsave(out, filename = "figures/signal.pdf", width = 7, height = 3.5)
  return(out)
}
