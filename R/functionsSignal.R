# custom functions - signal

# load covariance matrices
loadCovMat <- function(file, log) {
  # load distance matrix
  out <- 
    read_excel(file, na = "") %>%
    dplyr::select(-ISO) %>%
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
    dplyr::select(c("Country", "1990", "2000", "2010", 
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

# load gdp data
loadGDPData <- function(fileGDP, iso) {
  read_csv(fileGDP, skip = 4) %>%
    # join iso2 data
    left_join(dplyr::select(iso, iso3, iso2), by = c("Country Code" = "iso3")) %>%
    # drop countries without iso2 codes
    drop_na(iso2) %>%
    # relevant columns
    dplyr::select(-`Indicator Name`, -`Indicator Code`, -`...67`, -`Country Code`) %>%
    # pivot longer
    pivot_longer(
      cols = !c(`Country Name`, iso2),
      names_to = "year",
      values_to = "GDPPerCapita"
    ) %>%
    # drop empty cases
    drop_na(GDPPerCapita) %>%
    mutate(
      # log gdp per capita
      logGDPPerCapita = log(GDPPerCapita),
      # standardise
      logGDPPerCapitaStd = as.numeric(scale(logGDPPerCapita))
      )
}

# load gdp growth data
loadGDPGrowthData <- function(fileGDPGrowth, iso) {
  read_csv(fileGDPGrowth, skip = 4) %>%
    # join iso2 data
    left_join(dplyr::select(iso, iso3, iso2), by = c("Country Code" = "iso3")) %>%
    # drop countries without iso2 codes
    drop_na(iso2) %>%
    # relevant columns
    dplyr::select(-`Indicator Name`, -`Indicator Code`, -`...67`, -`Country Code`) %>%
    # pivot longer
    pivot_longer(
      cols = !c(`Country Name`, iso2),
      names_to = "year",
      values_to = "gdpPerCapitaGrowth"
    ) %>%
    # drop empty cases
    drop_na(gdpPerCapitaGrowth) %>%
    # standardised
    mutate(gdpPerCapitaGrowthStd = as.numeric(scale(gdpPerCapitaGrowth)))
}

# load gini data
loadGiniData <- function(fileGini, iso) {
  read_csv(fileGini, skip = 4) %>%
    # join iso2 data
    left_join(dplyr::select(iso, iso3, iso2), by = c("Country Code" = "iso3")) %>%
    # drop countries without iso2 codes
    drop_na(iso2) %>%
    # relevant columns
    dplyr::select(-`Indicator Name`, -`Indicator Code`, -`...67`, -`Country Code`) %>%
    # pivot longer
    pivot_longer(
      cols = !c(`Country Name`, iso2),
      names_to = "year",
      values_to = "gini"
    ) %>%
    # drop empty cases
    drop_na(gini) %>%
    # gini between 0 and 1
    mutate(gini = gini / 100) -> d
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

# signal for gdp
fitGDPSignal <- function(gdp, geoCov, linCov) {
  # add iso variables for modelling
  gdp$isoGeo <- gdp$iso2
  gdp$isoLin <- gdp$iso2
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% gdp$iso2],
                   colnames(geoCov)[colnames(geoCov) %in% gdp$iso2]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% gdp$iso2],
                   colnames(linCov)[colnames(linCov) %in% gdp$iso2]]
  # fit model
  brm(logGDPPerCapitaStd ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)) + (1 | iso2),
      data = gdp, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(3), class = sd),
                prior(exponential(3), class = sigma)),
      seed = 2113, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.9))
}

# signal for gdp growth
fitGDPGrowthSignal <- function(gdpGrowth, geoCov, linCov) {
  # add iso variables for modelling
  gdpGrowth$isoGeo <- gdpGrowth$iso2
  gdpGrowth$isoLin <- gdpGrowth$iso2
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% gdpGrowth$iso2],
                   colnames(geoCov)[colnames(geoCov) %in% gdpGrowth$iso2]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% gdpGrowth$iso2],
                   colnames(linCov)[colnames(linCov) %in% gdpGrowth$iso2]]
  # fit model
  brm(gdpPerCapitaGrowthStd ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)) + (1 | iso2),
      data = gdpGrowth, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0, 0.5), class = b),
                prior(exponential(3), class = sd),
                prior(exponential(3), class = sigma)),
      seed = 2113, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.9))
}

# signal for gini
fitGiniSignal <- function(gini, geoCov, linCov) {
  # add iso variables for modelling
  gini$isoGeo <- gini$iso2
  gini$isoLin <- gini$iso2
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% gini$iso2],
                   colnames(geoCov)[colnames(geoCov) %in% gini$iso2]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% gini$iso2],
                   colnames(linCov)[colnames(linCov) %in% gini$iso2]]
  # fit model
  brm(gini ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)) + (1 | iso2),
      data = gini, data2 = list(geoCov = geoCov, linCov = linCov),
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
    dplyr::select(
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

# signal for tightness
fitTightnessSignal <- function(tightness, geoCov, linCov) {
  # add iso variables for modelling
  tightness$isoGeo <- tightness$iso
  tightness$isoLin <- tightness$iso
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% tightness$iso],
                   colnames(geoCov)[colnames(geoCov) %in% tightness$iso]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% tightness$iso],
                   colnames(linCov)[colnames(linCov) %in% tightness$iso]]
  # fit model
  brm(Tightness ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)),
      data = tightness, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(9), class = sd),
                prior(exponential(9), class = sigma)),
      iter = 6000, seed = 2113, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.9999, max_treedepth = 15))
}

# signal for individualism
fitIndividualismSignal <- function(fincherData, geoCov, linCov) {
  # add iso variables for modelling
  fincherData$isoGeo <- fincherData$iso
  fincherData$isoLin <- fincherData$iso
  # subset matrices
  geoCov <- geoCov[rownames(geoCov)[rownames(geoCov) %in% fincherData$iso],
                   colnames(geoCov)[colnames(geoCov) %in% fincherData$iso]]
  linCov <- linCov[rownames(linCov)[rownames(linCov) %in% fincherData$iso],
                   colnames(linCov)[colnames(linCov) %in% fincherData$iso]]
  # outcome variable between 0 and 1
  fincherData <- mutate(fincherData, individualismHofstede = individualismHofstede / 100)
  # fit model
  brm(individualismHofstede ~ 0 + Intercept + (1 | gr(isoGeo, cov = geoCov)) + (1 | gr(isoLin, cov = linCov)),
      data = fincherData, data2 = list(geoCov = geoCov, linCov = linCov),
      prior = c(prior(normal(0.5, 0.1), class = b),
                prior(exponential(8), class = sd),
                prior(exponential(8), class = sigma)),
      iter = 6000, seed = 2113, cores = 4, sample_prior = "yes",
      control = list(adapt_delta = 0.9999, max_treedepth = 15))
}

# plot signal
plotGeoLinSignal <- function(geoHDI, geoGDP, geoGrow, geoGini,
                             geoTrad, geoSurv, geoTight, geoInd,
                             linHDI, linGDP, linGrow, linGini,
                             linTrad, linSurv, linTight, linInd) {
  # vector of outcomes
  outcomes <- c("HDI",
                "GDP per capita",
                "GDP per capita growth",
                "Gini index",
                "Traditional values",
                "Survival values",
                "Tightness",
                "Individualism")
  # plot
  out <-
    tibble(
      Geographic = list(geoHDI$samples$H1, geoGDP$samples$H1, 
                        geoGrow$samples$H1, geoGini$samples$H1,
                        geoTrad$samples$H1, geoSurv$samples$H1,
                        geoTight$samples$H1, geoInd$samples$H1),
      Cultural = list(linHDI$samples$H1, linGDP$samples$H1,
                      linGrow$samples$H1, linGini$samples$H1,
                      linTrad$samples$H1, linSurv$samples$H1,
                      linTight$samples$H1, linInd$samples$H1),
      Outcome = factor(outcomes, levels = outcomes)
    ) %>%
    unnest(c(Geographic, Cultural)) %>%
    mutate(Residual = 1 - (Geographic + Cultural)) %>%
    pivot_longer(
      cols = !Outcome,
      names_to = "par",
      values_to = "signal"
    ) %>%
    mutate(par = factor(par, levels = c("Residual", "Cultural", "Geographic"))) %>%
    ggplot(aes(x = signal, y = par)) +
    stat_halfeye(normalize = "xy") +
    facet_wrap(. ~ Outcome, nrow = 2) +
    labs(x = "Proportion of national-level variance explained", y = NULL) +
    xlim(c(0, 1)) +
    theme_classic() +
    theme(panel.spacing = unit(1.1, "lines"),
          axis.title.x = element_text(margin = margin(7, 0, 0, 0)))
  # save
  ggsave(out, filename = "figures/signal.pdf", width = 7, height = 5)
  return(out)
}

# make table of signal values
makeTableGeoLinSignal <- function(geoHDI, geoGDP, geoGrow, geoGini,
                                  geoTrad, geoSurv, geoTight, geoInd,
                                  linHDI, linGDP, linGrow, linGini,
                                  linTrad, linSurv, linTight, linInd) {
  # tidy function
  tidySignal <- function(signal) {
    paste0(
      format(round(signal$hypothesis$Estimate, 2), nsmall = 2), ", 95% CI [",
      format(round(signal$hypothesis$CI.Lower, 2), nsmall = 2), ", ",
      format(round(signal$hypothesis$CI.Upper, 2), nsmall = 2), "], BF ",
      ifelse(
        signal$hypothesis$Evid.Ratio <= 0, "> 100",
        ifelse(
          1 / signal$hypothesis$Evid.Ratio > 100, "> 100",
          paste0("= ", format(round(1 / signal$hypothesis$Evid.Ratio, 2), nsmall = 2)))
        )
      )
  }
  # make table
  tibble(
    Outcome = c("HDI", "GDPpc", "GDPpc growth", "Gini index",
                "Traditional values", "Survival values", "Tightness", "Individualism"),
    `Geographic signal` = c(
      tidySignal(geoHDI),
      tidySignal(geoGDP),
      tidySignal(geoGrow),
      tidySignal(geoGini),
      tidySignal(geoTrad),
      tidySignal(geoSurv),
      tidySignal(geoTight),
      tidySignal(geoInd)
    ),
    `Cultural phylogenetic signal` = c(
      tidySignal(linHDI),
      tidySignal(linGDP),
      tidySignal(linGrow),
      tidySignal(linGini),
      tidySignal(linTrad),
      tidySignal(linSurv),
      tidySignal(linTight),
      tidySignal(linInd)
    )
  )
}
