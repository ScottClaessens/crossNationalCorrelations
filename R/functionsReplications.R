# custom functions - replications

# load data for Adamczyk 2009
loadDataAdamczyk2009 <- function(wvs, iso) {
  wvs %>%
    # wvs only
    filter(studyno == 4001) %>%
    # wave 4
    filter(S002 == 4) %>%
    # keep relevant variables
    select(
      # country code
      S009,
      # survival vs. self-expression values
      Y002, A008, E025, A165, # does not include F118 as it is dv
      # disapproval of homosexuality (dependent variable)
      F118,
      # covariates
      A006, # importance of religion
      F025, # religious denomination
      X007, # marital status
      X011, # number of children
      X025, # education
      X001, # gender
      C006, # satisfaction with financial situation
      X002  # year of birth
    ) %>%
    # drop missings
    drop_na() %>%
    mutate(
      # create survival vs. self-expression values composite
      surv = as.numeric(scale(4-Y002) + scale(A008) + scale(E025) + scale(A165)),
      # reverse code dependent variable
      F118 = 11 - F118,
      # reverse code religion importance
      A006 = 5 - A006,
      # 1 = married, 0 = all else
      X007 = ifelse(X007 == 1, "Married", "Not married"),
      # birth cohorts
      X002 = ifelse(X002 %in% 1920:1929, "1920s",
                    ifelse(X002 %in% 1930:1939, "1930s",
                           ifelse(X002 %in% 1940:1949, "1940s",
                                  ifelse(X002 %in% 1950:1959, "1950s",
                                         ifelse(X002 %in% 1960:1969, "1960s",
                                                ifelse(X002 %in% 1970:1979, "1970s", "1980s")))))),
      # principal religious traditions for countries
      countryReligTrad = ifelse(S009 %in% c("CA","PR","MX","VE","ES",
                                            "AR","PE","CL","UG","PH"), "Catholic",
                                ifelse(S009 %in% c("IN"), "Hindu",
                                       ifelse(S009 %in% c("MD","RS","ME","MK"), "Orthodox Christian",
                                              ifelse(S009 %in% c("US","KR","ZA","ZW"), "Protestant",
                                                     ifelse(S009 %in% c("JP","VN"), "Buddhist", "Islam"))))),
      # as factors
      F025 = as.factor(F025),
      X001 = as.factor(X001),
      X002 = as.factor(X002),
      countryReligTrad = as.factor(countryReligTrad)
    ) %>%
    # country level survival value means
    group_by(S009) %>%
    mutate(countrySurv = mean(surv)) %>%
    ungroup() %>%
    # join lon lat
    left_join(iso, by = c("S009" = "iso2")) %>%
    # new var for modelling
    mutate(iso2lin = S009)
}

# fit model to Adamczyk 2009 data
fitModelAdamczyk2009 <- function(adamczykData, linCov, control = FALSE) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% adamczykData$iso2lin,
                   colnames(linCov) %in% adamczykData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.5), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(5), class = sd),
              prior(exponential(5), class = sigma))
  # model formula (Model 5 in paper - Table 4) and additional priors
  if (control == "none") {
    model <- bf(F118 ~ surv + X007 + X011 + X025 + X001 + C006 + 
                  X002 + F025 + A006 + countryReligTrad + countrySurv + 
                  (1 | S009))
  } else if (control == "spatial") {
    model <- bf(F118 ~ surv + X007 + X011 + X025 + X001 + C006 + 
                  X002 + F025 + A006 + countryReligTrad + countrySurv + 
                  gp(Longitude..average., Latitude..average., k = 5, c = 5/4, gr = TRUE) + 
                  (1 | S009))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(F118 ~ surv + X007 + X011 + X025 + X001 + C006 + 
                  X002 + F025 + A006 + countryReligTrad + countrySurv + 
                  (1 | gr(iso2lin, cov = linCov)) + (1 | S009))
  } else if (control == "both") {
    model <- bf(F118 ~ surv + X007 + X011 + X025 + X001 + C006 + 
                  X002 + F025 + A006 + countryReligTrad + countrySurv + 
                  gp(Longitude..average., Latitude..average., k = 5, c = 5/4, gr = TRUE) + 
                  (1 | gr(iso2lin, cov = linCov)) + (1 | S009))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  adamczykData %>%
    # remove labels for modelling
    haven::zap_labels() %>%
    # standardise numeric variables
    mutate(
      F118 = as.numeric(scale(F118)),
      surv = as.numeric(scale(surv)),
      X011 = as.numeric(scale(X011)),
      X025 = as.numeric(scale(X025)),
      C006 = as.numeric(scale(C006)),
      A006 = as.numeric(scale(A006)),
      countrySurv = as.numeric(scale(countrySurv))
    ) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000, 
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Alesina 2013
loadDataAlesina2013 <- function(fileAlesina, iso) {
  read_xlsx(fileAlesina) %>%
    left_join(iso, by = c("iso" = "iso3")) %>%
    mutate(iso2 = ifelse(is.na(iso2), "NA", iso2),
           # new var for modelling
           iso2lin = iso2)
}

# fit model to Alesina 2013 data
fitModelAlesina2013 <- function(alesinaData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% alesinaData$iso2lin,
                   colnames(linCov) %in% alesinaData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(FLFP2000 ~ tradPloughMean)
  } else if (control == "spatial") {
    model <- bf(FLFP2000 ~ tradPloughMean + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(FLFP2000 ~ tradPloughMean + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(FLFP2000 ~ tradPloughMean + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) +
                  (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sdgp))
    priors <- c(priors, prior(exponential(5), class = sd))
  }
  # begin pipe
  alesinaData %>%
    # drop NAs manually
    drop_na(FLFP2000, tradPloughMean) %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000, 
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Beck 2003
loadDataBeck2003 <- function(fileBeck1, iso) {
  read_xlsx(fileBeck1) %>%
    left_join(iso, by = c("Country code" = "iso3")) %>%
    mutate(
      # iso 2 code for Zaire
      iso2 = ifelse(is.na(iso2), "ZR", iso2),
      # log settler mortality
      logSettlerMortality = log(`Settler mortality`),
      # new var for modelling
      iso2lin = iso2)
}

# fit model to Beck 2003 data
fitModelBeck2003 <- function(beckData1, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% beckData1$iso2lin,
                   colnames(linCov) %in% beckData1$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(stockMarketDev ~ logSettlerMortality)
  } else if (control == "spatial") {
    model <- bf(stockMarketDev ~ logSettlerMortality + gp(Longitude..average., Latitude..average.))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(stockMarketDev ~ logSettlerMortality + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(stockMarketDev ~ logSettlerMortality + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  beckData1 %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # rename for modelling
    rename(stockMarketDev = `Stock market development`) %>%
    # remove ZR - not in linguistic matrix
    filter(iso2lin != "ZR") %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Beck 2005
loadDataBeck2005 <- function(fileBeck2, iso) {
  read_xlsx(fileBeck2) %>%
    left_join(iso, by = "iso2") %>%
    # new vars for modelling
    mutate(iso2lin = iso2)
}

# fit model to Beck 2005 data
fitModelBeck2005 <- function(beckData2, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% beckData2$iso2lin,
                   colnames(linCov) %in% beckData2$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(GDPpc ~ SME250)
  } else if (control == "spatial") {
    model <- bf(GDPpc ~ SME250 + gp(Longitude..average., Latitude..average.))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(GDPpc ~ SME250 + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(GDPpc ~ SME250 + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    # stronger priors to help model converge
    priors <- c(priors, prior(exponential(10), class = sd))
    priors <- c(priors, prior(exponential(10), class = sdgp))
  }
  # begin pipe
  beckData2 %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Bockstette 2002
loadDataBockstette2002 <- function(fileBockstette1, fileBockstette2, iso) {
  # load state hist data
  stateHist <- read_xlsx(fileBockstette1)
  # load gdp growth data
  gdpGrowth <- 
    read_xlsx(fileBockstette2, na = "..") %>%
    # keep countries only
    slice(1:217) %>%
    transmute(iso3 = `Country Code`,
              gdpGrowth = rowMeans(select(., starts_with("19")), na.rm = TRUE) / 100) %>%
    # remove missings
    drop_na()
  # join datasets
  stateHist %>%
    # link iso
    left_join(iso, by = "iso2") %>%
    # link gdp growth
    left_join(gdpGrowth, by = "iso3") %>%
    drop_na() %>%
    # add variable for modelling
    mutate(iso2lin = iso2)
}

# fit model to Bockstette 2002 data
fitModelBockstette2002 <- function(bockstetteData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% bockstetteData$iso2lin,
                   colnames(linCov) %in% bockstetteData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(gdpGrowth ~ Statehist5)
  } else if (control == "spatial") {
    model <- bf(gdpGrowth ~ Statehist5 + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(gdpGrowth ~ Statehist5 + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(gdpGrowth ~ Statehist5 + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  bockstetteData %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Easterly 2003
loadDataEasterly2003 <- function(fileEasterly1, fileEasterly2, iso) {
  # log gdp data from appendix
  lg <- read_xlsx(fileEasterly1)
  # institutions index (average of six wgi variables in 1998)
  wgi1 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 2) %>% transmute(iso3 = Code, wgi1 = Estimate...9)
  wgi2 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 3) %>% transmute(wgi2 = Estimate...9)
  wgi3 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 4) %>% transmute(wgi3 = Estimate...9)
  wgi4 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 5) %>% transmute(wgi4 = Estimate...9)
  wgi5 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 6) %>% transmute(wgi5 = Estimate...9)
  wgi6 <- read_xlsx(fileEasterly2, skip = 14, na = "#N/A", sheet = 7) %>% transmute(wgi6 = Estimate...9)
  instit <- 
    cbind(wgi1, wgi2, wgi3, wgi4, wgi5, wgi6) %>%
    transmute(iso3 = iso3,
              institIndex = rowMeans(select(., -iso3), na.rm = TRUE)) %>%
    drop_na()
  # put together
  lg %>%
    left_join(instit, by = c("isocode" = "iso3")) %>%
    left_join(iso, by = c("isocode" = "iso3")) %>%
    mutate(
      # fix problem with NA iso
      iso2 = ifelse(is.na(iso2), "NA", iso2),
      # new var for modelling
      iso2lin = iso2
      )
}

# fit model to Easterly 2003 data
fitModelEasterly2003 <- function(easterlyData1, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% easterlyData1$iso2lin,
                   colnames(linCov) %in% easterlyData1$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(logGDPpc1995 ~ institIndex)
  } else if (control == "spatial") {
    model <- bf(logGDPpc1995 ~ institIndex + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(logGDPpc1995 ~ institIndex + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(logGDPpc1995 ~ institIndex + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  easterlyData1 %>%
    # drop NAs manually
    drop_na(logGDPpc1995, institIndex) %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Easterly 2007
loadDataEasterly2007 <- function(fileEasterly2, fileEasterly3, iso) {
  # individual files
  lwheatsugar <- read_xlsx(fileEasterly2)
  gini <-
    read_xls(fileEasterly3) %>%
    # last three decades (from 2007)
    filter(Year < 1987) %>%
    group_by(Country3) %>%
    summarise(gini = mean(Gini))
  # put together
  gini %>%
    left_join(iso, by = c("Country3" = "iso3")) %>%
    left_join(select(lwheatsugar, -Country), by = "iso2") %>%
    mutate(iso2lin = iso2)
}

# fit model to Easterly 2007 data
fitModelEasterly2007 <- function(easterlyData2, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% easterlyData2$iso2lin,
                   colnames(linCov) %in% easterlyData2$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(gini ~ LWHEATSUGAR)
  } else if (control == "spatial") {
    model <- bf(gini ~ LWHEATSUGAR + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(gini ~ LWHEATSUGAR + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(gini ~ LWHEATSUGAR + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  easterlyData2 %>%
    # drop NAs manually
    drop_na(gini, LWHEATSUGAR) %>%
    # standardise variables
    mutate_if(is.numeric, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Fincher 2008
loadDataFincher2008 <- function(fileFincher, iso) {
  read_xlsx(fileFincher) %>%
    # combine england and northern ireland into GB
    group_by(iso) %>%
    summarise_if(is.double, mean) %>%
    # join iso data
    left_join(iso, by = c("iso" = "iso2")) %>%
    # new var for modelling
    mutate(iso2lin = iso)
}

# fit model to Fincher 2008 data
fitModelFincher2008 <- function(fincherData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% fincherData$iso2lin,
                   colnames(linCov) %in% fincherData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(individualismHofstede ~ pathPrevHistorical)
  } else if (control == "spatial") {
    model <- bf(individualismHofstede ~ pathPrevHistorical + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(individualismHofstede ~ pathPrevHistorical + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(individualismHofstede ~ pathPrevHistorical + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  fincherData %>%
    # drop NAs manually
    drop_na(individualismHofstede, pathPrevHistorical) %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Gelfand 2011
loadDataGelfand2011 <- function(fileGelfand1, fileGelfand2, fileGelfand3, 
                                fileGelfand4, fileGelfand5, iso) {
  # individual files
  pathogens <- 
    read_xlsx(fileGelfand1) %>%
    group_by(iso) %>%
    summarise(histDiseasePrevalence9 = mean(histDiseasePrevalence9),
              histDiseasePrevalence7 = mean(histDiseasePrevalence7))
  tightness <- 
    read_xlsx(fileGelfand2) %>%
    group_by(iso) %>%
    summarise(tightness = mean(tightness))
  esi <- 
    read_xls(fileGelfand3, sheet = "Raw_Data_All_Countries", skip = 1, na = "..") %>%
    select(...1, GR2050, DISCAS) %>%
    rename(iso = ...1)
  gni <-
    read_csv(fileGelfand4, skip = 4) %>%
    transmute(iso = `Country Code`,
              logGNI2000 = log(`2000`))
  tuberculosis <-
    read_csv(fileGelfand5) %>%
    transmute(iso = SpatialDimValueCode,
              logTuberculosis = log(FactValueNumeric))
  # put together
  tightness %>%
    left_join(pathogens, by = "iso") %>%
    left_join(esi, by = "iso") %>%
    left_join(gni, by = "iso") %>%
    left_join(tuberculosis, by = "iso") %>%
    # join iso2 codes
    left_join(iso, by = c("iso" = "iso3")) %>%
    # new var for modelling
    mutate(iso2lin = iso2)
}

# fit model to Gelfand 2011 data
fitModelGelfand2011 <- function(gelfandData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% gelfandData$iso2lin,
                   colnames(linCov) %in% gelfandData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(tightness ~ DISCAS + logGNI2000)
  } else if (control == "spatial") {
    model <- bf(tightness ~ DISCAS + logGNI2000 + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(tightness ~ DISCAS + logGNI2000 + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(tightness ~ DISCAS + logGNI2000 + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  gelfandData %>%
    # remove Venzuela outlier
    mutate(DISCAS = ifelse(iso2 == "VE", NA, DISCAS)) %>%
    # drop NAs manually
    drop_na(tightness, DISCAS, logGNI2000) %>%
    # standardise variables
    mutate_if(is.double, function(x) as.numeric(scale(x))) %>%
    # fit model (partial correlation)
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 4000,
        prior = priors, control = list(adapt_delta = 0.9999))
}

# load data for Inglehart 2000
loadDataInglehart2000 <- function(fileInglehart, wvs, iso) {
  # load world bank data
  wb <- 
    read_csv(fileInglehart, skip = 4) %>%
    # need to take 1991 data, earliest on record
    transmute(iso3 = `Country Code`, indust = `1991`)
  # world values survey
  wvs <-
    wvs %>%
    filter(
      # some countries from wvs wave 3
      (S009 %in% c("US","AU","NZ","CN","JP","TW","KR",
                   "TR","BD","IN","PK","PH","AM","AZ",
                   "GE","GB","DE","CH","NO","SE","FI",
                   "ES","RU","UA","BY","EE","LV","LT",
                   "MD","PL","BG","BA","SI","HR","MK",
                   "NG","ZA","GH","AR","BR","CL","CO",
                   "DO","MX","PE","PR","UY","VE") & studyno == 4001 & S002 == 3) |
      # other countries from evs wave 2
      (S009 %in% c("CA","FR","IT","PT","NL","BE","DK",
                   "IS","IE","AT","HU","CZ","SK","RO") & studyno == 7503 & S002EVS == 2)
      ) %>%
    # keep only trad-values variables and occupation information
    select(
      # country code
      S009,
      # traditional vs. secular-rational values
      F063, F120, G006, E018, Y003,
      # survival self-expression values
      A008, A165, E025, F118, Y002,
      # occupation
      X036
    ) %>%
    # drop missings
    drop_na() %>%
    filter(Y003 != -5) %>%
    # summarise at country level
    group_by(S009) %>%
    summarise_if(is.numeric, mean) %>%
    ungroup()
  # get traditional vs. secular-rational values
  m <- pca(select(wvs, -S009, -X036), nfactors = 2, rotate = "varimax")
  # summarise at country level
  wvs %>%
    mutate(trad = as.numeric(scale(m$scores[,1]))) %>%
    # link iso data
    left_join(iso, by = c("S009" = "iso2")) %>%
    # link world bank data
    left_join(wb, by = "iso3") %>%
    # new var for modelling
    mutate(iso2lin = S009) %>%
    # drop missings
    drop_na(trad, indust)
}

# fit model to Inglehart 2000 data
fitModelInglehart2000 <- function(inglehartData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% inglehartData$iso2lin,
                   colnames(linCov) %in% inglehartData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(trad ~ indust)
  } else if (control == "spatial") {
    model <- bf(trad ~ indust + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(trad ~ indust + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(trad ~ indust + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  inglehartData %>%
    # standardise variables
    mutate_if(is.numeric, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Knack 1997
loadDataKnack1997 <- function(fileKnack, iso) {
  read_xlsx(fileKnack) %>%
    # sweden's ISO code is SE not SW
    mutate(ISO = ifelse(ISO == "SW", "SE", ISO)) %>%
    # link iso data
    left_join(iso, by = c("ISO" = "iso2")) %>%
    # new var for modelling
    mutate(iso2lin = ISO)
}

# fit model to Knack 1997 data
fitModelKnack1997 <- function(knackData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% knackData$iso2lin,
                   colnames(linCov) %in% knackData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(conf ~ Trust)
  } else if (control == "spatial") {
    model <- bf(conf ~ Trust + gp(Longitude..average., Latitude..average., k = 5, c = 5/4))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(conf ~ Trust + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(conf ~ Trust + gp(Longitude..average., Latitude..average., k = 5, c = 5/4) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  knackData %>%
    # rename variable
    rename(conf = `Confidence in government`) %>%
    # drop NAs manually
    drop_na(Trust, conf) %>%
    # standardise variables
    mutate_if(is.numeric, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# load data for Skidmore 2002
loadDataSkidmore2002 <- function(fileSkidmore1, fileSkidmore2, fileSkidmore3, fileSkidmore4, iso) {
  # individual files
  disasters <-
    read_xlsx(fileSkidmore1, skip = 6) %>%
    mutate(disaster = `Disaster Type` %in% c("Landslide", "Earthquake", "Volcanic activity", "Flood", "Storm")) %>%
    group_by(ISO) %>%
    summarise(totalDisasters = sum(disaster))
  econGrowth <-
    read_xls(fileSkidmore2) %>%
    filter(Year %in% 1960:1990) %>%
    group_by(Country) %>%
    mutate(lagRGDPL = lag(RGDPL),
           growthRateGDP = (RGDPL / lagRGDPL) - 1) %>%
    summarise(growthRateGDP = median(growthRateGDP, na.rm = TRUE))
  link <- 
    read.csv(fileSkidmore3)
  landArea <-
    read_csv(fileSkidmore4, skip = 4) %>%
    transmute(ISO = `Country Code`, landArea = `1990`)
  # put together
  econGrowth %>%
    left_join(link, by = "Country") %>%
    # east / west germany repeat
    group_by(ISO) %>%
    summarise(growthRateGDP = mean(growthRateGDP)) %>%
    # join disaster data
    left_join(disasters, by = "ISO") %>%
    # join land area data
    left_join(landArea, by = "ISO") %>%
    # fill in taiwan's missing land area: https://en.wikipedia.org/wiki/Geography_of_Taiwan
    mutate(landArea = ifelse(ISO == "TWN", 32260, landArea)) %>%
    # log variables
    mutate(
      logDisasters = log(totalDisasters + 1),
      logDisastersPerLandArea = log((totalDisasters / (landArea / 1e+06)) + 1)
    ) %>%
    # remove countries not in paper
    filter(!(ISO %in% c("AGO", "BHS", "BHR", "BLZ", "BEN", "BTN", "BGR", "BFA", 
                        "BDI", "CPV", "TCD", "CHN", "COM", "CZE", "DJI", "DMA", 
                        "EGY", "ETH", "GAB", "GMB", "GRD", "GIN", "GNB", "HUN", 
                        "CIV", "KWT", "LAO", "LUX", "MDG", "MLT", "MRT", "MNG", 
                        "MAR", "MMR", "NAM", "NGA", "OMN", "POL", "PRI", "QAT", 
                        "REU", "ROU", "RWA", "SAU", "SYC", "SLE", "SLB", "SOM", 
                        "KNA", "LCA", "VCT", "SDN", "SUR", "TZA", "TON", "RUS", 
                        "ARE", "VUT", "WSM", "YEM", "YUG", "ZAR"))) %>%
    # join iso codes
    left_join(iso, by = c("ISO" = "iso3")) %>%
    mutate(iso2lin = iso2)
}

# fit model to Skidmore 2002 data
fitModelSkidmore2002 <- function(skidmoreData, linCov, control) {
  # subset covariance matrix
  linCov <- linCov[rownames(linCov) %in% skidmoreData$iso2lin,
                   colnames(linCov) %in% skidmoreData$iso2lin]
  # initial priors
  priors <- c(prior(normal(0, 0.4), class = Intercept),
              prior(normal(0, 0.4), class = b),
              prior(exponential(5), class = sigma))
  # model formula and additional priors
  if (control == "none") {
    model <- bf(growthRateGDP ~ logDisastersPerLandArea)
  } else if (control == "spatial") {
    model <- bf(growthRateGDP ~ logDisastersPerLandArea + gp(Longitude..average., Latitude..average.))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  } else if (control == "cultural") {
    model <- bf(growthRateGDP ~ logDisastersPerLandArea + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
  } else if (control == "both") {
    model <- bf(growthRateGDP ~ logDisastersPerLandArea + gp(Longitude..average., Latitude..average.) + (1 | gr(iso2lin, cov = linCov)))
    priors <- c(priors, prior(exponential(5), class = sd))
    priors <- c(priors, prior(exponential(5), class = sdgp))
  }
  # begin pipe
  skidmoreData %>%
    # drop NAs manually
    drop_na(growthRateGDP, logDisastersPerLandArea) %>%
    # standardise variables
    mutate_if(is.numeric, function(x) as.numeric(scale(x))) %>%
    # fit model
    brm(model, data = ., data2 = list(linCov = linCov),
        cores = 4, seed = 2113, iter = 5000,
        prior = priors, control = list(adapt_delta = 0.999))
}

# plot replication results - densities
plotReplicationResults1 <- function(slopeList) {
  # type vector
  types <- c("No control", "Spatial control", "Cultural control", "Both controls")
  # plot
  out <-
    tibble(
      paper = rep(c("*Adamczyk and Pitt (2009)*<br>Homosexuality disapproval ~ Survival values",
                    "*Alesina et al. (2013)*<br>FLFP ~ Traditional plough use",
                    "*Beck et al. (2003)*<br>Stock market development ~ Settler mortality",
                    "*Beck et al. (2005)*<br>GDP ~ SME sector share",
                    "*Bockstette et al. (2002)*<br>GDP growth ~ State antiquity",
                    "*Easterly and Levine (2003)*<br>Log GDP ~ Institutional development",
                    "*Easterly (2007)*<br>Gini ~ Log wheat sugar ratio",
                    "*Fincher et al. (2008)*<br>Individualism ~ Hist. pathogen prevalence",
                    "*Gelfand et al. (2011)*<br>Tightness ~ Nature disaster vulnerability",
                    "*Inglehart and Baker (2000)*<br>Traditional values ~ % industrial sector",
                    "*Knack and Keefer (1997)*<br>Confidence in institutions ~ % trusting",
                    "*Skidmore and Toya (2002)*<br>GDP growth ~ Log number of natural disasters"), each = 4),
      review = rep(c(rep("Cultural values", 2), rep("Economic development", 5),
                     rep("Cultural values", 4), "Economic development"), each = 4),
      type = factor(rep(types, times = 12), levels = types),
      slope = slopeList,
    ) %>%
    mutate(review = factor(review, levels = c("Economic development", "Cultural values"))) %>%
    unnest(c(slope)) %>%
    ggplot(aes(x = slope, y = fct_rev(paper), fill = type)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(size = 1, position = position_dodge(width = -0.8)) +
    facet_grid(review ~ ., switch = "y", scales = "free_y") +
    scale_x_continuous(name = "Cross-cultural correlation", limits = c(-1.3, 1.3)) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_markdown(),
          legend.title = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", margin = margin(0, 10, 0, 0)))
  # save plot
  ggsave(out, filename = "figures/replications1.pdf", width = 8, height = 6.5)
  return(out)
}

# plot replication results - predictions
plotReplicationResults2 <- function(dataList, condList) {
  # plotting function
  plotFun <- function(data, cond1, cond2, cond3, cond4, x, y, xlab, ylab, title) {
    ggplot() +
      geom_text(data = data, aes(x = !!sym(x), y = !!sym(y), label = iso2lin), size = 2) +
      geom_ribbon(data = cond1, aes(x = effect1__, ymin = lower__, ymax = upper__), fill = "#F8766D", alpha = 0.25) +
      geom_line(data = cond1, aes(x = effect1__, y = estimate__), colour = "#F8766D") +
      geom_ribbon(data = cond2, aes(x = effect1__, ymin = lower__, ymax = upper__), fill = "#7CAE00", alpha = 0.25) +
      geom_line(data = cond2, aes(x = effect1__, y = estimate__), colour = "#7CAE00") +
      geom_ribbon(data = cond3, aes(x = effect1__, ymin = lower__, ymax = upper__), fill = "#00BFC4", alpha = 0.25) +
      geom_line(data = cond3, aes(x = effect1__, y = estimate__), colour = "#00BFC4") +
      geom_ribbon(data = cond4, aes(x = effect1__, ymin = lower__, ymax = upper__), fill = "#C77CFF", alpha = 0.25) +
      geom_line(data = cond4, aes(x = effect1__, y = estimate__), colour = "#C77CFF") +
      scale_x_continuous(name = xlab, limits = c(-3.5, 3.5), breaks = c(-2, 0, 2)) +
      scale_y_continuous(name = ylab, limits = c(-3.5, 3.5), breaks = c(-2, 0, 2)) +
      ggtitle(title) +
      theme_classic() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8))
  }
  # plots
  pA <- plotFun(dataList[[1]] %>% group_by(iso2lin) %>% summarise(countrySurv = mean(countrySurv), F118 = mean(F118)),
                condList[[1]], condList[[2]], condList[[3]], condList[[4]], "countrySurv", "F118", "Survival\nvalues (avg)", 
                "Homosexuality\ndisapproval (avg)", "Adamczyk and Pitt (2009)")
  pB <- plotFun(dataList[[2]], condList[[5]], condList[[6]], condList[[7]], condList[[8]],
                "tradPloughMean", "FLFP2000", "Traditional\nplough use", 
                "Female labour\nforce participation", "Alesina et al. (2013)")
  pC <- plotFun(dataList[[3]], condList[[9]], condList[[10]], condList[[11]], condList[[12]],
                "logSettlerMortality", "stockMarketDev", "Log settler\nmortality", 
                "Stock market\ndevelopment", "Beck et al. (2003)")
  pD <- plotFun(dataList[[4]], condList[[13]], condList[[14]], condList[[15]], condList[[16]],
                "SME250", "GDPpc", "SME sector\nshare", 
                "GDP per\ncapita", "Beck et al. (2005)")
  pE <- plotFun(dataList[[5]], condList[[17]], condList[[18]], condList[[19]], condList[[20]],
                "Statehist5", "gdpGrowth", "State\nantiquity", 
                "GDP per\ncapita growth", "Bockstette et al. (2002)")
  pF <- plotFun(dataList[[6]], condList[[21]], condList[[22]], condList[[23]], condList[[24]],
                "institIndex", "logGDPpc1995", "Institutional\ndevelopment", 
                "Log GDP\nper capita", "Easterly and Levine (2003)")
  pG <- plotFun(dataList[[7]], condList[[25]], condList[[26]], condList[[27]], condList[[28]],
                "LWHEATSUGAR", "gini", "Log wheat\nsugar ratio", 
                "Gini\ncoefficient", "Easterly (2007)")
  pH <- plotFun(dataList[[8]], condList[[29]], condList[[30]], condList[[31]], condList[[32]],
                "pathPrevHistorical", "individualismHofstede", "Historical pathogen\nprevalence", 
                " \nIndividualism", "Fincher et al. (2008)")
  pI <- plotFun(dataList[[9]], condList[[33]], condList[[34]], condList[[35]], condList[[36]],
                "DISCAS", "tightness", "Natural disaster\nvulnerability", 
                " \nTightness", "Gelfand et al. (2011)")
  pJ <- plotFun(dataList[[10]], condList[[37]], condList[[38]], condList[[39]], condList[[40]],
                "indust", "trad", "% employed in\nindustrial sector", 
                "Traditional\nvalues", "Inglehart and Baker (2000)")
  pK <- plotFun(dataList[[11]], condList[[41]], condList[[42]], condList[[43]], condList[[44]],
                "Trust", "conf", "% trusting\n ", 
                "Confidence in\ninstitutions", "Knack and Keefer (1997)")
  pL <- plotFun(dataList[[12]], condList[[45]], condList[[46]], condList[[47]], condList[[48]],
                "logDisastersPerLandArea", "growthRateGDP", "Log number of\nnatural disasters", 
                "GDP\ngrowth", "Skidmore and Toya (2002)")
  # put together
  top <- plot_grid(pC, pD, pE, pG, pF, pL, nrow = 2)
  bot <- plot_grid(pA, pB, pH, pI, pJ, pK, nrow = 2)
  out <- plot_grid(top, NULL, bot, ncol = 1, rel_heights = c(1, 0.1, 1),
                   labels = c("a", "", "b"))
  # add legend
  p <- 
    tibble(x = 1:4, y = 1:4, l = c("No control", "Spatial control", "Cultural control", "Both controls")) %>%
    mutate(l = factor(l, levels = .$l)) %>%
    ggplot(aes(x, y, colour = l)) + 
    geom_line() +
    theme_classic() +
    theme(legend.title = element_blank())
  out <- plot_grid(out, get_legend(p), rel_widths = c(1, 0.2))
  # save
  ggsave(out, filename = "figures/replications2.pdf", height = 8, width = 8)
  return(out)
}

# plot replication results - spatial autocorrelation plots
plotReplicationResults3 <- function(modelList, dataList) {
  # list of papers
  papers <- c("Adamczyk and Pitt (2009)", "Alesina et al. (2013)", 
              "Beck et al. (2003)", "Beck et al. (2005)", "Bockstette et al. (2002)", 
              "Easterly and Levine (2003)", "Easterly (2007)", "Fincher et al. (2008)", 
              "Gelfand et al. (2011)", "Inglehart and Baker (2000)", 
              "Knack and Keefer (1997)", "Skidmore and Toya (2002)")
  # plotting function
  plotFun <- function(model, data, title) {
    post <- posterior_samples(model, pars = c("lscale_gpLongitude..average.Latitude..average."))
    data <- data %>% group_by(iso2lin) %>% summarise_if(is.numeric, mean)
    maxDist <- max(as.matrix(distm(data[data$iso2lin %in% model$data$iso2lin, c("Longitude..average.","Latitude..average.")])))
    tibble(x = seq(50*1000, 50000*1000, length.out = 1000)) %>%
      mutate(
        median = map(x, function(x) median(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (x / maxDist)^2))),
        lower50 = map(x, function(x) quantile(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (x / maxDist)^2), 0.25)),
        upper50 = map(x, function(x) quantile(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (x / maxDist)^2), 0.75)),
        lower95 = map(x, function(x) quantile(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (x / maxDist)^2), 0.025)),
        upper95 = map(x, function(x) quantile(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (x / maxDist)^2), 0.975))
      ) %>%
      unnest(c(median, lower50, lower95, upper50, upper95)) %>%
      ggplot(aes(x = x, y = median)) +
      geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "grey") +
      geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = "darkgrey") +
      geom_line() +
      scale_x_log10(labels = function(x) x/1000) +
      theme_classic() +
      ggtitle(title) +
      labs(x = "", y = "")
  }
  # individual plots
  plotList <- list()
  for (i in 1:12) {
    plotList[[i]] <- plotFun(modelList[[i]], dataList[[i]], papers[i])
  }
  # put together
  out <-
    plot_grid(
      plot_grid(plotList[[3]] + ylab("Correlation"), plotList[[4]], plotList[[5]],
                plotList[[6]] + ylab("Correlation"), plotList[[7]] + xlab("Spatial distance (km)"), 
                plotList[[12]], nrow = 2),
      NULL,
      plot_grid(plotList[[1]] + ylab("Correlation"), plotList[[2]], plotList[[8]],
                plotList[[9]] + ylab("Correlation"), plotList[[10]] + xlab("Spatial distance (km)"), 
                plotList[[11]], nrow = 2),
      nrow = 3, rel_heights = c(1, 0.1, 1), labels = c("a", "", "b")
    )
  # save
  ggsave(out, filename = "figures/replications3.pdf", width = 8.2, height = 6)
  return(out)
}

# plot replication results - cultural autocorrelation plots
plotReplicationResults4 <- function(modelList) {
  # list of papers
  papers <- c("Adamczyk and Pitt (2009)", "Alesina et al. (2013)", 
              "Beck et al. (2003)", "Beck et al. (2005)", "Bockstette et al. (2002)", 
              "Easterly and Levine (2003)", "Easterly (2007)", "Fincher et al. (2008)", 
              "Gelfand et al. (2011)", "Inglehart and Baker (2000)", 
              "Knack and Keefer (1997)", "Skidmore and Toya (2002)")
  # calculate cultural phylogenetic signal
  getSignal <- function(model) {
    post <- posterior_samples(model, pars = c("sd_", "sigma"))
    total <- post$sd_iso2lin__Intercept^2 + post$sigma^2
    if ("sd_S009__Intercept" %in% names(post)) total <- post$sd_iso2lin__Intercept^2 + post$sd_S009__Intercept^2
    signal <- post$sd_iso2lin__Intercept^2 / total
    return(signal)
  }
  # plot
  out <-
    tibble(model = modelList, paper = papers,
           review = c(rep("Cultural values", 2), rep("Economic development", 5),
                      rep("Cultural values", 4), "Economic development")) %>%
    mutate(review = factor(review, levels = c("Economic development", "Cultural values"))) %>%
    mutate(signal = map(model, getSignal)) %>%
    select(-model) %>%
    unnest(c(signal)) %>%
    ggplot(aes(x = signal, y = fct_rev(paper))) +
    stat_halfeye(scale = 0.6, normalize = "groups") +
    facet_grid(review ~ ., switch = "y", scales = "free_y") +
    theme_classic() +
    xlab("Cultural phylogenetic signal") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", margin = margin(0, 10, 0, 0)))
  # save
  ggsave(out, filename = "figures/replications4.pdf", height = 5, width = 5)
  return(out)
}

# plot replication results - autocorrelation at attenuated effect sizes
plotReplicationResults5 <- function(modelListD, modelListA, dataList) {
  # list of papers
  papers <- c("Adamczyk and Pitt (2009)", "Alesina et al. (2013)", 
              "Beck et al. (2003)", "Beck et al. (2005)", "Bockstette et al. (2002)", 
              "Easterly and Levine (2003)", "Easterly (2007)", "Fincher et al. (2008)", 
              "Gelfand et al. (2011)", "Inglehart and Baker (2000)", 
              "Knack and Keefer (1997)", "Skidmore and Toya (2002)")
  # calculate spatial autocorrelation at 1000km distance
  calculateSpatial1000 <- function(model, data) {
    post <- posterior_samples(model, pars = c("lscale_gpLongitude..average.Latitude..average."))
    data <- data %>% group_by(iso2lin) %>% summarise_if(is.numeric, mean)
    maxDist <- max(as.matrix(distm(data[data$iso2lin %in% model$data$iso2lin, c("Longitude..average.","Latitude..average.")])))
    out <- median(exp(-(1 / (2 * post$lscale_gpLongitude..average.Latitude..average.^2)) * (1000*1000 / maxDist)^2))
    return(out)
  }
  # calculate cultural phylogenetic signal
  calculateSignal <- function(model, data) {
    post <- posterior_samples(model, pars = c("sd_", "sigma"))
    total <- post$sd_iso2lin__Intercept^2 + post$sigma^2
    if ("sd_S009__Intercept" %in% names(post)) total <- total + post$sd_S009__Intercept^2
    signal <- median(post$sd_iso2lin__Intercept^2 / total)
    return(signal)
  }
  # ratio of original to "control" effect size
  calculateRatioEffectSize <- function(modelA, modelD) {
    if (nrow(fixef(modelA)) == 29) {
      original <- fixef(modelA)[29,1]
      control  <- fixef(modelD)[29,1]
    } else {
      original <- fixef(modelA)[2,1]
      control  <- fixef(modelD)[2,1]
    }
    ratio <- control / original
    return(ratio)
  }
  # data
  d <-
    tibble(modelD = modelListD, modelA = modelListA, data = dataList, paper = papers) %>%
    mutate(corAt1000km = map2(modelD, data, calculateSpatial1000),
           culturalSignal = map2(modelD, data, calculateSignal),
           effRatio = map2(modelA, modelD, calculateRatioEffectSize)) %>%
    unnest(c(corAt1000km, culturalSignal, effRatio)) %>%
    select(-starts_with("model"), -data)
  # plots
  pA <-
    ggplot(d, aes(x = corAt1000km, y = effRatio, label = paper)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", alpha = 0.3) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_text_repel(size = 2.4) +
    labs(x = "Spatial autocorrelation at 1000km",
         y = "Ratio of effect size with control\nto original effect size") +
    ylim(c(0, 1.5)) +
    theme_classic()
  pB <-
    ggplot(d, aes(x = culturalSignal, y = effRatio, label = paper)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", alpha = 0.3) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_text_repel(size = 2.4) +
    labs(x = "Cultural phylogenetic signal", y = " \n ") +
    ylim(c(0, 1.5)) +
    theme_classic()
  # put together
  out <- plot_grid(pA, pB, nrow = 1, labels = c("a","b"))
  # save
  ggsave(out, filename = "figures/replications5.pdf", width = 8, height = 4)
  return(out)
}

getSpatialAutocorrelation1000km <- function(data, model) {
  post <- posterior_samples(model, pars = c("lscale_gpLongitude..average.Latitude..average."))
  maxDist <- max(as.matrix(distm(data[,c("Longitude..average.","Latitude..average.")])))
  sa <- exp(-(1 / (2 * post[,1]^2)) * (1000*1000 / maxDist)^2)
  return(sa)
}

getCulturalPhylogeneticSignal <- function(model) {
  post <- posterior_samples(model)
  signalG <- post$sd_iso2lin__Intercept^2 / (post$sd_iso2lin__Intercept^2 + post$sigma^2)
}
