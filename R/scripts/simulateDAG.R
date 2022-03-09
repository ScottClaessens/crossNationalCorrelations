library(brms)
library(tidyverse)

# simulate DAG to show that method works

# DAG:
# X -> Y
# U -> X
# U -> X
# D -> U
# where X = exposure, Y = outcome, U = unobserved confound
# and D is a distance matrix

# set seed
set.seed(1)

# sample size
n <- 20^2

# coordinates
lon <- rep(seq(0, 1, length.out = sqrt(n)), each = sqrt(n))
lat <- rep(seq(0, 1, length.out = sqrt(n)), times = sqrt(n))

# distance matrix
dist <- as.matrix(dist(data.frame(lon = lon, lat = lat)))
dist <- dist / max(dist)

# strength of autocorrelation
lambda <- 0.95

# generate spatially autocorrelated unobserved confound
U <- ((1 - lambda*dist)^-1 %*% rnorm(n)) + rnorm(n)
U <- as.numeric(scale(U)) # standardise

# U positively predicts X and negatively predicts Y
X <- as.numeric(scale(rnorm(n, 1*U)))
Y <- as.numeric(scale(rnorm(n, -1*U)))

# creates spurious bivariate relationship between X and Y
m1 <- lm(Y ~ X)
summary(m1)

# residuals are spatially autocorrelated
r1 <- resid(m1)
tibble(lat, lon, r1) %>%
  ggplot(aes(lat, lon, colour = r1)) +
  geom_point(size = 5) +
  scale_color_gradient2(high = "red", low = "blue") +
  theme_classic()

# control for confound U removes Y ~ X relationship
summary(lm(Y ~ X + U))

# but we can't control for U as it is unobserved!
# instead, control for spatial autocorrelation
# with approximate gaussian process
m2 <- brm(Y ~ X + gp(lat, lon, k = 5, c = 5/4),
          data = tibble(Y, X, lat, lon),
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          cores = 4, seed = 1)
summary(m2)

# removes association between X and Y by accounting for
# the covariation between X and Y induced by the
# spatially autocorrelated confound U
