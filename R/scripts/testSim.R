library(brms)

n <- 200

# create proximity matrix
distMat <- as.matrix(dist(rnorm(n)))
distMat <- distMat / max(distMat)
rownames(distMat) <- colnames(distMat)
proxMat <- 1 - distMat

# mock data
d <- data.frame(y = rnorm(n), x = rnorm(n), id = 1:n)

# signal for each variable
signalY <- 0.9
signalX <- 0.9
# true effect
r <- 0.8

# calculate sd and sigma values, assuming total variance sd^2 + sigma^2 = 1
# signal = sd^2 / sd^2 + sigma^2
# signal = sd^2 / 1
# sd = sqrt(signal)
sdY <- sqrt(signalY)
sdX <- sqrt(signalX)
# sd^2 + sigma^2 = 1
# sigma = sqrt(1 - sd^2)
sigmaY <- sqrt(1 - sdY^2)
sigmaX <- sqrt(1 - sdX^2)

# "data"
d <- data.frame(y = rnorm(n), x = rnorm(n), id = 1:n)

# simulation
mSim <- brm(bf(y ~ 0 + (1 | gr(id, cov = proxMat))) +
            bf(x ~ 0 + (1 | gr(id, cov = proxMat))) + set_rescor(TRUE),
            data = d, data2 = list(proxMat = proxMat),
            prior = c(prior_string(paste0("constant(", sdY, ")"), class = "sd", resp = "y"),
                      prior_string(paste0("constant(", sigmaY, ")"), class = "sigma", resp = "y"),
                      prior_string(paste0("constant(", sdX, ")"), class = "sd", resp = "x"),
                      prior_string(paste0("constant(", sigmaX, ")"), class = "sigma", resp = "x"),
                      prior_string(paste0("constant(cholesky_decompose([[1,", r, "],[", r, ",1]]))"), 
                                   class = "rescor")),
            sample_prior = "only", iter = 1, chains = 1)

# simulate data
out <- posterior_predict(mSim)

# standardise y and x
y <- as.numeric(scale(out[,,"y"]))
x <- as.numeric(scale(out[,,"x"]))

# data frame
d <- data.frame(y, x, id = 1:n)

# raw correlation
cor(d$y, d$x)

# get signal back
m1 <- brm(y ~ 1 + x + (1 | gr(id, cov = proxMat)),
          data = d, data2 = list(proxMat = proxMat),
          prior = c(prior(normal(0, 0.5), class = Intercept),
                    prior(normal(0, 0.5), class = b),
                    prior(exponential(1), class = sd),
                    prior(exponential(1), class = sigma)),
          cores = 4)
m2 <- brm(x ~ 1 + (1 | gr(id, cov = proxMat)),
          data = d, data2 = list(proxMat = proxMat),
          prior = c(prior(normal(0, 0.5), class = Intercept),
                    prior(exponential(1), class = sd),
                    prior(exponential(1), class = sigma)),
          cores = 4)

# signal estimates
post1 <- posterior_samples(m1)
quantile(post1$sd_id__Intercept^2 / (post1$sd_id__Intercept^2 + post1$sigma^2),
         c(0.025, 0.5, 0.975))
post2 <- posterior_samples(m2)
quantile(post2$sd_id__Intercept^2 / (post2$sd_id__Intercept^2 + post2$sigma^2),
         c(0.025, 0.5, 0.975))
