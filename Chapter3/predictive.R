library(readr); library(coda); library(rjags)

dat <- read.csv("callers.csv")

model_string <- "
model {
  for (i in 1:N) {
    calls[i] ~ dpois(mu[i])
    log(mu[i]) <- b0 + b1*age[i] + b2*isgroup2[i] + log(days_active[i])  # offset
  }
  b0 ~ dnorm(0, 1.0E-6)
  b1 ~ dnorm(0, 1.0E-6)
  b2 ~ dnorm(0, 1.0E-6)
}
"

data_jags <- list(
  calls = dat$calls,
  age = dat$age,
  isgroup2 = dat$isgroup2,
  days_active = dat$days_active,
  N = nrow(dat)
)

set.seed(1)
jm <- jags.model(textConnection(model_string), data = data_jags,
                 n.chains = 3, n.adapt = 1000)
update(jm, 2000)

samp <- coda.samples(jm, c("b0","b1","b2"), n.iter = 10000, thin = 5)

# Posterior predictive for age=29, group2=1, 30 days active:
A   <- as.matrix(samp)
eta <- A[, "b0"] + A[, "b1"]*29 + A[, "b2"]*1 + log(30)
mu  <- exp(eta)
ypp <- rpois(length(mu), mu)
mean(ypp >= 3)
