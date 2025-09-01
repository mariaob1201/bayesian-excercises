### Please visit the resources of JAGS, an introduction
# https://d3c33hcgiwev3.cloudfront.net/_4ec003e4706226af504d524aafb2c527_JAGSintro.html?Expires=1756771200&Signature=OdYkCeDjRpfgsmb-hrB~ae7gjZTlq5-eqo8M39pt9B9BTSlGxfyPzCygkRfd5-dbe0T5O3ZF-4ZKRujTtrjndGgpOQr-fSpeamvQ49JhdKP4h7F-jVBInqsGL2ARoNjAyQQYVZvKW82Rw98Q2fOrS3ksBWLqVaN~DzSOI-pAQrg_&Key-Pair-Id=APKAJLTNE6QMUY6HBC5A

library("rjags")

mod_string = " model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu, 1.0/sig2)
  }
  mu ~ dt(0.0, 1.0/1.0, 1.0) # location, inverse scale, degrees of freedom
  sig2 = 1.0
} "

## setting up the model
set.seed(50)
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
n = length(y)

data_jags = list(y=y, n=n)
params = c("mu")

inits = function() {
  inits = list("mu"=0.0)
} # optional (and fixed)

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits)

## run mcmc sampler
update(mod, 500) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=1000)

## post processing
summary(mod_sim)
library("coda")
plot(mod_sim)


############################### Example from callers
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
