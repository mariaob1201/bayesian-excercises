# Load required libraries
library("car")
library("MCMCpack")

# Load the Anscombe data
data("Anscombe")

# Fit the Bayesian linear regression model using MCMCpack
# Use a large number of iterations for reliable DIC calculation
set.seed(123)  # for reproducibility
mcmc_fit <- MCMCregress(
  education ~ income + young + urban, 
  data = Anscombe,
  b0 = 0,           # prior mean for coefficients
  B0 = 1e-6,        # prior precision for coefficients (diffuse prior)
  c0 = 0.001,       # prior shape parameter for error variance
  d0 = 0.001,       # prior scale parameter for error variance
  mcmc = 100000,    # large number of MCMC samples
  burnin = 10000,   # burn-in samples
  thin = 1          # no thinning
)

# MCMCpack doesn't have a built-in DIC function, so we'll calculate it manually
# DIC = D_bar + p_D, where D_bar is posterior mean deviance and p_D is effective number of parameters

# Extract the posterior samples
beta_samples <- mcmc_fit[, 1:4]  # coefficients (intercept + 3 predictors)
sigma2_samples <- mcmc_fit[, "sigma2"]  # error variance

# Prepare data
y <- Anscombe$education
X <- model.matrix(~ income + young + urban, data = Anscombe)
n <- length(y)

# Calculate deviance for each MCMC sample
deviances <- numeric(nrow(mcmc_fit))

for (i in 1:nrow(mcmc_fit)) {
  # Predicted values for this sample
  mu_pred <- X %*% as.numeric(beta_samples[i, ])
  
  # Log-likelihood for this sample
  log_lik <- sum(dnorm(y, mean = mu_pred, sd = sqrt(sigma2_samples[i]), log = TRUE))
  
  # Deviance = -2 * log-likelihood
  deviances[i] <- -2 * log_lik
}

# Calculate DIC components
D_bar <- mean(deviances)  # Posterior mean deviance

# Calculate deviance at posterior means
beta_hat <- colMeans(beta_samples)
sigma2_hat <- mean(sigma2_samples)
mu_hat <- X %*% beta_hat
log_lik_hat <- sum(dnorm(y, mean = mu_hat, sd = sqrt(sigma2_hat), log = TRUE))
D_hat <- -2 * log_lik_hat

# Effective number of parameters
p_D <- D_bar - D_hat

# DIC
DIC <- D_bar + p_D

# Print results
cat("Posterior mean deviance (D_bar):", round(D_bar, 2), "\n")
cat("Deviance at posterior means (D_hat):", round(D_hat, 2), "\n")
cat("Effective number of parameters (p_D):", round(p_D, 2), "\n")
cat("DIC:", round(DIC, 2), "\n")
cat("DIC (rounded to nearest whole number):", round(DIC), "\n")

# Summary of the model fit
# Print results in JAGS-style format
cat("Mean deviance:", round(D_bar, 2), "\n")
cat("penalty:", round(p_D, 2), "\n")
cat("Penalized deviance:", round(DIC, 2), "\n")
cat("\nDIC (Penalized deviance rounded to nearest whole number):", round(DIC), "\n")

# Summary of the model fit
summary(mcmc_fit)

############## posterior distribution
# Load required libraries
library("car")
library("MCMCpack")

# Load the Anscombe data
data("Anscombe")

# Fit the Bayesian linear regression model using MCMCpack
# (This is the original model from your code)
set.seed(123)  # for reproducibility
mcmc_fit <- MCMCregress(
  education ~ income + young + urban, 
  data = Anscombe,
  b0 = 0,           # prior mean for coefficients
  B0 = 1e-6,        # prior precision for coefficients (diffuse prior)
  c0 = 0.001,       # prior shape parameter for error variance
  d0 = 0.001,       # prior scale parameter for error variance
  mcmc = 100000,    # large number of MCMC samples
  burnin = 10000,   # burn-in samples
  thin = 1          # no thinning
)

# Extract the posterior samples for the income coefficient
income_coeff_samples <- mcmc_fit[, "income"]

# Calculate the posterior probability that income coefficient > 0
prob_positive <- mean(income_coeff_samples > 0)

# Print results
cat("Summary statistics for income coefficient:\n")
cat("Posterior mean:", round(mean(income_coeff_samples), 4), "\n")
cat("Posterior SD:", round(sd(income_coeff_samples), 4), "\n")
cat("95% Credible interval:", round(quantile(income_coeff_samples, c(0.025, 0.975)), 4), "\n")
cat("\nNumber of posterior samples:", length(income_coeff_samples), "\n")
cat("Number of samples where income coeff > 0:", sum(income_coeff_samples > 0), "\n")
cat("Posterior probability that income coefficient > 0:", round(prob_positive, 4), "\n")
cat("Answer (rounded to 2 decimal places):", round(prob_positive, 2), "\n")

# Optional: Create a histogram to visualize the posterior distribution
hist(income_coeff_samples, 
     main = "Posterior Distribution of Income Coefficient",
     xlab = "Income Coefficient",
     breaks = 50,
     col = "lightblue",
     border = "white")
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = mean(income_coeff_samples), col = "blue", lwd = 2)
legend("topright", 
       c("Zero line", "Posterior mean"), 
       col = c("red", "blue"), 
       lty = c(2, 1), 
       lwd = 2)