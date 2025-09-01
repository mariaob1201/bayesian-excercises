library("car")  # load the 'car' package
data("Anscombe")  # load the data set
?Anscombe  # read a description of the data
head(Anscombe)  # look at the first few lines of the data
pairs(Anscombe)  # scatter plots for each pair of variables

#######Bayesian linear regression model q4
library(LearnBayes)

y <- Anscombe$education
X <- model.matrix(~ income + young + urban, data = Anscombe)
fit <- blinreg(y, X, m = 100000)
fit$betadraw
post_mean_intercept <- local({
  if (!is.null(fit$betadraw) && is.matrix(fit$betadraw)) {
    colMeans(fit$betadraw)[1]
  } else if (!is.null(fit$beta) && is.matrix(fit$beta)) {
    colMeans(fit$beta)[1]
  } else if (!is.null(fit$coeff)) {
    fit$coeff[1]
  } else {
    # Fallback to OLS (correct for reference prior)
    coef(lm(education ~ income + young + urban, data = Anscombe))[1]
  }
})

round(unname(post_mean_intercept), 1)

##no independence apparently
plot(lm(education ~ income + young + urban, data = Anscombe))

############### q6
install.packages("MCMCpack")
library(MCMCpack)
fit <- MCMCregress(education ~ income + young + urban, data = Anscombe,
                   b0 = 0, B0 = 1e-6, c0 = 0.001, d0 = 0.001)
round(mean(fit[, "(Intercept)"]), 1)  # posterior mean intercept