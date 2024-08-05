# Functions for modelling
# Custom function to calculate AIC for multivariate models
multivariate_AIC <- function(model) {
  n <- nrow(model$fitted.values)  # Number of observations
  p <- ncol(model$fitted.values)  # Number of response variables
  
  # Calculate residual covariance matrix
  residuals_matrix <- as.matrix(residuals(model))
  sigma_hat <- cov(residuals_matrix)
  
  # Calculate log-likelihood for multivariate normal distribution
  logL <- -0.5 * n * (p * log(2 * pi) + log(det(sigma_hat)) + p)
  
  # Calculate AIC
  aic <- -2 * logL + 2 * length(coef(model))
  return(aic)
}
