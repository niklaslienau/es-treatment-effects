

# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)


simulate_dgp <- function(n = 1000, beta = 1, gamma = 0.5, pi = 1, rho = 0.7, 
                         scale = TRUE, cont = TRUE, linear = FALSE) {
  # Simulate instrument
  Z <- rnorm(n)
  
  # Simulate correlated errors (u, e)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  
  # Generate endogenous treatment
  D_latent <- pi * Z + e
  
  # Handle discrete vs continuous treatment
  if (!cont) {
    thresholds <- quantile(D_latent, probs = seq(0, 1, length.out = 17))
    D <- cut(D_latent, breaks = thresholds, labels = 0:15, include.lowest = TRUE)
    D <- as.numeric(as.character(D))
  } else {
    D <- D_latent
  }
  
  # Outcome equation
  if (scale) {
    if (linear) {
      # Linear location scale: sigma(D) = gamma * D
      Y <- beta * D + gamma * D * u
    } else {
      # Nonlinear scale case (e.g. sigma(D) = exp(gamma * D))
      Y <- beta * D + exp(gamma * D) * u
    }
  } else {
    # Location shift only (homoskedastic errors)
    Y <- beta * D + u
  }
  
  return(data.frame(Y = Y, D = D, Z = Z, u = u, e = e))
}



