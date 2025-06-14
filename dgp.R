

# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)


simulate_dgp <- function(n = 1000, beta = 1, pi = 1, gamma = 0.5, rho = 0.7, 
                         scale = TRUE, cont = TRUE) {
  # Instrument Z
  Z <- rnorm(n)
  
  # Jointly normal errors (u, e) with correlation rho
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  
  # Latent treatment variable (always continuous)
  D_latent <- pi * Z + e
  
  # Convert to treatment D
  if (cont) {
    D <- D_latent
  } else {
    # Discretize D_latent into 16 bins labeled 0 to 15 (schooling style)
    thresholds <- quantile(D_latent, probs = seq(0, 1, length.out = 17))
    D <- cut(D_latent, breaks = thresholds, labels = 0:15, include.lowest = TRUE)
    D <- as.numeric(as.character(D))  # Convert factor to numeric
  }
  
  # Outcome equation
  if (scale) {
    Y <- beta * D + exp(gamma * D) * u
  } else {
    Y <- beta * D + u
  }
  
  return(data.frame(Y = Y, D = D, Z = Z))
}



