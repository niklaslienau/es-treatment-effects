

# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)

# Function to simulate  from location SCALE 
simulate_loc_scale <- function(n = 1000, beta = 1, pi = 1, gamma = 0.5, rho = 0.7) {
  Z <- rnorm(n)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  D <- pi * Z + e
  Y <- beta * D + exp(gamma * D) * u
  data.frame(Y = Y, D = D, Z = Z)
}


# Function to simulate from a location shift model (homoskedastic errors)
simulate_loc_shift <- function(n = 1000, beta = 1, pi = 1, rho = 0.7) {
  Z <- rnorm(n)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  D <- pi * Z + e
  Y <- beta * D + u
  data.frame(Y = Y, D = D, Z = Z)
}







