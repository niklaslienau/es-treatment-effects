

# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)
library(quadprog)
library(esreg)
library(dplyr)
#----------
simulate_dgp_lin <- function(n = 1000, beta = 1, gamma = 0.5, pi = 1, rho = 0.7, cont = TRUE,  error_var= 1) {
  # Simulate instrument
  Z <- runif(n, min=0, max=2)
  
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
  
  # Additional independent noise
  epsilon <- rnorm(n, mean = 0, sd = error_var)
  
  #structural error
  v =  (gamma * abs(D)) * u + epsilon
  
  # Linear location scale: sigma(D) = gamma * D
  Y <- beta * D + v
  
    
  return(data.frame(Y = Y, D = D, Z = Z, u = u, e = e, v=v, epsilon=epsilon ))
}
##### no hetero  Option #####
sim_dgp_homo = function(n= 1000, beta =1 , pi=1, rho=0.5){
  
  # Simulate correlated errors (u, e)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  
  #Simulate Z from a runiform dist U(0,2)
  Z= runif(n, 0,2)
  
  #Generate D
  D = pi*Z + e
  #Generate Y
  Y= beta*D + u
  
  return(data.frame(Y = Y, D = D, Z = Z, u = u, e = e ))
  
}

#-----------


#### random coefficient model #######

#sim data
simulate_dgp = function(n=1000, pi=1, rho=0.5){
  ## DGP
  # Simulate correlated errors (u, e)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  
  #random coefficent
  v = rnorm(n, 0, 1)
  
  #Simulate Z from a normal dist
  Z= rnorm(n, 0,1)
  
  #Generate D
  D = 1 + pi*Z + e
  #Generate Y with random coefficent beta(v) = v
  Y= v*D + u
  
  data=data.frame(Y = Y, D = D, Z = Z, u = u, e = e )
  
  return(data)
  
  ## True Quantile Effects
  
}
#calc true qte (analytical closed form not available)
compute_qte_tau <- function(tau, pi = 1, n = 100000) {
  # Draw D ~ N(0, pi^2 + 1)
  D <- rnorm(n, mean = 1, sd = sqrt(pi^2 + 1))
  
  # Compute QTE(tau, d) pointwise
  QTE_d <- D / sqrt(D^2 + 1)
  
  # Multiply by Phi^-1(tau) to get QTE(tau)
  qte_tau <- qnorm(tau) * mean(QTE_d)
  
  return(qte_tau)
}


#Plot true QTES
taus <- seq(0.01, 0.99, by = 0.01)
qtes <- sapply(taus, compute_qte_tau)
plot(taus, qtes, type = "l", main = "True QTE(τ)", xlab = "τ", ylab = "QTE(τ)")


