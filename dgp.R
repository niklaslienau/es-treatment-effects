

# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)
library(quadprog)
library(esreg)


# gamma = heteroskedasticity coefficient
# pi = Reduced form coefficent
# rho = Error correlation that induces endogeneity
#cont = D is continous (default) as opposed to multivalued discrete
#error_var = variance of second normal error 

simulate_dgp <- function(n = 1000, beta = 1, gamma = 0.5, pi = 1, rho = 0.7, cont = TRUE,  error_var= 1) {
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
  
  # Additional independent noise
  epsilon <- rnorm(n, mean = 0, sd = error_var)
  
  #structural error
  v =  (gamma * abs(D)) * u + epsilon
  
  # Linear location scale: sigma(D) = gamma * D
  Y <- beta * D + v
  
    
  return(data.frame(Y = Y, D = D, Z = Z, u = u, e = e, v=v, epsilon=epsilon ))
}


##### EASY Option #####

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
