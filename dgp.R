
# Load libraries
library(MASS)
library(ggplot2)
library(quantreg)
library(quadprog)
library(esreg)
library(dplyr)
library(np)



#### random coefficient model #######

#sim data
simulate_dgp_rand = function(n=1000, pi=1, rho=0.5){
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

#compute true qte 
compute_qte_tau <- function(tau, pi = 1) {
  g <- function(d) (d / sqrt(d^2 + 1)) * dnorm(d, mean = 1, sd = sqrt(pi^2 + 1))
  ex <- integrate(g, lower = -Inf, upper = Inf, rel.tol = 1e-9)$value
  
  # Compute QTE(tau, d) pointwise
  qte_tau <- qnorm(tau) * ex
  

  return(qte_tau)
}

#Plot true QTES
taus <- seq(0.001, 0.999, by = 0.01)
qtes <- sapply(taus, compute_qte_tau)
plot(taus, qtes, type = "l", main = "True QTE(τ)", xlab = "τ", ylab = "QTE(τ)")

#compute the true CTATE (closed form available, could also integrate over QTE's)


 compute_ctate_tau<- function(tau, pi = 1) {
   z <- qnorm(tau)
   g <- function(d) (d / sqrt(d^2 + 1)) * dnorm(d, mean = 1, sd = sqrt(pi^2 + 1))
   ex <- integrate(g, lower = -Inf, upper = Inf, rel.tol = 1e-9)$value
   -(dnorm(z)/tau) * ex
 }
 taus <- seq(0.001, 0.999, by = 0.01)
 qtes <- sapply(taus, compute_ctate_tau)
 plot(taus, qtes, type = "l", main = "True CTATE(τ)", xlab = "τ", ylab = "CTATE(τ)")

 
########################################################################################
 #---------- Location Scale
 #sim data
 simulate_dgp_locscale <- function(n = 1000, beta = 1, gamma = 0.5, pi = 1, rho = 0.5) {
   # Simulate instrument
   Z <- runif(n, min=0, max=2)
   
   # Simulate correlated errors (u, e)
   Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
   errors <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
   u <- errors[, 1]
   e <- errors[, 2]
   
   # Generate endogenous treatment
   D <- pi * Z + e
   
   
   # Additional independent noise
   # epsilon <- rnorm(n, mean = 0, sd = error_var)
   
   #structural error
   v =  (gamma * abs(D)) * u 
   
   # Linear location scale: sigma(D) = gamma * D
   Y <- beta * D + v
   
   
   return(data.frame(Y = Y, D = D, Z = Z, u = u, e = e, v=v ))
 }
 
 
 #compute true qte
 comp_qte_locscale = function(tau, beta=1, gamma=0.5){
   qte = beta + gamma* qnorm(tau)
   return(qte)
 }
 #Plot true QTES
 taus <- seq(0.01, 0.99, by = 0.01)
 qtes <- sapply(taus, comp_qte_locscale)
 plot(taus, qtes, type = "l", main = "True QTE(τ)", xlab = "τ", ylab = "QTE(τ)")
 
 








