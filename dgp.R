# Load packages
library(MASS)
library(ggplot2)

# Set seed
#set.seed(42)

# Parameters
n <- 1000
beta <- 1
pi <- 1
gamma <- 0.5
rho <- 0.7

# Simulate Z ~ N(0, 1)
Z <- rnorm(n)

# Simulate bivariate normal errors (u, e) with correlation rho
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
errors <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
u <- errors[, 1]
e <- errors[, 2]

# Generate X and Y
X <- pi * Z + e
Y <- beta * X + exp(gamma * X) * u

# Put in data frame
df <- data.frame(X = X, Y = Y)

# Create grid of X values for plotting quantile functions
x_grid <- seq(min(X), max(X), length.out = 200)

# True quantile function values
qnorm_0.1 <- qnorm(0.1)
qnorm_0.5 <- qnorm(0.5)

q_0.1 <- beta * x_grid + exp(gamma * x_grid) * qnorm_0.1
q_0.5 <- beta * x_grid + exp(gamma * x_grid) * qnorm_0.5

# Quantile function dataframe
q_df <- data.frame(
  X = rep(x_grid, 2),
  Quantile = c(q_0.1, q_0.5),
  Tau = factor(rep(c("τ = 0.1", "τ = 0.5"), each = length(x_grid)))
)

# Plot: simulated data + quantile curves
ggplot(df, aes(x = X, y = Y)) +
  geom_point(alpha = 0.4, color = "grey40") +
  geom_line(data = q_df, aes(x = X, y = Quantile, color = Tau), size = 1) +
  labs(title = "Simulated Data with True Conditional Quantile Functions",
       x = "X", y = "Y", color = "Quantile Level") +
  theme_minimal()


########################
# Load packages
library(MASS)

# Set parameters
set.seed(123)
n <- 1000
R <- 50  # number of repetitions

# DGP parameters
beta <- 1
pi <- 1
gamma <- 0.5
rho <- 0.7
q_0.1_val <- qnorm(0.1)

# Initialize storage
below_counts <- numeric(R)

for (r in 1:R) {
  # Simulate instrument Z
  Z <- rnorm(n)
  
  # Simulate errors (u, e)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  
  # Generate X and Y
  X <- pi * Z + e
  Y <- beta * X + exp(gamma * X) * u
  
  # Compute true 10% quantile of Y|X for each X
  Q_0.1 <- beta * X + exp(gamma * X) * q_0.1_val
  
  # Check how many Y_i fall below Q_{Y|X_i}(0.1)
  below_counts[r] <- mean(Y < Q_0.1)
}

# Compute average across all repetitions
avg_below <- mean(below_counts)

# Print result
cat("Average percent of points below true 10% quantile line:", round(avg_below * 100, 2), "%\n")

