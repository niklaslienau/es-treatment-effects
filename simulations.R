
### MONTE CARLO: Naive QR vs CF Method across quantile values


# Simulate from location shift model
simulate_loc_shift <- function(n = 1000, beta = 1, pi = 1, rho = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z <- rnorm(n)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  errors <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  u <- errors[, 1]
  e <- errors[, 2]
  D <- pi * Z + e
  Y <- beta * D + u
  data.frame(Y = Y, D = D, Z = Z)
}

# Control function quantile estimator
cf_qr_estimate <- function(Y, D, Z, tau = 0.25, degree = 3) {
  first_stage <- rq(D ~ Z, tau = tau)
  e_hat <- resid(first_stage)
  series_terms <- sapply(1:degree, function(p) e_hat^p)
  colnames(series_terms) <- paste0("resid_", 1:degree)
  series_df <- as.data.frame(series_terms)
  design_df <- cbind(D = D, series_df)
  formula_str <- paste("Y ~ D +", paste(colnames(series_df), collapse = " + "))
  second_stage <- rq(as.formula(formula_str), data = cbind(Y = Y, design_df), tau = tau)
  return(second_stage$coefficients["D"])
}







# Initialize output data frame
mse_results <- data.frame(
  tau = quantiles,
  naive_mse = NA, naive_bias = NA, naive_var = NA,
  cf_mse = NA, cf_bias = NA, cf_var = NA
)

for (i in seq_along(quantiles)) {
  tau <- quantiles[i]
  naive_estimates <- numeric(R)
  cf_estimates <- numeric(R)
  
  for (r in 1:R) {
    df <- simulate_loc_shift(seed = r)
    
    # Naive QR with error check
    naive_est <- try(rq(Y ~ D, data = df, tau = tau), silent = TRUE)
    naive_estimates[r] <- if (inherits(naive_est, "try-error")) NA else naive_est$coefficients["D"]
    
    # Control function QR with error check
    cf_try <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau, degree = 3), silent = TRUE)
    cf_estimates[r] <- if (inherits(cf_try, "try-error") || is.na(cf_try)) NA else cf_try
  }
  
  # Naive metrics
  mse_results$naive_mse[i]  <- mean((naive_estimates - true_beta)^2, na.rm = TRUE)
  mse_results$naive_bias[i] <- mean(naive_estimates - true_beta, na.rm = TRUE)
  mse_results$naive_var[i]  <- var(naive_estimates, na.rm = TRUE)
  
  # Control function metrics
  mse_results$cf_mse[i]  <- mean((cf_estimates - true_beta)^2, na.rm = TRUE)
  mse_results$cf_bias[i] <- mean(cf_estimates - true_beta, na.rm = TRUE)
  mse_results$cf_var[i]  <- var(cf_estimates, na.rm = TRUE)
}


# Display result
print(mse_results)