
### MONTE CARLO for QTE: Naive QR vs CF Method on MSE, Bias & Variance

evaluate_qte_estimators <- function(
    quantiles = c(0.01, 0.05, 0.25, 0.5),
    R = 100,
    n = 1000,
    beta_true = 1,
    degree = 3,
    seed = 123
) {
  set.seed(seed)
  
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
      # Location shift DGP
      df <- simulate_loc_shift(n = n, seed = r)
      
      # Naive QR
      naive_fit <- try(rq(Y ~ D, data = df, tau = tau), silent = TRUE)
      naive_estimates[r] <- if (inherits(naive_fit, "try-error")) NA else naive_fit$coefficients["D"]
      
      # CF QR
      cf_try <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau, degree = degree), silent = TRUE)
      cf_estimates[r] <- if (inherits(cf_try, "try-error") || !is.numeric(cf_try) || length(cf_try) != 1 || is.na(cf_try)) NA else cf_try
    }
    
    # Store performance metrics
    mse_results$naive_mse[i]  <- mean((naive_estimates - beta_true)^2, na.rm = TRUE)
    mse_results$naive_bias[i] <- mean(naive_estimates - beta_true, na.rm = TRUE)
    mse_results$naive_var[i]  <- var(naive_estimates, na.rm = TRUE)
    
    mse_results$cf_mse[i]  <- mean((cf_estimates - beta_true)^2, na.rm = TRUE)
    mse_results$cf_bias[i] <- mean(cf_estimates - beta_true, na.rm = TRUE)
    mse_results$cf_var[i]  <- var(cf_estimates, na.rm = TRUE)
  }
  
  return(mse_results)
}



# For location-scale DGP
evaluate_qte_estimators(quantiles = c(0.01,0.1,0.25,0.5), n=250, R=500)








