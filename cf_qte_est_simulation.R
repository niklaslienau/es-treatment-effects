
### MONTE CARLO for QTE: Naive QR vs CF Method on MSE, Bias & Variance


#####1#####
#MSE Var and Bias Table for multiple quantiles
evaluate_cf_qte <- function(
    tau = 0.25,
    R = 1000,
    n = 1000,
    beta_true = 1,
    seed = 123
) {
  set.seed(seed)
  
  cf_estimates <- numeric(R)
  
  for (r in 1:R) {
    df <- simulate_loc_shift(n = n)
    
    cf_try <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
    cf_estimates[r] <- if (inherits(cf_try, "try-error") || !is.numeric(cf_try) || length(cf_try) != 1 || is.na(cf_try)) NA else cf_try
  }
  
  # Clean and compute MSE
  valid_estimates <- cf_estimates[!is.na(cf_estimates)]
  mse <- mean((valid_estimates - beta_true)^2)
  
  return(list(
    estimates = valid_estimates,
    mse = mse
  ))
}




res <- evaluate_cf_qte(tau = 0.25, R = 1000, n = 1000)
res$mse     # Mean squared error
hist(res$estimates)  # Vector of QTE estimates



#####1#####
#Monte Carlo Performance as a function of SAMPLE SIZE plot

### QTE estimator performance by sample size
evaluate_qte_performance_by_sample_size <- function(
    quantiles = c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5),
    sample_sizes = c(30, 50, 100, 150, 200, 300, 500, 1000),
    R = 100,
    beta_true = 1,
    seed = 123
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  results <- expand.grid(
    n = sample_sizes,
    tau = quantiles,
    stringsAsFactors = FALSE
  ) %>% 
    mutate(cf_mse = NA, cf_abs_bias = NA, cf_var = NA)
  
  for (i in 1:nrow(results)) {
    n <- results$n[i]
    tau <- results$tau[i]
    
    cat("Running for n =", n, ", tau =", tau, "\n")
    
    res <- evaluate_qte_estimators(
      quantiles = tau,
      R = R,
      n = n,
      beta_true = beta_true,
      seed = seed
    )
    
    results$cf_mse[i] <- res$cf_mse[1]
    results$cf_abs_bias[i] <- abs(res$cf_bias[1])
    results$cf_var[i] <- res$cf_var[1]
  }
  
  # Convert to long format for plotting
  results_long <- results %>%
    pivot_longer(cols = c("cf_mse", "cf_abs_bias", "cf_var"),
                 names_to = "metric", values_to = "value")
  
  # Plot
  plot <- ggplot(results_long, aes(x = n, y = value, color = as.factor(tau))) +
    geom_line() +
    geom_point() +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Performance of CF QTE Estimator vs. Sample Size",
      x = "Sample Size",
      y = "Value",
      color = "Quantile τ"
    ) +
    theme_minimal()
  
  return(list(results_table = results, plot = plot))
}

res <- evaluate_qte_performance_by_sample_size(quantiles = c( 0.01,0.05,0.1),
                                               sample_sizes = c(250, 500, 1000),
                                                                R = 1000, seed=123)
res$plot





