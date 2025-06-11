
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
      df <- simulate_loc_shift(n = n)
      
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



# Evaluate
evaluate_qte_estimators(quantiles = c(0.001,0.01,0.025,0.1, 0.25,0.5), n=100, R=200, seed=123)



### QTE estimator performance by sample size
evaluate_qte_performance_by_sample_size <- function(
    quantiles = c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5),
    sample_sizes = c(30, 50, 100, 150, 200, 300, 500, 1000),
    R = 100,
    beta_true = 1,
    degree = 3,
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
      degree = degree,
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

res <- evaluate_qte_performance_by_sample_size(quantiles = c(0.01, 0.05,0.1,0.25,0.5),
                                               sample_sizes = c(100, 200, 350, 500, 1000,10000),
                                                                R = 1000, seed=123)
res$plot






