
### MONTE CARLO for QTE: Naive QR vs CF Method on MSE, Bias & Variance


############## Var & Bias of QTE Estimates across tau values #########################
#Only for linear hetero DGP with continous D
evaluate_qte_bias_variance <- function(
    tau_grid = seq(0.05, 0.5, length.out = 10),
    R = 100,
    n = 500,
    beta = 1,
    gamma = 0.5,
    seed= 123
) {
  set.seed(seed)
  
  results <- data.frame(
    tau = tau_grid,
    abs_bias = NA,
    variance = NA
  )
  
  qte_matrix <- matrix(NA, nrow = R, ncol = length(tau_grid))
  
  for (r in 1:R) {
    df <- simulate_dgp(n = n, beta = beta, gamma = gamma)
    
    for (i in seq_along(tau_grid)) {
      tau <- tau_grid[i]
      est <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      
      if (!inherits(est, "try-error") && is.numeric(est) && !is.na(est)) {
        qte_matrix[r, i] <- est
      }
    }
  }
  
  # True QTE for each quantile level (constant in D)
  qte_true <- beta + gamma * qnorm(tau_grid)
  
  for (i in seq_along(tau_grid)) {
    estimates <- qte_matrix[, i]
    valid <- !is.na(estimates)
    results$abs_bias[i] <- mean(abs(estimates[valid] - qte_true[i]))
    results$variance[i] <- var(estimates[valid])
  }
  
  return(results)
}


#EXAMPLE USE
res <- evaluate_qte_bias_variance()

ggplot(res, aes(x = tau)) +
  #geom_line(aes(y = abs_bias, color = "Absolute Bias")) +
  geom_line(aes(y = variance, color = "Variance")) +
  labs(x = "Quantile Level (τ)", y = "Value", color = "Metric",
       title = "Bias and Variance of CF-QTE across Quantiles") +
  theme_minimal()






########## Emprical Distribution vs Naive Estimator ###########
library(quantreg)
library(ggplot2)

plot_mc_distribution_qte_estimators <- function(tau = 0.5, R = 100, n = 500, seed = 123, beta_true=1) {
  set.seed(seed)
  
  naive_estimates <- numeric(R)
  cf_estimates <- numeric(R)
  
  for (r in 1:R) {
    df <- sim_dgp_homo(
      n = n
    )
    
    # Naive QR
    naive_fit <- try(rq(Y ~ D, tau = tau, data = df), silent = TRUE)
    naive_estimates[r] <- if (!inherits(naive_fit, "try-error")) coef(naive_fit)["D"] else NA
    
    # Control Function QR
    cf_fit <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
    cf_estimates[r] <- if (!inherits(cf_fit, "try-error") && is.numeric(cf_fit)) cf_fit else NA
  }
  
  # True QTE
  qte_true <- beta_true
  
  # Combine estimates
  df_plot <- data.frame(
    estimate = c(naive_estimates, cf_estimates),
    method = factor(rep(c("Naive QR", "Control Function QR"), each = R))
  )
  
  df_plot <- df_plot[!is.na(df_plot$estimate), ]
  
  # Plot
  ggplot(df_plot, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(aes(xintercept = qte_true, color = "True QTE", linetype = "True QTE"), linewidth = 1) +
    scale_color_manual(name = "", values = c("True QTE" = "black"), labels = expression(beta(tau))) +
    scale_linetype_manual(name = "", values = c("True QTE" = "dashed"), labels = expression(beta(tau))) +
    labs(
      title = paste0("MC Distribution of QTE Estimators (τ = ", tau, ")"),
      x = "Estimate", y = "Density", fill = "Estimator"
    ) +
    theme_minimal()
  
}

#Example Use
quartz()
plot_mc_distribution_qte_estimators(tau = 0.25, R = 500, n=500, seed=123)














