
### MONTE CARLO for QTE: Naive QR vs CF Method 




plot_mc_distribution_qte_estimators <- function(tau = 0.5, R = 1000, n = 500, seed = 123, rho=0.5) {
  set.seed(seed)
  
  naive_estimates <- numeric(R)
  cf_estimates <- numeric(R)
  
  for (r in 1:R) {
    df <- simulate_dgp(
      n = n, rho = rho
    )
    
    # Naive QR
    naive_fit <- try(rq(Y ~ D, tau = tau, data = df), silent = TRUE)
    naive_estimates[r] <- if (!inherits(naive_fit, "try-error")) coef(naive_fit)["D"] else NA
    
    # Control Function QR
    cf_fit <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
    cf_estimates[r] <- if (!inherits(cf_fit, "try-error") && is.numeric(cf_fit)) cf_fit else NA
  }
  
  # True QTE
  qte_true <- compute_qte_tau(tau)
  
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
plot_mc_distribution_qte_estimators(tau = 0.05, R = 1000, n=1000, rho = 0.5)



########
library(quantreg)
library(ggplot2)

evaluate_mse_across_tau <- function(taus = seq(0.1, 0.9, by = 0.05),
                                    R = 500,
                                    n = 500,
                                    pi = 1,
                                    rho = 0.5,
                                    seed = 123) {
  set.seed(seed)
  
  mse_cf <- numeric(length(taus))
  mse_naive <- numeric(length(taus))
  
  for (i in seq_along(taus)) {
    tau <- taus[i]
    
    naive_estimates <- numeric(R)
    cf_estimates <- numeric(R)
    
    for (r in 1:R) {
      df <- simulate_dgp(n = n, rho = rho, pi = pi)
      
      # Naive QR
      naive_fit <- try(rq(Y ~ D, tau = tau, data = df), silent = TRUE)
      naive_estimates[r] <- if (!inherits(naive_fit, "try-error")) coef(naive_fit)["D"] else NA
      
      # CF QR
      cf_fit <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      cf_estimates[r] <- if (!inherits(cf_fit, "try-error") && is.numeric(cf_fit)) cf_fit else NA
    }
    
    # True QTE
    qte_true <- compute_qte_tau(tau)
    
    # Compute MSE (ignoring NAs)
    mse_naive[i] <- mean((naive_estimates - qte_true)^2, na.rm = TRUE)
    mse_cf[i] <- mean((cf_estimates - qte_true)^2, na.rm = TRUE)
  }
  
  # Plot
  df_plot <- data.frame(
    tau = rep(taus, 2),
    mse = c(mse_naive, mse_cf),
    method = rep(c("Naive QR", "Control Function QR"), each = length(taus))
  )
  
  ggplot(df_plot, aes(x = tau, y = mse, color = method)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(
      title = "MSE of QTE Estimators across Quantiles",
      x = expression(tau), y = "Mean Squared Error", color = "Estimator"
    ) +
    theme_minimal()
}





evaluate_mse_across_tau()
