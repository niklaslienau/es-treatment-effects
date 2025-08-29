
## Monte Carlo Simulations for weighted estimator
### Histogram Plots and Summary Stats (MSE, Bias, Variance)

### MC Performance comparison CTATE vs Status Quo ES estimators

plot_mc_distribution_es_estimators <- function(tau = 0.5, R = 100, n = 1000, seed = 123) {
  set.seed(seed)
  
  naive_joint_estimates <- numeric(R)   
  naive_twostep_estimates <- numeric(R)   
  lienau_estimates <- numeric(R)      
  for (r in 1:R) {   
    cat("Iteration", r, "of", R)
    df <- simulate_dgp_rand(n = n)          
    
    # Joint ES     
    naive_joint_fit <- esreg(df$Y ~  df$D, alpha=0.25)     
    naive_joint_estimates[r] <- naive_joint_fit$coefficients[4]          
    
    # Two Step ES    
    naive_twostep_estimates[r] <- two_step_es_est(df$Y, df$D,0.25)     
    
    # Lienau Estimator     
    lienau_fit = efficient_es_est(df$Y, df$D, df$Z, tau = tau)    
    lienau_estimates[r] <- lienau_fit$estimate   }      
  
  # True CTATE   
  es_true <- compute_ctate_tau(tau=tau)      
  # Combine estimates  
  df_plot <- data.frame(
    estimate = c(naive_joint_estimates, naive_twostep_estimates,lienau_estimates),     
    method = factor(rep(c("Naive Joint ES","Naive Two Step ES" , "Lienau ES"), each = R))   )     
  df_plot <- df_plot[!is.na(df_plot$estimate), ]    
  
  # Plot   
  ggplot(df_plot, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(aes(xintercept = es_true, linetype = "True CTATE"), color = "black", linewidth = 1) +
    scale_linetype_manual(name = "", values = c("True CTATE" = "dashed")) +
    labs(
      title = paste0("MC Distribution of CTATE Estimators (τ = ", tau, ")"),
      x = "Estimate", y = "Density", fill = "Estimator"
    ) +
    theme_minimal()    
}


quartz()
plot_mc_distribution_es_estimators(tau = 0.25, R = 200, n=2000, seed=123)


