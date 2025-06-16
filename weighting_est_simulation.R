## Monte Carlo Simulations for weighted estimator

############################################### UNIFORM ################################
#### 1#####
#UNIFORM MSE Output Table

evaluate_grid_averaged_estimator <- function(
    tau_max = 0.5,
    grid_points = 10,
    R = 1000,
    n = 500,
    beta_true = 1,
    seed = 123
) {
  set.seed(seed)
  
  grid_estimates <- numeric(R)
  
  for (r in 1:R) {
    # Draw sample from DGP
    df <- simulate_loc_shift(n = n)
    
    # Compute grid-averaged estimate
    grid_try <- try(uniform_es_est(df$Y, df$D, df$Z, 
                                        tau_max = tau_max, 
                                        grid_points = grid_points
    ), 
    silent = TRUE)
    
    grid_estimates[r] <- if (inherits(grid_try, "try-error")) {
      NA
    } else if (!is.list(grid_try) || is.null(grid_try$estimate)) {
      NA  
    } else if (!is.numeric(grid_try$estimate) || length(grid_try$estimate) != 1 || is.na(grid_try$estimate)) {
      NA
    } else {
      grid_try$estimate
    }
  }
  
  # Remove NAs for performance calculation
  valid_estimates <- grid_estimates[!is.na(grid_estimates)]
  
  # Calculate performance metrics
  mse <- mean((valid_estimates - beta_true)^2)
  bias <- mean(valid_estimates - beta_true)
  variance <- var(valid_estimates)
  
  # Return results
  results <- list(
    estimates = grid_estimates,
    valid_estimates = valid_estimates,
    n_valid = length(valid_estimates),
    mse = mse,
    bias = bias,
    variance = variance,
    tau_max = tau_max,
    grid_points = grid_points,
    R = R
  )
  
  return(results)
}

#MC Eval

mc_res = evaluate_grid_averaged_estimator(tau_max = 0.25, grid_points = 10, n=500, R=1000)


####2#####
# Uniform Perfermance vs #Grid Points

## Bias Variance Trade of in number of grid points
grid_number_evaluation <- function(
    tau_max = 0.25,
    grid_points_vec = c(3, 5, 10, 25, 50, 100),
    n = 250,
    R = 200,
    beta_true = 1
) {
  library(ggplot2)
  library(tidyr)
  
  results <- data.frame(
    grid_points = grid_points_vec,
    mse = NA,
    abs_bias = NA,
    variance = NA
  )
  
  for (i in seq_along(grid_points_vec)) {
    gp <- grid_points_vec[i]
    cat("Running for grid_points =", gp, "...\n")
    
    res <- evaluate_grid_averaged_estimator(
      tau_max = tau_max,
      grid_points = gp,
      n = n,
      R = R,
      beta_true = beta_true
    )
    
    results$mse[i] <- res$mse
    results$abs_bias[i] <- abs(res$bias)
    results$variance[i] <- res$variance
  }
  
  # Reshape for plotting
  results_long <- pivot_longer(results, cols = c("mse", "abs_bias", "variance"),
                               names_to = "metric", values_to = "value")
  
  # Plot
  plot <- ggplot(results_long, aes(x = grid_points, y = value, color = metric)) +
    geom_line() +
    geom_point() +
    labs(
      title = paste0("Bias Variance Trade off from MC Sim"),
      x = "Number of Grid Points",
      y = "Value",
      color = "Metric"
    ) +
    theme_minimal()
  
  return(list(results_table = results, plot = plot))
}

#Example Use 
#NOTE: Here is no Bias Variance Trade off because each individual grid estimator in unbiased for target paramter
# Also: for extreme tails in small samples QR are known to be biased -> Bias increases in grid size
res_plot_obj <- grid_number_evaluation(tau_max = 0.25, n = 500, R= 1000,grid_points_vec = c(3,5,10,25,50))
res_plot_obj$plot  # to display the ggplot
res_plot_obj$results_table  # to inspect the raw numbers


############################################### EFFICIENT WEIGHTING ################################


# Performance as a function of regularization parameter c
monte_carlo_es_performance <- function(c_values, R = 1000, n = 250, beta_true = 1, tau_max = 0.25, grid_points = 10) {
  results <- data.frame(
    c = numeric(),
    mse = numeric(),
    variance = numeric()
  )
  
  for (c_val in c_values) {
    estimates <- numeric(R)
    
    for (r in 1:R) {
      df <- simulate_loc_shift(n = n)
      est <- try(estimate_weighted_es_cf(df$Y, df$D, df$Z, tau_max = tau_max, grid_points = grid_points, c = c_val), silent = TRUE)
      
      estimates[r] <- if (inherits(est, "try-error") || is.null(est$estimate) || is.na(est$estimate)) NA else est$estimate
    }
    
    estimates <- na.omit(estimates)
    
    results <- rbind(results, data.frame(
      c = c_val,
      mse = mean((estimates - beta_true)^2),
      variance = var(estimates)
    ))
  }
  
  return(results)
}


# Run it
c_values <- c(1e-6,1e-4, 1e-2,1e0,1e2,1e4 ,1e6)
mc_results <- monte_carlo_es_performance(c_values, R = 1000, n = 250)

# Plot
ggplot(mc_results, aes(x = log10(c))) +
  geom_line(aes(y = mse, color = "MSE"), size = 1.2) +
  geom_point(aes(y = mse, color = "MSE")) +
  geom_line(aes(y = variance, color = "Variance"), size = 1.2, linetype = "dashed") +
  geom_point(aes(y = variance, color = "Variance")) +
  labs(x = expression(log[10](c)), y = "Value", color = "Metric",
       title = "Monte Carlo Performance as a function of Regularization Parameter c") +
  theme_minimal()

seq(1e-6, 1e5, length.out = 100)








######################### Bootstrap vs MC Variance Estimation ###########


### against quantile level
compare_variance_bootstrap_mc_by_tau <- function(tau_vec = c(0.01, 0.05, 0.25, 0.5), 
                                                 n = 500, R = 200, B = 100) {
  results <- data.frame(
    tau = numeric(),
    bootstrap_var = numeric(),
    mc_var = numeric()
  )
  
  for (tau in tau_vec) {
    # Store QTEs from each MC iteration
    qte_estimates <- numeric(R)
    bootstrap_vars <- numeric(R)
    
    for (r in 1:R) {
      df <- simulate_dgp(n = n, scale = FALSE)
      
      # QTE estimation
      qte <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      qte_estimates[r] <- if (!inherits(qte, "try-error") && is.numeric(qte)) qte else NA
      
      # Bootstrap variance estimation
      vboot <- try(bootstrap_qte_variance_single(df$Y, df$D, df$Z, tau = tau, B = B), silent = TRUE)
      bootstrap_vars[r] <- if (!inherits(vboot, "try-error") && is.numeric(vboot)) vboot else NA
    }
    
    # Final summary for this tau
    results <- rbind(results, data.frame(
      tau = tau,
      bootstrap_var = mean(bootstrap_vars, na.rm = TRUE),
      mc_var = var(na.omit(qte_estimates))
    ))
  }
  
  return(results)
}


#Example Use
res <- compare_variance_bootstrap_mc_by_tau(tau_vec =c(0.01,0.05,0.1,0.2, 0.25, 0.3,0.4,0.5), n = 500, R = 1000, B = 250)

ggplot(res, aes(x = tau)) +
  geom_line(aes(y = mc_var, color = "Monte Carlo")) +
  geom_line(aes(y = bootstrap_var, color = "Avg Nonparametric Bootstrap")) +
  labs(
    x = "Quantile Level (τ)",
    y = "Variance Estimate",
    title = "Variance for Lienau QTE estimator across τ",
    color = "Estimator"
  ) +
  theme_minimal()


######## Effect of Regularization on weights #############
plot_weight_sensitivity <- function(n = 500, tau_max = 0.25, grid_points = 10,
                                    c_values = c(1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6),
                                    B = 200) {
  library(ggplot2)
  library(quadprog)
  library(viridis)
  
  # 1. Simulate data once
  df <- simulate_dgp(n = n, scale = FALSE)
  
  # 2. Bootstrap covariance matrix
  boot <- bootstrap_qte_variances(Y = df$Y, D = df$D, Z = df$Z,
                                  tau_max = tau_max, grid_points = grid_points, B = B)
  tau_grid <- boot$tau_grid
  cov_matrix <- boot$cov_matrix
  G <- length(tau_grid)
  u <- rep(1 / G, G)
  
  # 3. Container for all weights
  weight_df <- data.frame()
  
  # 4. Loop over c values
  for (c in c_values) {
    Dmat <- cov_matrix + c * diag(G)
    dvec <- c * u
    Amat <- matrix(1, nrow = G, ncol = 1)
    bvec <- 1
    
    sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
    if (!inherits(sol, "try-error")) {
      weights <- sol$solution
      weight_df <- rbind(weight_df, data.frame(
        tau = tau_grid,
        weight = weights,
        c_label = paste0("c = ", format(c, scientific = TRUE)),
        c_value = c
      ))
    }
  }
  
  # 5. Add uniform weights
  uniform_df <- data.frame(
    tau = tau_grid,
    weight = rep(1 / G, G),
    c_label = "Uniform",
    c_value = NA
  )
  
  weight_df <- rbind(weight_df, uniform_df)
  
  # 6. Define color scale
  ordered_labels <- c(paste0("c = ", format(c_values, scientific = TRUE)), "Uniform")
  weight_df$c_label <- factor(weight_df$c_label, levels = ordered_labels)
  
  # 7. Plot
  ggplot(weight_df, aes(x = tau, y = weight, color = c_label, linetype = c_label)) +
    geom_line(size = 1) +
    scale_color_manual(
      values = c(
        viridis(length(c_values), direction = -1),  # darker to lighter
        "black"                                    # uniform
      )
    ) +
    scale_linetype_manual(
      values = c(rep("solid", length(c_values)), "dashed")
    ) +
    labs(
      title = "Optimal Weights vs. Regularization Parameter c",
      x = "Quantile Level (τ)",
      y = "Weight"
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

plot_weight_sensitivity(n = 500, tau_max = 0.25, grid_points = 10,c_values = c(1e-6, 1e-5,1e-4, 1e-3,1e-2, 1e-1,1e0, 1e2,1e3), B = 100)

