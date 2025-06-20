## Monte Carlo simulations for linear hetero location scale DGP




# Performance as a function of regularization parameter c
# Function to run a full Monte Carlo simulation of the efficient ES estimator for different values of c
#Visualize bias variance tradeof
monte_carlo_weighted_es <- function(
    c_values,
    R = 1000,
    n = 500,
    tau_max = 0.25,
    grid_points = 10,
    B = 200,
    beta = 1,
    gamma = 0.5
) {
  start_time <- Sys.time()
  results_list <- list()
  
  # Compute true ES-TE for comparison
  tau_grid_fine <- seq(0, tau_max, length.out = 1000)[-1]  # remove 0
  integrand <- qnorm(tau_grid_fine)
  es_true <- beta + gamma * mean(integrand)
  
  for (r in 1:R) {
    cat("Monte Carlo iteration", r, "of", R, "\n")
    
    # Step 1: Simulate from linear location-scale DGP
    df <- simulate_dgp(n = n, beta = beta, gamma = gamma, scale = TRUE, linear = TRUE)
    
    # Step 2: Estimate QTEs
    tau_grid <- seq(0, tau_max, length.out = grid_points)
    tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
    G <- length(tau_grid)
    qte_estimates <- numeric(G)
    
    for (g in seq_along(tau_grid)) {
      tau <- tau_grid[g]
      model <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      qte_estimates[g] <- if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) model else NA
    }
    
    valid <- which(!is.na(qte_estimates))
    if (length(valid) < 2) next
    
    qte_estimates <- qte_estimates[valid]
    tau_grid <- tau_grid[valid]
    G <- length(valid)
    
    # Step 3: Bootstrap covariance matrix
    boot <- bootstrap_qte_variances(Y = df$Y, D = df$D, Z = df$Z, tau_max = tau_max, grid_points = grid_points, B = B)
    cov_matrix <- boot$cov_matrix[valid, valid, drop = FALSE]
    u <- rep(1 / G, G)
    
    # Step 4: Weighted ES estimates for each c
    estimates <- numeric(length(c_values))
    for (i in seq_along(c_values)) {
      c_val <- c_values[i]
      Dmat <- cov_matrix + c_val * diag(G)
      dvec <- c_val * u
      Amat <- matrix(1, nrow = G, ncol = 1)
      bvec <- 1
      
      sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
      if (!inherits(sol, "try-error")) {
        weights <- sol$solution
        estimates[i] <- sum(weights * qte_estimates)
      } else {
        estimates[i] <- NA
      }
    }
    
    results_list[[r]] <- estimates
  }
  
  # Combine MC results
  result_matrix <- do.call(rbind, results_list)
  colnames(result_matrix) <- paste0("c_", c_values)
  
  # Summary statistics vs true ES treatment effect
  summary_stats <- data.frame(
    c = c_values,
    mse = apply(result_matrix, 2, function(x) mean((x - es_true)^2, na.rm = TRUE)),
    variance = apply(result_matrix, 2, var, na.rm = TRUE),
    bias = apply(result_matrix, 2, function(x) mean(x - es_true, na.rm = TRUE)),
    abs_bias = apply(result_matrix, 2, function(x) mean(abs(x - es_true), na.rm = TRUE)),
    squared_bias = apply(result_matrix, 2, function(x) mean(x - es_true, na.rm = TRUE)^2)
  )
  
  total_time <- Sys.time() - start_time
  cat("Total computation time:", total_time, "\n")
  
  return(list(
    estimates = result_matrix,
    summary = summary_stats,
    runtime = total_time,
    true_es = es_true
  ))
}

c_values <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)
mc_plots <- monte_carlo_weighted_es(c_values = c_values, R = 1000, n = 500, tau_max = 0.5)

ggplot(mc_plots$summary, aes(x = log10(c))) +
  geom_line(aes(y = mse, color = "MSE"), size = 1.2) +
  geom_point(aes(y = mse, color = "MSE")) +
  labs(
    x = expression(log[10](c)),
    y = "Value",
    color = "Metric",
    title = "Monte Carlo Performance vs Regularization Parameter c"
  ) +
  theme_minimal()




###Compare bootstrap against MC estimates
## See how well bootstraping performs on average to estimate the variance
compare_variance_bootstrap_mc_by_tau <- function(tau_vec = seq(0.01,0.5, length.out=20), 
                                                 n = 500, R = 1000, B = 200) {
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
      df <- simulate_dgp(n = n, scale = FALSE, linear=TRUE)
      
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
res_1 <- compare_variance_bootstrap_mc_by_tau(tau_vec =seq(0.05, 0.5, length.out = 10), n = 500, R = 1000, B = 200)

ggplot(res_1, aes(x = tau)) +
  geom_line(aes(y = mc_var, color = "Monte Carlo")) +
  geom_line(aes(y = bootstrap_var, color = "Avg Nonparametric Bootstrap")) +
  labs(
    x = "Quantile Level (τ)",
    y = "Variance Estimate",
    title = "Variance for CF-QTE estimator across τ",
    color = "Estimator"
  ) +
  theme_minimal()




######## Effect of Regularization on weights #############

library(quadprog)
library(ggplot2)

weights_vs_regul_plot <- function(R = 100,
                                  n = 500,
                                  tau_max = 0.25,
                                  grid_points_vec = c(5, 10, 15, 20),
                                  c_values = c(1e-6, 1e-3, 1e0, 1e3),
                                  B = 200) {
  start_time <- Sys.time()
  results_list <- list()
  
  for (grid_points in grid_points_vec) {
    cat("Grid Points:", grid_points, "\n")
    G_weights <- list()
    
    for (r in 1:R) {
      cat("  Iteration", r, "of", R, "\n")
      df <- simulate_dgp(n = n, scale = FALSE, linear = TRUE)
      
      boot <- bootstrap_qte_variances(Y = df$Y, D = df$D, Z = df$Z,
                                      tau_max = tau_max, grid_points = grid_points, B = B)
      tau_grid <- boot$tau_grid
      cov_matrix <- boot$cov_matrix
      G <- length(tau_grid)
      u <- rep(1 / G, G)
      
      weights_mat <- matrix(NA, nrow = G, ncol = length(c_values))
      colnames(weights_mat) <- paste0("c=", c_values)
      
      for (j in seq_along(c_values)) {
        c_val <- c_values[j]
        Dmat <- cov_matrix + c_val * diag(G)
        dvec <- c_val * u
        Amat <- matrix(1, nrow = G, ncol = 1)
        bvec <- 1
        
        sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
        if (!inherits(sol, "try-error")) {
          weights_mat[, j] <- sol$solution
        }
      }
      
      G_weights[[r]] <- weights_mat
    }
    
    avg_weights <- Reduce("+", G_weights) / R
    df_avg <- data.frame(
      tau = rep(tau_grid, times = length(c_values)),
      weight = as.vector(avg_weights),
      c = rep(paste0("c=", c_values), each = length(tau_grid)),
      grid_points = grid_points
    )
    results_list[[as.character(grid_points)]] <- df_avg
  }
  
  all_results <- do.call(rbind, results_list)
  
  plot_obj <- ggplot(all_results, aes(x = tau, y = weight, color = c)) +
    geom_line(size = 1) +
    facet_wrap(~ grid_points, scales = "free_y") +
    geom_hline(aes(yintercept = 1 / grid_points), color = "black", linetype = "dashed", inherit.aes = FALSE) +
    labs(title = "Average Efficient Weights across Monte Carlo Simulations",
         x = "Quantile Level (tau)", y = "Average Weight") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  runtime <- Sys.time() - start_time
  cat("Total runtime:", runtime, "\n")
  
  return(list(
    plot = plot_obj,
    weights = all_results,
    runtime = runtime
  ))
}

#use
quartz()
res_2 <- weights_vs_regul_plot(R = 1000, n = 500, grid_points_vec = c(3, 5, 10, 15))
print(res_2$plot)


########### Effect of Grid Size on MSE #################
# calibrate c with 10^-2 
library(quadprog)
library(ggplot2)
evaluate_grid_fineness <- function(
    grid_points_vec = c(3, 5, 7, 9, 10),
    R = 1000,
    n = 500,
    tau_max = 0.25,
    B = 200,
    beta = 1,
    gamma = 1,
    c_val = 1e-2
) {
  start_time <- Sys.time()
  results_list <- list()
  
  for (r in 1:R) {
    cat("Monte Carlo iteration", r, "of", R, "\n")
    df <- simulate_dgp(n = n, beta = beta, gamma = gamma, scale = TRUE, linear = TRUE)  # << changed
    estimates_per_grid <- numeric(length(grid_points_vec))
    
    for (i in seq_along(grid_points_vec)) {
      grid_points <- grid_points_vec[i]
      tau_grid <- seq(0, tau_max, length.out = grid_points)
      tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
      G <- length(tau_grid)
      
      # True ES for current grid
      es_true <- beta + gamma * mean(qnorm(tau_grid))  # << added
      
      # Estimate QTEs
      qte_estimates <- numeric(G)
      for (g in seq_along(tau_grid)) {
        tau <- tau_grid[g]
        model <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
        qte_estimates[g] <- if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) model else NA
      }
      
      valid <- which(!is.na(qte_estimates))
      if (length(valid) < 2) {
        estimates_per_grid[i] <- NA
        next
      }
      
      qte_estimates <- qte_estimates[valid]
      tau_grid_valid <- tau_grid[valid]
      G <- length(valid)
      
      # Bootstrap variance-covariance matrix
      boot <- bootstrap_qte_variances(df$Y, df$D, df$Z, tau_max = tau_max, grid_points = grid_points, B = B)
      cov_matrix <- boot$cov_matrix[valid, valid, drop = FALSE]
      u <- rep(1 / G, G)
      
      # Compute optimal weights and ES
      Dmat <- cov_matrix + c_val * diag(G)
      dvec <- c_val * u
      Amat <- matrix(1, nrow = G, ncol = 1)
      bvec <- 1
      
      sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
      if (!inherits(sol, "try-error")) {
        weights <- sol$solution
        estimates_per_grid[i] <- sum(weights * qte_estimates)
      } else {
        estimates_per_grid[i] <- NA
      }
    }
    
    results_list[[r]] <- estimates_per_grid
  }
  
  result_matrix <- do.call(rbind, results_list)
  colnames(result_matrix) <- paste0("G_", grid_points_vec)
  
  # True values for each grid
  true_vals <- sapply(grid_points_vec, function(G) {
    taus <- seq(0, tau_max, length.out = G)
    taus <- taus[taus > 0 & taus < 1]
    beta + gamma * mean(qnorm(taus))
  })
  
  summary_df <- data.frame(
    grid_points = grid_points_vec,
    grid_fineness = grid_points_vec / tau_max,
    mse = mapply(function(i, truth) mean((result_matrix[, i] - truth)^2, na.rm = TRUE),
                 seq_along(grid_points_vec), true_vals),
    variance = mapply(function(i) var(result_matrix[, i], na.rm = TRUE), seq_along(grid_points_vec)),
    bias = mapply(function(i, truth) mean(result_matrix[, i] - truth, na.rm = TRUE),
                  seq_along(grid_points_vec), true_vals),
    squared_bias = mapply(function(i, truth) mean(result_matrix[, i] - truth, na.rm = TRUE)^2,
                          seq_along(grid_points_vec), true_vals)
  )
  
  runtime <- Sys.time() - start_time
  cat("Total computation time:", runtime, "\n")
  
  return(list(
    estimates = result_matrix,
    summary = summary_df,
    runtime = runtime
  ))
}

res_3 <- evaluate_grid_fineness()

ggplot(res_3$summary, aes(x = grid_fineness)) +
  geom_line(aes(y = mse, color = "MSE")) +
  labs(x = "Grid Fineness (Grid Points / Tau Max)", y = "Metric Value", color = "Metric",
       title = "Effect of Grid Fineness on Estimator Performance") +
  theme_minimal()

