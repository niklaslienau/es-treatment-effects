### Simulations on the asymptotically efficient weighting estimator ###




# Performance as a function of regularization parameter c
# Function to run a full Monte Carlo simulation of the efficient ES estimator for different values of c
#Visualize bias variance tradeof
monte_carlo_weighted_es <- function(
    c_values, R = 1000, n = 500, tau_max = 0.25, grid_points = 10, B = 200
){
  start_time <- Sys.time()
  results_list <- list()
  
  for (r in 1:R) {
    cat("Monte Carlo iteration", r, "of", R, "\n")
    
    # 1) simulate
    df <- simulate_dgp_rand(n = n)
    
    # 2) fixed tau-grid for this rep
    tau_grid <- seq(0, tau_max, length.out = grid_points)
    tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
    G <- length(tau_grid)
    
    # QTE estimates 
    qte_estimates <- rep(NA_real_, G)
    for (g in seq_along(tau_grid)) {
      tau <- tau_grid[g]
      est <- try(cf_qr_series_avg(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      if (!inherits(est, "try-error") && is.finite(est)) qte_estimates[g] <- est
    }
    
    valid <- which(is.finite(qte_estimates))
    if (length(valid) < 2L) next
    
    qte_estimates <- qte_estimates[valid]
    tau_grid <- tau_grid[valid]
    G <- length(valid)
    u <- rep(1/G, G)
    
    # 3) bootstrap covariance
    boot <- try(bootstrap_qte_variances(Y = df$Y, D = df$D, Z = df$Z,
                                        B = B, qte_est = qte_estimates, grid = tau_grid),
                silent = TRUE)
    if (inherits(boot, "try-error")) next
    cov_matrix <- boot$cov_estimate
    
    # 4) ES for each c (correct Dmat/dvec)
    estimates <- rep(NA_real_, length(c_values))
    for (i in seq_along(c_values)) {
      c_val <- c_values[i]
      A <- cov_matrix + c_val * diag(G)
      Dmat <- 2 * A
      dvec <- 2 * c_val * u
      Amat <- matrix(1, nrow = G, ncol = 1)
      bvec <- 1
      sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
      if (!inherits(sol, "try-error")) {
        w <- sol$solution
        estimates[i] <- sum(w * qte_estimates)
      }
    }
    
    if (any(is.finite(estimates))) results_list[[length(results_list)+1L]] <- estimates
  }
  
  if (!length(results_list))
    stop("No successful MC replications — check that QTE and bootstrap succeed for your DGP.")
  
  result_matrix <- do.call(rbind, results_list)
  if (is.null(dim(result_matrix))) result_matrix <- matrix(result_matrix, nrow = 1)
  colnames(result_matrix) <- paste0("c_", c_values)
  
  theta_true <- compute_ctate_tau(tau_max)  # true ES at tau_max
  
  summary_stats <- data.frame(
    c = c_values,
    mse = apply(result_matrix, 2, function(x) mean((x - theta_true)^2, na.rm = TRUE)),
    variance = apply(result_matrix, 2, var, na.rm = TRUE),
    bias = apply(result_matrix, 2, function(x) mean(x - theta_true, na.rm = TRUE)),
    abs_bias = apply(result_matrix, 2, function(x) mean(abs(x - theta_true), na.rm = TRUE)),
    squared_bias = apply(result_matrix, 2, function(x) mean(x - theta_true, na.rm = TRUE)^2)
  )
  
  total_time <- Sys.time() - start_time
  cat("Total computation time:", total_time, "\n")
  
  list(estimates = result_matrix, summary = summary_stats, runtime = total_time)
}

# Run it
c_values <- c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2)
mc_plots <- monte_carlo_weighted_es(c= c_values, R = 500, n = 5000, tau_max =0.5, grid_points = 15)


# Plot
ggplot(mc_plots$summary, aes(x = log10(c))) +
  geom_line(aes(y = variance, color = "variance"), size = 1.2) +
  geom_point(aes(y = variance, color = "variance")) +
  labs(
    x = expression(log[10](c)),
    y = "Value",
    color = "Metric",
    title = "Monte Carlo Performance as a Function of Regularization Parameter c"
  ) +
  theme_minimal()
















######## Effect of Regularization on weights #############

library(quadprog)
library(ggplot2)

weights_vs_regul_plot <- function(R = 100,
                                        n = 500,
                                        tau_max = 0.5,
                                        grid_points_vec = c(5, 10, 15, 20),
                                        c_values = c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2),
                                        B = 200) {
  start_time <- Sys.time()
  results_list <- list()
  
  for (grid_points in grid_points_vec) {
    cat("Grid Points:", grid_points, "\n")
    G_weights <- list()
    
    for (r in 1:R) {
      cat("  Iteration", r, "of", R, "\n")
      df <- simulate_dgp(n = n, scale = FALSE)
      
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
    grid_points_vec = c(3, 5,7,9, 10),
    R = 1000,
    n = 500,
    tau_max = 0.5,
    B = 200,
    beta_true = 1,
    c_val = 1e-2
) {
  start_time <- Sys.time()
  results_list <- list()
  
  for (r in 1:R) {
    cat("Monte Carlo iteration", r, "of", R, "\n")
    df <- simulate_dgp(n = n, scale = FALSE)
    estimates_per_grid <- numeric(length(grid_points_vec))
    
    for (i in seq_along(grid_points_vec)) {
      grid_points <- grid_points_vec[i]
      tau_grid <- seq(0, tau_max, length.out = grid_points)
      tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
      G <- length(tau_grid)
      
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
      tau_grid <- tau_grid[valid]
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
  
  summary_df <- data.frame(
    grid_points = grid_points_vec,
    grid_fineness = grid_points_vec / tau_max,
    mse = apply(result_matrix, 2, function(x) mean((x - beta_true)^2, na.rm = TRUE)),
    variance = apply(result_matrix, 2, var, na.rm = TRUE),
    bias = apply(result_matrix, 2, function(x) mean(x - beta_true, na.rm = TRUE)),
    squared_bias = apply(result_matrix, 2, function(x) mean(x - beta_true, na.rm = TRUE)^2)
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









######



# weights over many taus for different regularization values
plot_weight_sensitivity <- function(n = 500, tau_max = 0.5, grid_points = 10,
                                    c_values = c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2),
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
      title = "The effect of c on the optimal weights",
      x = "Quantile Level (τ)",
      y = "Weight"
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

library(gridExtra)

# Generate 8 plots and store them in a list
plot_list <- vector("list", length = 8)
for (i in 1:8) {
  plot_list[[i]] <- plot_weight_sensitivity(
  )
}

# Display the plots in a grid (e.g. 2 rows × 4 columns)
quartz()  # macOS only
grid.arrange(grobs = plot_list, nrow = 2, ncol = 4)
















