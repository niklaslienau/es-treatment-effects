## Monte Carlo Simulations for weighted estimator

############################################### (1) EFFICIENT WEIGHTING ################################


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
    beta_true = 1
) {
  start_time <- Sys.time()
  results_list <- list()
  
  for (r in 1:R) {
    cat("Monte Carlo iteration", r, "of", R, "\n")
    
    # Step 1: Simulate data
    df <- simulate_dgp(n = n, scale = FALSE)
    
    # Step 2: Estimate QTEs across tau grid
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
    
    # Step 3: Bootstrap variance-covariance matrix
    boot <- bootstrap_qte_variances(Y = df$Y, D = df$D, Z = df$Z, tau_max = tau_max, grid_points = grid_points, B = B)
    cov_matrix <- boot$cov_matrix[valid, valid, drop = FALSE]
    u <- rep(1 / G, G)
    
    # Step 4: Calculate weights and compute ES estimates for each c
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
  
  # Combine into matrix
  result_matrix <- do.call(rbind, results_list)
  colnames(result_matrix) <- paste0("c_", c_values)
  
  # Compute summary statistics
  summary_stats <- data.frame(
    c = c_values,
    mse = apply(result_matrix, 2, function(x) mean((x - beta_true)^2, na.rm = TRUE)),
    variance = apply(result_matrix, 2, var, na.rm = TRUE),
    bias = apply(result_matrix, 2, function(x) mean(x - beta_true, na.rm = TRUE)),
    abs_bias = apply(result_matrix, 2, function(x) mean(abs(x - beta_true), na.rm = TRUE)),
    squared_bias = apply(result_matrix, 2, function(x) mean(x - beta_true, na.rm = TRUE)^2)
  )
  
  total_time <- Sys.time() - start_time
  cat("Total computation time:", total_time, "\n")
  
  return(list(
    estimates = result_matrix,
    summary = summary_stats,
    runtime = total_time
  ))
}



# Run it
c_values <- c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4 )
mc_plots <- monte_carlo_weighted_es(c= c_values, R = 1000, n = 500, tau_max =0.5)


# Plot
ggplot(mc_plots$summary, aes(x = log10(c))) +
  geom_line(aes(y = mse, color = "MSE"), size = 1.2) +
  geom_point(aes(y = mse, color = "MSE")) +
  labs(
    x = expression(log[10](c)),
    y = "Value",
    color = "Metric",
    title = "Monte Carlo Performance as a Function of Regularization Parameter c"
  ) +
  theme_minimal()










######################### Bootstrap diagnostics ###########


# Generate 50 bootstrap variance curves from 50 different random samples
# Analyze variance in bootstrap estimates and how estimates behave across taus
#library(ggplot2)
#library(gridExtra)

#plot_bootstrap_variance_curves <- function(n = 500, B = 200, grid_points = 10, tau_max = 0.25, reps = 50) {
  # Helper to run one replication for a given sample size
  run_for_n <- function(n_single) {
    all_variances <- list()
    tau_grid <- NULL
    
    for (i in 1:reps) {
      df <- simulate_dgp(n = n_single, scale = FALSE)
      
      boot <- bootstrap_qte_variances(df$Y, df$D, df$Z,
                                      tau_max = tau_max,
                                      grid_points = grid_points,
                                      B = B)
      
      if (is.null(tau_grid)) tau_grid <- boot$tau_grid
      
      curve_df <- data.frame(
        tau = boot$tau_grid,
        variance = diag(boot$cov_matrix),
        sim = paste0("sim_", i)
      )
      
      all_variances[[i]] <- curve_df
    }
    
    df_all <- do.call(rbind, all_variances)
    
    ggplot(df_all, aes(x = tau, y = variance, group = sim)) +
      geom_line(alpha = 0.5, color = "steelblue") +
      labs(title = paste0("Bootstrap Variance Estimates (n = ", n_single, ")"),
           x = "Quantile Level (tau)", y = "Bootstrap Variance") +
      theme_minimal()
  }
  
  if (length(n) == 1) {
    # Single sample size
    return(run_for_n(n))
  } else {
    # Multiple sample sizes → plot grid
    plot_list <- lapply(n, run_for_n)
    
    quartz()  # macOS only; skip or change on other OS
    grid.arrange(grobs = plot_list, ncol = 2)
  }
}
#plot_bootstrap_variance_curves(n=c(250,500,1000,2000), reps=30) #Run




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
    tau_max = 0.25,
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



############### SAVE SIMULATION RESULTS  ############### 
save(mc_plots, res_1,res_2,res_3, file = "sim_res.RData")








############### STOP SIMULATING HERE !!!!!!!!!!!!!!!!!  ###############


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

















####################################### (1) Uniform WEIGHTING ################################ 


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

mc_res = evaluate_grid_averaged_estimator(tau_max = 0.25, grid_points = 10, n=500, R=500)


####2#####
# Uniform Performance vs #Grid Points

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
res_plot_obj <- grid_number_evaluation(tau_max = 0.25, n = 500, R= 5000,grid_points_vec = c(3,5,10,25,50))
res_plot_obj$plot  # to display the ggplot
res_plot_obj$results_table  # to inspect the raw numbers

