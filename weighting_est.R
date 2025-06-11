library(ggplot2)
library(tidyr)
### Grab variance estimate for QR coefficent (for now from regression, maybe later bootstrap ?)
data= simulate_loc_shift()
cf_qr_res=cf_qr_estimate(return_model = TRUE, tau=0.25, Y= data$Y, D= data$D, Z= data$Z)
summary(cf_qr_res, se= "nid")$coefficients["D",2] # grab sd from regression output


### First Uniform weighting
#u = 1/ grid points (if evenly spaced)


# Weighted Quantile Estimator (UNIFORM)
cf_qr_grid_estimate <- function(Y, D,Z,  tau_max, grid_points, degree = 3) {
  # Create evenly spaced grid from 0 to tau_max
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  
  # Remove tau = 0 if it's included (quantile regression undefined at 0)
  tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
  
  estimates <- numeric(length(tau_grid))
  
  # Estimate at each grid point
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    cf_try <- try(cf_qr_estimate(Y, D, Z, tau = tau, degree = degree, return_model = TRUE), silent = TRUE)
    
    # Extract the D coefficient from the regression output
    if (inherits(cf_try, "try-error")) {
      estimates[i] <- NA
    } else {
      # Assuming cf_qr_estimate returns a regression object with coefficients
      coef_try <- try(cf_try$coefficients["D"], silent = TRUE)
      estimates[i] <- if (inherits(coef_try, "try-error") || !is.numeric(coef_try) || 
                          length(coef_try) != 1 || is.na(coef_try)) NA else coef_try
    }
  }
  
  # Return both the averaged estimate and the grid information
  return(list(
    estimate = mean(estimates, na.rm = TRUE),
    tau_grid = tau_grid,
    grid_estimates = estimates
  ))
}

# Monte Carlo evaluation function for grid-averaged estimator
evaluate_grid_averaged_estimator <- function(
    tau_max = 0.5,
    grid_points = 10,
    R = 100,
    n = 100,
    beta_true = 1,
    degree = 3,
    seed = 123
) {
  set.seed(seed)
  
  grid_estimates <- numeric(R)
  
  for (r in 1:R) {
    # Draw sample from DGP
    df <- simulate_loc_shift(n = n)
    
    # Compute grid-averaged estimate
    grid_try <- try(cf_qr_grid_estimate(df$Y, df$D, df$Z, 
                                        tau_max = tau_max, 
                                        grid_points = grid_points, 
                                        degree = degree), 
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

mc_res = evaluate_grid_averaged_estimator(tau_max = 0.25, grid_points = 10, n=250, R=200)


## Bias Variance Trade of in number of grid points
grid_number_evaluation <- function(
    tau_max = 0.25,
    grid_points_vec = c(3, 5, 10, 25, 50, 100),
    n = 250,
    R = 200,
    degree = 3,
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
      beta_true = beta_true,
      degree = degree
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
# Also: for extreme tails in small samples QR are known to be biased 
res_plot_obj <- grid_number_evaluation(grid_points_vec = c(3,5,10,25,50,100))
res_plot_obj$plot  # to display the ggplot
res_plot_obj$results_table  # to inspect the raw numbers






### Then asymptotically efficient weighting (with variance estimates)
#A = identity matrix -> sum squared entries 
#c = see table 1 of Leorato et al (2012) paper
#u = 1/ grid points (if evenly spaced)


