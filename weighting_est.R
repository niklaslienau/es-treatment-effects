library(quadprog)
############################################### UNIFORM ################################

# Weighted Quantile Estimator (UNIFORM)
uniform_es_est <- function(Y, D,Z,  tau_max, grid_points) {
  # Create evenly spaced grid from 0 to tau_max
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  
  # Remove tau = 0 if it's included (quantile regression undefined at 0)
  tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
  
  estimates <- numeric(length(tau_grid))
  
  # Estimate at each grid point
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    cf_try <- try(cf_qr_estimate(Y, D, Z, tau = tau, return_model = TRUE), silent = TRUE)
    
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





############################################### Efficient Weighting ################################


### asymptotically efficient weighting (with  bootstrapped variance)
#A = identity matrix 
#c = calibration similar to table 1 of Leorato et al (2012) paper
#u = 1/ grid points (because evenly spaced)






#efficient weighting estimator no bootstrap
#->variance downward biased and covariance ignored
estimate_weighted_es_cf <- function(Y, D, Z, tau_max = 0.25, grid_points = 10, c = 1e-2) {
  # Step 1: Create tau grid
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]
  G <- length(tau_grid)
  
  # Step 2: Run cf_qr_estimate and store QTEs + variances
  qte_estimates <- numeric(G)
  variances <- numeric(G)
  
  for (g in seq_along(tau_grid)) {
    tau <- tau_grid[g]
    model <- try(cf_qr_estimate(Y, D, Z, tau = tau, return_model = TRUE), silent = TRUE)
    
    if (inherits(model, "try-error")) {
      qte_estimates[g] <- NA
      variances[g] <- NA
    } else {
      qte <- try(coef(model)["D"], silent = TRUE)
      se <- try(summary(model, se = "nid")$coefficients["D", "Std. Error"], silent = TRUE)
      
      qte_estimates[g] <- if (!inherits(qte, "try-error")) qte else NA
      variances[g] <- if (!inherits(se, "try-error") && !is.na(se)) se^2 else NA
    }
  }
  
  # Step 3: Remove failed fits
  valid <- which(!is.na(qte_estimates) & !is.na(variances))
  qte_estimates <- qte_estimates[valid]
  variances <- variances[valid]
  tau_grid <- tau_grid[valid]
  G <- length(valid)
  u <- rep(1 / G, G)
  
  # Step 4: Solve QP for optimal weights
  Dmat <- diag(variances + c)
  dvec <- c * u
  Amat <- matrix(1, nrow = G, ncol = 1)
  bvec <- 1
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  weights <- sol$solution
  
  # Step 5: Weighted average to estimate ES
  ES_estimate <- sum(weights * qte_estimates)
  
  return(list(
    estimate = ES_estimate,
    weights = weights,
    tau_grid = tau_grid,
    qte_estimates = qte_estimates,
    variances = variances
  ))
}






# standard nonparametric bootstrap for covariance matrix
bootstrap_qte_variances <- function(Y, D, Z, tau_max = 0.25, grid_points = 10, B = 200) {
  n <- length(Y)
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]  # avoid endpoints
  G <- length(tau_grid)
  
  qte_matrix <- matrix(NA, nrow = B, ncol = G)
  
  for (b in 1:B) {
    idx <- sample(1:n, size = n, replace = TRUE)
    Yb <- Y[idx]; Db <- D[idx]; Zb <- Z[idx]
    
    for (g in seq_along(tau_grid)) {
      tau <- tau_grid[g]
      model <- try(cf_qr_estimate(Yb, Db, Zb, tau = tau), silent = TRUE)
      
      if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) {
        qte_matrix[b, g] <- model
      }
    }
  }
  
  # Remove rows with too many NAs
  valid_rows <- complete.cases(qte_matrix)
  qte_matrix <- qte_matrix[valid_rows, , drop = FALSE]
  
  cov_matrix <- cov(qte_matrix, use = "pairwise.complete.obs")
  
  return(list(
    tau_grid = tau_grid,
    cov_matrix = cov_matrix,
    qte_matrix = qte_matrix
  ))
}


#Efficient weighting estimator with boostraped cov matrix
efficient_es_est <- function(Y, D, Z, tau_max = 0.25, grid_points = 10, B = 200, c = 1e-2) {
  # Step 1: Bootstrap to estimate variances and covariance matrix
  boot <- bootstrap_qte_variances(Y, D, Z, tau_max = tau_max, grid_points = grid_points, B = B)
  tau_grid <- boot$tau_grid
  cov_matrix <- boot$cov_matrix
  G <- length(tau_grid)
  
  # Step 2: Estimate QTEs at each tau
  qte_estimates <- numeric(G)
  for (g in seq_along(tau_grid)) {
    tau <- tau_grid[g]
    model <- try(cf_qr_estimate(Y, D, Z, tau = tau), silent = TRUE)
    qte_estimates[g] <- if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) model else NA
  }
  
  # Step 3: Remove invalid entries (NA)
  valid <- which(!is.na(qte_estimates))
  qte_estimates <- qte_estimates[valid]
  tau_grid <- tau_grid[valid]
  cov_matrix <- cov_matrix[valid, valid, drop = FALSE]
  G <- length(valid)
  u <- rep(1 / G, G)
  
  # Step 4: Solve QP for efficient weights
  Dmat <- cov_matrix + c * diag(G)
  dvec <- c * u
  Amat <- matrix(1, nrow = G, ncol = 1)
  bvec <- 1
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  weights <- sol$solution
  
  # Step 5: Compute final weighted estimate
  ES_estimate <- sum(weights * qte_estimates)
  
  return(list(
    estimate = ES_estimate,
    weights = weights,
    tau_grid = tau_grid,
    qte_estimates = qte_estimates,
    cov_matrix = cov_matrix
  ))
}




