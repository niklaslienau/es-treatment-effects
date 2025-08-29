


##### CTATE ESTIMATORS ######

##### 1- Asymptotically Efficient Weights #########
# standard nonparametric bootstrap for covariance matrix
bootstrap_qte_variances <- function(Y, D, Z,  B = 200, qte_est, grid ) {
  n <- length(Y)
  tau_grid = grid
  J <- length(tau_grid)
  
  qte_matrix <- matrix(NA, nrow = B, ncol = J)
  
  for (b in 1:B) {
    idx <- sample(1:n, size = n, replace = TRUE)
    Yb <- Y[idx]; Db <- D[idx]; Zb <- Z[idx]
    
    for (j in seq_along(tau_grid)) {
      tau <- tau_grid[j]
      model <- try(cf_qr_series_avg(Yb, Db, Zb, tau = tau, degree = 3), silent = TRUE)
      
      if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) {
        qte_matrix[b, j] <- model
      }
    }
  }
  
  # Remove rows with too many NAs
  valid_rows <- complete.cases(qte_matrix)
  qte_matrix <- qte_matrix[valid_rows, , drop = FALSE]
  
  #Demean each vector of QTE estimates in every bootstrap iteration  by vector of initial sample estimates
  qte_demeaned <- sweep(qte_matrix, MARGIN = 2, STATS = qte_est, FUN = "-")
  
  # Initialize the sum
  V_hat <- matrix(0, nrow = J, ncol = J)
  
  # Loop over rows and sum outer products
  for (i in 1:B) {
    v_i <- qte_demeaned[i, ]
    V_hat <- V_hat + tcrossprod(v_i)  # tcrossprod is v_i %*% t(v_i)
  }
  
  # Divide by B
  V_hat <- V_hat / B
  
  return(list(
    tau_grid = tau_grid,
    cov_estimate = V_hat,
    bootstrap_estimates = qte_matrix
  ))
}

#Efficient weighting estimator with boostraped cov matrix
efficient_es_est <- function(Y, D, Z, tau_max = 0.25, grid_points = 8, B = 200, c = 1e-2, degree=3) {
 
  # Step 0 : Set Up Tau Grid
  n <- length(Y)
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  tau_grid <- tau_grid[tau_grid > 0 & tau_grid < 1]  # avoid endpoints
  J <- length(tau_grid)

  
  # Step 1: Estimate QTEs at each tau
  qte_estimates <- numeric(J)
  for (g in seq_along(tau_grid)) {
    tau <- tau_grid[g]
    model <- try(cf_qr_series_avg(Y, D, Z, tau = tau, degree = degree), silent = TRUE)
    qte_estimates[g] <- if (!inherits(model, "try-error") && is.numeric(model) && !is.na(model)) model else NA
  }
  
  # Step 2: Bootstrap to estimate variances and covariance matrix
  boot <- bootstrap_qte_variances(Y, D, Z, B = B, qte_est = qte_estimates, grid= tau_grid)
  cov_matrix <- boot$cov_estimate
  
  
  
  # Step 3: Remove invalid entries (NA)
  if (anyNA(qte_estimates)) {
    message("Found some NA values in estimated qte vector ")
  } # print message
  valid <- which(!is.na(qte_estimates))
  qte_estimates <- qte_estimates[valid]
  tau_grid <- tau_grid[valid]
  cov_matrix <- cov_matrix[valid, valid, drop = FALSE]
  G <- length(valid)
  u <- rep(1 / G, G)
  
  # Step 4: Solve QP for efficient weights
  A = cov_matrix + c * G
  
  # Define the quadratic and linear terms for quadprog
  Dmat <- 2 * A
  dvec <- 2 * c * u
  
  # Constraint: sum(w) = 1 --> A_eq w = b_eq
  Amat <- matrix(1, nrow = G, ncol = 1)  # 1^T w = 1
  bvec <- 1
  #Solve
  solution <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  # Extract solution
  weights <- solution$solution
  
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


##### 2- Deterministic Weights #########


##### 3- Status Quo Two Step Estimator  #########

two_step_es_est <- function(Y, D, tau) {
  # First Step
  first_stage <- rq(Y ~ D, tau = tau)
  pred_quant <- first_stage$fitted.values
  
  # Logical mask: which obs are below the predicted quantile
  idx <- Y < pred_quant
  
  # Subset both Y and D using the mask
  Y_sub <- Y[idx]
  D_sub <- D[idx]
  
  # Run OLS on the subset
  second_stage <- lm(Y_sub ~ D_sub)
  es_estimate <- coef(second_stage)[2]
  
  return(es_estimate)
}


