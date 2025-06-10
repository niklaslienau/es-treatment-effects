# Control function quantile regression with series expansion for residuals
cf_qr_estimate <- function(Y, D, Z, tau = 0.25, degree = 3) {
  
  # Step 1: First stage - quantile regression of D on Z
  first_stage <- rq(D ~ Z, tau = tau)
  e_hat <- resid(first_stage)
  
  # Step 2: Build polynomial series of residuals
  series_terms <- sapply(1:degree, function(p) e_hat^p)
  colnames(series_terms) <- paste0("resid_", 1:degree)
  series_df <- as.data.frame(series_terms)
  
  # Step 3: Combine D and series into a design matrix
  design_df <- cbind(D = D, series_df)
  
  # Step 4: Second stage - quantile regression of Y on D and residual polynomials
  formula_str <- paste("Y ~ D +", paste(colnames(series_df), collapse = " + "))
  second_stage <- rq(as.formula(formula_str), data = cbind(Y = Y, design_df), tau = tau)
  
  return(list(
    first_stage = first_stage,
    second_stage = second_stage,
    e_hat = e_hat,
    design_matrix = design_df
  ))
}


