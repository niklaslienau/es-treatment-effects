
#CTATE estimator with uniform weights
#Returns avg qte over empirical distribution of D

#tau_max = CTATE tau level
#tau_min = truncation level
#degree= polynomial degree for D and control function
ctate_estimator_avg <- function(Y, D, Z,
                                tau_max = 0.25,
                                tau_min = 0.01,
                                grid_points = 50,
                                degree = 3) {
  # 1) Tau grid from 0 to tau_max
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  
  # 2) Estimate QTE at tau_min once (for imputation)
  qte_tau_min <- cf_qr_series_avg(Y, D, Z, tau = tau_min, degree = degree)
  
  # 3) Fill QTEs on the grid
  qte_grid <- vapply(tau_grid, function(tau) {
    if (tau < tau_min) {
      qte_tau_min
    } else {
      cf_qr_series_avg(Y, D, Z, tau = tau, degree = degree)
    }
  }, numeric(1))
  
  # 4) CTATE ≈ average of grid QTEs
  ctate_estimate <- mean(qte_grid)
  
  return(list(
    estimate    = ctate_estimate,
    tau_grid    = tau_grid,
    qte_grid    = qte_grid
  ))
}






