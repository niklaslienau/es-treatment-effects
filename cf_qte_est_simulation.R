





################################################################
# Unified Monte Carlo framework for QTE estimators
### Histogram Plots and Summary Stats (MSE, Bias, Variance)

# Choose estimators from:  c("naiveQR", "cf_linear", "cf_series", "cf_kernel")
#Choose DGP from : c("random_coeff", "loc_scale")

run_mc_qte <- function(
    tau        = 0.5,
    R          = 1000,
    n          = 2000,
    seed       = 123,
    rho        = 0.5,
    dgp        = c("random_coeff", "loc_scale"),
    estimators = c("naiveQR", "cf_linear", "cf_series", "cf_kernel"),
    # controls for CF-series
    series_degree = 3,
    # controls for CF-kernel
    np_trim       = c(0.01, 0.99),
    np_grid_len   = 30,
    np_bw         = NULL,
    np_deg_V      = 3,
    # misc
    quiet = FALSE
) {
  stopifnot(length(tau) == 1, tau > 0, tau < 1)
  dgp        <- match.arg(dgp)
  estimators <- match.arg(estimators,
                          choices = c("naiveQR", "cf_linear", "cf_series", "cf_kernel"),
                          several.ok = TRUE)
  if (!quiet) message("Estimators: ", paste(estimators, collapse = ", "),
                      " | DGP: ", dgp, " | R = ", R, " | n = ", n,
                      " | tau = ", tau)
  
  set.seed(seed)
  
  # --- DGP + true QTE hooks ---
  gen_fun <- switch(dgp,
                    random_coeff = function(n, rho) simulate_dgp_rand(n = n, rho = rho),
                    loc_scale    = function(n, rho) simulate_dgp_locscale(n = n, rho = rho)
  )
  true_qte_fun <- switch(dgp,
                         random_coeff = function(tau) compute_qte_tau(tau),
                         loc_scale    = function(tau) comp_qte_locscale(tau)
  )
  qte_true <- true_qte_fun(tau)
  
  # --- Estimator wrappers (must return numeric scalar or NA) ---
  est_wrappers <- list(
    naiveQR = function(df) {
      fit <- try(quantreg::rq(Y ~ D, tau = tau, data = df), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        as.numeric(coef(fit)["D"])
      } else NA_real_
    },
    cf_linear = function(df) {
      val <- try(cf_qr_estimate(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      if (!inherits(val, "try-error") && is.numeric(val)) as.numeric(val) else NA_real_
    },
    cf_series = function(df) {
      val <- try(cf_qr_series_avg(df$Y, df$D, df$Z, tau = tau, degree = series_degree), silent = TRUE)
      if (!inherits(val, "try-error") && is.numeric(val)) as.numeric(val) else NA_real_
    },
    cf_kernel = function(df) {
      # average local slopes over a grid of d0 (your qte_avg_cf)
      val <- try(qte_avg_cf(
        Y = df$Y, D = df$D, Z = df$Z, tau = tau,
        trim = np_trim, grid_len = np_grid_len,
        bw = np_bw, deg_V = np_deg_V
      ), silent = TRUE)
      if (!inherits(val, "try-error") && is.numeric(val)) as.numeric(val) else NA_real_
    }
  )
  
  # keep only requested wrappers
  est_wrappers <- est_wrappers[estimators]
  
  # --- Monte Carlo loop ---
  # pre-allocate list of numeric vectors, one per estimator
  out_list <- lapply(est_wrappers, function(...) rep(NA_real_, R))
  
  for (r in seq_len(R)) {
    if (!quiet && (r %% max(1, floor(R/10)) == 0)) {
      message(sprintf("MC iteration %d / %d", r, R))
    }
    df <- gen_fun(n, rho)
    
    # compute each estimator on this replication
    for (m in names(est_wrappers)) {
      out_list[[m]][r] <- est_wrappers[[m]](df)
    }
  }
  
  # --- Stack to long data.frame for plotting/summaries ---
  df_long <- do.call(rbind, lapply(names(out_list), function(m) {
    data.frame(method = m, estimate = out_list[[m]], stringsAsFactors = FALSE)
  }))
  df_long <- df_long[is.finite(df_long$estimate), , drop = FALSE]
  
  # --- Summaries per estimator ---
  summarize_one <- function(x) {
    c(
      mean   = mean(x, na.rm = TRUE),
      sd     = sd(x, na.rm = TRUE),
      bias   = mean(x - qte_true, na.rm = TRUE),
      var    = var(x, na.rm = TRUE),
      mse    = mean((x - qte_true)^2, na.rm = TRUE),
      n_ok   = sum(is.finite(x))
    )
  }
  summary_df <- do.call(rbind, lapply(split(df_long$estimate, df_long$method), summarize_one))
  summary_df <- data.frame(method = rownames(summary_df), summary_df, row.names = NULL)
  
  # --- Plot distributions + true QTE line ---
  plt <- ggplot2::ggplot(df_long, ggplot2::aes(x = estimate, fill = method)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = qte_true, color = "True QTE", linetype = "True QTE"),
                        linewidth = 1) +
    ggplot2::scale_color_manual(name = "", values = c("True QTE" = "black"),
                                labels = expression(beta(tau))) +
    ggplot2::scale_linetype_manual(name = "", values = c("True QTE" = "dashed"),
                                   labels = expression(beta(tau))) +
    ggplot2::labs(
      title = sprintf("MC Distribution of QTE Estimators (τ = %.2f, DGP = %s, R=%d, n=%d)",
                      tau, dgp, R, n),
      x = "Estimate", y = "Density", fill = "Estimator"
    ) +
    ggplot2::theme_minimal()
  
  list(
    data      = df_long,   # long data frame with columns: method, estimate
    summary   = summary_df,
    true_qte  = qte_true,
    plot      = plt,
    settings  = list(
      tau = tau, R = R, n = n, seed = seed, rho = rho, dgp = dgp,
      estimators = estimators,
      series_degree = series_degree,
      np_trim = np_trim, np_grid_len = np_grid_len, np_bw = np_bw, np_deg_V = np_deg_V
    )
  )
}


# Example 1: random-coefficient DGP, all estimators
res1 <- run_mc_qte(
  tau = 0.1, R = 500, n = 5000, dgp = "random_coeff",
  estimators =c("naiveQR", "cf_linear", "cf_series", "cf_kernel"),
  series_degree = 3, np_trim = c(0.01, 0.99)
)
res1$summary
print(res1$plot)










