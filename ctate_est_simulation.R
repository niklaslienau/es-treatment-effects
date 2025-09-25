
## Monte Carlo Simulations for weighted estimator
### Histogram Plots 

#----------------------
# Flexible MC plotting for any subset of ES/CTATE estimators
# estimators: character vector among c("naive_joint","two_step","efficient","ctate_unif")
plot_mc_distribution_es_estimators <- function(
    tau = 0.5, R = 100, n = 1000, seed = 123,
    tau_min = 0.01, grid_points = 10,
    estimators = c("naive_joint", "two_step", "efficient", "ctate_unif")) {
  
  set.seed(seed)
  
  # --- Registry: how to compute each estimator from a simulated df ---
  # Each entry returns a single numeric estimate.
  estimator_funs <- list(
    naive_joint = function(df) {
      fit <- esreg(df$Y ~ df$D, alpha = tau)
      as.numeric(fit$coefficients[4])
    },
    two_step = function(df) {
      as.numeric(two_step_es_est(df$Y, df$D, tau = tau))
    },
    efficient = function(df) {
      # Try common return patterns: numeric OR list with $estimate
      res <- try(efficient_es_est(df$Y, df$D, df$Z, tau = tau), silent = TRUE)
      if (inherits(res, "try-error")) {
        as.numeric(NA)
      } else if (is.list(res) && !is.null(res$estimate)) {
        as.numeric(res$estimate)
      } else {
        as.numeric(res)
      }
    },
    ctate_unif = function(df) {
      # ctate_estimator_avg implements uniform weights + flat-tail imputation
      res <- ctate_estimator_avg(
        df$Y, df$D, df$Z,
        tau_max = tau, tau_min = tau_min, grid_points = grid_points
      )
      as.numeric(res$estimate)
    }
  )
  
  # Pretty labels for the legend
  estimator_labels <- c(
    naive_joint = "Naive Joint ES",
    two_step    = "Naive Two Step ES",
    efficient   = "Efficient ES",
    ctate_unif  = "CTATE (uniform)"
  )
  
  # Keep only valid/known estimators
  estimators <- intersect(estimators, names(estimator_funs))
  if (length(estimators) == 0) stop("No valid estimators selected.")
  
  # Storage
  est_mat <- matrix(NA_real_, nrow = R, ncol = length(estimators))
  colnames(est_mat) <- estimators
  
  # --- Monte Carlo loop ---
  for (r in seq_len(R)) {
    cat("Iteration", r, "of", R, "\n")
    df <- simulate_dgp_rand(n = n)
    for (k in seq_along(estimators)) {
      est_fun <- estimator_funs[[ estimators[k] ]]
      est_mat[r, k] <- est_fun(df)
    }
  }
  
  # True CTATE
  es_true <- compute_ctate_tau(tau = tau)
  
  # Long data frame for ggplot
  df_plot <- data.frame(
    estimate = as.vector(est_mat),
    method   = rep(estimators, each = R),
    stringsAsFactors = FALSE
  )
  df_plot$method <- factor(estimator_labels[df_plot$method], levels = estimator_labels[estimators])
  df_plot <- df_plot[!is.na(df_plot$estimate), ]
  
  # Plot
  ggplot(df_plot, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(aes(xintercept = es_true, linetype = "True CTATE"),
               color = "black", linewidth = 1) +
    scale_linetype_manual(name = "", values = c("True CTATE" = "dashed")) +
    labs(
      title = paste0("MC Distribution of CTATE/ES Estimators (τ = ", tau, ")"),
      x = "Estimate", y = "Density", fill = "Estimator"
    ) +
    theme_minimal()
}

quartz()
plot_mc_distribution_es_estimators(tau = 0.25, R = 200, n=1000, seed=123,tau_min = 0.01, grid_points = 25 ,estimators = c("naive_joint", "two_step", "ctate_unif") )

