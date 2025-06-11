source("/Users/niklaslienau/Dropbox/PhD/Research/ES_cont_treatment/dgp.R")
#data visualization

###Plot Counterfactual function vs Quantile Regression (Endogeniety)

plot_quantile_shift <- function(tau = 0.25, beta = 1, n = 1000, rho = 0.7, pi = 1, seed = 123) {
  set.seed(123)
  # Simulate data
  df <- simulate_loc_shift(n = n, beta = beta, pi = pi, rho = rho)
  
  # Fit linear quantile regression
  qreg <- rq(Y ~ D, tau = tau, data = df)
  
  # Grid for plotting lines
  D_seq <- seq(min(df$D), max(df$D), length.out = 100)
  
  # True quantile line (counterfactual)
  true_line <- beta * D_seq + qnorm(tau)
  
  # Estimated quantile line
  est_line <- predict(qreg, newdata = data.frame(D = D_seq))
  
  # Plot
  plot(df$D, df$Y, pch = 16, col = rgb(0, 0, 0, 0.2),
       xlab = "D (Treatment)", ylab = "Y (Outcome)",
       main = paste0("Quantile τ = ", tau))
  
  lines(D_seq, true_line, col = "blue", lwd = 2, lty = 2)
  lines(D_seq, est_line, col = "red", lwd = 2)
  
  legend("topleft", legend = c("counterfactual quantile function", "Estimated quantile regression"),
         col = c("blue", "red"), lty = c(2,1), lwd = 2, bty = "n")
}

#Plot it
# Open a quartz plotting device (on macOS)
quartz(width = 8, height = 8)

# 2x2 layout
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Adjust margins

# Plot for different quantiles
plot_quantile_shift(tau = 0.05)
plot_quantile_shift(tau = 0.10)
plot_quantile_shift(tau = 0.25)
plot_quantile_shift(tau = 0.50)

