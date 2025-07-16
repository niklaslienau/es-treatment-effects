source("/Users/niklaslienau/Dropbox/PhD/Research/ES_cont_treatment/dgp.R")
#data visualization

###Plot Counterfactual function vs Quantile Regression (Endogeniety)

plot_endogeniety <- function(tau = 0.25, beta = 1, gamma = 0.5, n = 1000 ) {
  
  # Simulate data
  df <- sim_dgp_homo()
  
  # Fit quantile regression
  qreg <- rq(Y ~ D, tau = tau, data = df)
  
  #Grid for plotting
  D_vals= seq(min(df$D), max(df$D), length.out = 100)
  # True structural quantile function 
 
  qnorm_tau <- qnorm(tau)
  # Y(d) = βd + u  → Q = βd + Φ⁻¹(τ)
  
  true_line <- beta * D_vals + qnorm_tau
    
  
  # Fitted quantile regression line
  est_line <- predict(qreg, newdata = data.frame(D = D_vals))
  
  # Plot
  plot(df$D, df$Y, pch = 16, col = rgb(0, 0, 0, 0.2),
       xlab = "D (Treatment)", ylab = "Y (Outcome)",
       main = paste0("Quantile τ = ", tau
                     ))
  
  lines(D_vals, true_line, col = "blue", lwd = 2, lty = 2)
  lines(D_vals, est_line, col = "red", lwd = 2)
  
  legend("topleft",
         legend = c("Counterfactual Quantile Function", "Estimated Quantile Regression"),
         col = c("blue", "red"), lty = c(2, 1), lwd = 2, bty = "n")
}




#Plot it
# Open a quartz plotting device (on macOS)
quartz(width = 8, height = 8)

# 2x2 layout
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Adjust margins

#
plot_endogeniety(tau = 0.05)
plot_endogeniety(tau = 0.10)
plot_endogeniety(tau = 0.25)
plot_endogeniety(tau = 0.50)


