source("/Users/niklaslienau/Dropbox/PhD/Research/ES_cont_treatment/dgp.R")
#data visualization

###Plot Counterfactual function vs Quantile Regression (Endogeniety)


plot_endogeniety<- function(tau = 0.25, beta = 1, gamma = 0.5, n = 1000,
                            seed = 123, scale = FALSE, cont = TRUE) {
  set.seed(seed)
  
  # Simulate data
  df <- simulate_dgp(n = n, beta = beta, gamma = gamma, scale = scale, cont = cont)
  
  # Fit quantile regression
  qreg <- rq(Y ~ D, tau = tau, data = df)
  
  # Grid for plotting
  if (cont) {
    D_vals <- seq(min(df$D), max(df$D), length.out = 100)
  } else {
    D_vals <- 0:15
  }
  
  # True counterfactual quantile function
  qnorm_tau <- qnorm(tau)
  if (scale) {
    true_line <- beta * D_vals + exp(gamma * D_vals) * qnorm_tau
  } else {
    true_line <- beta * D_vals + qnorm_tau
  }
  
  # Predicted quantile regression line
  est_line <- predict(qreg, newdata = data.frame(D = D_vals))
  
  # Plot
  plot(df$D, df$Y, pch = 16, col = rgb(0, 0, 0, 0.2),
       xlab = "D (Treatment)", ylab = "Y (Outcome)",
       main = paste0("Quantile τ = ", tau,
                     ifelse(scale, ", Location-Scale", ", Location-Shift"),
                     ifelse(cont, ", Continuous D", ", Discrete D")))
  
  lines(D_vals, true_line, col = "blue", lwd = 2, lty = 2)
  lines(D_vals, est_line, col = "red", lwd = 2)
  
  legend("topleft",
         legend = c("Counterfactual Quantile Function", "Estimated Quantile Regression"),
         col = c("blue", "red"), lty = c(2,1), lwd = 2, bty = "n")
}




#Plot it
# Open a quartz plotting device (on macOS)
quartz(width = 8, height = 8)

# 2x2 layout
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Adjust margins

# LOCATION SHIFT & CONTINIOUS D
plot_endogeniety(tau = 0.05, scale= FALSE, cont= TRUE)
plot_endogeniety(tau = 0.10, scale= FALSE, cont= TRUE)
plot_endogeniety(tau = 0.25, scale= FALSE, cont= TRUE)
plot_endogeniety(tau = 0.50, scale= FALSE, cont= TRUE)

# LOCATION SCALE & Continous D
plot_endogeniety(tau = 0.05, scale= TRUE, cont= TRUE)
plot_endogeniety(tau = 0.10, scale= TRUE, cont= TRUE)
plot_endogeniety(tau = 0.25, scale= TRUE, cont= TRUE)
plot_endogeniety(tau = 0.50, scale= TRUE, cont= TRUE)

# LOCATION SHIFT & Discrete D
plot_endogeniety(tau = 0.05, scale= FALSE, cont= FALSE)
plot_endogeniety(tau = 0.10, scale= FALSE, cont= FALSE)
plot_endogeniety(tau = 0.25, scale= FALSE, cont= FALSE)
plot_endogeniety(tau = 0.50, scale= FALSE, cont= FALSE)

# LOCATION SCALE & Discrete D
plot_endogeniety(tau = 0.05, scale= TRUE, cont= FALSE)
plot_endogeniety(tau = 0.10, scale= TRUE, cont= FALSE)
plot_endogeniety(tau = 0.25, scale= TRUE, cont= FALSE)
plot_endogeniety(tau = 0.50, scale= TRUE, cont= FALSE)

