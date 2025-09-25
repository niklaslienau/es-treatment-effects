# V(tau) shape under Normal errors:
# V(tau) ∝ tau*(1 - tau) / dnorm(qnorm(tau))^2

# grid (avoid exactly 0 and 1 to prevent Inf)
tau <- seq(1e-2, 1 - 1e-2, length.out = 2000)

z    <- qnorm(tau)
phi  <- dnorm(z)
Vtau <- tau * (1 - tau) / (phi^2)
V_tau_inv = 1/Vtau

# (optional) scale to mean 1 so the y-axis is easier to read
Vtau_scaled <- Vtau / mean(Vtau)
Vtau_inv_scaled = V_tau_inv / mean(V_tau_inv)

# Linear-scale plot
plot(tau, Vtau_scaled, type = "l", lwd = 2,
     xlab = expression(tau),
     ylab = expression(paste("Scaled  ", V(tau))),
     main = expression(paste("Asymptotic variance shape  V(tau)  (Normal errors)")))


# Linear-scale plot
plot(tau, Vtau_inv_scaled, type = "l", lwd = 2,
     xlab = expression(tau),
     ylab = expression(paste("Scaled Inverted ", V(tau))),
     main = expression(paste("Inverted Asymptotic variance shape  V(tau) ")))

# (optional) log-scale version—useful because tails explode
plot(tau, Vtau_scaled, type = "l", lwd = 2, log = "y",
     xlab = expression(tau),
     ylab = expression(paste("Scaled  ", V(tau), " (log y)")),
     main = expression(paste("V(tau) on log scale")))
abline(v = c(0.1, 0.9), col = "gray80", lty = 3)
