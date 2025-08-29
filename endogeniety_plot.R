source("/Users/niklaslienau/Dropbox/PhD/Research/ES_cont_treatment/dgp.R")
#data visualization

###Plot Counterfactual quantile functions 



set.seed(1)
df <- simulate_dgp_rand(n = 2000)   # your function

taus <- c(0.10, 0.25, 0.50, 0.75)

# replicate the scatter data across panels (one copy per tau)
df_rep <- do.call(rbind, lapply(taus, function(t) {
  tmp <- df
  tmp$tau <- factor(t, levels = taus,
                    labels = sprintf("τ = %.2f", taus))
  tmp
}))

# make smooth d-grid for the structural quantile curves
d_grid <- seq(min(df$D), max(df$D), length.out = 400)
curve_df <- do.call(rbind, lapply(taus, function(t) {
  data.frame(
    d   = d_grid,
    q   = qnorm(t) * sqrt(d_grid^2 + 1),
    tau = factor(t, levels = taus,
                 labels = sprintf("τ = %.2f", taus))
  )
}))

# (optional) open a separate window on macOS
# quartz(width = 9, height = 7)

ggplot() +
  geom_point(data = df_rep, aes(x = D, y = Y),
             alpha = 0.25, size = 0.6) +
  geom_line(data = curve_df, aes(x = d, y = q),
            linewidth = 1, color= "blue") +
  facet_wrap(~ tau, ncol = 2) +
  labs(
    title = "D vs Y with Structural Quantile Curves",
    x = "D", y = "Y"
  ) +
  theme_minimal()
