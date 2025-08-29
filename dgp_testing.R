# This file explores potential problem with the Location Scale DGP

### ANALYZE ENDOGENIETY; HETEROSKEDASTICITY AND EXLCUSION RESTRICTION




#### ENDGENIETY ######
# Assuming simulate_dgp() returns a data frame with D and v
df_ne <- simulate_dgp(rho=0)
df_e  <- simulate_dgp(rho=0.5)

# Add labels to each dataset
df_ne$type <- "No Endogeneity (ρ = 0)"
df_e$type  <- "Endogeneity (ρ = 0.7)"

# Combine into one dataframe
df_all <- bind_rows(df_ne, df_e)

# Plot
ggplot(df_all, aes(x = D, y = u)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~type) +
  labs(
    title = "",
    x = "D",
    y = "u"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  )





#### Exlcusion restricition depends on endogeniety ######
library(ggplot2)
library(dplyr)

# Simulate data
df_ne <- simulate_dgp(rho=0)
df_e  <- simulate_dgp(rho=0.5)

# Label each
df_ne$type <- "No Endogeneity (ρ = 0)"
df_e$type  <- "Endogeneity (ρ = 0.7)"

# Combine
df_all <- bind_rows(df_ne, df_e)

# Compute covariance for annotation
cov_labels <- df_all %>%
  group_by(type) %>%
  summarise(
    cov_zv = cov(Z, u),
    xpos = min(Z) + 0.1 * (max(Z) - min(Z)),
    ypos = max(u) - 0.1 * (max(u) - min(u)),
    label = paste0("Cov(Z, v) = ", round(cov_zv, 3))
  )

# Plot
ggplot(df_all, aes(x = Z, y = u)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~type) +
  geom_text(
    data = cov_labels,
    aes(x = xpos, y = ypos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "italic"
  ) +
  labs(
    title = "Instrument Exogeneity Depends on ρ",
    x = "Instrument Z",
    y = "Structural Error u"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold")
  )



#### ONLY FOR LOCATION SCALE ###############

######### linear Heteroskedasticity with NORMAL ERROR  kills Endogeniety because D is dist symetrically around 0 ############
# -> Go uniform 

#why sigma needs abs(D) not D
#df= simulate_dgp()
df = simulate_dgp(error_var = 0.2)
# Plot 1: u vs. D*u
p1 <- ggplot(df, aes(x = D, y = D*u+ epsilon)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "v = D*u", x = "D", y = "v") +
  theme_minimal()

# Plot 2: u vs. |D|*u
p2 <- ggplot(df, aes(x = D, y = abs(D)*u + epsilon)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "v= abs(D) *u", x = "D", y = "v") +
  theme_minimal()

# Show both plots side by side
gridExtra::grid.arrange(p1, p2, ncol = 2)





