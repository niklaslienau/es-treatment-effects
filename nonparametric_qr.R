

# Local QTE at d0 using:
#  - control function residuals Vhat from median QR of D on Z
#  - bandwidth h chosen by silvermann plug in rule,  unless provided
#  - kernel weights from npksum (use kw)
#  - weighted QR: Y ~ (D - d0) + poly(Vhat, degree = deg_V, raw = TRUE)
qr_cf_local<- function(Y, D, Z,
                                tau, d0,
                                deg_V   = 3,
                                tau_fs  = 0.5,
                                ckertype = "gaussian",
                                ckerorder = 2,
                                bw = NULL) {
  stopifnot(length(Y) == length(D), length(D) == length(Z))
  # 1) control function residuals (first stage)
  fs <- rq(D ~ Z, tau = tau_fs)
  Vhat <- resid(fs)
  
  # 2) bandwidth selection for now with silvermann/LL plug in Rule 
  if (is.null(bw)) {
    h= 1.06 * sd(D) * length(D)^(-1/5)
  } else {
    h <- bw
  }
  
  # 3) kernel weights at evaluation point d0 from npksum
  wobj <- npksum(bws = h,
                 txdat = data.frame(D),
                 exdat = data.frame(d0),
                 ckertype = ckertype,
                 ckerorder = ckerorder, return.kernel.weights=TRUE)
  w <- as.numeric(wobj$kw)    # weights for each observation
  
  # 4) local linear in D (centered) + polynomial in Vhat up to deg_V
  Dc <- D - d0 # center around d0
  # build polynomial columns for Vhat as a matrix, degree 0..deg_V
  make_poly_V <- function(Vhat, deg_V = 3) {
    n <- length(Vhat)
    if (deg_V == 0) {
      M <- matrix(1, nrow = n, ncol = 1)
      colnames(M) <- "V0"
      return(M)
    }
    M <- cbind(1, outer(Vhat, 1:deg_V, `^`))  # V^1, V^2, ..., V^deg_V
    colnames(M) <- paste0("V", 0:deg_V)
    M
  }
  PV <- make_poly_V(Vhat, deg_V)
  PV <- PV[, -1, drop = FALSE]  #drop constant because qr adds one by default
  X <- data.frame(Dc = Dc, PV)
  
  # 5) weighted quantile regression; slope on Dc is ∂Q/∂d at d0
  fit <- rq(Y ~ ., tau = tau, weights = w,data = data.frame(Y = Y, X))
  qte_d=unname(coef(fit)["Dc"])
  return(list(qte = qte_d, h = h, n_eff = sum(w > 1e-6)))
}




##Global averaged QTE

qte_avg_cf <- function(Y, D, Z, tau,
                       d_grid = NULL,
                       trim = c(0.1, 0.9),
                       grid_len = 15,
                       bw = NULL,
                       deg_V = 3) {
  # pick evaluation points
  if (is.null(d_grid)) {
    qs <- quantile(D, probs = seq(trim[1], trim[2], length.out = grid_len), names = FALSE)
    d_grid <- as.numeric(qs)
  }
  
  # choose/reuse bandwidth (Silverman/LL if not supplied)
  if (is.null(bw)) {
    bw <- 1.06 * sd(D) * length(D)^(-1/5)
  }
  
  # compute local slopes on the grid
  slopes <- sapply(d_grid, function(d0) {
    res <- qr_cf_local(Y = Y, D = D, Z = Z,
                       tau = tau, d0 = d0,
                       deg_V = deg_V,
                       bw = bw)
    res$qte   # extract only the slope
  })
  
  # simple average over grid 
  return(mean(slopes, na.rm = TRUE))
}


#Example Use
df <- simulate_dgp_rand(n = 10000)
qte_tau25 <- qte_avg_cf(df$Y, df$D, df$Z, tau = 0.25)
qte_tau50 <- qte_avg_cf(df$Y, df$D, df$Z, tau = 0.50)
qte_tau75 <- qte_avg_cf(df$Y, df$D, df$Z, tau = 0.75)
qte_tau25


#### Series Estimator


# Build a raw polynomial basis (always returns an n x degree matrix)
build_poly <- function(x, degree) {
  if (degree <= 0) stop("degree must be >= 1")
  M <- outer(x, 1:degree, `^`)   # columns: x^1, x^2, ..., x^degree
  colnames(M) <- paste0(deparse(substitute(x)), "_", 1:degree)
  M
}

# Control-function QR with series in D and in residuals (no interactions)
cf_qr_series_local <- function(Y, D, Z, tau = 0.5, degree = 3, return_model = FALSE) {
  # 1) First stage (median)
  fs   <- rq(D ~ Z, tau = 0.5)
  ehat <- resid(fs)
  
  # 2) Series bases (force proper matrices)
  Dbasis <- build_poly(D,    degree)   # D_1, D_2, ...
  Vbasis <- build_poly(ehat, degree)   # ehat_1, ehat_2, ...
  
  design <- data.frame(Dbasis, Vbasis)
  
  # 3) Second stage QR
  fit <- rq(Y ~ ., tau = tau, data = data.frame(Y = Y, design))
  
  if (return_model) return(fit)
  
  # Extract coefficients on the D-polynomial only
  beta_D <- coef(fit)[colnames(Dbasis)]           
  
  # Return coefficients + a handy evaluator for the slope at any d0
  list(
    coef_D    = as.numeric(beta_D),
    qte_at    = function(d0) {
      p <- seq_len(degree)
      sum(p * as.numeric(beta_D) * d0^(p - 1))
    },
    model = fit
  )
}



# Average QTE using the series-based CF-QR (nonparametric in D)
# Assumes cf_qr_np() from before is in scope.
cf_qr_series_avg  <- function(Y, D, Z, tau, degree = 3,
                              d_vals = NULL,          # where to evaluate; default = all D
                              return_details = FALSE  # return vector of local QTEs too?
) {
  # 1) Fit the series CF-QR once
  fit_obj <- cf_qr_series_local(Y = Y, D = D, Z = Z, tau = tau, degree = degree)
  
  # 2) Choose evaluation points (default: every observed D)
  if (is.null(d_vals)) d_vals <- D
  
  # 3) Evaluate local QTE at each d0
  qte_vec <- vapply(d_vals, fit_obj$qte_at, numeric(1))
  
  # 4) Average
  qte_avg <- mean(qte_vec, na.rm = TRUE)
  
  if (!return_details) {
    return(qte_avg)
  } else {
    return(list(
      qte_avg = qte_avg,
      qte_vec = qte_vec,
      d_vals  = d_vals
    ))
  }
}


