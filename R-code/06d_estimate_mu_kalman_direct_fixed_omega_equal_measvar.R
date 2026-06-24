# ============================================================
# 06d_estimate_mu_kalman_direct_fixed_omega_equal_measvar.R
#
# Purpose:
# Estimate the CES consumption system with a state-space
# preference process using a direct-system Kalman approach.
#
# This version uses the original two measurement equations directly:
#   1. relative demand equation
#   2. CES aggregation equation
#
# Difference from 06c:
#   In 06c, var(epsilon_1) and var(epsilon_2) were estimated freely.
#   This allowed one measurement equation to dominate the likelihood.
#
#   In this version, the two measurement errors are restricted to have
#   the same variance:
#
#     var(epsilon_1) = var(epsilon_2)
#
#   This forces the relative demand equation and the CES aggregation
#   equation to receive equal measurement-error weight in the likelihood.
#
# State equations:
#   mu_t  = tau_t + c_t
#   tau_t = tau_{t-1} + eta_t
#   c_t   = rho * c_{t-1} + nu_t
#
# Fixed parameter:
#   omega = calibrated Paasche omega
#
# Estimated parameters:
#   sigma
#   rho
#   var(eta_t)
#   var(nu_t)
#   common var(epsilon_t)
#
# Estimation:
#   Maximum likelihood using an extended/direct Kalman filter.
#
# Main outputs:
#   Data/results_mu_kalman_direct_equal_measvar_1983_2017_paasche.csv
#   Data/residuals_mu_kalman_direct_equal_measvar_1983_2017_paasche.csv
#   Data/states_mu_kalman_direct_equal_measvar_1983_2017_paasche.csv
#
# Note:
# This version tests whether imposing equal measurement variance across
# the two equations improves identification of sigma relative to 06c.
# ============================================================


rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(minpack.lm)
library(gt)
library(webshot2)
library(gridExtra)
library(ggplot2)
library(grid)


setwd("P:/leawin/Forbrugselasticitet")

# ------------------------------------------------------------
# 1. Vælg samme periode som omega-filen
# ------------------------------------------------------------

start_year <- 1983
end_year   <- 2017

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

reg_file <- paste0(
  "Data/regression_data_",
  start_year, "_", end_year,
  "_paasche.csv"
)

reg_data <- read_csv(reg_file) %>%
  mutate(
    nest = tolower(nest),
    good1 = tolower(good1),
    good2 = tolower(good2),
    t = as.integer(t),
    C_obs = as.numeric(C_obs),
    P_obs = as.numeric(P_obs),
    V_obs = as.numeric(V_obs),
    C1 = as.numeric(C1),
    C2 = as.numeric(C2),
    P1 = as.numeric(P1),
    P2 = as.numeric(P2),
    V1 = as.numeric(V1),
    V2 = as.numeric(V2),
    omega = as.numeric(omega),
    sigma_old = as.numeric(sigma_old)
  )


# ------------------------------------------------------------
# 3. Hjælpefunktioner
# ------------------------------------------------------------

logit <- function(x) {
  log(x / (1 - x))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

safe_solve <- function(A) {
  
  # Sørg for symmetri
  A <- (A + t(A)) / 2
  
  jitters <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
  
  for (j in jitters) {
    out <- tryCatch(
      solve(A + diag(j, nrow(A))),
      error = function(e) NULL
    )
    
    if (!is.null(out) && all(is.finite(out))) {
      return(out)
    }
  }
  
  # Sidste fallback: pseudo-invers via eigen-decomposition
  eig <- eigen(A, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  vals_adj <- ifelse(abs(vals) < 1e-6, 1e-6, vals)
  
  out <- vecs %*% diag(1 / vals_adj, length(vals_adj)) %*% t(vecs)
  
  return(out)
}

safe_logdet <- function(A) {
  
  A <- (A + t(A)) / 2
  
  jitters <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
  
  for (j in jitters) {
    out <- tryCatch(
      determinant(A + diag(j, nrow(A)), logarithm = TRUE),
      error = function(e) NULL
    )
    
    if (!is.null(out)) {
      val <- as.numeric(out$modulus[1])
      sign_val <- out$sign
      
      if (is.finite(val) && sign_val > 0) {
        return(val)
      }
    }
  }
  
  return(NA_real_)
}

omega_from_raw <- function(omega_raw, omega_center, width = 0.15) {
  
  lower <- max(0.001, omega_center - width)
  upper <- min(0.999, omega_center + width)
  
  lower + (upper - lower) * inv_logit(omega_raw)
}

ces_log_and_derivative <- function(mu, C1, C2, omega, rho) {
  
  if (abs(rho) < 1e-6) {
    
    # Cobb-Douglas limit
    log_c_model <- omega * (mu + log(C1)) + (1 - omega) * log(C2)
    d_logc_d_mu <- omega
    
  } else {
    
    term1 <- omega * (exp(mu) * C1)^rho
    term2 <- (1 - omega) * C2^rho
    A <- term1 + term2
    
    if (!is.finite(A) || A <= 0) {
      return(
        list(
          log_c_model = NA_real_,
          d_logc_d_mu = NA_real_
        )
      )
    }
    
    log_c_model <- (1 / rho) * log(A)
    
    # Derivative of log(C_model) with respect to mu
    d_logc_d_mu <- term1 / A
  }
  
  list(
    log_c_model = log_c_model,
    d_logc_d_mu = d_logc_d_mu
  )
}

# ------------------------------------------------------------
# 4. Kalman filter likelihood: direct system
# ------------------------------------------------------------

kalman_loglik_one_nest <- function(theta, d) {
  
  # ----------------------------------------------------------
  # Transform parameters to valid ranges
  # ----------------------------------------------------------
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma <- sigma_min + (sigma_max - sigma_min) * inv_logit(theta["sigma_raw"])
  omega <- as.numeric(unique(d$omega)[1])
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e   <- sd_min + exp(theta["log_sd_e"])
  
  rho <- (sigma - 1) / sigma
  
  n <- nrow(d)
  
  # ----------------------------------------------------------
  # State system
  # tau_t = tau_{t-1} + eta_t
  # c_t   = rho_c * c_{t-1} + nu_t
  # ----------------------------------------------------------
  
  Tmat <- matrix(
    c(1, 0,
      0, rho_c),
    nrow = 2,
    byrow = TRUE
  )
  
  Q <- matrix(
    c(sd_eta^2, 0,
      0, sd_nu^2),
    nrow = 2,
    byrow = TRUE
  )
  
  # Equal measurement variance
  R <- matrix(
    c(sd_e^2, 0,
      0, sd_e^2),
    nrow = 2,
    byrow = TRUE
  )
  
  # Initial state
  a <- matrix(c(0, 0), nrow = 2)
  P <- diag(c(100, 100))
  
  loglik <- 0
  
  for (i in seq_len(n)) {
    
    mu_pred <- as.numeric(a[1] + a[2])
    
    # --------------------------------------------------------
    # Measurement equation 1: relative demand
    # --------------------------------------------------------
    
    y1_obs <- log((d$P1[i] * d$C1[i]) / (d$P2[i] * d$C2[i]))
    
    y1_model <- log(omega / (1 - omega)) +
      rho * (mu_pred + log(d$C1[i] / d$C2[i]))
    
    # --------------------------------------------------------
    # Measurement equation 2: CES aggregator directly
    # --------------------------------------------------------
    
    y2_obs <- log(d$C_obs[i])
    
    ces_i <- ces_log_and_derivative(
      mu = mu_pred,
      C1 = d$C1[i],
      C2 = d$C2[i],
      omega = omega,
      rho = rho
    )
    
    y2_model <- ces_i$log_c_model
    d_y2_d_mu <- ces_i$d_logc_d_mu
    
    if (
      !is.finite(y1_obs) ||
      !is.finite(y2_obs) ||
      !is.finite(y1_model) ||
      !is.finite(y2_model) ||
      !is.finite(d_y2_d_mu)
    ) {
      return(1e10)
    }
    
    y <- matrix(c(y1_obs, y2_obs), nrow = 2)
    h <- matrix(c(y1_model, y2_model), nrow = 2)
    
    # Extended Kalman observation matrix
    # derivatives wrt tau and c, since mu = tau + c
    H <- matrix(
      c(rho, rho,
        d_y2_d_mu, d_y2_d_mu),
      nrow = 2,
      byrow = TRUE
    )
    
    v <- y - h
    Fmat <- H %*% P %*% t(H) + R
    
    logdetF <- safe_logdet(Fmat)
    Finv <- safe_solve(Fmat)
    
    if (any(!is.finite(Finv)) || !is.finite(logdetF)) {
      return(1e10)
    }
    
    quad <- as.numeric(t(v) %*% Finv %*% v)
    
    if (!is.finite(quad) || quad < 0) {
      return(1e10)
    }
    
    loglik <- loglik - 0.5 * (
      logdetF +
        quad +
        2 * log(2 * pi)
    )
    
    # Update
    K <- P %*% t(H) %*% Finv
    a <- a + K %*% v
    P <- P - K %*% H %*% P
    P <- (P + t(P)) / 2 + diag(1e-8, nrow(P))
    
    # Predict next
    a <- Tmat %*% a
    P <- Tmat %*% P %*% t(Tmat) + Q
    P <- (P + t(P)) / 2 + diag(1e-8, nrow(P))
  }
  
  return(-as.numeric(loglik))
}
# ------------------------------------------------------------
# 5. Estimér ét nest med Kalman ML
# ------------------------------------------------------------

estimate_kalman_one_nest <- function(data, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t)
  
  sigma_start <- as.numeric(unique(d$sigma_old)[1])
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma_start <- min(max(sigma_start, sigma_min + 1e-4), sigma_max - 1e-4)
  
  sigma_raw_start <- logit(
    (sigma_start - sigma_min) / (sigma_max - sigma_min)
  )
  
  
  theta_start <- c(
    sigma_raw = sigma_raw_start,
    rho_c_raw = 0,
    log_sd_eta = log(0.05),
    log_sd_nu  = log(0.05),
    log_sd_e   = log(0.05)
  )
  
  fit <- optim(
    par = theta_start,
    fn = kalman_loglik_one_nest,
    d = d,
    method = "BFGS",
    control = list(
      maxit = 2000,
      reltol = 1e-9
    )
  )
  
  theta <- fit$par
  
  sigma <- sigma_min + (sigma_max - sigma_min) * inv_logit(theta["sigma_raw"])
  omega <- as.numeric(unique(d$omega)[1])
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e   <- sd_min + exp(theta["log_sd_e"])
  
  summary <- tibble(
    nest = nest_name,
    sigma = sigma,
    omega = omega,
    rho_c = rho_c,
    var_eta = sd_eta^2,
    var_nu = sd_nu^2,
    var_e = sd_e^2,
    neg_loglik = fit$value,
    converged = fit$convergence == 0
  )
  
  list(
    summary = summary,
    theta = theta,
    data = d
  )
}

# ------------------------------------------------------------
# 6. Filtrer states og beregn residualer
# ------------------------------------------------------------

filter_states_kalman_one_nest <- function(fit_object) {
  
  d <- fit_object$data
  theta <- fit_object$theta
  nest_name <- unique(d$nest)[1]
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma <- sigma_min + (sigma_max - sigma_min) * inv_logit(theta["sigma_raw"])
  omega <- as.numeric(unique(d$omega)[1])
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e   <- sd_min + exp(theta["log_sd_e"])
  
  rho <- (sigma - 1) / sigma
  
  Tmat <- matrix(
    c(1, 0,
      0, rho_c),
    nrow = 2,
    byrow = TRUE
  )
  
  Q <- matrix(
    c(sd_eta^2, 0,
      0, sd_nu^2),
    nrow = 2,
    byrow = TRUE
  )
  
  R <- matrix(
    c(sd_e^2, 0,
      0, sd_e^2),
    nrow = 2,
    byrow = TRUE
  )
  
  n <- nrow(d)
  
  a <- matrix(c(0, 0), nrow = 2)
  P <- diag(c(100, 100))
  
  out <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    mu_pred <- as.numeric(a[1] + a[2])
    
    y1_obs <- log((d$P1[i] * d$C1[i]) / (d$P2[i] * d$C2[i]))
    
    y1_model <- log(omega / (1 - omega)) +
      rho * (mu_pred + log(d$C1[i] / d$C2[i]))
    
    y2_obs <- log(d$C_obs[i])
    
    ces_i <- ces_log_and_derivative(
      mu = mu_pred,
      C1 = d$C1[i],
      C2 = d$C2[i],
      omega = omega,
      rho = rho
    )
    
    y2_model <- ces_i$log_c_model
    d_y2_d_mu <- ces_i$d_logc_d_mu
    
    y <- matrix(c(y1_obs, y2_obs), nrow = 2)
    h <- matrix(c(y1_model, y2_model), nrow = 2)
    
    H <- matrix(
      c(rho, rho,
        d_y2_d_mu, d_y2_d_mu),
      nrow = 2,
      byrow = TRUE
    )
    
    v <- y - h
    Fmat <- H %*% P %*% t(H) + R
    Finv <- safe_solve(Fmat)
    
    K <- P %*% t(H) %*% Finv
    a <- a + K %*% v
    P <- P - K %*% H %*% P
    P <- (P + t(P)) / 2 + diag(1e-8, nrow(P))
    
    tau <- as.numeric(a[1])
    c_comp <- as.numeric(a[2])
    mu <- tau + c_comp
    
    # Recompute model values using updated state
    foc_model_updated <- log(omega / (1 - omega)) +
      rho * (mu + log(d$C1[i] / d$C2[i]))
    
    ces_updated <- ces_log_and_derivative(
      mu = mu,
      C1 = d$C1[i],
      C2 = d$C2[i],
      omega = omega,
      rho = rho
    )
    
    ces_model_updated <- ces_updated$log_c_model
    
    out[[i]] <- tibble(
      nest = nest_name,
      t = d$t[i],
      sigma = sigma,
      omega = omega,
      rho_c = rho_c,
      tau = tau,
      c = c_comp,
      mu = mu,
      foc_obs = y1_obs,
      foc_model = foc_model_updated,
      resid_foc = y1_obs - foc_model_updated,
      ces_obs = y2_obs,
      ces_model = ces_model_updated,
      resid_ces = y2_obs - ces_model_updated
    )
    
    # Predict next
    a <- Tmat %*% a
    P <- Tmat %*% P %*% t(Tmat) + Q
    P <- (P + t(P)) / 2 + diag(1e-8, nrow(P))
  }
  
  residuals <- bind_rows(out)
  
  states <- residuals %>%
    select(nest, t, tau, c, mu)
  
  list(
    residuals = residuals,
    states = states
  )
}

# ------------------------------------------------------------
# 7. Kør Kalman ML for alle nests
# ------------------------------------------------------------

kalman_fits <- lapply(unique(reg_data$nest), function(n) {
  estimate_kalman_one_nest(reg_data, n)
})

results_kalman <- bind_rows(
  lapply(kalman_fits, function(x) x$summary)
)

kalman_filtered <- lapply(kalman_fits, filter_states_kalman_one_nest)

residuals_kalman <- bind_rows(
  lapply(kalman_filtered, function(x) x$residuals)
)

states_kalman <- bind_rows(
  lapply(kalman_filtered, function(x) x$states)
)

results_kalman
residuals_kalman
states_kalman

# ------------------------------------------------------------
# 8. Gem resultater
# ------------------------------------------------------------

write_csv(
  results_kalman,
  paste0(
    "Data/results_mu_kalman_direct_equal_measvar_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

write_csv(
  residuals_kalman,
  paste0(
    "Data/residuals_mu_kalman_direct_equal_measvar_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

write_csv(
  states_kalman,
  paste0(
    "Data/states_mu_kalman_direct_equal_measvar_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

# ------------------------------------------------------------
# 9. Plots
# ------------------------------------------------------------

ggplot(states_kalman, aes(x = t, y = mu)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Kalman-estimeret præferencekomponent mu_t",
    x = "År",
    y = "mu_t"
  )

ggplot(states_kalman, aes(x = t, y = tau)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Kalman-estimeret trendkomponent tau_t",
    x = "År",
    y = "tau_t"
  )

ggplot(states_kalman, aes(x = t, y = c)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Kalman-estimeret cyklisk komponent c_t",
    x = "År",
    y = "c_t"
  )

ggplot(residuals_kalman, aes(x = t, y = resid_foc)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra relativ demand equation, Kalman system",
    x = "År",
    y = "FOC residual"
  )

ggplot(residuals_kalman, aes(x = t, y = resid_ces)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra CES-aggregator, Kalman system",
    x = "År",
    y = "CES residual"
  )
