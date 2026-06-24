# ============================================================
# 06b_estimate_mu_kalman_free_omega_bounded.R
#
# Purpose:
# Estimate the full CES consumption system with a state-space
# preference process using Kalman maximum likelihood.
#
# This version estimates omega jointly with sigma and the state-space
# variance parameters, but restricts omega to remain within a bounded
# interval around the calibrated Paasche value.
#
# The model combines:
#   1. relative demand equation
#   2. CES aggregation equation
#   3. state-space preference process
#
# State equations:
#   mu_t  = tau_t + c_t
#   tau_t = tau_{t-1} + eta_t
#   c_t   = rho * c_{t-1} + nu_t
#
# Estimated parameters:
#   sigma
#   omega, bounded around calibrated Paasche omega
#   rho
#   var(eta_t)
#   var(nu_t)
#   var(epsilon_1)
#   var(epsilon_2)
#
# Estimation:
#   Maximum likelihood using Kalman filter.
#
# Main outputs:
#   Data/results_mu_kalman_free_omega_bounded_1983_2017_paasche.csv
#   Data/residuals_mu_kalman_free_omega_bounded_1983_2017_paasche.csv
#   Data/states_mu_kalman_free_omega_bounded_1983_2017_paasche.csv
#
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

# ------------------------------------------------------------
# 4. Kalman filter likelihood
# ------------------------------------------------------------

kalman_loglik_one_nest <- function(theta, d) {
  
  # ----------------------------------------------------------
  # Transform parameters to valid ranges
  # ----------------------------------------------------------
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma <- sigma_min + (sigma_max - sigma_min) * inv_logit(theta["sigma_raw"])
  omega_center <- as.numeric(unique(d$omega)[1])
  omega <- omega_from_raw(
    omega_raw = theta["omega_raw"],
    omega_center = omega_center,
    width = 0.15
  )
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e1  <- sd_min + exp(theta["log_sd_e1"])
  sd_e2  <- sd_min + exp(theta["log_sd_e2"])
  
  rho <- (sigma - 1) / sigma
  
  # Avoid numerical problems close to Cobb-Douglas
  if (abs(rho) < 1e-5) {
    return(1e10)
  }
  
  C  <- d$C_obs
  C1 <- d$C1
  C2 <- d$C2
  P1 <- d$P1
  P2 <- d$P2
  
  # ----------------------------------------------------------
  # Measurement equation 1: relative demand / FOC
  #
  # log((P1*C1)/(P2*C2)) =
  #   log(omega/(1-omega)) + rho * (mu_t + log(C1/C2)) + e1
  #
  # Rearranged:
  # y1_t = rho * mu_t + e1
  # ----------------------------------------------------------
  
  y1 <- log((P1 * C1) / (P2 * C2)) -
    log(omega / (1 - omega)) -
    rho * log(C1 / C2)
  
  # ----------------------------------------------------------
  # Measurement equation 2: CES aggregator
  #
  # C_t^rho = omega * (exp(mu_t)*C1_t)^rho
  #           + (1-omega) * C2_t^rho
  #
  # Rearranged to implied mu:
  # y2_t = mu_t + e2
  # ----------------------------------------------------------
  
  inside <- (C^rho - (1 - omega) * C2^rho) / (omega * C1^rho)
  
  if (any(!is.finite(inside)) || any(inside <= 0)) {
    return(1e10)
  }
  
  y2 <- log(inside) / rho
  
  if (any(!is.finite(y1)) || any(!is.finite(y2))) {
    return(1e10)
  }
  
  # ----------------------------------------------------------
  # State system:
  #
  # mu_t = tau_t + c_t
  # tau_t = tau_{t-1} + eta_t
  # c_t = rho_c * c_{t-1} + nu_t
  # ----------------------------------------------------------
  
  n <- nrow(d)
  
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
    c(sd_e1^2, 0,
      0, sd_e2^2),
    nrow = 2,
    byrow = TRUE
  )
  
  # Observation matrix
  # y1_t = rho * (tau_t + c_t) + e1
  # y2_t =       (tau_t + c_t) + e2
  
  H <- matrix(
    c(rho, rho,
      1,   1),
    nrow = 2,
    byrow = TRUE
  )
  
  # Initial state
  a <- matrix(c(mean(y2, na.rm = TRUE), 0), nrow = 2)
  P <- diag(c(100, 100))
  
  loglik <- 0
  
  for (i in seq_len(n)) {
    
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    
    v <- y - H %*% a
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
  omega_start <- as.numeric(unique(d$omega)[1])
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma_start <- min(max(sigma_start, sigma_min + 1e-4), sigma_max - 1e-4)
  
  sigma_raw_start <- logit(
    (sigma_start - sigma_min) / (sigma_max - sigma_min)
  )
  
  omega_width <- 0.15
  omega_lower <- max(0.001, omega_start - omega_width)
  omega_upper <- min(0.999, omega_start + omega_width)
  
  omega_scaled_start <- (omega_start - omega_lower) / (omega_upper - omega_lower)
  omega_scaled_start <- min(max(omega_scaled_start, 0.001), 0.999)
  
  omega_raw_start <- logit(omega_scaled_start)
  
  theta_start <- c(
    sigma_raw = sigma_raw_start,
    omega_raw = omega_raw_start,
    rho_c_raw = 0,
    log_sd_eta = log(0.05),
    log_sd_nu  = log(0.05),
    log_sd_e1  = log(0.05),
    log_sd_e2  = log(0.05)
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
  
  sigma_min <- 0.10
  sigma_max <- 2.50
  
  sigma <- sigma_min + (sigma_max - sigma_min) * inv_logit(theta["sigma_raw"])
  omega_center <- as.numeric(unique(d$omega)[1])
  omega <- omega_from_raw(
    omega_raw = theta["omega_raw"],
    omega_center = omega_center,
    width = 0.15
  )
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e1  <- sd_min + exp(theta["log_sd_e1"])
  sd_e2  <- sd_min + exp(theta["log_sd_e2"])
  
  summary <- tibble(
    nest = nest_name,
    sigma = sigma,
    omega = omega,
    rho_c = rho_c,
    var_eta = sd_eta^2,
    var_nu = sd_nu^2,
    var_e1 = sd_e1^2,
    var_e2 = sd_e2^2,
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
  omega_center <- as.numeric(unique(d$omega)[1])
  omega <- omega_from_raw(
    omega_raw = theta["omega_raw"],
    omega_center = omega_center,
    width = 0.15
  )
  rho_c <- 0.98 * tanh(theta["rho_c_raw"])
  
  sd_min <- 1e-4
  
  sd_eta <- sd_min + exp(theta["log_sd_eta"])
  sd_nu  <- sd_min + exp(theta["log_sd_nu"])
  sd_e1  <- sd_min + exp(theta["log_sd_e1"])
  sd_e2  <- sd_min + exp(theta["log_sd_e2"])
  
  rho <- (sigma - 1) / sigma
  
  C  <- d$C_obs
  C1 <- d$C1
  C2 <- d$C2
  P1 <- d$P1
  P2 <- d$P2
  
  y1 <- log((P1 * C1) / (P2 * C2)) -
    log(omega / (1 - omega)) -
    rho * log(C1 / C2)
  
  inside <- (C^rho - (1 - omega) * C2^rho) / (omega * C1^rho)
  y2 <- log(inside) / rho
  
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
    c(sd_e1^2, 0,
      0, sd_e2^2),
    nrow = 2,
    byrow = TRUE
  )
  
  H <- matrix(
    c(rho, rho,
      1,   1),
    nrow = 2,
    byrow = TRUE
  )
  
  n <- nrow(d)
  
  a <- matrix(c(mean(y2, na.rm = TRUE), 0), nrow = 2)
  P <- diag(c(100, 100))
  
  tau_vec <- numeric(n)
  c_vec <- numeric(n)
  mu_vec <- numeric(n)
  
  foc_obs_vec <- numeric(n)
  foc_model_vec <- numeric(n)
  ces_obs_vec <- numeric(n)
  ces_model_vec <- numeric(n)
  
  for (i in seq_len(n)) {
    
    y <- matrix(c(y1[i], y2[i]), nrow = 2)
    
    v <- y - H %*% a
    Fmat <- H %*% P %*% t(H) + R
    Finv <- safe_solve(Fmat)
    
    K <- P %*% t(H) %*% Finv
    a <- a + K %*% v
    P <- P - K %*% H %*% P
    
    tau <- as.numeric(a[1])
    c_comp <- as.numeric(a[2])
    mu <- tau + c_comp
    
    tau_vec[i] <- tau
    c_vec[i] <- c_comp
    mu_vec[i] <- mu
    
    foc_obs <- log((P1[i] * C1[i]) / (P2[i] * C2[i]))
    
    foc_model <- log(omega / (1 - omega)) +
      rho * (mu + log(C1[i] / C2[i]))
    
    C_model <- (
      omega * (exp(mu) * C1[i])^rho +
        (1 - omega) * C2[i]^rho
    )^(1 / rho)
    
    foc_obs_vec[i] <- foc_obs
    foc_model_vec[i] <- foc_model
    ces_obs_vec[i] <- log(C[i])
    ces_model_vec[i] <- log(C_model)
    
    a <- Tmat %*% a
    P <- Tmat %*% P %*% t(Tmat) + Q
  }
  
  residuals <- tibble(
    nest = nest_name,
    t = d$t,
    sigma = sigma,
    omega = omega,
    rho_c = rho_c,
    tau = tau_vec,
    c = c_vec,
    mu = mu_vec,
    foc_obs = foc_obs_vec,
    foc_model = foc_model_vec,
    resid_foc = foc_obs_vec - foc_model_vec,
    ces_obs = ces_obs_vec,
    ces_model = ces_model_vec,
    resid_ces = ces_obs_vec - ces_model_vec
  )
  
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
    "Data/results_mu_kalman_free_omega_bounded_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

write_csv(
  residuals_kalman,
  paste0(
    "Data/residuals_mu_kalman_free_omega_bounded_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

write_csv(
  states_kalman,
  paste0(
    "Data/states_mu_kalman_free_omega_bounded_",
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
