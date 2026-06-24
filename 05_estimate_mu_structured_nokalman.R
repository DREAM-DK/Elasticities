# ============================================================
# 05_estimate_mu_structured_nokalman.R
#
# Purpose:
# Estimate a CES system with a structured preference process for mu_t
# without using a Kalman filter.
#
# Model:
#   mu_t   = tau_t + c_t
#   tau_t  = tau_{t-1} + eta_t
#   c_t    = rho * c_{t-1} + v_t
#
# The model is estimated with nls.lm by combining:
#   1. relative budget-share FOC residuals
#   2. CES aggregator residuals
#   3. smoothness/state residuals for tau_t
#   4. AR(1) residuals for c_t
#
# Main task:
#   - Load regression data
#   - Estimate sigma jointly with structured preference components
#   - Keep omega fixed at calibrated Paasche values
#   - Compare fit against the baseline and flexible-mu models
#   - Plot tau_t, c_t, mu_t, and residuals
#
# Main outputs:
#   Data/results_structured_mu_nokalman_1983_2017_paasche.csv
#   Data/residuals_structured_mu_nokalman_1983_2017_paasche.csv
#   Data/states_structured_mu_nokalman_1983_2017_paasche.csv
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
# 3. Estimation med struktureret mu_t uden Kalman
# ------------------------------------------------------------

estimate_structured_mu_one_nest <- function(data, nest_name,
                                            lambda_tau = 0.10,
                                            lambda_c = 0.10) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      C = C_obs
    )
  
  omega_fixed <- as.numeric(unique(d$omega)[1])
  sigma_start <- as.numeric(unique(d$sigma_old)[1])
  
  n_years <- nrow(d)
  
  # Startværdier
  par_start <- c(
    sig = sigma_start,
    rho_c_raw = 0,
    rep(0, n_years),  # tau_t
    rep(0, n_years)   # c_t
  )
  
  names(par_start)[3:(2 + n_years)] <- paste0("tau_", d$t)
  names(par_start)[(3 + n_years):(2 + 2 * n_years)] <- paste0("c_", d$t)
  
  lower <- c(
    sig = 0.05,
    rho_c_raw = -Inf,
    rep(-Inf, n_years),
    rep(-Inf, n_years)
  )
  
  upper <- c(
    sig = 5,
    rho_c_raw = Inf,
    rep(Inf, n_years),
    rep(Inf, n_years)
  )
  
  residual_fun <- function(par) {
    
    sig <- as.numeric(par["sig"])
    omg <- omega_fixed
    
    # Transformerer rho_c_raw, så rho_c altid ligger mellem -0.98 og 0.98
    rho_c <- 0.98 * tanh(as.numeric(par["rho_c_raw"]))
    
    rho <- (sig - 1) / sig
    
    tau <- as.numeric(par[paste0("tau_", d$t)])
    c_comp <- as.numeric(par[paste0("c_", d$t)])
    mu <- tau + c_comp
    
    # --------------------------------------------------------
    # 1. Relativ FOC / budgetandelsligning
    # --------------------------------------------------------
    
    foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
    
    foc_model <- log(omg / (1 - omg)) +
      rho * (mu + log(d$C1 / d$C2))
    
    res_foc <- foc_obs - foc_model
    
    # --------------------------------------------------------
    # 2. CES aggregator
    # --------------------------------------------------------
    
    C_model <- (
      omg * (exp(mu) * d$C1)^rho +
        (1 - omg) * d$C2^rho
    )^(1 / rho)
    
    res_ces <- log(d$C) - log(C_model)
    
    # --------------------------------------------------------
    # 3. State residual for tau_t
    #    tau_t = tau_{t-1} + eta_t
    # --------------------------------------------------------
    
    res_tau <- diff(tau)
    
    # --------------------------------------------------------
    # 4. AR(1)-residual for c_t
    #    c_t = rho_c * c_{t-1} + v_t
    # --------------------------------------------------------
    
    res_c <- c_comp[-1] - rho_c * c_comp[-n_years]
    
    # --------------------------------------------------------
    # Normalisering / identifikation
    # --------------------------------------------------------
    # Uden denne kan tau og c flytte niveau mellem sig.
    
    res_mean_c <- mean(c_comp)
    
    c(
      res_foc,
      res_ces,
      sqrt(lambda_tau) * res_tau,
      sqrt(lambda_c) * res_c,
      10 * res_mean_c
    )
  }
  
  fit <- nls.lm(
    par = par_start,
    fn = residual_fun,
    lower = lower,
    upper = upper,
    control = nls.lm.control(maxiter = 2000)
  )
  
  par <- coef(fit)
  
  sig <- as.numeric(par["sig"])
  omg <- omega_fixed
  rho_c <- 0.98 * tanh(as.numeric(par["rho_c_raw"]))
  rho <- (sig - 1) / sig
  
  tau <- as.numeric(par[paste0("tau_", d$t)])
  c_comp <- as.numeric(par[paste0("c_", d$t)])
  mu <- tau + c_comp
  
  foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
  
  foc_model <- log(omg / (1 - omg)) +
    rho * (mu + log(d$C1 / d$C2))
  
  C_model <- (
    omg * (exp(mu) * d$C1)^rho +
      (1 - omg) * d$C2^rho
  )^(1 / rho)
  
  residuals <- tibble(
    nest = nest_name,
    t = d$t,
    sigma = sig,
    omega = omg,
    rho_c = rho_c,
    tau = tau,
    c = c_comp,
    mu = mu,
    foc_obs = foc_obs,
    foc_model = foc_model,
    resid_foc = foc_obs - foc_model,
    ces_obs = log(d$C),
    ces_model = log(C_model),
    resid_ces = log(d$C) - log(C_model)
  )
  
  states <- tibble(
    nest = nest_name,
    t = d$t,
    tau = tau,
    c = c_comp,
    mu = mu
  )
  
  summary <- tibble(
    nest = nest_name,
    sigma = sig,
    omega = omg,
    rho_c = rho_c,
    ssr = sum(residual_fun(par)^2),
    converged = fit$info %in% c(1, 2, 3, 4)
  )
  
  list(
    summary = summary,
    residuals = residuals,
    states = states
  )
}

# ------------------------------------------------------------
# 4. Kør struktureret mu_t for alle nests
# ------------------------------------------------------------

structured_mu_fits <- lapply(unique(reg_data$nest), function(n) {
  estimate_structured_mu_one_nest(
    data = reg_data,
    nest_name = n,
    lambda_tau = 0.10,
    lambda_c = 0.10
  )
})

results_structured_mu <- bind_rows(
  lapply(structured_mu_fits, function(x) x$summary)
)

residuals_structured_mu <- bind_rows(
  lapply(structured_mu_fits, function(x) x$residuals)
)

states_structured_mu <- bind_rows(
  lapply(structured_mu_fits, function(x) x$states)
)

results_structured_mu
residuals_structured_mu
states_structured_mu

# ------------------------------------------------------------
# 5. Gem resultater
# ------------------------------------------------------------

results_file <- paste0(
  "Data/results_structured_mu_nokalman_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(results_structured_mu, results_file)

residuals_file <- paste0(
  "Data/residuals_structured_mu_nokalman_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(residuals_structured_mu, residuals_file)

states_file <- paste0(
  "Data/states_structured_mu_nokalman_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(states_structured_mu, states_file)


# ------------------------------------------------------------
# 6. Plot states
# ------------------------------------------------------------

ggplot(states_structured_mu, aes(x = t, y = mu)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Struktureret præferencekomponent mu_t",
    x = "År",
    y = "mu_t"
  )

ggplot(states_structured_mu, aes(x = t, y = tau)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Trendkomponent tau_t",
    x = "År",
    y = "tau_t"
  )

ggplot(states_structured_mu, aes(x = t, y = c)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Cyklisk komponent c_t",
    x = "År",
    y = "c_t"
  )

# ------------------------------------------------------------
# 7. Residualplots
# ------------------------------------------------------------

ggplot(residuals_structured_mu, aes(x = t, y = resid_foc)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra relativ FOC, struktureret mu_t",
    x = "År",
    y = "FOC residual"
  )

ggplot(residuals_structured_mu, aes(x = t, y = resid_ces)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra CES-aggregator, struktureret mu_t",
    x = "År",
    y = "CES residual"
  )