# ============================================================
# 04_estimate_mu_flexible.R
#
# Purpose:
# Test a more flexible preference component, mu_t, before imposing
# a structured state process.
#
# Model:
#   mu_t is estimated freely year by year
#
# The model is estimated with nls.lm by combining:
#   1. relative budget-share FOC residuals
#   2. CES aggregator residuals
#
# Main task:
#   - Load regression data
#   - Estimate sigma and one free mu_t for each year and nest
#   - Inspect whether a flexible preference component improves fit
#   - Plot estimated mu_t and residuals
#
# Main outputs:
#   Data/results_flexible_mu_1983_2017_paasche.csv
#   Data/residuals_flexible_mu_1983_2017_paasche.csv
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
# 3. Estimation med fleksibel mu_t
# ------------------------------------------------------------

estimate_flexible_mu_one_nest <- function(data, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      C = C_obs
    )
  
  omega_fixed <- as.numeric(unique(d$omega)[1])
  sigma_start <- as.numeric(unique(d$sigma_old)[1])
  
  n_years <- nrow(d)
  
  par_start <- c(
    sig = sigma_start,
    rep(0, n_years)
  )
  
  names(par_start)[-1] <- paste0("mu_", d$t)
  
  lower <- c(
    sig = 0.05,
    rep(-Inf, n_years)
  )
  
  upper <- c(
    sig = 5,
    rep(Inf, n_years)
  )
  
  residual_fun <- function(par) {
    
    sig <- as.numeric(par["sig"])
    omg <- omega_fixed
    rho <- (sig - 1) / sig
    
    mu <- as.numeric(par[paste0("mu_", d$t)])
    
    # Relativ FOC / budgetandelsligning
    foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
    
    foc_model <- log(omg / (1 - omg)) +
      rho * (mu + log(d$C1 / d$C2))
    
    res_foc <- foc_obs - foc_model
    
    # CES aggregator
    C_model <- (
      omg * (exp(mu) * d$C1)^rho +
        (1 - omg) * d$C2^rho
    )^(1 / rho)
    
    res_ces <- log(d$C) - log(C_model)
    
    c(res_foc, res_ces)
  }
  
  fit <- nls.lm(
    par = par_start,
    fn = residual_fun,
    lower = lower,
    upper = upper,
    control = nls.lm.control(maxiter = 1000)
  )
  
  par <- coef(fit)
  
  sig <- as.numeric(par["sig"])
  omg <- omega_fixed
  rho <- (sig - 1) / sig
  
  mu <- as.numeric(par[paste0("mu_", d$t)])
  
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
    mu = mu,
    foc_obs = foc_obs,
    foc_model = foc_model,
    resid_foc = foc_obs - foc_model,
    ces_obs = log(d$C),
    ces_model = log(C_model),
    resid_ces = log(d$C) - log(C_model)
  )
  
  summary <- tibble(
    nest = nest_name,
    sigma = sig,
    omega = omg,
    ssr = sum(residual_fun(par)^2),
    converged = fit$info %in% c(1, 2, 3, 4)
  )
  
  list(
    summary = summary,
    residuals = residuals
  )
}

# ------------------------------------------------------------
# 4. Kør fleksibel mu_t for alle nests
# ------------------------------------------------------------

flexible_mu_fits <- lapply(unique(reg_data$nest), function(n) {
  estimate_flexible_mu_one_nest(reg_data, n)
})

results_flexible_mu <- bind_rows(
  lapply(flexible_mu_fits, function(x) x$summary)
)

residuals_flexible_mu <- bind_rows(
  lapply(flexible_mu_fits, function(x) x$residuals)
)

results_flexible_mu
residuals_flexible_mu

# ------------------------------------------------------------
# 5. Gem resultater
# ------------------------------------------------------------

results_file <- paste0(
  "Data/results_flexible_mu_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(results_flexible_mu, results_file)

residuals_file <- paste0(
  "Data/residuals_flexible_mu_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(residuals_flexible_mu, residuals_file)

# ------------------------------------------------------------
# 6. Plot mu_t
# ------------------------------------------------------------

ggplot(residuals_flexible_mu, aes(x = t, y = mu)) +
  geom_line() +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Estimeret fleksibel præferencekomponent mu_t",
    x = "År",
    y = "mu_t"
  )

ggplot(residuals_flexible_mu, aes(x = t, y = resid_foc)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra relativ FOC, fleksibel mu_t",
    x = "År",
    y = "FOC residual"
  )

ggplot(residuals_flexible_mu, aes(x = t, y = resid_ces)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra CES-aggregator, fleksibel mu_t",
    x = "År",
    y = "CES residual"
  )
