# ============================================================
# 03_estimate_baseline_fixed_omega.R
#
# Purpose:
# Estimate a baseline CES system with fixed calibrated omega.
#
# Model:
#   mu_t = a + b * trend_t
#
# The model is estimated with nls.lm by combining:
#   1. relative budget-share FOC residuals
#   2. CES aggregator residuals
#
# Main task:
#   - Load regression data
#   - Estimate sigma for each CES nest
#   - Keep omega fixed at calibrated Paasche values
#   - Compare new sigma estimates with old paper estimates
#   - Save residual diagnostics and result tables
#
# Main outputs:
#   Data/results_fixed_omega_1983_2017_paasche.csv
#   Data/residuals_fixed_omega_1983_2017_paasche.csv
#   Tables/comparison_sigma_fixed_omega_1983_2017_paasche.pdf
#   Plots/residuals_foc_fixed_omega_1983_2017_paasche.pdf
#   Plots/residuals_ces_fixed_omega_1983_2017_paasche.pdf
# ============================================================


rm(list = ls())
#install.packages("minpack.lm")
#install.packages("webshot2")
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

start_year <- as.integer(Sys.getenv("START_YEAR", unset = "1983"))
end_year <- as.integer(Sys.getenv("END_YEAR", unset = "2017"))

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
    C1 = as.numeric(C1),
    C2 = as.numeric(C2),
    P1 = as.numeric(P1),
    P2 = as.numeric(P2),
    omega = as.numeric(omega),
    sigma_old = as.numeric(sigma_old)
  )

# ------------------------------------------------------------
# 3. Estimating
# ------------------------------------------------------------

estimate_system_one_nest <- function(data, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      trend = t - mean(t),
      C = C_obs
    )
  
  omega_start <- unique(d$omega)[1]
  sigma_start <- unique(d$sigma_old)[1]
  
  # Residualfunktion til nls.lm
  residual_fun <- function(par) {
    
    sig <- par["sig"]
    omg <- par["omg"]
    a   <- par["a"]
    b   <- par["b"]
    
    mu <- a + b * d$trend
    
    # Relative FOC
    res_foc <- log((d$P1 * d$C1) / (d$P2 * d$C2)) -
      (
        log(omg / (1 - omg)) +
          ((sig - 1) / sig) * (mu + log(d$C1 / d$C2))
      )
    
    # CES aggregator
    rho <- (sig - 1) / sig
    
    C_model <- (
      omg * (exp(mu) * d$C1)^rho +
        (1 - omg) * d$C2^rho
    )^(1 / rho)
    
    res_ces <- log(d$C) - log(C_model)
    
    c(res_foc, res_ces)
  }
  
  fit <- nls.lm(
    par = c(
      sig = sigma_start,
      omg = omega_start,
      a = 0,
      b = 0
    ),
    fn = residual_fun,
    lower = c(
      sig = 0.05,
      omg = 0.001,
      a = -Inf,
      b = -Inf
    ),
    upper = c(
      sig = 5,
      omg = 0.999,
      a = Inf,
      b = Inf
    ),
    control = nls.lm.control(maxiter = 1000)
  )
  
  par <- coef(fit)
  
  tibble(
    nest = nest_name,
    sigma = par["sig"],
    omega = par["omg"],
    a = par["a"],
    b = par["b"],
    ssr = sum(residual_fun(par)^2),
    converged = fit$info %in% c(1, 2, 3, 4)
  )
}

#########################################################
#med fixed omega
#########################################################

estimate_system_one_nest_fixed_omega <- function(data, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      trend = t - mean(t),
      C = C_obs
    )
  
  omega_fixed <- unique(d$omega)[1]
  sigma_start <- unique(d$sigma_old)[1]
  
  residual_fun <- function(par) {
    
    sig <- par["sig"]
    a   <- par["a"]
    b   <- par["b"]
    
    omg <- omega_fixed
    mu <- a + b * d$trend
    rho <- (sig - 1) / sig
    
    res_foc <- log((d$P1 * d$C1) / (d$P2 * d$C2)) -
      (
        log(omg / (1 - omg)) +
          rho * (mu + log(d$C1 / d$C2))
      )
    
    C_model <- (
      omg * (exp(mu) * d$C1)^rho +
        (1 - omg) * d$C2^rho
    )^(1 / rho)
    
    res_ces <- log(d$C) - log(C_model)
    
    c(res_foc, res_ces)
  }
  
  fit <- nls.lm(
    par = c(
      sig = sigma_start,
      a = 0,
      b = 0
    ),
    fn = residual_fun,
    lower = c(
      sig = 0.05,
      a = -Inf,
      b = -Inf
    ),
    upper = c(
      sig = 5,
      a = Inf,
      b = Inf
    ),
    control = nls.lm.control(maxiter = 1000)
  )
  
  par <- coef(fit)
  
  tibble(
    nest = nest_name,
    sigma = par["sig"],
    omega = omega_fixed,
    a = par["a"],
    b = par["b"],
    ssr = sum(residual_fun(par)^2),
    converged = fit$info %in% c(1, 2, 3, 4)
  )
}

#########################################################
#med fixed omega normaliseret
#########################################################
estimate_system_one_nest_fixed_omega_normalized_v2 <- function(data, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      trend = t - mean(t),
      C = C_obs
    )
  
  omega_fixed <- unique(d$omega)[1]
  sigma_start <- unique(d$sigma_old)[1]
  
  C_bar  <- exp(mean(log(d$C),  na.rm = TRUE))
  C1_bar <- exp(mean(log(d$C1), na.rm = TRUE))
  C2_bar <- exp(mean(log(d$C2), na.rm = TRUE))
  
  residual_fun <- function(par) {
    
    sig <- par["sig"]
    g1  <- par["g1"]
    g2  <- par["g2"]
    xi  <- par["xi"]
    
    omg <- omega_fixed
    rho <- (sig - 1) / sig
    
    alpha1 <- exp(g1 * d$trend)
    alpha2 <- exp(g2 * d$trend)
    
    # Relative FOC, normalized like in the paper
    log_rel_q_norm <- log((d$C1 / C1_bar) / (d$C2 / C2_bar))
    
    foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
    
    foc_model <- log(omg / (1 - omg)) +
      rho * ((g1 - g2) * d$trend + log_rel_q_norm)
    
    res_foc <- foc_obs - foc_model
    
    # Normalized CES aggregator
    C_model <- exp(xi) * C_bar * (
      omg * (alpha1 * d$C1 / C1_bar)^rho +
        (1 - omg) * (alpha2 * d$C2 / C2_bar)^rho
    )^(1 / rho)
    
    res_ces <- log(d$C) - log(C_model)
    
    c(res_foc, res_ces)
  }
  
  fit <- nls.lm(
    par = c(
      sig = sigma_start,
      g1 = 0,
      g2 = 0,
      xi = 0
    ),
    fn = residual_fun,
    lower = c(
      sig = 0.05,
      g1 = -Inf,
      g2 = -Inf,
      xi = -Inf
    ),
    upper = c(
      sig = 5,
      g1 = Inf,
      g2 = Inf,
      xi = Inf
    ),
    control = nls.lm.control(maxiter = 1000)
  )
  
  par <- coef(fit)
  
  tibble(
    nest = nest_name,
    sigma = par["sig"],
    omega = omega_fixed,
    gamma1 = par["g1"],
    gamma2 = par["g2"],
    xi = par["xi"],
    ssr = sum(residual_fun(par)^2),
    converged = fit$info %in% c(1, 2, 3, 4)
  )
}

# ------------------------------------------------------------
# 4. Residualer for fixed omega-estimationen 
# ------------------------------------------------------------

get_residuals_fixed_omega <- function(data, results, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      trend = t - mean(t),
      C = C_obs
    )
  
  par <- results %>%
    filter(nest == nest_name)
  
  sig <- par$sigma[1]
  omg <- par$omega[1]
  a   <- par$a[1]
  b   <- par$b[1]
  
  mu <- a + b * d$trend
  rho <- (sig - 1) / sig
  
  foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
  
  foc_model <- log(omg / (1 - omg)) +
    rho * (mu + log(d$C1 / d$C2))
  
  C_model <- (
    omg * (exp(mu) * d$C1)^rho +
      (1 - omg) * d$C2^rho
  )^(1 / rho)
  
  ces_obs <- log(d$C)
  ces_model <- log(C_model)
  
  tibble(
    nest = nest_name,
    t = d$t,
    mu = mu,
    foc_obs = foc_obs,
    foc_model = foc_model,
    resid_foc = foc_obs - foc_model,
    ces_obs = ces_obs,
    ces_model = ces_model,
    resid_ces = ces_obs - ces_model
  )
}


get_residuals_fixed_omega_normalized_v2 <- function(data, results, nest_name) {
  
  d <- data %>%
    filter(nest == nest_name) %>%
    arrange(t) %>%
    mutate(
      trend = t - mean(t),
      C = C_obs
    )
  
  par <- results %>%
    filter(nest == nest_name)
  
  sig <- par$sigma[1]
  omg <- par$omega[1]
  g1  <- par$gamma1[1]
  g2  <- par$gamma2[1]
  xi  <- par$xi[1]
  
  rho <- (sig - 1) / sig
  
  C_bar  <- exp(mean(log(d$C),  na.rm = TRUE))
  C1_bar <- exp(mean(log(d$C1), na.rm = TRUE))
  C2_bar <- exp(mean(log(d$C2), na.rm = TRUE))
  
  alpha1 <- exp(g1 * d$trend)
  alpha2 <- exp(g2 * d$trend)
  
  log_rel_q_norm <- log((d$C1 / C1_bar) / (d$C2 / C2_bar))
  
  foc_obs <- log((d$P1 * d$C1) / (d$P2 * d$C2))
  
  foc_model <- log(omg / (1 - omg)) +
    rho * ((g1 - g2) * d$trend + log_rel_q_norm)
  
  C_model <- exp(xi) * C_bar * (
    omg * (alpha1 * d$C1 / C1_bar)^rho +
      (1 - omg) * (alpha2 * d$C2 / C2_bar)^rho
  )^(1 / rho)
  
  ces_obs <- log(d$C)
  ces_model <- log(C_model)
  
  tibble(
    nest = nest_name,
    t = d$t,
    gamma1 = g1,
    gamma2 = g2,
    xi = xi,
    foc_obs = foc_obs,
    foc_model = foc_model,
    resid_foc = foc_obs - foc_model,
    ces_obs = ces_obs,
    ces_model = ces_model,
    resid_ces = ces_obs - ces_model
  )
}

# ------------------------------------------------------------
# 5. Kører alle nests og residualer
# ------------------------------------------------------------

results_system <- bind_rows(
  lapply(
    unique(reg_data$nest),
    function(n) {
      estimate_system_one_nest(reg_data, n)
    }
  )
)

results_system

# omega=fixed
results_fixed_omega <- bind_rows(
  lapply(
    unique(reg_data$nest),
    function(n) {
      estimate_system_one_nest_fixed_omega(reg_data, n)
    }
  )
)

results_fixed_omega


results_fixed_omega_normalized_v2 <- bind_rows(
  lapply(unique(reg_data$nest), function(n) {
    estimate_system_one_nest_fixed_omega_normalized_v2(reg_data, n)
  })
)

results_fixed_omega_normalized_v2

residuals_fixed_omega <- bind_rows(
  lapply(unique(reg_data$nest), function(n) {
    get_residuals_fixed_omega(
      data = reg_data,
      results = results_fixed_omega,
      nest_name = n
    )
  })
)

residuals_fixed_omega

residuals_fixed_omega_normalized_v2 <- bind_rows(
  lapply(unique(reg_data$nest), function(n) {
    get_residuals_fixed_omega_normalized_v2(
      data = reg_data,
      results = results_fixed_omega_normalized_v2,
      nest_name = n
    )
  })
)


# ------------------------------------------------------------
# 6. saving
# ------------------------------------------------------------

results_file <- paste0(
  "Data/results_fixed_omega_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(results_fixed_omega, results_file)


residuals_file <- paste0(
  "Data/residuals_fixed_omega_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(residuals_fixed_omega, residuals_file)


results_v2_file <- paste0(
  "Data/results_fixed_omega_normalized_v2_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(results_fixed_omega_normalized_v2, results_v2_file)


residuals_v2_file <- paste0(
  "Data/residuals_fixed_omega_normalized_v2_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(residuals_fixed_omega_normalized_v2, residuals_v2_file)

# ------------------------------------------------------------
# 7. Sammenligning med gammelt paper
# ------------------------------------------------------------

comparison_table <- results_fixed_omega %>%
  left_join(
    reg_data %>%
      distinct(nest, sigma_old),
    by = "nest"
  ) %>%
  mutate(
    sigma_new = sigma,
    omega_calibrated = omega
  ) %>%
  select(
    nest,
    sigma_old,
    sigma_new,
    omega_calibrated,
    a,
    b,
    ssr
  )

# ------------------------------------------------------------
# 8. Gem tabel som PDF
# ------------------------------------------------------------

tables_dir <- "P:/leawin/Forbrugselasticitet/Tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

table_pdf <- file.path(
  tables_dir,
  paste0(
    "comparison_sigma_",
    start_year, "_", end_year,
    "_paasche.pdf"
  )
)

comparison_table_pdf <- comparison_table %>%
  mutate(
    nest = recode(
      nest,
      "cturctje" = "Tur vs Tje",
      "cturctjecvar" = "TurTje vs Var",
      "cturctjecvarcene" = "TurTjeVar vs Ene",
      "cikkebol" = "IkkeBol vs Bil"
    ),
    sigma_old = sprintf("%.3f", sigma_old),
    sigma_new = sprintf("%.3f", sigma_new),
    omega_calibrated = sprintf("%.3f", omega_calibrated),
    a = sprintf("%.3f", a),
    b = sprintf("%.3f", b),
    ssr = sprintf("%.3f", ssr)
  ) %>%
  rename(
    "Nest" = nest,
    "Old sigma" = sigma_old,
    "New sigma" = sigma_new,
    "Omega" = omega_calibrated,
    "a" = a,
    "b" = b,
    "SSR" = ssr
  )
draw_table <- function(df, title, note) {
  
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  
  x_left <- 0.08
  x_right <- 0.92
  table_width <- x_right - x_left
  
  y_top <- 0.80
  row_h <- 0.075
  
  # Bredder summerer til 1 og skaleres til table_width
  col_widths_raw <- c(0.26, 0.13, 0.13, 0.12, 0.10, 0.10, 0.10)
  col_widths <- col_widths_raw / sum(col_widths_raw) * table_width
  
  col_x_left <- c(x_left, x_left + cumsum(col_widths)[-length(col_widths)])
  col_x_mid <- col_x_left + col_widths / 2
  
  grid.text(
    title,
    x = 0.5,
    y = 0.93,
    gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "serif")
  )
  
  # Top rule
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_top, y_top),
    gp = gpar(lwd = 1.2)
  )
  
  # Header
  header_y <- y_top - row_h / 2
  
  for (j in seq_len(n_cols)) {
    grid.text(
      names(df)[j],
      x = if (j == 1) col_x_left[j] else col_x_mid[j],
      y = header_y,
      just = if (j == 1) "left" else "center",
      gp = gpar(fontsize = 9.5, fontface = "bold", fontfamily = "serif")
    )
  }
  
  # Mid rule
  y_midrule <- y_top - row_h
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_midrule, y_midrule),
    gp = gpar(lwd = 0.7)
  )
  
  # Body
  for (i in seq_len(n_rows)) {
    row_y <- y_midrule - row_h * (i - 0.5)
    
    for (j in seq_len(n_cols)) {
      grid.text(
        as.character(df[i, j][[1]]),
        x = if (j == 1) col_x_left[j] else col_x_mid[j],
        y = row_y,
        just = if (j == 1) "left" else "center",
        gp = gpar(fontsize = 9.8, fontfamily = "serif")
      )
    }
  }
  
  # Bottom rule
  y_bottom <- y_midrule - row_h * n_rows
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_bottom, y_bottom),
    gp = gpar(lwd = 1.0)
  )
  
  # Note
  grid.text(
    note,
    x = x_left,
    y = y_bottom - 0.065,
    just = "left",
    gp = gpar(fontsize = 8.2, fontfamily = "serif")
  )
}

pdf(table_pdf, width = 10.5, height = 5.8)

grid.newpage()

draw_table(
  df = comparison_table_pdf,
  title = paste0("Comparison of CES substitution elasticities, ", start_year, "-", end_year),
  note = "Note: New estimates use relative FOC and CES aggregator system. Omega is fixed at calibrated Paasche value."
)

dev.off()

file.exists(table_pdf)
table_pdf

# ------------------------------------------------------------
# 9. Gem residualplots for baseline fixed omega
# ------------------------------------------------------------

# ------------------------------------------------------------
# 9. Gem residualplots for baseline fixed omega
# ------------------------------------------------------------

plots_dir <- "P:/leawin/Forbrugselasticitet/Plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# FOC residualer
plot_foc_fixed_omega <- ggplot(residuals_fixed_omega, aes(x = t, y = resid_foc)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra relativ FOC, fixed omega baseline",
    x = "År",
    y = "FOC residual"
  )

ggsave(
  filename = file.path(
    plots_dir,
    paste0("residuals_foc_fixed_omega_", start_year, "_", end_year, "_paasche.pdf")
  ),
  plot = plot_foc_fixed_omega,
  width = 8,
  height = 6
)


# CES residualer
plot_ces_fixed_omega <- ggplot(residuals_fixed_omega, aes(x = t, y = resid_ces)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ nest, scales = "free_y") +
  labs(
    title = "Residualer fra CES-aggregator, fixed omega baseline",
    x = "År",
    y = "CES residual"
  )

ggsave(
  filename = file.path(
    plots_dir,
    paste0("residuals_ces_fixed_omega_", start_year, "_", end_year, "_paasche.pdf")
  ),
  plot = plot_ces_fixed_omega,
  width = 8,
  height = 6
)
