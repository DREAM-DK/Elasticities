# ============================================================
# 02_make_regression_data.R
#
# Purpose:
# Build the regression dataset used for CES system estimation.
#
# Main task:
#   - Load Paasche nest data and calibrated omega values
#   - Match each nest with its two input components
#   - Add observed aggregate consumption for each nest
#   - Construct relative quantity, price, and value variables
#   - Export one estimation-ready dataset
#
# Main output:
#   Data/regression_data_1983_2017_paasche.csv
# ============================================================


rm(list=ls())

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

setwd("P:/leawin/Forbrugselasticitet")


# ------------------------------------------------------------
# 1. Vælg samme periode som omega-filen
# ------------------------------------------------------------

start_year <- 1983
end_year   <- 2017

# ------------------------------------------------------------
# 2.+3. Load data + oprydning af datatyper
# ------------------------------------------------------------

omega_file <- paste0(
  "Data/omega_nests_",
  start_year, "_", end_year,
  "_paasche.csv"
)

nests_data_file <- paste0(
  "Data/forbrug_nests_",
  start_year, "_", end_year,
  "_paasche.csv"
)

omega_nests <- read_csv(omega_file) %>%
  mutate(
    nest = tolower(nest),
    good1 = tolower(good1),
    good2 = tolower(good2),
    sigma_old = as.numeric(sigma),
    omega = as.numeric(omega)
  )

df_nests <- read_csv(nests_data_file) %>%
  mutate(
    uni = tolower(uni),
    t = as.integer(t),
    qC = as.numeric(qC),
    pC = as.numeric(pC),
    value = as.numeric(value)
  )
# ------------------------------------------------------------
# 4. Regresionsdatasættet
# ------------------------------------------------------------

make_regression_data <- function(data, nest_name, good1_name, good2_name) {
  
  inputs_wide <- data %>%
    filter(uni %in% c(good1_name, good2_name)) %>%
    select(uni, t, qC, pC, value) %>%
    pivot_wider(
      names_from = uni,
      values_from = c(qC, pC, value)
    )
  
  nest_obs <- data %>%
    filter(uni == nest_name) %>%
    select(
      t,
      C_obs = qC,
      P_obs = pC,
      V_obs = value
    )
  
  inputs_wide %>%
    left_join(nest_obs, by = "t") %>%
    transmute(
      nest = nest_name,
      t = t,
      
      C_obs = C_obs,
      P_obs = P_obs,
      V_obs = V_obs,
      
      C1 = .data[[paste0("qC_", good1_name)]],
      C2 = .data[[paste0("qC_", good2_name)]],
      P1 = .data[[paste0("pC_", good1_name)]],
      P2 = .data[[paste0("pC_", good2_name)]],
      
      V1 = P1 * C1,
      V2 = P2 * C2,
      
      log_rel_q = log(C1 / C2),
      log_rel_p = log(P1 / P2),
      log_rel_value = log((P1 * C1) / (P2 * C2)),
      log_C_obs = log(C_obs)
    )
}

reg_data <- bind_rows(
  lapply(
    seq_len(nrow(omega_nests)),
    function(i) {
      make_regression_data(
        data = df_nests,
        nest_name = omega_nests$nest[i],
        good1_name = omega_nests$good1[i],
        good2_name = omega_nests$good2[i]
      )
    }
  )
) %>%
  left_join(
    omega_nests %>%
      select(nest, good1, good2, sigma_old, omega),
    by = "nest"
  )

reg_data

# ------------------------------------------------------------
# 7. Gem regressionsdatasæt
# ------------------------------------------------------------

reg_file <- paste0(
  "Data/regression_data_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(reg_data, reg_file)
