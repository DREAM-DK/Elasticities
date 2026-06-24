# ============================================================
# 01_make_omega_and_nests.R
#
# Purpose:
# Construct Paasche-based CES nest aggregates and calibrate fixed
# CES share parameters, omega, for each consumption nest.
#
# Main task:
#   - Load raw consumption data
#   - Create Paasche price and quantity aggregates for nested goods
#   - Define the CES nest structure
#   - Calibrate omega using old sigma estimates from the paper
#   - Export nest data and calibrated omega values
#
# Main outputs:
#   Data/omega_nests_1983_2017_paasche.csv
#   Data/forbrug_nests_1983_2017_paasche.csv
# ============================================================


rm(list=ls())

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

komponenter <- c("ctur", "ctje", "cvar", "cene", "cbil", "cbol")

setwd("P:/leawin/Forbrugselasticitet")

df <- read_csv("Data/forbrug_final.csv")

df <- df %>%
  mutate(
    t = as.integer(t),
    qC = as.numeric(qC),
    pC = as.numeric(pC),
    value = pC * qC
  )

# ------------------------------------------------------------
# 2. Vælg estimationsperiode
# ------------------------------------------------------------

start_year <- 1983
end_year   <- 2017

df_old <- df %>%
  filter(t >= start_year, t <= end_year)

# ------------------------------------------------------------
# 3. Funktion til at calibrere omega fra FOC
# ------------------------------------------------------------

calibrate_omega <- function(data, good1, good2, sigma, mu = 0) {
  
  nest_data <- data %>%
    filter(uni %in% c(good1, good2)) %>%
    select(uni, t, qC, pC, value) %>%
    tidyr::pivot_wider(
      names_from = uni,
      values_from = c(qC, pC, value)
    )
  
  # Gennemsnitsdata
  C1 <- mean(nest_data[[paste0("qC_", good1)]], na.rm = TRUE)
  C2 <- mean(nest_data[[paste0("qC_", good2)]], na.rm = TRUE)
  
  P1 <- mean(nest_data[[paste0("pC_", good1)]], na.rm = TRUE)
  P2 <- mean(nest_data[[paste0("pC_", good2)]], na.rm = TRUE)
  
  # FOC-komponenter
  R <- (P1 * C1) / (P2 * C2)
  
  Z <- ((exp(mu) * C1) / C2)^((sigma - 1) / sigma)
  
  omega <- R / (R + Z)
  
  return(
    tibble(
      good1 = good1,
      good2 = good2,
      sigma = sigma,
      mu = mu,
      R = R,
      Z = Z,
      omega = omega
    )
  )
}


# ------------------------------------------------------------
# 4. Lav aggregat-funktionen
# ------------------------------------------------------------
make_paasche_aggregate <- function(data, components, new_name) {
  
  tmp <- data %>%
    filter(uni %in% components) %>%
    mutate(value = pC * qC)
  
  years <- sort(unique(tmp$t))
  
  result <- vector("list", length(years))
  p_index_prev <- NA_real_
  
  for (i in seq_along(years)) {
    
    year_i <- years[i]
    
    data_i <- tmp %>%
      filter(t == year_i) %>%
      select(uni, qC, pC, value)
    
    value_i <- sum(data_i$value, na.rm = TRUE)
    
    if (i == 1) {
      p_index <- 1
    } else {
      
      year_prev <- years[i - 1]
      
      data_prev <- tmp %>%
        filter(t == year_prev) %>%
        select(uni, pC_prev = pC)
      
      joined <- data_i %>%
        left_join(data_prev, by = "uni")
      
      paasche_growth <- sum(joined$pC * joined$qC, na.rm = TRUE) /
        sum(joined$pC_prev * joined$qC, na.rm = TRUE)
      
      p_index <- p_index_prev * paasche_growth
    }
    
    q_index <- value_i / p_index
    
    result[[i]] <- tibble(
      uni = new_name,
      t = year_i,
      qC = q_index,
      pC = p_index,
      value = value_i
    )
    
    p_index_prev <- p_index
  }
  
  bind_rows(result)
}
# ------------------------------------------------------------
# 5. Lav de sammensatte forbrugsaggregater
# ------------------------------------------------------------
df_old2 <- df_old %>%
  mutate(value = pC * qC)

agg_cturctje <- make_paasche_aggregate(
  data = df_old2,
  components = c("ctur", "ctje"),
  new_name = "cturctje"
)

agg_cturctjecvar <- make_paasche_aggregate(
  data = df_old2,
  components = c("ctur", "ctje", "cvar"),
  new_name = "cturctjecvar"
)

agg_cturctjecvarcene <- make_paasche_aggregate(
  data = df_old2,
  components = c("ctur", "ctje", "cvar", "cene"),
  new_name = "cturctjecvarcene"
)

agg_cikkebol <- make_paasche_aggregate(
  data = df_old2,
  components = c("ctur", "ctje", "cvar", "cene", "cbil"),
  new_name = "cikkebol"
)

# ------------------------------------------------------------
# 6. Saml data
# ------------------------------------------------------------

df_nests <- bind_rows(
  df_old2,
  agg_cturctje,
  agg_cturctjecvar,
  agg_cturctjecvarcene,
  agg_cikkebol
)

# ------------------------------------------------------------
# 7. Definér nests med gamle sigmaer
# ------------------------------------------------------------

nests <- tibble(
  nest  = c("cturctje", "cturctjecvar", "cturctjecvarcene", "cikkebol"),
  good1 = c("ctur", "cturctje", "cturctjecvar", "cturctjecvarcene"),
  good2 = c("ctje", "cvar", "cene", "cbil"),
  sigma = c(1.25, 0.94, 0.26, 1.04),
  mu    = c(0, 0, 0, 0)
)
# ------------------------------------------------------------
# 8. Beregn omega pr. nest
# ------------------------------------------------------------

omega_nests <- bind_rows(
  lapply(
    seq_len(nrow(nests)),
    function(i) {
      calibrate_omega(
        data = df_nests,
        good1 = nests$good1[i],
        good2 = nests$good2[i],
        sigma = nests$sigma[i],
        mu = nests$mu[i]
      ) %>%
        mutate(nest = nests$nest[i], .before = 1)
    }
  )
)

omega_nests

# ------------------------------------------------------------
# 9. Gemmer Omegas
# ------------------------------------------------------------

# omega
omega_file <- paste0(
  "Data/omega_nests_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(omega_nests, omega_file)

# Paasche-aggregater over tid 
nests_data_file <- paste0(
  "Data/forbrug_nests_",
  start_year, "_", end_year,
  "_paasche.csv"
)

write_csv(df_nests, nests_data_file)
