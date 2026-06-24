library(dplyr)
library(readr)
library(knitr)
library(kableExtra)
#install.packages("kableExtra")

rm(list = ls())

setwd("P:/leawin/Forbrugselasticitet")

start_year <- as.integer(Sys.getenv("START_YEAR", unset = "1983"))
end_year <- as.integer(Sys.getenv("END_YEAR", unset = "2017"))

safe_read_first_existing <- function(files, label) {
  existing_file <- files[file.exists(files)][1]
  if (is.na(existing_file)) {
    stop("No file found for ", label, ". Tried: ", paste(files, collapse = ", "))
  }
  message("Using ", label, " file: ", existing_file)
  read_csv(existing_file, show_col_types = FALSE)
}

baseline <- safe_read_first_existing(
  c(
    paste0("Data/results_fixed_omega_", start_year, "_", end_year, "_paasche.csv"),
    paste0("Data/results_fixed_omega_normalized_v2_", start_year, "_", end_year, "_paasche.csv")
  ),
  "baseline"
)

structured <- read_csv(
  paste0("Data/results_structured_mu_nokalman_", start_year, "_", end_year, "_paasche.csv"),
  show_col_types = FALSE
)

kalman <- read_csv(
  paste0("Data/results_mu_kalman_direct_equal_measvar_", start_year, "_", end_year, "_paasche.csv"),
  show_col_types = FALSE
)

estimation_table <- baseline %>%
  select(nest, sigma_baseline = sigma) %>%
  left_join(structured %>% select(nest, sigma_structured = sigma), by = "nest") %>%
  left_join(kalman %>% select(nest, sigma_kalman = sigma), by = "nest")

write_csv(
  estimation_table,
  paste0("Data/estimation_table_", start_year, "_", end_year, "_paasche.csv")
)

print(estimation_table)
