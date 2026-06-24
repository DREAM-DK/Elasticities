# ============================================================
# 07_compare_all_sigma_results.R
#
# Purpose:
# Compare estimated CES substitution elasticities across all model
# versions.
#
# Models compared:
#   - old paper estimates
#   - 03 baseline fixed omega
#   - 05 structured mu without Kalman
#   - 06a Kalman implied-mu, fixed omega
#   - 06b Kalman implied-mu, bounded free omega
#   - 06c Kalman direct system, fixed omega, free measurement variances
#   - 06d Kalman direct system, fixed omega, equal measurement variance
#
# Main outputs:
#   Data/comparison_all_sigma_results_1983_2017_paasche.csv
#   Tables/comparison_all_sigma_results_1983_2017_paasche.pdf
# ============================================================

rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(grid)

setwd("P:/leawin/Forbrugselasticitet")

start_year <- as.integer(Sys.getenv("START_YEAR", unset = "1983"))
end_year <- as.integer(Sys.getenv("END_YEAR", unset = "2017"))

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

reg_file <- paste0(
  "Data/regression_data_",
  start_year, "_", end_year,
  "_paasche.csv"
)

reg_data <- read_csv(reg_file) %>%
  mutate(
    nest = tolower(nest),
    sigma_old = as.numeric(sigma_old)
  )

old_sigma <- reg_data %>%
  distinct(nest, sigma_old)


safe_read <- function(file) {
  if (file.exists(file)) {
    read_csv(file, show_col_types = FALSE)
  } else {
    message("File not found: ", file)
    NULL
  }
}

safe_read_first_existing <- function(files, label) {
  existing_file <- files[file.exists(files)][1]
  if (is.na(existing_file)) {
    message("No file found for ", label, ". Tried: ", paste(files, collapse = ", "))
    return(NULL)
  }
  message("Using ", label, " file: ", existing_file)
  read_csv(existing_file, show_col_types = FALSE)
}
# ------------------------------------------------------------
# 2. Load model results
# ------------------------------------------------------------

baseline_files <- c(
  paste0(
    "Data/results_fixed_omega_",
    start_year, "_", end_year,
    "_paasche.csv"
  ),
  paste0(
    "Data/results_fixed_omega_normalized_v2_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

structured_file <- paste0(
  "Data/results_structured_mu_nokalman_",
  start_year, "_", end_year,
  "_paasche.csv"
)

kalman_6a_file <- paste0(
  "Data/results_mu_kalman_fixed_omega_",
  start_year, "_", end_year,
  "_paasche.csv"
)

kalman_6b_file <- paste0(
  "Data/results_mu_kalman_free_omega_bounded_",
  start_year, "_", end_year,
  "_paasche.csv"
)

kalman_6c_file <- paste0(
  "Data/results_mu_kalman_direct_fixed_omega_",
  start_year, "_", end_year,
  "_paasche.csv"
)

kalman_6d_file <- paste0(
  "Data/results_mu_kalman_direct_equal_measvar_",
  start_year, "_", end_year,
  "_paasche.csv"
)


results_baseline <- safe_read_first_existing(baseline_files, "baseline")
results_structured <- safe_read(structured_file)
results_6a <- safe_read(kalman_6a_file)
results_6b <- safe_read(kalman_6b_file)
results_6c <- safe_read(kalman_6c_file)
results_6d <- safe_read(kalman_6d_file)

# ------------------------------------------------------------
# 3. Build comparison table
# ------------------------------------------------------------

comparison <- old_sigma

if (!is.null(results_baseline)) {
  comparison <- comparison %>%
    left_join(
      results_baseline %>%
        select(nest, sigma_baseline = sigma),
      by = "nest"
    )
}

if (!is.null(results_structured)) {
  comparison <- comparison %>%
    left_join(
      results_structured %>%
        select(nest, sigma_structured_nokalman = sigma),
      by = "nest"
    )
}

if (!is.null(results_6a)) {
  comparison <- comparison %>%
    left_join(
      results_6a %>%
        select(nest, sigma_6a = sigma),
      by = "nest"
    )
}

if (!is.null(results_6b)) {
  comparison <- comparison %>%
    left_join(
      results_6b %>%
        select(nest, sigma_6b = sigma),
      by = "nest"
    )
}

if (!is.null(results_6c)) {
  comparison <- comparison %>%
    left_join(
      results_6c %>%
        select(nest, sigma_6c = sigma),
      by = "nest"
    )
}

if (!is.null(results_6d)) {
  comparison <- comparison %>%
    left_join(
      results_6d %>%
        select(nest, sigma_6d = sigma),
      by = "nest"
    )
}

expected_sigma_cols <- c(
  "sigma_baseline",
  "sigma_structured_nokalman",
  "sigma_6a",
  "sigma_6b",
  "sigma_6c",
  "sigma_6d"
)

for (col_name in expected_sigma_cols) {
  if (!col_name %in% names(comparison)) {
    comparison[[col_name]] <- NA_real_
  }
}

comparison <- comparison %>%
  mutate(
    delta_baseline_vs_old = sigma_baseline - sigma_old,
    delta_6d_vs_old = sigma_6d - sigma_old,
    delta_6d_vs_baseline = sigma_6d - sigma_baseline
  )

comparison

# ------------------------------------------------------------
# 4. Save CSV
# ------------------------------------------------------------

write_csv(
  comparison,
  paste0(
    "Data/comparison_all_sigma_results_",
    start_year, "_", end_year,
    "_paasche.csv"
  )
)

# ------------------------------------------------------------
# 5. Prepare PDF table
# ------------------------------------------------------------

tables_dir <- "P:/leawin/Forbrugselasticitet/Tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

table_pdf <- file.path(
  tables_dir,
  paste0(
    "comparison_all_sigma_results_",
    start_year, "_", end_year,
    "_paasche.pdf"
  )
)

comparison_pdf <- comparison %>%
  mutate(
    nest = recode(
      nest,
      "cturctje" = "Tur vs Tje",
      "cturctjecvar" = "TurTje vs Var",
      "cturctjecvarcene" = "TurTjeVar vs Ene",
      "cikkebol" = "IkkeBol vs Bil"
    ),
    sigma_old = sprintf("%.3f", sigma_old),
    sigma_baseline = sprintf("%.3f", sigma_baseline),
    sigma_structured_nokalman = sprintf("%.3f", sigma_structured_nokalman),
    sigma_6a = sprintf("%.3f", sigma_6a),
    sigma_6b = sprintf("%.3f", sigma_6b),
    sigma_6c = sprintf("%.3f", sigma_6c),
    sigma_6d = sprintf("%.3f", sigma_6d)
  ) %>%
  select(
    nest,
    sigma_old,
    sigma_baseline,
    sigma_structured_nokalman,
    sigma_6a,
    sigma_6b,
    sigma_6c,
    sigma_6d
  ) %>%
  rename(
    "Nest" = nest,
    "Old" = sigma_old,
    "03 baseline" = sigma_baseline,
    "05 no Kalman" = sigma_structured_nokalman,
    "06a" = sigma_6a,
    "06b" = sigma_6b,
    "06c" = sigma_6c,
    "06d" = sigma_6d
  )

# ------------------------------------------------------------
# 6. Draw PDF table
# ------------------------------------------------------------

draw_table <- function(df, title, note) {
  
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  
  x_left <- 0.05
  x_right <- 0.95
  table_width <- x_right - x_left
  
  y_top <- 0.82
  row_h <- 0.075
  
  col_widths_raw <- c(
    0.22,  # Nest
    0.10,  # Old
    0.13,  # 03 baseline
    0.13,  # 05 no Kalman
    0.09,  # 06a
    0.09,  # 06b
    0.09,  # 06c
    0.09   # 06d
  )
  
  col_widths <- col_widths_raw / sum(col_widths_raw) * table_width
  
  col_x_left <- c(x_left, x_left + cumsum(col_widths)[-length(col_widths)])
  col_x_mid <- col_x_left + col_widths / 2
  
  grid.text(
    title,
    x = 0.5,
    y = 0.94,
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "serif")
  )
  
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_top, y_top),
    gp = gpar(lwd = 1.2)
  )
  
  header_y <- y_top - row_h / 2
  
  for (j in seq_len(n_cols)) {
    grid.text(
      names(df)[j],
      x = if (j == 1) col_x_left[j] else col_x_mid[j],
      y = header_y,
      just = if (j == 1) "left" else "center",
      gp = gpar(fontsize = 8.5, fontface = "bold", fontfamily = "serif")
    )
  }
  
  y_midrule <- y_top - row_h
  
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_midrule, y_midrule),
    gp = gpar(lwd = 0.7)
  )
  
  for (i in seq_len(n_rows)) {
    row_y <- y_midrule - row_h * (i - 0.5)
    
    for (j in seq_len(n_cols)) {
      grid.text(
        as.character(df[i, j][[1]]),
        x = if (j == 1) col_x_left[j] else col_x_mid[j],
        y = row_y,
        just = if (j == 1) "left" else "center",
        gp = gpar(fontsize = 9, fontfamily = "serif")
      )
    }
  }
  
  y_bottom <- y_midrule - row_h * n_rows
  
  grid.lines(
    x = c(x_left, x_right),
    y = c(y_bottom, y_bottom),
    gp = gpar(lwd = 1.0)
  )
  
  grid.text(
    note,
    x = x_left,
    y = y_bottom - 0.06,
    just = "left",
    gp = gpar(fontsize = 7.8, fontfamily = "serif")
  )
}


pdf(table_pdf, width = 12.5, height = 5.8)

grid.newpage()

draw_table(
  df = comparison_pdf,
  title = paste0("Comparison of CES substitution elasticities, ", start_year, "-", end_year),
  note = "Note: 06a-06d are Kalman system variants. Boundary estimates indicate weak identification in that model version."
)

dev.off()

file.exists(table_pdf)
table_pdf

