# ============================================================
# 00_data_fra_makro.R
#
# Purpose:
# Extract consumption quantities and prices from MAKRO/GDX data and
# create the raw consumption dataset used in the estimation project.
#
# Main task:
#   - Read qC and pC from MAKRO data
#   - Merge quantities and prices by consumption component and year
#   - Keep relevant consumption components
#   - Export cleaned data to CSV
#
# Main output:
#   Data/forbrug_final.csv
# ============================================================


rm(list = ls())

library(gamstransfer)
library(dplyr)

# Sæt working directory til mappen
setwd("P:/leawin/Forbrugselasticitet/Data")

komponenter <- c("ctur", "ctje", "cvar", "cene", "cbil", "cbol")

# Læs GDX-filen
m <- Container$new("makrobk.gdx")

# Læser data
qC <- m$data$get("qC")$records
pC <- m$data$get("pC")$records

pC_clean <- pC %>%
  select(uni, t, level) %>%
  rename(pC = level)

qC_clean <- qC %>%
  select(uni, t, level) %>%
  rename(qC = level)

# Definerer forbrug
forbrug <- left_join(qC_clean, pC_clean, by = c("uni", "t"))

forbrug_final <- forbrug[forbrug$uni %in% komponenter, ]

# Udskriver til CSV
write.csv(
  forbrug_final,
  "P:/leawin/Forbrugselasticitet/Data/forbrug_final.csv",
  row.names = FALSE
)
