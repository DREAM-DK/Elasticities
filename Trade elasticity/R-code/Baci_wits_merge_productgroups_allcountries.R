rm(list = ls())

library(data.table)
library(dplyr)
library(readxl)

# ---------------------------------------------
# Settings
# ---------------------------------------------
base_path <- "P:/leawin/Elasticiteter2503/Estimating LP Boehm et al"
baci_dir <- file.path(base_path, "Baci")
tariffs_dir <- file.path(base_path, "Tariffs_WITS")
pref_dir <- file.path(tariffs_dir, "Prefferential_WITS")
output_dir <- file.path(base_path, "Merge_all_countries", "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
years <- 1995:2023
exclude_usa_as_exporter <- FALSE

test_years <- trimws(Sys.getenv("TEST_YEARS", ""))
if (nzchar(test_years)) {
  years_parsed <- suppressWarnings(as.integer(unlist(strsplit(test_years, ","))))
  years_parsed <- years_parsed[is.finite(years_parsed)]
  if (length(years_parsed) > 0) {
    years <- sort(unique(years_parsed))
    message("TEST_YEARS active. Running only years: ", paste(years, collapse = ", "))
  } else {
    warning("TEST_YEARS was set but could not be parsed. Using default years 1995:2023.")
  }
}

test_pref_max_files <- suppressWarnings(as.integer(Sys.getenv("TEST_PREF_MAX_FILES", "")))
if (!is.finite(test_pref_max_files) || test_pref_max_files <= 0) {
  test_pref_max_files <- NA_integer_
}

stop_after_pref <- tolower(trimws(Sys.getenv("STOP_AFTER_PREF", ""))) %in% c("1", "true", "yes", "y")

# Only used when preferential data has reporter_code = 918 (EUN aggregate reporter),
# which must be expanded to EU member-country reporters.
eu_members_for_eun_mapping <- c(
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
  "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
  "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"
)

country_codes_file <- file.path(baci_dir, "country_codes_V202601.csv")
baci_files <- file.path(baci_dir, paste0("BACI_", years, ".csv"))

wits_h6_file <- file.path(output_dir, paste0("wits_h6_allcountries_", run_tag, ".csv"))
wits_pref_h6_file <- file.path(output_dir, paste0("wits_h6_preferential_allcountries_", run_tag, ".csv"))
baci_h6_file <- file.path(output_dir, paste0("baci_h6_allcountries_", run_tag, ".csv"))
merged_file <- file.path(output_dir, paste0("baci_wits_h6_merged_allcountries_", run_tag, ".csv"))

skip_step1_wits <- tolower(trimws(Sys.getenv("SKIP_STEP1_WITS", ""))) %in% c("1", "true", "yes", "y")
skip_step2_baci <- tolower(trimws(Sys.getenv("SKIP_STEP2_BACI", ""))) %in% c("1", "true", "yes", "y")
existing_wits_h6_file <- trimws(Sys.getenv("EXISTING_WITS_H6_FILE", ""))
existing_baci_h6_file <- trimws(Sys.getenv("EXISTING_BACI_H6_FILE", ""))

if (skip_step1_wits) {
  if (!nzchar(existing_wits_h6_file)) {
    stop("SKIP_STEP1_WITS is active, but EXISTING_WITS_H6_FILE is empty.")
  }
  if (!file.exists(existing_wits_h6_file)) {
    stop("EXISTING_WITS_H6_FILE does not exist: ", existing_wits_h6_file)
  }
  wits_h6_file <- existing_wits_h6_file
  message("SKIP_STEP1_WITS active. Using existing WITS file: ", wits_h6_file)
}

if (skip_step2_baci) {
  if (!nzchar(existing_baci_h6_file)) {
    stop("SKIP_STEP2_BACI is active, but EXISTING_BACI_H6_FILE is empty.")
  }
  if (!file.exists(existing_baci_h6_file)) {
    stop("EXISTING_BACI_H6_FILE does not exist: ", existing_baci_h6_file)
  }
  baci_h6_file <- existing_baci_h6_file
  message("SKIP_STEP2_BACI active. Using existing BACI file: ", baci_h6_file)
}

mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

to_hs6 <- function(x) {
  x_num <- suppressWarnings(as.integer(x))
  out <- sprintf("%06d", x_num)
  out[is.na(x_num)] <- NA_character_
  out
}

# ---------------------------------------------
# Shared metadata
# ---------------------------------------------
country_codes <- fread(country_codes_file, colClasses = c(country_code = "integer", country_iso3 = "character"))
country_codes <- unique(country_codes[, .(country_code, country_iso3)])
valid_codes <- unique(country_codes$country_code)
valid_iso3 <- unique(country_codes$country_iso3)
eu_member_codes_for_eun <- unique(country_codes[country_iso3 %in% eu_members_for_eun_mapping, country_code])

# ---------------------------------------------
# 1) Build WITS weighted tariffs from Tariffs_WITS (HS6)
# ---------------------------------------------
message("=== Step 1: Build WITS HS6 weighted tariffs from Tariffs_WITS ===")
if (skip_step1_wits) {
  message("Step 1 skipped.")
} else {
  wits_year_files <- file.path(tariffs_dir, paste0(years, ".csv"))
  existing_wits_year_files <- wits_year_files[file.exists(wits_year_files)]
  if (length(existing_wits_year_files) == 0) stop("No yearly files found in Tariffs_WITS.")

  first_wits_write <- TRUE
  total_wits_rows <- 0L

  for (file in existing_wits_year_files) {
    this_year <- suppressWarnings(as.integer(gsub("\\.csv$", "", basename(file), ignore.case = TRUE)))
    message("Reading WITS yearly file: ", basename(file))

    dt <- fread(
      file = file,
      select = c("Reporter", "Partner", "Product", "Tariff Year", "DutyType", "Weighted Average"),
      showProgress = TRUE
    )

    setDT(dt)
    setnames(dt, old = c("Tariff Year", "DutyType", "Weighted Average"), new = c("year", "duty_type", "weighted_average"))
    dt[, reporter_code := suppressWarnings(as.integer(Reporter))]
    dt[, partner_code := suppressWarnings(as.integer(Partner))]
    dt[, hs6 := to_hs6(Product)]

    dt <- dt[
      year == this_year &
        reporter_code %in% valid_codes &
        partner_code %in% valid_codes &
        duty_type %in% c("MFN", "AHS") &
        !is.na(hs6)
    ]

    dt[, tariff_value := mean_or_na(weighted_average), by = .(reporter_code, partner_code, year, hs6, duty_type)]
    dt <- unique(dt[, .(reporter_code, partner_code, year, hs6, duty_type, tariff_value)])
    dt[, indicator_clean := fifelse(duty_type == "MFN", "mfn_weighted_avg", "ahs_weighted_avg")]

    dt <- dcast(
      dt,
      reporter_code + partner_code + year + hs6 ~ indicator_clean,
      value.var = "tariff_value",
      fill = NA_real_
    )

    dt <- merge(dt, country_codes, by.x = "reporter_code", by.y = "country_code", all.x = TRUE)
    setnames(dt, "country_iso3", "importer")
    dt <- merge(dt, country_codes, by.x = "partner_code", by.y = "country_code", all.x = TRUE)
    setnames(dt, "country_iso3", "exporter")

    dt <- dt[importer %in% valid_iso3 & exporter %in% valid_iso3]
    if (exclude_usa_as_exporter) dt <- dt[exporter != "USA"]

    out <- dt[, .(importer, exporter, year, hs6, mfn_weighted_avg, ahs_weighted_avg)]
    n_out <- nrow(out)
    total_wits_rows <- total_wits_rows + n_out

    if (n_out > 0) {
      fwrite(out, wits_h6_file, append = !first_wits_write, col.names = first_wits_write)
      first_wits_write <- FALSE
    }
  }

  if (!file.exists(wits_h6_file)) stop("No weighted WITS HS6 rows were written.")
  message("WITS HS6 rows written: ", total_wits_rows)
  message("WITS HS6 file: ", wits_h6_file)
}

# ---------------------------------------------
# 1b) Build preferential tariffs from Prefferential_WITS (HS6)
# ---------------------------------------------
message("=== Step 1b: Build preferential tariffs (HS6) ===")

if (!dir.exists(pref_dir)) {
  stop("Prefferential_WITS folder not found: ", pref_dir)
}

mapping_candidates <- unique(c(
  list.files(
    pref_dir,
    pattern = "TRAINS.*Preference.*Ben(e)?f.*\\.(csv|xls|xlsx)$",
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    tariffs_dir,
    pattern = "TRAINS.*Preference.*Ben(e)?f.*\\.(csv|xls|xlsx)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
))
mapping_candidates <- mapping_candidates[order(
  !grepl("^trainspreferencebenficiaries\\.(xls|xlsx|csv)$", basename(mapping_candidates), ignore.case = TRUE),
  !grepl("^trainspreferencebeneficiaries\\.(xls|xlsx|csv)$", basename(mapping_candidates), ignore.case = TRUE)
)]

benef_map <- data.table(region_code = character(0), partner_code = integer(0))
if (length(mapping_candidates) > 0) {
  mapping_file <- mapping_candidates[[1]]
  message("Using TRAINS mapping file: ", basename(mapping_file))

  ext <- tolower(tools::file_ext(mapping_file))
  benef_raw <- if (ext %in% c("xls", "xlsx")) {
    as.data.table(read_excel(mapping_file, sheet = 1))
  } else {
    as.data.table(fread(mapping_file, showProgress = FALSE))
  }

  region_col <- names(benef_raw)[tolower(names(benef_raw)) %in% c("regioncode", "region_code")]
  partner_col <- names(benef_raw)[tolower(names(benef_raw)) %in% c("partner", "partner_code")]

  if (length(region_col) == 0 || length(partner_col) == 0) {
    warning("Could not identify RegionCode/Partner columns in mapping file. Group partner codes will be dropped.")
  } else {
    benef_map <- benef_raw[
      ,
      .(
        region_code = toupper(trimws(as.character(get(region_col[[1]])))),
        partner_code = suppressWarnings(as.integer(get(partner_col[[1]])))
      )
    ]
    benef_map <- unique(benef_map[!is.na(region_code) & !is.na(partner_code), .(region_code, partner_code)])
  }
} else {
  warning("No TRAINS preference beneficiary mapping file found. Non-numeric Partner region codes will be dropped.")
}

# Preferential extracts in this project are mostly named Pref_H5_* for
# historical years and Pref_H6_* for recent years. Both share the same
# schema we use below (Reporter_ISO_N, Year, ProductCode, Partner,
# SimpleAverage), so include all Pref_H* files.
pref_files <- list.files(pref_dir, pattern = "Pref_H[0-9]+_.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (length(pref_files) == 0) {
  stop("No Pref_H6 files found in Prefferential_WITS.")
}
pref_files <- sort(pref_files)
if (!is.na(test_pref_max_files)) {
  pref_files <- head(pref_files, test_pref_max_files)
  message("TEST_PREF_MAX_FILES active. Using first ", length(pref_files), " preferential files.")
}
message("Found preferential H6 files: ", length(pref_files))

pref_list <- vector("list", length(pref_files))
pref_total_rows <- 0L
pref_rows_numeric_partner <- 0L
pref_rows_group_partner <- 0L
pref_rows_group_partner_mapped <- 0L
pref_rows_group_partner_unmapped <- 0L

for (idx in seq_along(pref_files)) {
  file <- pref_files[[idx]]
  message(sprintf("Reading preferential file %d/%d: %s", idx, length(pref_files), basename(file)))

  dtp <- fread(
    file = file,
    select = c("Reporter_ISO_N", "Partner", "Year", "ProductCode", "SimpleAverage"),
    showProgress = FALSE
  )

  setDT(dtp)
  setnames(
    dtp,
    old = c("Reporter_ISO_N", "Partner", "Year", "ProductCode", "SimpleAverage"),
    new = c("reporter_code", "partner_raw", "year", "hs6_raw", "pref_simple_avg")
  )

  dtp[, reporter_code := suppressWarnings(as.integer(reporter_code))]
  dtp[, year := suppressWarnings(as.integer(year))]
  dtp[, partner_raw := toupper(trimws(as.character(partner_raw)))]
  dtp <- dtp[year %in% years & reporter_code %in% c(valid_codes, 918L)]
  pref_total_rows <- pref_total_rows + nrow(dtp)

  dt_num <- dtp[
    grepl("^[0-9]+$", partner_raw),
    .(
      reporter_code,
      partner_code = suppressWarnings(as.integer(partner_raw)),
      year,
      hs6_raw,
      pref_simple_avg
    )
  ]
  pref_rows_numeric_partner <- pref_rows_numeric_partner + nrow(dt_num)

  dt_grp <- dtp[!grepl("^[0-9]+$", partner_raw), .(reporter_code, partner_raw, year, hs6_raw, pref_simple_avg)]
  pref_rows_group_partner <- pref_rows_group_partner + nrow(dt_grp)
  if (nrow(dt_grp) > 0 && nrow(benef_map) > 0) {
    pref_rows_group_partner_mapped <- pref_rows_group_partner_mapped + nrow(dt_grp[partner_raw %chin% benef_map$region_code])
    pref_rows_group_partner_unmapped <- pref_rows_group_partner_unmapped + nrow(dt_grp[!partner_raw %chin% benef_map$region_code])
  } else if (nrow(dt_grp) > 0) {
    pref_rows_group_partner_unmapped <- pref_rows_group_partner_unmapped + nrow(dt_grp)
  }

  if (nrow(dt_grp) > 0 && nrow(benef_map) > 0) {
    dt_grp <- merge(
      dt_grp,
      benef_map,
      by.x = "partner_raw",
      by.y = "region_code",
      all = FALSE,
      allow.cartesian = TRUE
    )[
      ,
      .(reporter_code, partner_code, year, hs6_raw, pref_simple_avg)
    ]
  } else {
    dt_grp <- data.table(
      reporter_code = integer(0),
      partner_code = integer(0),
      year = integer(0),
      hs6_raw = character(0),
      pref_simple_avg = numeric(0)
    )
  }

  dt_expanded <- rbindlist(list(dt_num, dt_grp), use.names = TRUE, fill = TRUE)

  # Reporter 918 (EUN) is EU-level reporter: replicate to each EU member reporter code.
  dt_eun <- dt_expanded[reporter_code == 918L]
  if (nrow(dt_eun) > 0) {
    dt_eun <- dt_eun[
      ,
      .(reporter_code = eu_member_codes_for_eun),
      by = .(partner_code, year, hs6_raw, pref_simple_avg)
    ]
    dt_expanded <- rbindlist(
      list(dt_expanded[reporter_code != 918L], dt_eun),
      use.names = TRUE,
      fill = TRUE
    )
  }

  if (nrow(dt_expanded) == 0) {
    pref_list[[idx]] <- data.table(
      reporter_code = integer(0),
      partner_code = integer(0),
      year = integer(0),
      hs6 = character(0),
      pref_simple_avg = numeric(0)
    )
    next
  }

  dt_expanded[, hs6 := to_hs6(hs6_raw)]
  dt_expanded[, pref_simple_avg := suppressWarnings(as.numeric(pref_simple_avg))]
  dt_expanded <- dt_expanded[
    partner_code %in% valid_codes &
      !is.na(hs6)
  ]

  dt_expanded <- dt_expanded[
    ,
    .(pref_simple_avg = mean_or_na(pref_simple_avg)),
    by = .(reporter_code, partner_code, year, hs6)
  ]
  pref_list[[idx]] <- dt_expanded
}

pref_dt <- rbindlist(pref_list, use.names = TRUE, fill = TRUE)
pref_dt <- pref_dt[, .(pref_simple_avg = mean_or_na(pref_simple_avg)), by = .(reporter_code, partner_code, year, hs6)]

pref_dt <- merge(pref_dt, country_codes, by.x = "reporter_code", by.y = "country_code", all.x = TRUE)
setnames(pref_dt, "country_iso3", "importer")
pref_dt <- merge(pref_dt, country_codes, by.x = "partner_code", by.y = "country_code", all.x = TRUE)
setnames(pref_dt, "country_iso3", "exporter")

pref_dt <- pref_dt[importer %in% valid_iso3 & exporter %in% valid_iso3]
if (exclude_usa_as_exporter) pref_dt <- pref_dt[exporter != "USA"]

pref_out <- pref_dt[, .(importer, exporter, year, hs6, pref_simple_avg)]
fwrite(pref_out, wits_pref_h6_file)
message("Preferential HS6 rows written: ", nrow(pref_out))
if (nrow(pref_out) > 0) {
  pref_year_min <- pref_out[, min(year, na.rm = TRUE)]
  pref_year_max <- pref_out[, max(year, na.rm = TRUE)]
  pref_rows_9518 <- pref_out[year >= 1995L & year <= 2018L, .N]
  message("Preferential year range in output: ", pref_year_min, "-", pref_year_max)
  message("Preferential rows in paper window (1995-2018): ", pref_rows_9518)
  if (pref_rows_9518 == 0L) {
    warning(
      "No preferential observations overlap 1995-2018. ",
      "Any replication using sample years 1995-2018 will have preferential_simple_avg = NA throughout."
    )
  }
}
message(
  "Preferential diagnostics | input rows: ", pref_total_rows,
  " | numeric partner rows: ", pref_rows_numeric_partner,
  " | non-numeric partner rows: ", pref_rows_group_partner,
  " | mappable non-numeric rows: ", pref_rows_group_partner_mapped,
  " | unmapped non-numeric rows: ", pref_rows_group_partner_unmapped
)
message("Preferential HS6 file: ", wits_pref_h6_file)

if (stop_after_pref) {
  message("STOP_AFTER_PREF active. Test run finished after preferential step.")
  quit(save = "no", status = 0)
}

# ---------------------------------------------
# 2) Aggregate BACI to HS6 (all countries)
# ---------------------------------------------
message("=== Step 2: BACI aggregation (HS6, all countries) ===")
if (skip_step2_baci) {
  message("Step 2 skipped.")
} else {
  baci_agg_list <- vector("list", length(baci_files))

  for (idx in seq_along(baci_files)) {
    file <- baci_files[[idx]]
    if (!file.exists(file)) next
    message(sprintf("Reading %s (%d/%d)", basename(file), idx, length(baci_files)))

    dt <- fread(
      file = file,
      select = c("t", "i", "j", "k", "v"),
      colClasses = list(integer = c("t", "i", "j"), character = "k", numeric = "v"),
      showProgress = FALSE
    )

    dt <- dt[t %in% years]
    # BACI convention (see Baci/Readme.txt): i = EXPORTER, j = IMPORTER.
    # WITS tariffs are keyed by Reporter = importer, Partner = exporter, so we
    # must map i -> exporter and j -> importer for the trade flow to line up
    # with the duty the importer charges on goods from the exporter.
    dt <- merge(dt, country_codes, by.x = "i", by.y = "country_code", all.x = TRUE)
    setnames(dt, "country_iso3", "exporter")
    dt <- merge(dt, country_codes, by.x = "j", by.y = "country_code", all.x = TRUE)
    setnames(dt, "country_iso3", "importer")

    dt <- dt[importer %in% valid_iso3]
    dt <- dt[exporter %in% valid_iso3]
    if (exclude_usa_as_exporter) dt <- dt[exporter != "USA"]

    dt[, hs6 := to_hs6(k)]
    dt <- dt[!is.na(hs6)]

    dt_agg <- dt[, .(baci_trade_value = sum(v, na.rm = TRUE)), by = .(importer, exporter, year = t, hs6)]
    baci_agg_list[[idx]] <- dt_agg
  }

  baci_h6 <- rbindlist(baci_agg_list, use.names = TRUE, fill = TRUE)
  fwrite(baci_h6, baci_h6_file)
  message("BACI HS6 rows: ", nrow(baci_h6))
  message("BACI HS6 file: ", baci_h6_file)
}

# ---------------------------------------------
# 3) Merge BACI + WITS + preferential (HS6)
# ---------------------------------------------
message("=== Step 3: Merge BACI + WITS + preferential (HS6) ===")

wits_h6 <- fread(
  file = wits_h6_file,
  colClasses = list(character = c("importer", "exporter", "hs6"), integer = "year")
)
baci_h6 <- fread(
  file = baci_h6_file,
  colClasses = list(character = c("importer", "exporter", "hs6"), integer = "year")
)
pref_h6 <- fread(
  file = wits_pref_h6_file,
  colClasses = list(character = c("importer", "exporter", "hs6"), integer = "year")
)

if (nrow(pref_h6) > 0) {
  setDT(pref_h6)
  setnames(pref_h6, "pref_simple_avg", "preferential_simple_avg")
  setkey(pref_h6, importer, exporter, year, hs6)
  setDT(wits_h6)
  setkey(wits_h6, importer, exporter, year, hs6)
  wits_h6 <- pref_h6[wits_h6]
  pref_nonmissing_wits <- wits_h6[!is.na(preferential_simple_avg), .N]
  message(
    "Post-merge (WITS+preferential): ",
    pref_nonmissing_wits, "/", nrow(wits_h6),
    " rows have non-missing preferential_simple_avg (",
    sprintf("%.2f", 100 * pref_nonmissing_wits / nrow(wits_h6)), "%)."
  )
}

setDT(baci_h6)
setkey(baci_h6, importer, exporter, year, hs6)
if (!is.data.table(wits_h6)) setDT(wits_h6)
setkey(wits_h6, importer, exporter, year, hs6)

merged <- wits_h6[baci_h6]

fwrite(merged, merged_file, na = "")

message("Merged rows: ", nrow(merged))
message("Merged file: ", merged_file)
message("Done.")
