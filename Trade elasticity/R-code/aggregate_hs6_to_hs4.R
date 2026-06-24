rm(list = ls())

library(data.table)

base_path <- "P:/leawin/Elasticiteter2503/Estimating LP Boehm et al"
hs6_dir <- file.path(base_path, "Merge_all_countries", "outputs")
hs4_dir <- file.path(base_path, "Merge_all_countries", "output_HS4")
dir.create(hs4_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1 && nzchar(args[[1]])) {
  hs6_file <- args[[1]]
} else {
  hs6_files <- list.files(
    hs6_dir,
    pattern = "^baci_wits_h6_merged_allcountries_.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hs6_files) == 0) {
    stop("No HS6 merged files found in: ", hs6_dir)
  }
  hs6_file <- hs6_files[which.max(file.info(hs6_files)$mtime)]
}

if (!file.exists(hs6_file)) {
  stop("HS6 file not found: ", hs6_file)
}

run_tag <- gsub("^baci_wits_h6_merged_allcountries_|\\.csv$", "", basename(hs6_file))
hs4_file <- file.path(hs4_dir, paste0("baci_wits_hs4_merged_allcountries_", run_tag, ".csv"))

message("Reading HS6 merged file: ", hs6_file)
dt <- fread(hs6_file)

required_cols <- c("importer", "exporter", "year", "hs6", "baci_trade_value")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("HS6 file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

for (nm in c("mfn_weighted_avg", "ahs_weighted_avg", "preferential_simple_avg")) {
  if (!(nm %in% names(dt))) {
    dt[, (nm) := NA_real_]
  }
}

dt[, hs4 := substr(sprintf("%06s", hs6), 1, 4)]
dt <- dt[!is.na(hs4) & nchar(hs4) == 4]

wmean_or_na <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

message("Aggregating to HS4 ...")
hs4 <- dt[
  ,
  .(
    baci_trade_value = sum(baci_trade_value, na.rm = TRUE),
    mfn_weighted_avg = wmean_or_na(mfn_weighted_avg, baci_trade_value),
    ahs_weighted_avg = wmean_or_na(ahs_weighted_avg, baci_trade_value),
    preferential_simple_avg = wmean_or_na(preferential_simple_avg, baci_trade_value)
  ),
  by = .(importer, exporter, year, hs4)
]

fwrite(hs4, hs4_file)
message("HS4 rows: ", nrow(hs4))
message("HS4 output file: ", hs4_file)
message("Done.")
