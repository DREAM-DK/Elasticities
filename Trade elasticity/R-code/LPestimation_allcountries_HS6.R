library(data.table)
library(fixest)

# ===================================================================
# Replication of Boehm, Levchenko & Pandalai-Nayar (AER 2023):
# "The Long and Short (Run) of Trade Elasticities" -- Figures 1 & 2.
#
# This is the HS6 version of LPestimation_allcountries_HS4.R. It fixes
# the three issues that prevented the HS4 pipeline from replicating the
# paper:
#
#   (1) PRODUCT LEVEL = HS6 (the paper's unit of observation is an HS6
#       code of a given revision, tracked over time). Only the
#       fixed-effect / multilateral-resistance terms are at HS4. The
#       HS4 pipeline collapsed trade and tariffs to HS4 with
#       trade-weighted average tariffs, which washes out the discrete
#       MFN tariff shocks that identify the elasticity and biases the
#       estimates toward zero (more so at longer horizons).
#       => Here the panel is importer x exporter x HS6; FE use HS4.
#
#   (2) ENDOGENOUS TARIFF = the bilateral APPLIED (effectively applied /
#       AHS) tariff tau_{ijpt}, instrumented by the MFN-binding tariff
#       change (Eq. 2.5). The HS4 pipeline used the statutory MFN rate
#       as the regressor for every observation, including PTA/control
#       pairs that never face MFN changes.
#       => Here ln_tau = log1p(applied rate); the instrument uses the
#          MFN change only on observations where applied == MFN at t and
#          t-1.
#
#   (3) CALENDAR-YEAR-AWARE LEADS/LAGS. BACI only records positive
#       trade flows, so a country-pair-product panel has gaps. The HS4
#       pipeline shifted by row position, so a "one-year" change could
#       span a multi-year gap and a horizon-h lead could land on the
#       wrong calendar year. Here every lead/lag is resolved by an
#       explicit join on (importer, exporter, hs6, year + offset), so a
#       horizon-h value is used only if the observation at exactly
#       year + h exists.
#
# Specification (matches the paper):
#   FE: exporter^hs4^year + importer^hs4^year + exporter^importer^hs4
#   Instrument (2.5): 1{tau=tau^MFN}_t * 1{tau=tau^MFN}_{t-1} * dln tau^MFN
#   Sample restriction (2.6) for the baseline specifications.
#   Figure 1: h = -H_pre..H, with/without pretrend controls.
#             For h<0 the dependent variable is the one-year change.
#   Figure 2: baseline IV (2.6 filter) and "all data / all tariffs 2SLS"
#             (instrument dh ln tau with d0 ln tau, no 2.6 filter).
#   95% CIs, SE clustered at country-pair-product (= HS6 panel id).
# ===================================================================

base_path <- "P:/leawin/Elasticiteter2503/Estimating LP Boehm et al"

data_file <- file.path(
  base_path,
  "Merge_all_countries",
  "outputs",
  "baci_wits_h6_merged_allcountries_20260611_133326.csv"
)
out_dir   <- file.path(base_path, "Merge_all_countries", "LP_results_HS6_paper")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
diag_file <- file.path(out_dir, paste0("allcountries_hs6_diagnostics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))

if (!file.exists(data_file)) stop("HS6 merged file not found: ", data_file)
message("Using HS6 merged file: ", data_file)

# ---------------------- USER SETTINGS ------------------------------
H                   <- 10            # max forward horizon
H_pre               <- 6             # number of negative (placebo) horizons
ci_level            <- 0.95
top_n_major_partner <- 10
mfn_equal_tol       <- 1e-8
sample_start_year   <- 1995
sample_end_year     <- 2018

# Units of mfn_weighted_avg / ahs_weighted_avg in the merged file.
# WITS' "Weighted Average" is a PERCENTAGE (5.5 = 5.5%).
tariff_unit         <- "percent"     # "percent" or "fraction"

sample_mode         <- "eu"       # "paper", "eu", "eu_plus_dk_dummy", "dk"
target_country      <- "DNK"

# TEST MODE: set to an integer to keep only the N largest importers (by
# number of tariff observations) for a quick end-to-end validation run.
# Set to NULL for the full paper sample. Can also be set without editing
# this file via the environment variable LP_TEST_N (e.g. LP_TEST_N=15).
test_n_importers    <- NULL
.env_test_n <- suppressWarnings(as.integer(Sys.getenv("LP_TEST_N", "")))
if (is.finite(.env_test_n) && .env_test_n > 0L) test_n_importers <- .env_test_n

fixest_threads      <- NULL          # NULL = fixest default
verbose_summary     <- FALSE
# -------------------------------------------------------------------

z_crit <- qnorm((1 + ci_level) / 2)
setDTthreads(0L)
if (!is.null(fixest_threads)) setFixest_nthreads(fixest_threads)

eu_importers <- c(
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
  "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
  "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"
)



sample_mode    <- tolower(sample_mode)
valid_modes    <- c("paper", "eu", "eu_plus_dk_dummy", "dk")
if (!(sample_mode %in% valid_modes))
  stop("sample_mode must be one of: ", paste(valid_modes, collapse = ", "))
target_country <- toupper(target_country)
uses_dummy     <- identical(sample_mode, "eu_plus_dk_dummy")
single_importer <- identical(sample_mode, "dk")

mode_tag <- switch(
  sample_mode,
  paper            = "paper",
  eu               = "eu",
  eu_plus_dk_dummy = paste0("eu_dummy_", target_country),
  dk               = paste0("country_",  target_country)
)
sample_tag <- paste0(
  "y",
  if (is.null(sample_start_year)) "min" else as.character(sample_start_year),
  "_",
  if (is.null(sample_end_year)) "max" else as.character(sample_end_year)
)
file_tag <- paste0(mode_tag, "_", sample_tag)

# ===================================================================
# 1) Load & filter
# ===================================================================
keep_cols_in <- c("importer", "exporter", "hs6", "year",
                  "baci_trade_value", "mfn_weighted_avg", "ahs_weighted_avg",
                  "preferential_simple_avg")
dt <- fread(data_file, select = keep_cols_in)

importer_keep <- switch(
  sample_mode,
  paper            = NULL,
  eu               = eu_importers,
  eu_plus_dk_dummy = eu_importers,
  dk               = target_country
)
if (!is.null(importer_keep)) dt <- dt[importer %in% importer_keep]
if (!is.null(sample_start_year)) dt <- dt[year >= sample_start_year]
if (!is.null(sample_end_year))   dt <- dt[year <= sample_end_year]

# Normalise HS6 to a 6-character string, derive HS4 (used for the FE).
dt[, hs6 := gsub("\\s+", "", as.character(hs6))]
short_idx <- which(nchar(dt$hs6) > 0L & nchar(dt$hs6) < 6L)
if (length(short_idx))
  set(dt, short_idx, "hs6",
      paste0(strrep("0", 6L - nchar(dt$hs6[short_idx])), dt$hs6[short_idx]))
long_idx <- which(nchar(dt$hs6) > 6L)
if (length(long_idx))
  set(dt, long_idx, "hs6", substr(dt$hs6[long_idx], 1L, 6L))
dt <- dt[nchar(hs6) == 6L]
dt[, hs4 := substr(hs6, 1L, 4L)]

message("Rows after load/filter: ", format(nrow(dt), big.mark = ","))

# The HS6 normalisation above (zero-padding short codes and truncating codes
# longer than six digits to their first six) can map several distinct
# original rows onto the same (importer, exporter, hs6, year) cell. So can
# duplicates introduced upstream in the BACI/WITS merge. The calendar-aware
# lead/lag joins in Section 3 (and the major-partner lag join in Section 4)
# REQUIRE a unique panel key; a non-unique key makes the keyed join return
# more rows than dt and triggers the "Supplied N items to be assigned to M
# items" error. Collapse any duplicates now: trade values are summed and the
# tariff rates are aggregated as trade-weighted means (falling back to a
# simple mean over non-missing rates when there is no positive trade weight).
if (anyDuplicated(dt, by = c("importer", "exporter", "hs6", "year")) > 0L) {
  n_before <- nrow(dt)
  dt[, `:=`(
    .vw_mfn  = fifelse(is.na(mfn_weighted_avg),        NA_real_, baci_trade_value),
    .vw_ahs  = fifelse(is.na(ahs_weighted_avg),        NA_real_, baci_trade_value),
    .vw_pref = fifelse(is.na(preferential_simple_avg), NA_real_, baci_trade_value)
  )]
  dt <- dt[, .(
    hs4              = hs4[1L],
    baci_trade_value = sum(baci_trade_value, na.rm = TRUE),
    .num_mfn  = sum(mfn_weighted_avg        * .vw_mfn,  na.rm = TRUE),
    .den_mfn  = sum(.vw_mfn,  na.rm = TRUE),
    .cnt_mfn  = sum(!is.na(mfn_weighted_avg)),
    .sum_mfn  = sum(mfn_weighted_avg,        na.rm = TRUE),
    .num_ahs  = sum(ahs_weighted_avg        * .vw_ahs,  na.rm = TRUE),
    .den_ahs  = sum(.vw_ahs,  na.rm = TRUE),
    .cnt_ahs  = sum(!is.na(ahs_weighted_avg)),
    .sum_ahs  = sum(ahs_weighted_avg,        na.rm = TRUE),
    .num_pref = sum(preferential_simple_avg * .vw_pref, na.rm = TRUE),
    .den_pref = sum(.vw_pref, na.rm = TRUE),
    .cnt_pref = sum(!is.na(preferential_simple_avg)),
    .sum_pref = sum(preferential_simple_avg, na.rm = TRUE)
  ), by = .(importer, exporter, hs6, year)]
  dt[, `:=`(
    mfn_weighted_avg        = fifelse(.den_mfn  > 0, .num_mfn  / .den_mfn,
                               fifelse(.cnt_mfn  > 0L, .sum_mfn  / .cnt_mfn,  NA_real_)),
    ahs_weighted_avg        = fifelse(.den_ahs  > 0, .num_ahs  / .den_ahs,
                               fifelse(.cnt_ahs  > 0L, .sum_ahs  / .cnt_ahs,  NA_real_)),
    preferential_simple_avg = fifelse(.den_pref > 0, .num_pref / .den_pref,
                               fifelse(.cnt_pref > 0L, .sum_pref / .cnt_pref, NA_real_))
  )]
  dt[, c(".num_mfn", ".den_mfn", ".cnt_mfn", ".sum_mfn",
         ".num_ahs", ".den_ahs", ".cnt_ahs", ".sum_ahs",
         ".num_pref", ".den_pref", ".cnt_pref", ".sum_pref") := NULL]
  message("Collapsed duplicate (importer, exporter, hs6, year) cells: ",
          format(n_before, big.mark = ","), " -> ",
          format(nrow(dt), big.mark = ","), " rows.")
  invisible(gc(verbose = FALSE))
}

# ===================================================================
# 2) Tariff / log-variables
# ===================================================================
if (tariff_unit == "percent") {
  dt[, `:=`(mfn_rate = mfn_weighted_avg / 100,
            ahs_rate = ahs_weighted_avg / 100,
            preferential_rate = preferential_simple_avg / 100)]
} else if (tariff_unit == "fraction") {
  dt[, `:=`(mfn_rate = mfn_weighted_avg,
            ahs_rate = ahs_weighted_avg,
            preferential_rate = preferential_simple_avg)]
} else stop("tariff_unit must be 'percent' or 'fraction'.")

# The bilateral applied tariff faced by (i, j, p) is the effectively
# applied (AHS) rate. When no preference applies, AHS == MFN.
dt[, applied_rate := ahs_rate]

dt[, `:=`(
  v           = baci_trade_value,
  ln_x        = log1p(baci_trade_value),
  ln_tau      = log1p(applied_rate),     # endogenous regressor base (Eq. 2.4)
  ln_mfn      = log1p(mfn_rate),         # instrument base (Eq. 2.5)
  mfn_treated = as.integer(
    !is.na(mfn_rate) & !is.na(applied_rate) &
      abs(mfn_rate - applied_rate) <= mfn_equal_tol
  )
)]

# Diagnostics.
diag_dt <- dt[, .(
  n_obs = .N,
  n_with_applied = sum(!is.na(ln_tau)),
  n_with_mfn = sum(!is.na(ln_mfn)),
  share_mfn_binding = mean(mfn_treated == 1L, na.rm = TRUE),
  n_pref_nonmissing = sum(!is.na(preferential_rate)),
  share_pref_nonmissing = mean(!is.na(preferential_rate))
)]
fwrite(diag_dt, diag_file)
message("Saved diagnostics: ", diag_file)
if (diag_dt$n_pref_nonmissing[1] == 0L) {
  warning(
    "preferential_simple_avg is entirely missing after year/sample filtering (",
    if (is.null(sample_start_year)) "min" else sample_start_year, "-",
    if (is.null(sample_end_year)) "max" else sample_end_year, "). ",
    "This usually means the preferential source data does not overlap the estimation years."
  )
}

# Drop importers that have NO tariff data at all (they can never enter
# an IV regression and only inflate memory).
imp_with_tariff <- dt[!is.na(mfn_rate) | !is.na(applied_rate), unique(importer)]
if (length(imp_with_tariff) < uniqueN(dt$importer)) {
  n_before <- nrow(dt)
  dt <- dt[importer %in% imp_with_tariff]
  message("Pruned importers without tariff data: kept ",
          uniqueN(dt$importer), " importers, ",
          format(nrow(dt), big.mark = ","), "/",
          format(n_before, big.mark = ","), " rows.")
}

# Optional TEST subset: keep only the N largest importers by tariff obs.
if (!is.null(test_n_importers)) {
  imp_rank <- dt[!is.na(mfn_rate) | !is.na(applied_rate), .N, by = importer
                 ][order(-N)][seq_len(min(test_n_importers, .N)), importer]
  dt <- dt[importer %in% imp_rank]
  message("TEST MODE: kept ", length(imp_rank), " importers (",
          paste(imp_rank, collapse = ", "), "), ",
          format(nrow(dt), big.mark = ","), " rows.")
}

# ===================================================================
# 3) Calendar-year-aware leads / lags (panel = importer x exporter x HS6)
# ===================================================================
gcols <- c("importer", "exporter", "hs6")
setkeyv(dt, c(gcols, "year"))

# Hard guarantee that the panel key is unique. Every lead/lag below resolves
# values via a keyed join on (importer, exporter, hs6, year); a duplicate key
# would make the join return more rows than dt and break the assignments.
if (anyDuplicated(dt, by = c(gcols, "year")) > 0L)
  stop("Panel key (importer, exporter, hs6, year) is not unique after the ",
       "Section 1 collapse. Cannot build calendar-aware leads/lags.")

horizons_pos <- 0:H                # 0 ... H
horizons_neg <- -H_pre:-1L         # -H_pre ... -1
horizons_all <- c(horizons_neg, horizons_pos)

h_tag <- function(h) ifelse(h < 0, paste0("m", abs(h)), as.character(h))
y_cols   <- paste0("y_h",   h_tag(horizons_all))
tau_cols <- paste0("tau_h", h_tag(horizons_all))

# Lookup table holding the base series, keyed by (group, year). The
# panel key is unique per (importer, exporter, hs6, year).
lk <- dt[, c(gcols, "year", "ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
         with = FALSE]
setkeyv(lk, c(gcols, "year"))

# For each row of dt, fetch base columns observed at (group, year + off).
# Non-existent (group, year + off) cells resolve to NA -> the row drops
# from any regression that uses them. This is what makes the horizons
# calendar-correct in the presence of BACI gaps.
add_offset <- function(off, src_cols, dst_cols) {
  ktab <- dt[, c(gcols, "year"), with = FALSE]
  ktab[, year := year + off]
  res <- lk[ktab, on = c(gcols, "year")]
  for (k in seq_along(dst_cols)) dt[, (dst_cols[k]) := res[[src_cols[k]]]]
  invisible(NULL)
}

# Permanent lag columns (offsets -1 and -2).
add_offset(-1L, c("ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
                c("ln_x_l1", "ln_tau_l1", "ln_mfn_l1", "mfn_treated_l1"))
add_offset(-2L, c("ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
                c("ln_x_l2", "ln_tau_l2", "ln_mfn_l2", "mfn_treated_l2"))

# Impact / pretrend regressors (applied tariff) and instruments (MFN).
dt[, `:=`(
  d0_ln_tau   = ln_tau - ln_tau_l1,        # Delta_0 ln tau  (applied)
  dlag_ln_tau = ln_tau_l1 - ln_tau_l2,     # pretrend tariff change
  dlag_ln_x   = ln_x_l1 - ln_x_l2,         # pretrend trade change
  d0_ln_mfn   = ln_mfn - ln_mfn_l1,        # MFN impact change
  dlag_ln_mfn = ln_mfn_l1 - ln_mfn_l2      # MFN pretrend change
)]

# Baseline instrument (Eq. 2.5): MFN change for observations that trade
# on MFN terms at both t and t-1; zero for the control group.
dt[, d0_ln_tau_instr_25 := fifelse(
  mfn_treated == 1L & mfn_treated_l1 == 1L, d0_ln_mfn, 0)]
# Lag of the baseline instrument (used as the pretrend instrument).
dt[, dlag_d0_ln_tau_instr_25 := fifelse(
  mfn_treated_l1 == 1L & mfn_treated_l2 == 1L, dlag_ln_mfn, 0)]
# All-data / all-tariffs pretrend instrument = lag of the applied impact
# change (mathematically equal to dlag_ln_tau).
dt[, dlag_d0_ln_tau := dlag_ln_tau]

# Horizon-specific dependent variables (trade) and tariff changes.
#   h >= 0 : y_h = ln X_{t+h} - ln X_{t-1};  tau_h = ln tau_{t+h} - ln tau_{t-1}
#   h <  0 : y_h = ln X_{t+h} - ln X_{t+h-1} (one-year placebo change)
for (i in seq_along(horizons_all)) {
  h  <- horizons_all[i]
  yc <- y_cols[i]; tc <- tau_cols[i]
  if (h >= 0L) {
    add_offset(h, c("ln_x", "ln_tau"), c(".__Xh", ".__TAUh"))
    dt[, (yc) := .__Xh   - ln_x_l1]
    dt[, (tc) := .__TAUh - ln_tau_l1]
    dt[, c(".__Xh", ".__TAUh") := NULL]
  } else {
    add_offset(h,      c("ln_x", "ln_tau"), c(".__Xh",  ".__TAUh"))
    add_offset(h - 1L, c("ln_x", "ln_tau"), c(".__Xh1", ".__TAUh1"))
    dt[, (yc) := .__Xh   - .__Xh1]
    dt[, (tc) := .__TAUh - .__TAUh1]
    dt[, c(".__Xh", ".__TAUh", ".__Xh1", ".__TAUh1") := NULL]
  }
}

rm(lk); invisible(gc(verbose = FALSE))

# ===================================================================
# 4) Major-partner flags (paper: importer-HS4 level)
# ===================================================================
major_agg <- dt[, .(.v = sum(v, na.rm = TRUE)),
                by = .(importer, year, exporter)]
setorder(major_agg, importer, year, -`.v`)
major_agg[, major_agg_t := as.integer(seq_len(.N) <= top_n_major_partner),
          by = .(importer, year)]
major_agg[, .v := NULL]

major_prod <- dt[, .(.v = sum(v, na.rm = TRUE)),
                 by = .(importer, year, hs4, exporter)]
setorder(major_prod, importer, year, hs4, -`.v`)
major_prod[, major_prod_t := as.integer(seq_len(.N) <= top_n_major_partner),
           by = .(importer, year, hs4)]
major_prod[, .v := NULL]

dt[, `:=`(major_agg_t = 0L, major_prod_t = 0L)]
dt[major_agg,  on = c("importer", "year", "exporter"),
   major_agg_t  := i.major_agg_t]
dt[major_prod, on = c("importer", "year", "hs4", "exporter"),
   major_prod_t := i.major_prod_t]
dt[is.na(major_agg_t),  major_agg_t  := 0L]
dt[is.na(major_prod_t), major_prod_t := 0L]
rm(major_agg, major_prod); invisible(gc(verbose = FALSE))

setkeyv(dt, c(gcols, "year"))
# Calendar-aware t-1 of the major-partner flags.
lk2 <- dt[, c(gcols, "year", "major_agg_t", "major_prod_t"), with = FALSE]
setkeyv(lk2, c(gcols, "year"))
ktab2 <- dt[, c(gcols, "year"), with = FALSE]
ktab2[, year := year - 1L]
res2 <- lk2[ktab2, on = c(gcols, "year")]
dt[, major_agg_t1  := res2$major_agg_t]
dt[, major_prod_t1 := res2$major_prod_t]
rm(lk2, ktab2, res2); invisible(gc(verbose = FALSE))

dt[, major_partner_any := as.integer(
  (fcoalesce(major_agg_t1, 0L) + fcoalesce(major_prod_t1, 0L) +
     major_agg_t + major_prod_t) > 0L
)]
# Eq. (2.6): drop obs binding at t and t-1 AND with any major-partner flag.
dt[, eq26_drop := as.integer(
  mfn_treated == 1L & mfn_treated_l1 == 1L & major_partner_any == 1L
)]

# ===================================================================
# 5) Dummy interactions (only in eu_plus_dk_dummy mode)
# ===================================================================
if (uses_dummy) {
  dt[, dk := as.integer(importer == target_country)]
  dt[, `:=`(
    d0_ln_tau_dk                 = d0_ln_tau                 * dk,
    dlag_ln_tau_dk               = dlag_ln_tau               * dk,
    dlag_ln_x_dk                 = dlag_ln_x                 * dk,
    d0_ln_tau_instr_25_dk        = d0_ln_tau_instr_25        * dk,
    dlag_d0_ln_tau_instr_25_dk   = dlag_d0_ln_tau_instr_25   * dk,
    dlag_d0_ln_tau_dk            = dlag_d0_ln_tau            * dk
  )]
  dt[, (paste0(tau_cols, "_dk")) := lapply(.SD, function(x) x * dk),
     .SDcols = tau_cols]
}

# Cluster id = country-pair-product (HS6 panel).
dt[, pair_product_id := .GRP, by = gcols]

# Drop columns no longer needed for estimation.
drop_cols <- c("baci_trade_value", "mfn_weighted_avg", "ahs_weighted_avg",
               "preferential_simple_avg", "mfn_rate", "ahs_rate",
               "preferential_rate", "applied_rate",
               "v", "ln_x", "ln_tau", "ln_mfn",
               "ln_x_l2", "ln_tau_l2", "ln_mfn_l1", "ln_mfn_l2",
               "d0_ln_mfn", "dlag_ln_mfn",
               "mfn_treated_l2",
               "major_agg_t", "major_prod_t",
               "major_agg_t1", "major_prod_t1", "major_partner_any")
dt[, (intersect(drop_cols, names(dt))) := NULL]
invisible(gc(verbose = FALSE))

if (verbose_summary) print(skimr::skim(dt))

# ===================================================================
# 6) Estimation -- cluster-robust 2SLS via FE demeaning
# ===================================================================
# NOTE ON THE ESTIMATOR. We do NOT use fixest's IV syntax here. The paper's
# instrument equals the endogenous tariff change on the treated (MFN-binding)
# observations, and many control (PTA) observations have flat tariffs, so the
# first stage is often (near-)perfect. In that case fixest's automatic IV
# F-statistic crashes ("invert_posdef_mat"), aborting the whole estimation
# even though the coefficient is well defined. We therefore use fixest only
# to absorb the (high-dimensional) fixed effects via demeaning, and compute
# the 2SLS coefficient and the cluster-robust variance directly. This matches
# fixest's IV point estimates and standard errors to numerical precision on
# well-behaved models, and is robust to perfect first stages.

# Fixed-effect groupings (multilateral resistance terms at HS4).
fe_specs <- list(c("exporter", "hs4", "year"),
                 c("importer", "hs4", "year"),
                 c("exporter", "importer", "hs4"))

# Cluster-robust 2SLS on FE-demeaned data. Returns the full coefficient
# vector and variance matrix (named by regressor) so the caller can read a
# single coefficient or a linear combination (dummy mode: main + interaction).
fe_iv <- function(dat, yname, endo, exog, inst, cluster_col) {
  exog2 <- setdiff(exog, "1")
  vars  <- unique(c(yname, endo, exog2, inst))
  cc <- stats::complete.cases(dat[, ..vars])
  d  <- dat[cc]
  if (nrow(d) < 5L) return(NULL)
  fes <- lapply(fe_specs, function(cols) do.call(paste, c(as.list(d[, ..cols]), sep = "\r")))
  M  <- as.matrix(d[, ..vars])
  Md <- fixest::demean(M, fes)
  yv <- Md[, yname]
  Xe <- Md[, endo, drop = FALSE]
  Xx <- if (length(exog2)) Md[, exog2, drop = FALSE] else NULL
  Zi <- Md[, inst, drop = FALSE]
  X  <- cbind(Xe, Xx); Z <- cbind(Zi, Xx)
  ZtZ_inv <- tryCatch(solve(crossprod(Z)), error = function(e) NULL)
  if (is.null(ZtZ_inv)) return(NULL)
  PzX <- Z %*% (ZtZ_inv %*% crossprod(Z, X))
  A <- tryCatch(solve(crossprod(PzX, X)), error = function(e) NULL)
  if (is.null(A)) return(NULL)
  beta  <- A %*% crossprod(PzX, yv)
  resid <- as.vector(yv - X %*% beta)
  g <- d[[cluster_col]]
  S <- rowsum(PzX * resid, g)
  meat <- crossprod(S)
  N <- nrow(d); G <- nrow(S)
  K <- ncol(X) + sum(vapply(fes, function(f) length(unique(f)), 1L)) - (length(fes) - 1L)
  cc_adj <- (G / (G - 1)) * ((N - 1) / (N - K))
  V <- cc_adj * (A %*% meat %*% t(A))
  cn <- colnames(X)
  b <- as.vector(beta); names(b) <- cn
  dimnames(V) <- list(cn, cn)
  list(b = b, V = V, n = N, n_clusters = G)
}

# IRF point estimate / SE for a target (or sum of targets in dummy mode).
irf_from <- function(res, targets, h_val) {
  if (is.null(res) || !all(targets %in% names(res$b)) || any(is.na(res$b[targets])))
    return(data.table(h = h_val, beta = NA_real_, se = NA_real_))
  Vsub <- res$V[targets, targets, drop = FALSE]
  data.table(h = h_val, beta = sum(res$b[targets]), se = sqrt(sum(Vsub)))
}

# RHS (endogenous / exogenous / instrument / target) definitions.
fig1_rhs <- function(with_pretrend) {
  if (with_pretrend) {
    if (uses_dummy)
      list(endo = c("d0_ln_tau", "d0_ln_tau_dk", "dlag_ln_tau", "dlag_ln_tau_dk"),
           exog = c("dlag_ln_x", "dlag_ln_x_dk"),
           inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk",
                    "dlag_d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25_dk"),
           targets = c("d0_ln_tau", "d0_ln_tau_dk"))
    else
      list(endo = c("d0_ln_tau", "dlag_ln_tau"), exog = "dlag_ln_x",
           inst = c("d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25"),
           targets = "d0_ln_tau")
  } else {
    if (uses_dummy)
      list(endo = c("d0_ln_tau", "d0_ln_tau_dk"), exog = "1",
           inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk"),
           targets = c("d0_ln_tau", "d0_ln_tau_dk"))
    else
      list(endo = "d0_ln_tau", exog = "1", inst = "d0_ln_tau_instr_25",
           targets = "d0_ln_tau")
  }
}

# Figure 2 always includes the pretrend controls.
fig2_rhs <- function(h_val, spec) {
  tau    <- paste0("tau_h", h_tag(h_val))
  tau_dk <- paste0(tau, "_dk")
  if (spec == "baseline") {
    if (uses_dummy)
      list(endo = c(tau, tau_dk, "dlag_ln_tau", "dlag_ln_tau_dk"),
           exog = c("dlag_ln_x", "dlag_ln_x_dk"),
           inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk",
                    "dlag_d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25_dk"),
           targets = c(tau, tau_dk))
    else
      list(endo = c(tau, "dlag_ln_tau"), exog = "dlag_ln_x",
           inst = c("d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25"),
           targets = tau)
  } else {  # all_data: instrument dh ln tau with d0 ln tau (no 2.6 filter).
    if (uses_dummy)                            # at h = 0, d0 ln tau == tau_h0
      list(endo = c(tau, tau_dk),              # so 2SLS collapses to OLS (handled
           exog = c("dlag_ln_x", "dlag_ln_x_dk",   # automatically by fe_iv).
                    "dlag_ln_tau", "dlag_ln_tau_dk"),
           inst = c("d0_ln_tau", "d0_ln_tau_dk"),
           targets = c(tau, tau_dk))
    else
      list(endo = tau, exog = c("dlag_ln_x", "dlag_ln_tau"),
           inst = "d0_ln_tau", targets = tau)
  }
}

if (!"eq26_drop" %in% names(dt))
  stop("eq26_drop was not built -- Section 4 did not complete. Re-run from the top.")
dt_paper_sample <- dt[eq26_drop == 0L]

write_irf <- function(tab, out_path) {
  tab <- copy(tab)
  tab[, `:=`(ci_lo = beta - z_crit * se,
             ci_hi = beta + z_crit * se)]
  setorder(tab, h)
  fwrite(tab, out_path, na = "")
  message("Wrote: ", basename(out_path))
}

# ===================================================================
# 7) Figure 1: tariff (Eq. 2.3) and trade (Eq. 2.2), with/without pretrend
# ===================================================================
estimate_fig1 <- function(panel, with_pretrend) {
  cols <- if (panel == "tariff") tau_cols else y_cols
  spec <- fig1_rhs(with_pretrend)
  message("  Fig.1 ", panel, " | pretrend = ", with_pretrend,
          " | ", length(cols), " horizons")
  rbindlist(lapply(seq_along(horizons_all), function(i) {
    res <- fe_iv(dt_paper_sample, cols[i], spec$endo, spec$exog, spec$inst,
                 "pair_product_id")
    out <- irf_from(res, spec$targets, horizons_all[i])
    rm(res); invisible(gc(verbose = FALSE))
    out
  }))
}

message("=== Figure 1 [", sample_mode, "] ===")
for (panel in c("tariff", "trade")) {
  for (with_pret in c(FALSE, TRUE)) {
    tag <- if (with_pret) "with_pretrend" else "no_pretrend"
    res <- estimate_fig1(panel, with_pret)
    write_irf(res, file.path(out_dir,
              paste0("fig1_", panel, "_", file_tag, "_", tag, ".csv")))
  }
}

# ===================================================================
# 8) Figure 2: trade elasticity from Eq. (2.4), baseline vs all-data
# ===================================================================
estimate_fig2 <- function(spec) {
  data_in <- if (spec == "baseline") dt_paper_sample else dt
  message("  Fig.2 ", spec, " | n = ", format(nrow(data_in), big.mark = ","))
  rbindlist(lapply(horizons_pos, function(h_val) {
    rhs <- fig2_rhs(h_val, spec)
    res <- fe_iv(data_in, paste0("y_h", h_tag(h_val)),
                 rhs$endo, rhs$exog, rhs$inst, "pair_product_id")
    out <- irf_from(res, rhs$targets, h_val)
    rm(res); invisible(gc(verbose = FALSE))
    out
  }))
}

message("=== Figure 2 [", sample_mode, "] ===")
for (spec in c("baseline", "all_data")) {
  res <- estimate_fig2(spec)
  write_irf(res, file.path(out_dir,
            paste0("fig2_elast_", file_tag, "_", spec, ".csv")))
}

message("Done. Outputs in: ", out_dir)
