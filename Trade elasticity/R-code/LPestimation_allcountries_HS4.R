library(data.table)
library(fixest)

# ===================================================================
# Replication of Boehm, Levchenko & Pandalai-Nayar (AER 2023):
# "The Long and Short (Run) of Trade Elasticities" -- Figures 1 & 2.
#
# Specification details that match the paper:
#   * Panel is at HS6 (each HS6 code of a given revision, tracked
#     over time); fixed effects are at HS4.
#     FE: exporter^hs4^year + importer^hs4^year + exporter^importer^hs4
#   * Major / minor partner classification at importer-HS4 (top-10 by
#     trade value), in aggregate and at the product level.
#   * Pre-trend controls: lag of d ln tau (instrumented with the lag
#     of the baseline instrument) and lag of d ln X.
#   * Equations (2.2)/(2.3) and (2.4) are estimated by IV using the
#     baseline instrument (2.5). Sample restriction (2.6) is applied
#     for the "baseline" specifications.
#   * Figure 1: h = -H_pre .. H, two lines (with / without pretrend).
#     For h < 0 the dependent variable is the one-year change (placebo).
#   * Figure 2: two lines: baseline IV with (2.6) filter, and
#     "all data / all tariffs 2SLS" where Delta_h ln tau is
#     instrumented by Delta_0 ln tau on the full sample.
#   * 95% confidence intervals, SE clustered at country-pair-product.
#
# Sample modes (set via `sample_mode`):
#   "paper"            : all 183 economies (needs the merged file
#                        to NOT be EU-filtered -- see notes below)
#   "eu"               : EU importers only
#   "eu_plus_dk_dummy" : EU importers + dummy for target_country;
#                        every reported IRF is (beta_main + beta_dummy)
#   "dk"               : only target_country
#
# NOTE: Baci_wits_merge_H6.R as shipped filters BACI+WITS to
# `eu_importers`. To use sample_mode = "paper" you must re-build the
# merged file with that filter removed.
# ===================================================================

base_path <- "P:/leawin/Elasticiteter2503/Estimating LP Boehm et al"
data_file <- file.path(
  base_path,
  "Merge_all_countries",
  "output_HS4",
  "baci_wits_hs4_merged_allcountries_20260604_133857.csv.csv"
)
out_dir   <- file.path(base_path, "Merge_all_countries", "LP_results_HS4_paper")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
diag_file <- file.path(out_dir, paste0("allcountries_hs4_preferential_diagnostics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))

if (!file.exists(data_file)) stop("HS4 merged file not found: ", data_file)
message("Using HS4 merged file: ", data_file)

# ---------------------- USER SETTINGS ------------------------------
H                   <- 10            # max forward horizon
H_pre               <- 6             # number of negative (placebo) horizons
ci_level            <- 0.95
top_n_major_partner <- 10
mfn_equal_tol       <- 1e-8
sample_start_year   <- 1995
sample_end_year     <- 2018

# Units of mfn_weighted_avg / ahs_weighted_avg in the merged file.
# WITS' "Weighted Average" is typically a PERCENTAGE (5.5 = 5.5%);
# switch to "fraction" if your data is already in fractions (0.055).
tariff_unit         <- "percent"     # "percent" or "fraction"

sample_mode         <- "paper"          # "paper", "eu", "eu_plus_dk_dummy", "dk"
target_country      <- "DNK"

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
keep_cols_in <- c("importer", "exporter", "hs4", "year",
                  "baci_trade_value", "mfn_weighted_avg", "ahs_weighted_avg", "preferential_simple_avg")
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

if (sample_mode == "paper") {
  n_imp <- uniqueN(dt$importer)
  if (n_imp <= length(eu_importers) + 5L) {
    warning(
      "sample_mode = 'paper' but the merged file contains only ", n_imp,
      " importers (looks EU-filtered). Re-build the merge without the ",
      "EU importer filter for a true paper-baseline sample."
    )
  }
}

# Normalise HS4 from merged HS4 data.
dt[, hs4 := gsub("\\s+", "", as.character(hs4))]
short_idx <- which(nchar(dt$hs4) < 4L)
if (length(short_idx))
  set(dt, short_idx, "hs4",
      paste0(strrep("0", 4L - nchar(dt$hs4[short_idx])), dt$hs4[short_idx]))
long_idx <- which(nchar(dt$hs4) > 4L)
if (length(long_idx))
  set(dt, long_idx, "hs4", substr(dt$hs4[long_idx], 1L, 4L))

message("Rows after load/filter: ", format(nrow(dt), big.mark = ","))

# ===================================================================
# 2) Tariff / log-variables (HS4 version with HS6 method)
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

dt[, applied_rate := ahs_rate]

dt[, `:=`(
  v           = baci_trade_value,
  ln_x        = log1p(baci_trade_value),
  ln_tau      = log1p(applied_rate),
  ln_mfn      = log1p(mfn_rate),
  mfn_treated = as.integer(
    !is.na(mfn_rate) & !is.na(applied_rate) &
      abs(mfn_rate - applied_rate) <= mfn_equal_tol
  )
)]

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

imp_with_tariff <- dt[!is.na(mfn_rate) | !is.na(applied_rate), unique(importer)]
if (length(imp_with_tariff) < uniqueN(dt$importer)) {
  n_before <- nrow(dt)
  dt <- dt[importer %in% imp_with_tariff]
  message("Pruned importers without tariff data: kept ",
          uniqueN(dt$importer), " importers, ",
          format(nrow(dt), big.mark = ","), "/",
          format(n_before, big.mark = ","), " rows.")
}

# ===================================================================
# 3) Calendar-year-aware leads / lags (panel = importer x exporter x HS4)
# ===================================================================
gcols <- c("importer", "exporter", "hs4")
setkeyv(dt, c(gcols, "year"))

horizons_pos <- 0:H
horizons_neg <- -H_pre:-1L
horizons_all <- c(horizons_neg, horizons_pos)

h_tag <- function(h) ifelse(h < 0, paste0("m", abs(h)), as.character(h))
y_cols   <- paste0("y_h",   h_tag(horizons_all))
tau_cols <- paste0("tau_h", h_tag(horizons_all))

lk <- dt[, c(gcols, "year", "ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
         with = FALSE]
setkeyv(lk, c(gcols, "year"))

add_offset <- function(off, src_cols, dst_cols) {
  ktab <- dt[, c(gcols, "year"), with = FALSE]
  ktab[, year := year + off]
  res <- lk[ktab, on = c(gcols, "year")]
  for (k in seq_along(dst_cols)) dt[, (dst_cols[k]) := res[[src_cols[k]]]]
  invisible(NULL)
}

add_offset(-1L, c("ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
           c("ln_x_l1", "ln_tau_l1", "ln_mfn_l1", "mfn_treated_l1"))
add_offset(-2L, c("ln_x", "ln_tau", "ln_mfn", "mfn_treated"),
           c("ln_x_l2", "ln_tau_l2", "ln_mfn_l2", "mfn_treated_l2"))

dt[, `:=`(
  d0_ln_tau   = ln_tau - ln_tau_l1,
  dlag_ln_tau = ln_tau_l1 - ln_tau_l2,
  dlag_ln_x   = ln_x_l1 - ln_x_l2,
  d0_ln_mfn   = ln_mfn - ln_mfn_l1,
  dlag_ln_mfn = ln_mfn_l1 - ln_mfn_l2
)]

dt[, d0_ln_tau_instr_25 := fifelse(
  mfn_treated == 1L & mfn_treated_l1 == 1L, d0_ln_mfn, 0)]
dt[, dlag_d0_ln_tau_instr_25 := fifelse(
  mfn_treated_l1 == 1L & mfn_treated_l2 == 1L, dlag_ln_mfn, 0)]
dt[, dlag_d0_ln_tau := dlag_ln_tau]

for (i in seq_along(horizons_all)) {
  h  <- horizons_all[i]
  yc <- y_cols[i]
  tc <- tau_cols[i]
  if (h >= 0L) {
    add_offset(h, c("ln_x", "ln_tau"), c(".__Xh", ".__TAUh"))
    dt[, (yc) := .__Xh - ln_x_l1]
    dt[, (tc) := .__TAUh - ln_tau_l1]
    dt[, c(".__Xh", ".__TAUh") := NULL]
  } else {
    add_offset(h, c("ln_x", "ln_tau"), c(".__Xh", ".__TAUh"))
    add_offset(h - 1L, c("ln_x", "ln_tau"), c(".__Xh1", ".__TAUh1"))
    dt[, (yc) := .__Xh - .__Xh1]
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
dt[major_agg, on = c("importer", "year", "exporter"),
   major_agg_t := i.major_agg_t]
dt[major_prod, on = c("importer", "year", "hs4", "exporter"),
   major_prod_t := i.major_prod_t]
dt[is.na(major_agg_t), major_agg_t := 0L]
dt[is.na(major_prod_t), major_prod_t := 0L]
rm(major_agg, major_prod); invisible(gc(verbose = FALSE))

setkeyv(dt, c(gcols, "year"))
lk2 <- dt[, c(gcols, "year", "major_agg_t", "major_prod_t"), with = FALSE]
setkeyv(lk2, c(gcols, "year"))
ktab2 <- dt[, c(gcols, "year"), with = FALSE]
ktab2[, year := year - 1L]
res2 <- lk2[ktab2, on = c(gcols, "year")]
dt[, major_agg_t1 := res2$major_agg_t]
dt[, major_prod_t1 := res2$major_prod_t]
rm(lk2, ktab2, res2); invisible(gc(verbose = FALSE))

dt[, major_partner_any := as.integer(
  (fcoalesce(major_agg_t1, 0L) + fcoalesce(major_prod_t1, 0L) +
     major_agg_t + major_prod_t) > 0L
)]
dt[, eq26_drop := as.integer(
  mfn_treated == 1L & mfn_treated_l1 == 1L & major_partner_any == 1L
)]

# ===================================================================
# 5) Dummy interactions (only in eu_plus_dk_dummy mode)
# ===================================================================
if (uses_dummy) {
  dt[, dk := as.integer(importer == target_country)]
  dt[, `:=`(
    d0_ln_tau_dk                 = d0_ln_tau * dk,
    dlag_ln_tau_dk               = dlag_ln_tau * dk,
    dlag_ln_x_dk                 = dlag_ln_x * dk,
    d0_ln_tau_instr_25_dk        = d0_ln_tau_instr_25 * dk,
    dlag_d0_ln_tau_instr_25_dk   = dlag_d0_ln_tau_instr_25 * dk,
    dlag_d0_ln_tau_dk            = dlag_d0_ln_tau * dk
  )]
  dt[, (paste0(tau_cols, "_dk")) := lapply(.SD, function(x) x * dk),
     .SDcols = tau_cols]
}

dt[, pair_product_id := .GRP, by = gcols]

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
fe_specs <- list(c("exporter", "hs4", "year"),
                 c("importer", "hs4", "year"),
                 c("exporter", "importer", "hs4"))

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
  X  <- cbind(Xe, Xx)
  Z  <- cbind(Zi, Xx)
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
  N <- nrow(d)
  G <- nrow(S)
  K <- ncol(X) + sum(vapply(fes, function(f) length(unique(f)), 1L)) - (length(fes) - 1L)
  cc_adj <- (G / (G - 1)) * ((N - 1) / (N - K))
  V <- cc_adj * (A %*% meat %*% t(A))
  cn <- colnames(X)
  b <- as.vector(beta)
  names(b) <- cn
  dimnames(V) <- list(cn, cn)
  list(b = b, V = V, n = N, n_clusters = G)
}

irf_from <- function(res, targets, h_val) {
  if (is.null(res) || !all(targets %in% names(res$b)) || any(is.na(res$b[targets]))) {
    return(data.table(h = h_val, beta = NA_real_, se = NA_real_))
  }
  Vsub <- res$V[targets, targets, drop = FALSE]
  data.table(h = h_val, beta = sum(res$b[targets]), se = sqrt(sum(Vsub)))
}

fig1_rhs <- function(with_pretrend) {
  if (with_pretrend) {
    if (uses_dummy) {
      list(
        endo = c("d0_ln_tau", "d0_ln_tau_dk", "dlag_ln_tau", "dlag_ln_tau_dk"),
        exog = c("dlag_ln_x", "dlag_ln_x_dk"),
        inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk",
                 "dlag_d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25_dk"),
        targets = c("d0_ln_tau", "d0_ln_tau_dk")
      )
    } else {
      list(
        endo = c("d0_ln_tau", "dlag_ln_tau"),
        exog = "dlag_ln_x",
        inst = c("d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25"),
        targets = "d0_ln_tau"
      )
    }
  } else {
    if (uses_dummy) {
      list(
        endo = c("d0_ln_tau", "d0_ln_tau_dk"),
        exog = "1",
        inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk"),
        targets = c("d0_ln_tau", "d0_ln_tau_dk")
      )
    } else {
      list(endo = "d0_ln_tau", exog = "1", inst = "d0_ln_tau_instr_25", targets = "d0_ln_tau")
    }
  }
}

fig2_rhs <- function(h_val, spec) {
  tau <- paste0("tau_h", h_tag(h_val))
  tau_dk <- paste0(tau, "_dk")
  if (spec == "baseline") {
    if (uses_dummy) {
      list(
        endo = c(tau, tau_dk, "dlag_ln_tau", "dlag_ln_tau_dk"),
        exog = c("dlag_ln_x", "dlag_ln_x_dk"),
        inst = c("d0_ln_tau_instr_25", "d0_ln_tau_instr_25_dk",
                 "dlag_d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25_dk"),
        targets = c(tau, tau_dk)
      )
    } else {
      list(
        endo = c(tau, "dlag_ln_tau"),
        exog = "dlag_ln_x",
        inst = c("d0_ln_tau_instr_25", "dlag_d0_ln_tau_instr_25"),
        targets = tau
      )
    }
  } else {
    if (uses_dummy) {
      list(
        endo = c(tau, tau_dk),
        exog = c("dlag_ln_x", "dlag_ln_x_dk", "dlag_ln_tau", "dlag_ln_tau_dk"),
        inst = c("d0_ln_tau", "d0_ln_tau_dk"),
        targets = c(tau, tau_dk)
      )
    } else {
      list(
        endo = tau,
        exog = c("dlag_ln_x", "dlag_ln_tau"),
        inst = "d0_ln_tau",
        targets = tau
      )
    }
  }
}

if (!"eq26_drop" %in% names(dt)) {
  stop("eq26_drop was not built -- Section 4 did not complete. Re-run from the top.")
}
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
# 7) Figure 1
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
# 8) Figure 2
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
