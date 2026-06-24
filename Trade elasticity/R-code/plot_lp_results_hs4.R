library(data.table)
library(ggplot2)

base_path <- "P:/leawin/Elasticiteter2503/Estimating LP Boehm et al"
res_dir <- file.path(base_path, "Merge_all_countries", "LP_results_HS4_paper")

if (!dir.exists(res_dir)) {
  stop("Results directory not found: ", res_dir)
}

latest_file <- function(pattern) {
  files <- list.files(res_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No files found for pattern: ", pattern)
  }
  files[which.max(file.info(files)$mtime)]
}

fig1_tariff_no <- latest_file("^fig1_tariff_.*_no_pretrend\\.csv$")
fig1_tariff_yes <- latest_file("^fig1_tariff_.*_with_pretrend\\.csv$")
fig1_trade_no <- latest_file("^fig1_trade_.*_no_pretrend\\.csv$")
fig1_trade_yes <- latest_file("^fig1_trade_.*_with_pretrend\\.csv$")
fig2_base <- latest_file("^fig2_elast_.*_baseline\\.csv$")
fig2_all <- latest_file("^fig2_elast_.*_all_data\\.csv$")

message("Using Figure 1 tariff (no pretrend): ", basename(fig1_tariff_no))
message("Using Figure 1 tariff (with pretrend): ", basename(fig1_tariff_yes))
message("Using Figure 1 trade (no pretrend): ", basename(fig1_trade_no))
message("Using Figure 1 trade (with pretrend): ", basename(fig1_trade_yes))
message("Using Figure 2 baseline: ", basename(fig2_base))
message("Using Figure 2 all-data: ", basename(fig2_all))

spec_label <- c(
  no_pretrend = "No pretrend controls",
  with_pretrend = "Pretrend controls"
)

read_fig1 <- function(file, panel_name, spec_name) {
  dt <- fread(file)
  req <- c("h", "beta", "ci_lo", "ci_hi")
  missing <- setdiff(req, names(dt))
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(file), ": ", paste(missing, collapse = ", "))
  }
  dt[, panel := panel_name]
  dt[, spec := spec_label[[spec_name]]]
  dt[]
}

fig1_dt <- rbindlist(list(
  read_fig1(fig1_tariff_no, "Tariffs", "no_pretrend"),
  read_fig1(fig1_tariff_yes, "Tariffs", "with_pretrend"),
  read_fig1(fig1_trade_no, "Trade", "no_pretrend"),
  read_fig1(fig1_trade_yes, "Trade", "with_pretrend")
), use.names = TRUE, fill = TRUE)

setorder(fig1_dt, panel, spec, h)

fig1_png <- file.path(res_dir, "figure1_hs4_paper_tariffs_trade.png")
fig1_csv <- file.path(res_dir, "figure1_hs4_paper_plotdata.csv")

fwrite(fig1_dt, fig1_csv)

p1 <- ggplot(fig1_dt, aes(x = h, y = beta, color = spec, linetype = spec, group = interaction(spec, panel))) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(
    data = fig1_dt[panel == "Trade"],
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.15,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("No pretrend controls" = "#D55E00", "Pretrend controls" = "#0072B2")) +
  scale_linetype_manual(values = c("No pretrend controls" = "dashed", "Pretrend controls" = "solid")) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  labs(
    x = "Horizon (years)",
    y = "Estimates",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = "white", color = "grey70")
  )

ggsave(fig1_png, p1, width = 11, height = 5.5, dpi = 300)

fig2_base_dt <- fread(fig2_base)
fig2_all_dt <- fread(fig2_all)

for (nm in c("h", "beta", "ci_lo", "ci_hi")) {
  if (!(nm %in% names(fig2_base_dt))) stop("Missing column ", nm, " in ", basename(fig2_base))
  if (!(nm %in% names(fig2_all_dt))) stop("Missing column ", nm, " in ", basename(fig2_all))
}

fig2_base_dt[, series := "Baseline"]
fig2_all_dt[, series := "All data / all tariffs 2SLS"]
fig2_dt <- rbindlist(list(fig2_base_dt, fig2_all_dt), use.names = TRUE, fill = TRUE)
setorder(fig2_dt, series, h)
fig2_dt[, label := sprintf("%.2f", beta)]
fig2_dt[, label_y := fifelse(series == "Baseline", beta - 0.08, beta + 0.08)]

fig2_png <- file.path(res_dir, "figure2_hs4_paper_trade_elasticity.png")
fig2_csv <- file.path(res_dir, "figure2_hs4_paper_plotdata.csv")

fwrite(fig2_dt, fig2_csv)

p2 <- ggplot(fig2_dt, aes(x = h, y = beta, color = series, shape = series, linetype = series, group = series)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5, fill = "white", stroke = 1) +
  geom_text(aes(y = label_y, label = label), size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("Baseline" = "#0072B2", "All data / all tariffs 2SLS" = "#D55E00")) +
  scale_shape_manual(values = c("Baseline" = 21, "All data / all tariffs 2SLS" = 24)) +
  scale_linetype_manual(values = c("Baseline" = "solid", "All data / all tariffs 2SLS" = "solid")) +
  labs(
    x = "Horizon (years)",
    y = expression("Estimates of " * epsilon^h),
    color = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.73, 0.90))

ggsave(fig2_png, p2, width = 10, height = 7, dpi = 300)

message("Saved Figure 1 plot data: ", fig1_csv)
message("Saved Figure 1 PNG: ", fig1_png)
message("Saved Figure 2 plot data: ", fig2_csv)
message("Saved Figure 2 PNG: ", fig2_png)
