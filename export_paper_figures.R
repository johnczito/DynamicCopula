# ==============================================================================
# Producer for tables and figures needed for Section 6.2, 7.1 and 7.2
#
# figures/
#   fig_margin_coverage.pdf   Fig, 6.2
#   fig_sharpness.pdf         Fig, 7.2
#   fig_ablation.pdf          Fig, 7.1+7.2
# tables/
#   tab_margin_coverage.tex   6.2
#   tab_crime.tex             7.1
#   tab_macro.tex             7.2
#   tab_dm.tex                7.1+7.2
#   tab_ablation.tex          7.1+7.2
#
# Each .tex is a complete tabular, so the manuscript wraps it:
#   \begin{table}[!htbp]\centering\small
#   \input{tables/tab_macro}
#   \caption{...}\label{tab:macro}
#   \end{table}
#
# Structure:
#   1. Load the three run outputs plus the two ablation CSVs and, if present,
#      the per-origin crime CRPS cache
#   2. Figures first, unchanged from the previous version of this file
#   3. Tables. Crime and macro reproduce the published two-up block layout:
#      CRIME pairs split-half (series 1 with 6, 2 with 7, ...), as the crime
#      pipeline's cbind does; MACRO pairs adjacent (1 with 2, 3 with 4, ...),
#      as the macro pipeline's rbind of cbinds does. Do not swap them
#   4. Bolding follows the published captions, such as smallest MAE/MSE, SIZE 
#      and CRPS per series. Coverage is never bolded in the paper's own tables
#   5. Significance reuses nw_t() as results_summary.Rmd defines it
#
# Supersedes export_comparison_tables.R, which did items 4 and 5 only
#
# Run from the repo root
# ==============================================================================

SIM_FILE    <- "margincmp_R100_H10_T200_s8675309.rds"
MACRO_FILE  <- "macro_k14_s8675309_stbl/macro_metrics_k14_s8675309_stbl.rds"
CRIME_LAB   <- "crime_k10_s8675309/crime_fcasts_labeled_k10_s8675309.csv"
CRIME_CACHE <- "crime_k10_s8675309/crime_crps_per_origin_k10_s8675309.rds"

pick <- function(...) { p <- c(...); p <- p[file.exists(p)]; if (length(p)) p[1] else NA }
CRIME_ABL <- pick("crime_k10_s8675309/crime_ablation_k10_s8675309.csv",
                  "crime_ablation_k10_s8675309.csv")
MACRO_ABL <- pick("macro_k14_s8675309_stbl/macro_ablation_k14_s8675309_stbl.csv",
                  "macro_ablation_k14_s8675309_stbl.csv")

FIG_DIR <- "figures"; dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR <- "tables";  dir.create(TAB_DIR, showWarnings = FALSE)

W <- 6.5   # manuscript text width at 1in margins on letter

# ------------------------------------------------------------------------------
# shared look, matching results_summary.Rmd
# ------------------------------------------------------------------------------

COL  <- c(DGFC = "#0072B2", BVAR = "#E69F01",
          DLM  = "#56B4E9", `SV-BVAR` = "#CC79A7")
PAIR <- c("grey25", "grey70")
PCH  <- c(DGFC = 21, BVAR = 22, DLM = 24, `SV-BVAR` = 23)

ABL_COL <- c(COL[["DGFC"]], "grey70")
ABL_PCH <- c(21, 22)

MLAB <- c(gdpgrowth    = "GDP growth",      pcegrowth    = "PCE growth",
          bfigrowth    = "BFI growth",      resinvgrowth = "Resid. inv.",
          ipgrowth     = "IP growth",       capu         = "Capacity util.",
          employ       = "Employment",      hours        = "Hours",
          ur           = "Unemployment",    gdppigrowth  = "GDP PI infl.",
          pcepigrowth  = "PCE infl.",       ffr          = "Fed funds",
          spread       = "Yield spread",    realSPgrowth = "Real S&P 500")

# the longer names the manuscript uses in its table block headers
MFULL <- c(gdpgrowth    = "Real GDP growth",
           pcegrowth    = "PCE growth",
           bfigrowth    = "BFI growth",
           resinvgrowth = "Residential investment growth",
           ipgrowth     = "Industrial production growth",
           capu         = "Capacity utilization",
           employ       = "Payroll employment growth",
           hours        = "Growth in aggregate hours",
           ur           = "Unemployment rate",
           gdppigrowth  = "GDP PI inflation",
           pcepigrowth  = "PCE inflation",
           ffr          = "Federal funds rate",
           spread       = "Yield spread",
           realSPgrowth = "Real S\\&P 500 growth")

CFULL <- c(Murder                 = "Murder",
           Attempted              = "Attempted murder",
           Accessory              = "Accessory to murder",
           Manslaughter           = "Manslaughter",
           `Non-domestic assault` = "Non-domestic assault",
           `Officer assault`      = "Assaulting an officer",
           Abduction              = "Abduction",
           Blackmail              = "Blackmail",
           `Escape custody`       = "Escaping custody",
           `Resist officer`       = "Resisting an officer")

# Newey-West t, from results_summary.Rmd. Positive => DGFC is WORSE.
# Do not add an HLN correction or a t reference; the draft's counts came from
# this exact function at the 1.96 normal threshold
nw_t <- function(dd) {
  dd <- dd[is.finite(dd)]
  Tn <- length(dd)
  if (Tn < 10 || sd(dd) == 0) return(NA_real_)
  l <- floor(4 * (Tn / 100)^(2/9))
  g <- acf(dd, lag.max = l, type = "covariance", plot = FALSE)$acf
  v <- g[1] + 2 * sum((1 - seq_len(l) / (l + 1)) * g[-1])
  mean(dd) / sqrt(max(v, .Machine$double.eps) / Tn)
}
CRIT <- 1.96

wr <- function(lines, f, ncol = NULL) {
  if (!is.null(ncol)) check_cells(lines, ncol, f)
  writeLines(lines, file.path(TAB_DIR, f))
  cat("wrote ", file.path(TAB_DIR, f), "\n", sep = "")
}

# MLAB is written for plot axes, where & is an ordinary character. The one
# label containing it needs escaping before it reaches a .tex, or
# "Real S&P 500" opens a fifth column in a four-column table
# Spelled out rather than done with gsub
MTAB <- MLAB
MTAB[["realSPgrowth"]] <- "Real S\\&P 500"

# Catches this whole class of bug before a bad tabular reaches disk. Counts
# UNESCAPED ampersands per emitted row; only the pattern is passed to gsub, and
# a fixed pattern has no escaping ambiguity
check_cells <- function(lines, ncol, f) {
  body <- grep("\\\\\\\\$", lines, value = TRUE)
  body <- body[!grepl("multicolumn", body)]
  bare <- vapply(body, function(l)
    nchar(gsub("[^&]", "", gsub("\\&", "", l, fixed = TRUE))), 0L)
  bad <- which(bare != ncol - 1L)
  if (length(bad))
    stop(f, ": ", length(bad), " row(s) have the wrong cell count.\n  first: ",
         body[bad[1]], "\n  expected ", ncol - 1L, " separators, found ",
         bare[bad[1]], ". An unescaped & in a label is the usual cause.")
  
  stray <- grep("\\\\&", body, fixed = TRUE, value = TRUE)
  if (length(stray))
    stop(f, ": a row contains \\\\& , a line break before an ampersand.\n  first: ",
         stray[1], "\n  A label was escaped twice.")
  
  invisible(TRUE)
}

# ==============================================================================
# 1. Figure, Section 6.2: margin adjustment against a frozen ECDF
# ==============================================================================

SIM  <- readRDS(SIM_FILE)
fc   <- SIM$fc
DGPS <- c("varma", "varch", "varma_copula")
DLAB <- c("VARMA", "VARCH", "VARMA Copula")
TS   <- sort(unique(fc$T))

f1 <- aggregate(cov ~ dgp + T + arm, data = fc[fc$h == 1, ], FUN = mean)
YL <- range(f1$cov, 0.95)
YL <- c(floor(YL[1] * 50) / 50, ceiling(YL[2] * 50) / 50)

cov_panel <- function(x, ys, cols, pchs, ylab, main, ylim) {
  plot(x, ys[[1]], type = "n", log = "x", ylim = ylim, xaxt = "n",
       xlab = "sample size T", ylab = ylab, main = main)
  axis(1, at = TS, labels = TS)
  abline(h = 0.95, lty = 2, col = "grey40")
  for (i in seq_along(ys)) {
    lines(x, ys[[i]], col = cols[i], lwd = 2)
    points(x, ys[[i]], pch = pchs[i], bg = cols[i], col = "white",
           lwd = 0.6, cex = 1.3)
  }
}

pdf(file.path(FIG_DIR, "fig_margin_coverage.pdf"),
    width = W, height = 2.6, pointsize = 10)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.4, 0.8))
for (i in seq_along(DGPS)) {
  z  <- f1[f1$dgp == DGPS[i], ]
  zm <- z[z$arm == "ma",   ]; zm <- zm[order(zm$T), ]
  ze <- z[z$arm == "ecdf", ]; ze <- ze[order(ze$T), ]
  cov_panel(zm$T, list(zm$cov, ze$cov), ABL_COL, ABL_PCH,
            if (i == 1) "forecast interval coverage" else "", DLAB[i], YL)
  if (i == 1)
    legend("bottomright", c("margin adjustment", "frozen ECDF"),
           col = "white", pt.bg = ABL_COL, pch = ABL_PCH, pt.lwd = 0.6,
           bty = "n", cex = 0.85)
}
dev.off()
cat("wrote ", file.path(FIG_DIR, "fig_margin_coverage.pdf"), "\n", sep = "")

# ==============================================================================
# 2. Table, Section 6.2: pooled coverage and the exact McNemar test
#
# Paired because both arms come off one posterior fit, so the coverage
# indicators pair at the level of the individual forecast.
# ==============================================================================

key <- c("dgp", "rep", "T", "var", "h")
w <- merge(fc[fc$arm == "ma",   c(key, "cov")],
           fc[fc$arm == "ecdf", c(key, "cov")],
           by = key, suffixes = c(".ma", ".ec"))

pool <- aggregate(cov ~ T + arm, data = fc[fc$h == 1, ], FUN = mean)
pw   <- reshape(pool, idvar = "T", timevar = "arm", direction = "wide")
pw   <- pw[order(pw$T), ]

mcn <- function(T0) {
  z   <- w[w$h == 1 & w$T == T0, ]
  b01 <- sum(z$cov.ma == 1 & z$cov.ec == 0)
  b10 <- sum(z$cov.ma == 0 & z$cov.ec == 1)
  if (b01 + b10 == 0) return(NA_real_)
  binom.test(b01, b01 + b10, 0.5)$p.value
}

sgn <- function(x) sprintf(if (x < 0) "$-$%.3f" else "\\phantom{-}%.3f", abs(x))

wr(c("\\begin{tabular}{@{}rrrrrr@{}}",
     "\\toprule",
     "$T$ & frozen ECDF & margin adj. & ECDF shortfall & adj. shortfall & McNemar $p$ \\\\",
     "\\midrule",
     sprintf("%-3d & %.3f & %.3f & %s & %s & %.3f \\\\",
             pw$T, pw$cov.ecdf, pw$cov.ma,
             vapply(0.95 - pw$cov.ecdf, sgn, ""),
             vapply(0.95 - pw$cov.ma,   sgn, ""),
             vapply(pw$T, mcn, 0)),
     "\\bottomrule",
     "\\end{tabular}"), "tab_margin_coverage.tex", 6)

# ==============================================================================
# 3. Figure, Section 7.2: sharpness against calibration
# ==============================================================================

m <- readRDS(MACRO_FILE)
f <- m$finals

sz  <- cbind(f$dfc[, 3], f$var[, 3], f$dlm[, 3], f$svbvar[, 3])
cv  <- cbind(f$dfc[, 2], f$var[, 2], f$dlm[, 2], f$svbvar[, 2])
szn <- sz / sz[, 1]

pdf(file.path(FIG_DIR, "fig_sharpness.pdf"),
    width = W, height = 4.2, pointsize = 11)
par(mar = c(4.2, 4.2, 1.0, 7.5), xpd = NA)
plot(NA, xlim = range(szn) * c(0.8, 1.25), ylim = c(0.86, 1.01), log = "x",
     xlab = "interval width, relative to DGFC (log scale)",
     ylab = "empirical coverage", cex.lab = 0.9, cex.axis = 0.85)
par(xpd = FALSE)
abline(h = 0.95, lty = 2)
abline(v = 1, col = "grey85")
for (j in c(1, 3, 2, 4))
  points(szn[, j], cv[, j], pch = PCH[j], bg = COL[j],
         col = "white", lwd = 0.6, cex = 1.3)
par(xpd = NA)
legend(par("usr")[2] * 1.15, 1.0, names(COL), pch = PCH, pt.bg = COL,
       col = COL, bty = "o", cex = 0.85, xjust = -1, yjust = 2.5)
par(xpd = FALSE)
dev.off()
cat("wrote ", file.path(FIG_DIR, "fig_sharpness.pdf"), "\n", sep = "")

# ==============================================================================
# 4. Figure, Sections 7.1 and 7.2: ECDF ablation coverage
# ==============================================================================

ab_m <- data.frame(series = unname(MLAB[m$var_names]),
                   cov_ma = f$dfc[, 2], cov_ec = f$ecdf[, 2],
                   stringsAsFactors = FALSE)
ab_c <- read.csv(CRIME_ABL, stringsAsFactors = FALSE)

pdf(file.path(FIG_DIR, "fig_ablation.pdf"),
    width = W, height = 4.4, pointsize = 10)
par(mfrow = c(1, 2), mar = c(7.5, 4, 2.6, 0.8))
for (nm in c("Crime", "Macro")) {
  a <- if (nm == "Crime") ab_c else ab_m
  barplot(rbind(a$cov_ma, a$cov_ec), beside = TRUE, las = 2,
          names.arg = a$series, ylim = c(0.70, 1.0), xpd = FALSE,
          border = NA, col = ABL_COL, cex.names = 0.55,
          cex.axis = 0.8, ylab = "empirical coverage", cex.lab = 0.85,
          main = paste(nm, "coverage"), cex.main = 1)
  abline(h = 0.95, lty = 2)
  if (nm == "Crime")
    legend("topright", c("margin adjustment", "ECDF"), fill = ABL_COL,
           border = NA, bty = "n", cex = 0.75)
}
dev.off()
cat("wrote ", file.path(FIG_DIR, "fig_ablation.pdf"), "\n", sep = "")

# ==============================================================================
# 5. Tables 1 and 2: the published two-up block layout
# ==============================================================================

# bold the smallest, leave a column alone with bold = FALSE
fmt_col <- function(x, dg, bold = TRUE) {
  s <- formatC(round(x, dg), format = "f", digits = dg)
  if (bold) s[round(x, dg) == min(round(x, dg), na.rm = TRUE)] <-
      paste0("\\textbf{", s[round(x, dg) == min(round(x, dg), na.rm = TRUE)], "}")
  s
}

# one block of four model rows for one series, as a character matrix
block <- function(M, mods, dg, boldcol) {
  cells <- vapply(seq_len(ncol(M)),
                  function(j) fmt_col(M[, j], dg[j], boldcol[j]), character(nrow(M)))
  cbind(mods, cells)
}

two_up <- function(left, right, lab_l, lab_r, mods, dg, boldcol, colspec) {
  # left and right are LISTS of matrices, so colnames() on them is NULL and
  # the header comes out blank. Take the names off the first block.
  hdr <- paste(c("", colnames(left[[1]]), "", colnames(right[[1]])),
               collapse = " & ")
  # MFULL and CFULL are already LaTeX-safe as written; do not re-escape them.
  out <- c(sprintf("\\begin{tabular}{%s}", colspec), "\\toprule",
           paste0(hdr, " \\\\"), "\\midrule")
  for (b in seq_along(lab_l)) {
    L <- block(left[[b]],  mods, dg, boldcol)
    R <- block(right[[b]], mods, dg, boldcol)
    out <- c(out,
             sprintf("\\multicolumn{5}{c}{\\textit{%s}} & \\multicolumn{5}{c}{\\textit{%s}} \\\\",
                     lab_l[b], lab_r[b]),
             "\\cmidrule(r){1-5}\\cmidrule(l){6-10}",
             sprintf("%-9s & %s & %-9s & %s \\\\",
                     L[, 1], apply(L[, -1, drop = FALSE], 1, paste, collapse = " & "),
                     R[, 1], apply(R[, -1, drop = FALSE], 1, paste, collapse = " & ")),
             if (b < length(lab_l)) "\\addlinespace" else NULL)
  }
  c(out, "\\bottomrule", "\\end{tabular}")
}

# ---- Table 1, crime. Pairs are SPLIT-HALF: series i with series i+5. --------

if (file.exists(CRIME_LAB)) {
  cl <- read.csv(CRIME_LAB, check.names = FALSE, stringsAsFactors = FALSE)
  names(cl)[names(cl) == "INT-COV"]  <- "COV"
  names(cl)[names(cl) == "INT-SIZE"] <- "SIZE"
  cmods <- c("DGFC", "BVAR", "Pois-DGLM", "NB-DGLM")
  csrs  <- unique(cl$series)
  
  grab <- function(s) {
    b <- cl[cl$series == s, ]
    b <- b[match(cmods, b$model), c("MAE", "COV", "SIZE", "CRPS")]
    as.matrix(b)
  }
  half <- length(csrs) / 2
  wr(two_up(lapply(csrs[1:half],               grab),
            lapply(csrs[(half + 1):length(csrs)], grab),
            CFULL[csrs[1:half]], CFULL[csrs[(half + 1):length(csrs)]],
            cmods, c(2, 2, 2, 2), c(TRUE, FALSE, TRUE, TRUE),
            "@{}lrrrr@{\\qquad}lrrrr@{}"), "tab_crime.tex", 10)
} else {
  cat("skipped tab_crime: ", CRIME_LAB, " not found\n", sep = "")
}

# ---- Table 2, macro. Pairs are ADJACENT: series 1 with 2, 3 with 4, ... -----

mmods <- c("DGFC", "BVAR", "DLM", "SV-BVAR")
gm <- function(i) {
  M <- rbind(f$dfc[i, ], f$var[i, ], f$dlm[i, ], f$svbvar[i, ])[, c(1, 2, 3, 4)]
  colnames(M) <- c("MSE", "COV", "SIZE", "CRPS"); M
}
odd <- seq(1, length(m$var_names), by = 2)
wr(two_up(lapply(odd, gm), lapply(odd + 1, gm),
          MFULL[m$var_names[odd]], MFULL[m$var_names[odd + 1]],
          mmods, c(2, 3, 2, 3), c(TRUE, FALSE, TRUE, TRUE),
          "@{}lrrrr@{\\qquad}lrrrr@{}"), "tab_macro.tex", 10)

# ==============================================================================
# 6. Table 3: Newey-West t on the CRPS differential, both applications
# ==============================================================================

# braces are load-bearing: without them R ends the function at the first line
# and the dangling `else` is a syntax error
tfmt <- function(x) {
  if (is.na(x)) "--"
  else if (abs(x) <= CRIT) sprintf("\\textit{%.2f}$^{\\dagger}$", x)
  else sprintf("%.2f", x)
}

mcomp <- c("BVAR", "DLM", "SV-BVAR")
oth   <- list(BVAR = m$var[, 4, ], DLM = m$dlm[, 4, ], `SV-BVAR` = m$svbvar[, 4, ])
nsm   <- length(m$var_names)
tm <- vapply(mcomp, function(cm)
  vapply(seq_len(nsm), function(i) nw_t(m$dfc[, 4, i] - oth[[cm]][, i]), numeric(1)),
  numeric(nsm))

rows <- c("\\begin{tabular}{@{}lrrr@{}}", "\\toprule",
          "\\multicolumn{4}{@{}l}{\\textit{Macroeconomic aggregates}} \\\\",
          "series & vs BVAR & vs DLM & vs SV-BVAR \\\\", "\\midrule",
          sprintf("%-22s & %s & %s & %s \\\\", MTAB[m$var_names],
                  vapply(tm[, 1], tfmt, ""), vapply(tm[, 2], tfmt, ""),
                  vapply(tm[, 3], tfmt, "")))

tc <- NULL
if (file.exists(CRIME_CACHE)) {
  cc    <- readRDS(CRIME_CACHE)
  ccomp <- c("BVAR", "Pois-DGLM", "NB-DGLM")
  nsc   <- length(cc$labs)
  tc <- vapply(ccomp, function(cm)
    vapply(seq_len(nsc), function(i) nw_t(cc$crps[, "DGFC", i] - cc$crps[, cm, i]),
           numeric(1)), numeric(nsc))
  rows <- c(rows, "\\addlinespace", "\\midrule",
            "\\multicolumn{4}{@{}l}{\\textit{Crime counts}} \\\\",
            "series & vs BVAR & vs P-DGLM & vs NB-DGLM \\\\", "\\midrule",
            sprintf("%-22s & %s & %s & %s \\\\", CFULL[cc$labs],
                    vapply(tc[, 1], tfmt, ""), vapply(tc[, 2], tfmt, ""),
                    vapply(tc[, 3], tfmt, "")))
} else {
  cat("NOTE: ", CRIME_CACHE, " absent, so tab_dm.tex has the macro panel only.\n",
      "      Knit results_summary.Rmd once to build it.\n", sep = "")
}
wr(c(rows, "\\bottomrule", "\\end{tabular}"), "tab_dm.tex", 4)

# ==============================================================================
# 7. Table 4: ECDF ablation, pooled
# ==============================================================================

abm_full <- if (!is.na(MACRO_ABL)) read.csv(MACRO_ABL, stringsAsFactors = FALSE) else
  data.frame(cov_ma = f$dfc[, 2], cov_ec = f$ecdf[, 2],
             size_pct = 100 * (f$ecdf[, 3] / f$dfc[, 3] - 1),
             crps_pct = 100 * (f$ecdf[, 4] / f$dfc[, 4] - 1))

ab_row <- function(lab, a)
  sprintf("%-26s & %.3f & %.3f & %+.1f\\%% & %d of %d & %d of %d \\\\",
          lab, mean(a$cov_ma), mean(a$cov_ec), mean(a$size_pct),
          sum(a$size_pct < 0), nrow(a), sum(a$crps_pct > 0), nrow(a))

abc_full <- read.csv(CRIME_ABL, stringsAsFactors = FALSE)
wr(c("\\begin{tabular}{@{}lrrrrr@{}}", "\\toprule",
     paste("application & cov. (adj.) & cov. (ECDF) & width &",
           "ECDF narrower & CRPS prefers adj. \\\\"),
     "\\midrule",
     ab_row("Crime counts ($n = 10$)",      abc_full),
     ab_row("Macro aggregates ($n = 14$)",  abm_full),
     "\\bottomrule", "\\end{tabular}"), "tab_ablation.tex", 6)

# ==============================================================================
# every count the draft quotes
# ==============================================================================

cat("\n================ CHECK AGAINST THE DRAFT ================\n")
for (j in seq_along(mcomp))
  cat(sprintf("7.2 vs %-8s  DGFC better %2d, worse %2d, tied %2d  (of %d)\n",
              mcomp[j], sum(tm[, j] < -CRIT, na.rm = TRUE),
              sum(tm[, j] > CRIT, na.rm = TRUE),
              sum(abs(tm[, j]) <= CRIT, na.rm = TRUE), nsm))
if (!is.null(tc)) for (j in 1:3)
  cat(sprintf("7.1 vs %-8s  DGFC better %2d, worse %2d, tied %2d  (of %d)\n",
              c("BVAR", "P-DGLM", "NB-DGLM")[j],
              sum(tc[, j] < -CRIT, na.rm = TRUE), sum(tc[, j] > CRIT, na.rm = TRUE),
              sum(abs(tc[, j]) <= CRIT, na.rm = TRUE), nrow(tc)))

cat(sprintf("7.2 narrower than SV-BVAR on %d of %d; mean cov %.3f vs %.3f\n",
            sum(f$dfc[, 3] < f$svbvar[, 3]), nsm,
            mean(f$dfc[, 2]), mean(f$svbvar[, 2])))
cat(sprintf("7.2 MSE  vs SV-BVAR %d-%d,  vs DLM %d-%d\n",
            sum(f$dfc[, 1] < f$svbvar[, 1]), sum(f$svbvar[, 1] < f$dfc[, 1]),
            sum(f$dfc[, 1] < f$dlm[, 1]),    sum(f$dlm[, 1] < f$dfc[, 1])))
cat(sprintf("7.2 CRPS raw wins  vs BVAR %d of %d,  vs DLM %d of %d\n",
            sum(f$dfc[, 4] < f$var[, 4]), nsm, sum(f$dfc[, 4] < f$dlm[, 4]), nsm))
if (!is.null(m$svbvar_nstable))
  cat(sprintf("7.2 SV-BVAR stable draws  median %.0f, below half on %d of %d\n",
              median(m$svbvar_nstable), sum(m$svbvar_nstable < 2500),
              length(m$svbvar_nstable)))
for (nm in c("7.1 crime", "7.2 macro")) {
  a <- if (grepl("crime", nm)) abc_full else abm_full
  cat(sprintf("%s ablation  cov %.3f -> %.3f, width %+.1f%%, CRPS prefers adj. %d of %d\n",
              nm, mean(a$cov_ma), mean(a$cov_ec), mean(a$size_pct),
              sum(a$crps_pct > 0), nrow(a)))
}
if (!is.na(MACRO_ABL)) {
  lev <- abm_full$series %in% c("capu", "ur", "ffr", "spread")
  cat(sprintf("7.2 ablation  CRPS improves on %d of %d growth rates; capu %.0f%%, ur %.0f%%, ffr %.0f%%\n",
              sum(abm_full$crps_pct[!lev] < 0), sum(!lev),
              abm_full$crps_pct[abm_full$series == "capu"],
              abm_full$crps_pct[abm_full$series == "ur"],
              abm_full$crps_pct[abm_full$series == "ffr"]))
}

cat("\ndone\n")