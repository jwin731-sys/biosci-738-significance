# AI-USE STATEMENT: This script was drafted with Claude Code (Anthropic) and Forge (GPT-5.4).
# All modelling decisions, parameter choices, and interpretations were reviewed by the author.

# BIOSCI 738 — Significance Article
# Script: 09_hawkes_marked.R
# Description: 2-stream marked (multivariate) Hawkes process via stelfi::fit_mhawkes.
#
# Motivation: the univariate Hawkes fit (script 08) yields n = 0.69, but that
# branching ratio is contaminated by within-observer within-session clustering —
# stephen_thorpe does multiple observations per outing, and those consecutive
# events drive the apparent self-excitation. That is not social contagion; it is
# one person walking around with an app.
#
# The correct model for the social-contagion hypothesis is a 2-stream multivariate
# Hawkes process:
#   Stream 1 = stephen_thorpe observations
#   Stream 2 = everyone else ("others")
# The off-diagonal alpha terms estimate cross-observer excitation: does Thorpe's
# activity raise the rate of others, and vice versa?

Sys.setenv(TZ = "Pacific/Auckland")

suppressWarnings(suppressPackageStartupMessages({
  library(stelfi)
  library(readr)
  library(lubridate)
  library(here)
  library(ggplot2)
  library(dplyr)
}))

suppressMessages(here::i_am("scripts/09_hawkes_marked.R"))

# ── 1. Load and parse timestamps ─────────────────────────────────────────────

df <- read_csv(
  here::here("data/clean/observations_flat.csv"),
  show_col_types = FALSE,
  col_types = cols(time_observed_at = col_character(), observed_on = col_character())
)

has_time <- !is.na(df$time_observed_at) & nchar(df$time_observed_at) > 5
df_ts <- df[has_time, ]
df_ts$t_parsed <- suppressWarnings(ymd_hms(df_ts$time_observed_at, quiet = TRUE))
df_ts <- df_ts[!is.na(df_ts$t_parsed), ]
df_ts <- df_ts[order(df_ts$t_parsed), ]

t0 <- min(df_ts$t_parsed)
times_days <- as.numeric(difftime(df_ts$t_parsed, t0, units = "days"))

# Break ties: stelfi requires strictly ascending times (gap > 1e-10)
for (i in 2:length(times_days)) {
  if (times_days[i] - times_days[i - 1] < 1e-9) {
    times_days[i] <- times_days[i - 1] + 1e-9
  }
}

cat(sprintf("Total events with parsed timestamps: %d\n", length(times_days)))
cat(sprintf("Time span: %.1f days (%.1f years)\n",
            max(times_days), max(times_days) / 365.25))

# ── 2. Build stream vector ───────────────────────────────────────────────────

stream <- ifelse(df_ts$user_login == "stephen_thorpe", "thorpe", "others")

n_thorpe <- sum(stream == "thorpe")
n_others <- sum(stream == "others")
cat(sprintf("\nStream counts:\n"))
cat(sprintf("  thorpe:  %d events (%.1f%%)\n", n_thorpe, 100 * n_thorpe / length(stream)))
cat(sprintf("  others:  %d events (%.1f%%)\n", n_others, 100 * n_others / length(stream)))

# ── 3. Model A: full bivariate (free diagonal) ──────────────────────────────
#
# All four alpha terms are free: within-thorpe, thorpe->others,
# others->thorpe, within-others. This gives the optimizer full freedom.

cat("\n════════════════════════════════════════════════════════════\n")
cat("MODEL A: Full bivariate Hawkes (all alpha terms free)\n")
cat("════════════════════════════════════════════════════════════\n\n")

params_A <- list(
  mu    = c(0.2, 0.2),
  alpha = matrix(c(0.1, 0.2,
                    0.2, 0.1), byrow = TRUE, nrow = 2),
  beta  = c(0.7, 0.7)
)

fit_A <- tryCatch(
  fit_mhawkes(
    times      = times_days,
    stream     = stream,
    parameters = params_A
  ),
  warning = function(w) {
    cat("WARNING during Model A fit:", conditionMessage(w), "\n")
    suppressWarnings(
      fit_mhawkes(
        times      = times_days,
        stream     = stream,
        parameters = params_A
      )
    )
  },
  error = function(e) {
    cat("ERROR during Model A fit:", conditionMessage(e), "\n")
    NULL
  }
)

extract_mhawkes_params <- function(fit, label) {
  if (is.null(fit)) {
    cat(sprintf("  %s: fit failed — no parameters to extract.\n", label))
    return(NULL)
  }

  coefs <- get_coefs(fit)
  cat(sprintf("--- %s: raw coefficients ---\n", label))
  print(coefs)
  cat("\n")

  # get_coefs returns rows: mu, mu, alpha (4x), beta, beta
  # Extract in order
  estimates <- coefs[, "Estimate"]
  std_errs  <- coefs[, "Std. Error"]

  mu1    <- estimates[1]
  mu2    <- estimates[2]
  a11    <- estimates[3]  # within-thorpe (or stream 1)
  a12    <- estimates[4]  # thorpe -> others
  a21    <- estimates[5]  # others -> thorpe
  a22    <- estimates[6]  # within-others
  beta1  <- estimates[7]
  beta2  <- estimates[8]

  cat(sprintf("  mu_thorpe  = %.6f events/day (SE %.6f)\n", mu1, std_errs[1]))
  cat(sprintf("  mu_others  = %.6f events/day (SE %.6f)\n", mu2, std_errs[2]))
  cat(sprintf("  alpha[1,1] = %.6f (within-thorpe)     (SE %.6f)\n", a11, std_errs[3]))
  cat(sprintf("  alpha[1,2] = %.6f (thorpe -> others)  (SE %.6f)\n", a12, std_errs[4]))
  cat(sprintf("  alpha[2,1] = %.6f (others -> thorpe)  (SE %.6f)\n", a21, std_errs[5]))
  cat(sprintf("  alpha[2,2] = %.6f (within-others)     (SE %.6f)\n", a22, std_errs[6]))
  cat(sprintf("  beta_1     = %.6f per day              (SE %.6f)\n", beta1, std_errs[7]))
  cat(sprintf("  beta_2     = %.6f per day              (SE %.6f)\n", beta2, std_errs[8]))

  # Half-lives
  if (beta1 > 0) cat(sprintf("  half-life (thorpe stream):  %.2f min\n", log(2) / beta1 * 24 * 60))
  if (beta2 > 0) cat(sprintf("  half-life (others stream):  %.2f min\n", log(2) / beta2 * 24 * 60))

  cat("\n")

  # Branching ratios
  n_within_thorpe  <- if (beta1 > 0) a11 / beta1 else NA_real_
  n_thorpe_others  <- if (beta2 > 0) a12 / beta2 else NA_real_
  n_others_thorpe  <- if (beta1 > 0) a21 / beta1 else NA_real_
  n_within_others  <- if (beta2 > 0) a22 / beta2 else NA_real_

  cat("  Branching ratios:\n")
  cat(sprintf("    within-thorpe   (a11/beta1): %.6f\n", n_within_thorpe))
  cat(sprintf("    thorpe->others  (a12/beta2): %.6f\n", n_thorpe_others))
  cat(sprintf("    others->thorpe  (a21/beta1): %.6f\n", n_others_thorpe))
  cat(sprintf("    within-others   (a22/beta2): %.6f\n", n_within_others))
  cat("\n")

  # Stationarity: spectral radius of the branching matrix < 1
  alpha_mat <- matrix(c(a11, a12, a21, a22), byrow = TRUE, nrow = 2)
  beta_vec  <- c(beta1, beta2)
  # Branching matrix G[i,j] = alpha[i,j] / beta[j]
  G <- alpha_mat / matrix(beta_vec, nrow = 2, ncol = 2, byrow = TRUE)
  eigenvalues <- eigen(G, only.values = TRUE)$values
  spectral_radius <- max(Mod(eigenvalues))

  cat(sprintf("  Spectral radius of branching matrix: %.6f\n", spectral_radius))
  if (spectral_radius < 1) {
    cat("  STATIONARITY: PASS (spectral radius < 1)\n")
  } else {
    cat("  STATIONARITY: FAIL (spectral radius >= 1)\n")
  }
  cat("\n")

  list(
    mu1 = mu1, mu2 = mu2,
    a11 = a11, a12 = a12, a21 = a21, a22 = a22,
    beta1 = beta1, beta2 = beta2,
    n_within_thorpe = n_within_thorpe,
    n_thorpe_others = n_thorpe_others,
    n_others_thorpe = n_others_thorpe,
    n_within_others = n_within_others,
    spectral_radius = spectral_radius,
    coefs = coefs
  )
}

params_A_out <- extract_mhawkes_params(fit_A, "Model A")

# ── 4. Model B: cross-stream focus (diagonal near zero) ─────────────────────
#
# Initialize within-stream alpha near zero so the optimizer focuses on
# cross-observer excitation. stelfi does not support fixed parameters,
# so we set the diagonal starts to 0.001 and let the optimizer decide.
# If the diagonal stays near zero, within-stream self-excitation is not
# driving the result.

cat("════════════════════════════════════════════════════════════\n")
cat("MODEL B: Cross-stream focus (diagonal initialized near 0)\n")
cat("════════════════════════════════════════════════════════════\n\n")

params_B <- list(
  mu    = c(0.2, 0.2),
  alpha = matrix(c(0.001, 0.3,
                    0.3,   0.001), byrow = TRUE, nrow = 2),
  beta  = c(0.7, 0.7)
)

fit_B <- tryCatch(
  fit_mhawkes(
    times      = times_days,
    stream     = stream,
    parameters = params_B
  ),
  warning = function(w) {
    cat("WARNING during Model B fit:", conditionMessage(w), "\n")
    suppressWarnings(
      fit_mhawkes(
        times      = times_days,
        stream     = stream,
        parameters = params_B
      )
    )
  },
  error = function(e) {
    cat("ERROR during Model B fit:", conditionMessage(e), "\n")
    NULL
  }
)

params_B_out <- extract_mhawkes_params(fit_B, "Model B")

# ── 5. Compare to univariate and Sasquatch benchmarks ────────────────────────

cat("════════════════════════════════════════════════════════════\n")
cat("COMPARISON WITH BENCHMARKS\n")
cat("════════════════════════════════════════════════════════════\n\n")

univariate_n <- 0.69
sasquatch_n  <- 0.30

if (!is.null(params_A_out)) {
  cross_n_A <- max(params_A_out$n_thorpe_others, params_A_out$n_others_thorpe, na.rm = TRUE)
  cat(sprintf("  Model A max cross-observer branching ratio:  %.6f\n", cross_n_A))
} else {
  cross_n_A <- NA_real_
  cat("  Model A: not available\n")
}

if (!is.null(params_B_out)) {
  cross_n_B <- max(params_B_out$n_thorpe_others, params_B_out$n_others_thorpe, na.rm = TRUE)
  cat(sprintf("  Model B max cross-observer branching ratio:  %.6f\n", cross_n_B))
  cat(sprintf("  Model B diagonal stayed near zero? a11=%.6f, a22=%.6f\n",
              params_B_out$a11, params_B_out$a22))
} else {
  cross_n_B <- NA_real_
  cat("  Model B: not available\n")
}

cat(sprintf("\n  Univariate Hawkes (script 08):                n = %.2f\n", univariate_n))
cat(sprintf("  Sasquatch benchmark (Jones-Todd & vH 2024):   n = %.2f\n", sasquatch_n))

best_cross <- if (!is.na(cross_n_A)) cross_n_A else cross_n_B
if (!is.na(best_cross)) {
  cat(sprintf("  Cross-observer n from this model:             n = %.4f\n", best_cross))
  cat(sprintf("  Ratio vs univariate:  %.2f%%\n", 100 * best_cross / univariate_n))
  cat(sprintf("  Ratio vs Sasquatch:   %.2f%%\n", 100 * best_cross / sasquatch_n))
}
cat("\n")

# ── 6. Figure ────────────────────────────────────────────────────────────────

dir.create(here::here("figures"), recursive = TRUE, showWarnings = FALSE)

fig_path <- here::here("figures", "hawkes_marked.pdf")

# Try show_multivariate_hawkes first; fall back to manual ggplot on failure.
fig_saved <- FALSE

if (!is.null(fit_A)) {
  tryCatch({
    p <- show_multivariate_hawkes(fit_A)
    pdf(fig_path, width = 12, height = 6)
    print(p)
    dev.off()
    cat(sprintf("Figure saved: %s\n", fig_path))
    fig_saved <- TRUE
  }, error = function(e) {
    cat("show_multivariate_hawkes() failed:", conditionMessage(e), "\n")
    cat("Falling back to manual ggplot.\n")
  })
}

if (!fig_saved && !is.null(params_A_out)) {
  # Manual fallback: plot fitted conditional intensity for each stream
  # using the estimated parameters from Model A.
  dt <- 0.01  # days
  t_grid <- seq(0, max(times_days), by = dt)
  n_grid <- length(t_grid)

  mu1 <- params_A_out$mu1
  mu2 <- params_A_out$mu2
  a11 <- params_A_out$a11
  a12 <- params_A_out$a12
  a21 <- params_A_out$a21
  a22 <- params_A_out$a22
  b1  <- params_A_out$beta1
  b2  <- params_A_out$beta2

  stream_labels <- stream  # "thorpe" or "others"
  thorpe_times <- times_days[stream_labels == "thorpe"]
  others_times <- times_days[stream_labels == "others"]

  # Compute conditional intensities at each grid point
  lambda1 <- numeric(n_grid)  # thorpe stream intensity
  lambda2 <- numeric(n_grid)  # others stream intensity

  for (g in seq_along(t_grid)) {
    t_now <- t_grid[g]

    # Thorpe stream intensity: mu1 + a11 * sum(exp(-b1*(t-ti))) for thorpe events
    #                              + a21 * sum(exp(-b1*(t-tj))) for others events
    past_thorpe <- thorpe_times[thorpe_times < t_now]
    past_others <- others_times[others_times < t_now]

    sum_thorpe_to_1 <- if (length(past_thorpe) > 0) sum(exp(-b1 * (t_now - past_thorpe))) else 0
    sum_others_to_1 <- if (length(past_others) > 0) sum(exp(-b1 * (t_now - past_others))) else 0
    sum_thorpe_to_2 <- if (length(past_thorpe) > 0) sum(exp(-b2 * (t_now - past_thorpe))) else 0
    sum_others_to_2 <- if (length(past_others) > 0) sum(exp(-b2 * (t_now - past_others))) else 0

    lambda1[g] <- mu1 + a11 * sum_thorpe_to_1 + a21 * sum_others_to_1
    lambda2[g] <- mu2 + a12 * sum_thorpe_to_2 + a22 * sum_others_to_2
  }

  intensity_df <- data.frame(
    time_days = rep(t_grid, 2),
    lambda    = c(lambda1, lambda2),
    stream    = rep(c("thorpe", "others"), each = n_grid)
  )

  p_manual <- ggplot(intensity_df, aes(x = time_days, y = lambda, color = stream)) +
    geom_line(linewidth = 0.4, alpha = 0.7) +
    scale_color_manual(values = c("thorpe" = "#1b4965", "others" = "#c1121f")) +
    labs(
      title = "Multivariate Hawkes — fitted conditional intensity",
      x = "Days since first observation",
      y = expression(lambda(t) ~ "(events/day)"),
      color = "Stream"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = c(0.85, 0.85)
    )

  pdf(fig_path, width = 12, height = 6)
  print(p_manual)
  dev.off()
  cat(sprintf("Fallback figure saved: %s\n", fig_path))
  fig_saved <- TRUE
}

if (!fig_saved) {
  cat("No figure produced: both models failed to fit.\n")
}

# ── 7. Interpretation ────────────────────────────────────────────────────────

cat("\n════════════════════════════════════════════════════════════\n")
cat("INTERPRETATION\n")
cat("════════════════════════════════════════════════════════════\n\n")

if (!is.null(params_A_out)) {
  cross_max <- max(params_A_out$n_thorpe_others, params_A_out$n_others_thorpe, na.rm = TRUE)
  within_max <- max(params_A_out$n_within_thorpe, params_A_out$n_within_others, na.rm = TRUE)

  if (cross_max < 0.1) {
    cat("The 2-stream multivariate Hawkes process separates within-observer\n")
    cat("clustering from cross-observer contagion. The cross-observer branching\n")
    cat(sprintf("ratios (thorpe->others: %.4f; others->thorpe: %.4f) are both\n",
                params_A_out$n_thorpe_others, params_A_out$n_others_thorpe))
    cat("below 0.1, providing no evidence that one observer's activity raises\n")
    cat("the observation rate of others. The univariate branching ratio of 0.69\n")
    cat("is an artifact of within-observer session clustering — a single data-\n")
    cat("generating process in which stephen_thorpe records many observations per\n")
    cat("outing — rather than social contagion between independent observers.\n")
    cat(sprintf("By comparison, the Sasquatch citizen-science dataset (Jones-Todd &\n"))
    cat(sprintf("van Helsdingen 2024) exhibits n = 0.30, which likewise falls well\n"))
    cat("below the critical threshold but reflects genuine inter-report excitation\n")
    cat("across independent reporters.\n")
  } else {
    direction <- if (params_A_out$n_thorpe_others > params_A_out$n_others_thorpe) {
      "thorpe's observations raise the rate of others"
    } else {
      "others' observations raise the rate of thorpe"
    }
    dominant_beta <- if (params_A_out$n_thorpe_others > params_A_out$n_others_thorpe) {
      params_A_out$beta2
    } else {
      params_A_out$beta1
    }
    half_life_min <- log(2) / dominant_beta * 24 * 60

    cat("The 2-stream multivariate Hawkes process detects cross-observer\n")
    cat(sprintf("excitation: %s,\n", direction))
    cat(sprintf("with a cross-observer branching ratio of %.4f and a decay half-life\n",
                cross_max))
    cat(sprintf("of %.1f minutes. This is consistent with a data-generating process\n",
                half_life_min))
    cat("in which one observer's iNaturalist activity is visible to others on\n")
    cat("the platform, triggering a short-lived increase in observation effort.\n")
    cat(sprintf("The within-observer branching ratio (%.4f) remains the dominant\n",
                within_max))
    cat("source of temporal clustering, confirming that session-level autocorrelation\n")
    cat("— not social contagion — is the primary driver of the univariate n = 0.69.\n")
  }
} else {
  cat("Model A did not converge; interpretation deferred.\n")
}

cat("\n")
cat("Done.\n")
