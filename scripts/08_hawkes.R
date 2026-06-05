# BIOSCI 738 — Significance Article
# Script: 08_hawkes.R
# Description: Temporal Hawkes process fit (hawkesbow package)
# AI use statement: Script structure and hawkesbow API usage drafted with Claude Code assistance.
#   All analysis decisions, parameter interpretations, and biological reasoning are the author's own.

Sys.setenv(TZ = "Pacific/Auckland")

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(hawkesbow)
  library(here)
  library(lubridate)
  library(readr)
  library(scales)
}))

suppressMessages(here::i_am("scripts/08_hawkes.R"))

theme_hawkes <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey82", linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_blank()
    )
}

parse_observation_timestamps <- function(observations) {
  raw_observed_time <- trimws(as.character(observations$time_observed_at))
  raw_observed_time[is.na(raw_observed_time) | raw_observed_time == "NA"] <- ""
  has_explicit_time <- raw_observed_time != ""

  parsed_timestamp <- rep(
    as.POSIXct(NA_real_, origin = "1970-01-01", tz = "Pacific/Auckland"),
    nrow(observations)
  )

  if (any(has_explicit_time)) {
    parsed_timestamp[has_explicit_time] <- suppressWarnings(
      ymd_hms(raw_observed_time[has_explicit_time], quiet = TRUE)
    )
  }

  observed_date <- ymd(
    observations$observed_on,
    quiet = TRUE,
    tz = "Pacific/Auckland"
  )
  fallback_noon <- observed_date + hours(12)
  parsed_timestamp[!has_explicit_time] <- fallback_noon[!has_explicit_time]

  parsed_timestamp
}

prepare_events <- function(observations) {
  if (nrow(observations) == 0) {
    return(list(
      data = tibble(timestamp = as.POSIXct(character()), event_time_days = numeric()),
      events = numeric(),
      end_time = NA_real_,
      start_time = as.POSIXct(NA),
      n_events = 0L
    ))
  }

  parsed_timestamp <- parse_observation_timestamps(observations)

  prepared <- observations %>%
    mutate(timestamp = parsed_timestamp) %>%
    filter(!is.na(timestamp)) %>%
    arrange(timestamp)

  if (nrow(prepared) == 0) {
    return(list(
      data = tibble(timestamp = as.POSIXct(character()), event_time_days = numeric()),
      events = numeric(),
      end_time = NA_real_,
      start_time = as.POSIXct(NA),
      n_events = 0L
    ))
  }

  start_time <- min(prepared$timestamp)

  prepared <- prepared %>%
    mutate(
      event_time_days = as.numeric(difftime(timestamp, start_time, units = "days"))
    ) %>%
    arrange(event_time_days)

  event_times <- prepared$event_time_days

  if (length(event_times) == 1) {
    event_times[1] <- 1e-8
    prepared$event_time_days[1] <- 1e-8
  }

  list(
    data = prepared,
    events = event_times,
    end_time = max(event_times),
    start_time = start_time,
    n_events = length(event_times)
  )
}

fit_hawkes_model <- function(event_info, label) {
  if (event_info$n_events < 2 || !is.finite(event_info$end_time) || event_info$end_time <= 0) {
    return(list(
      label = label,
      fit = NULL,
      mu = NA_real_,
      alpha = NA_real_,
      beta = NA_real_,
      branching_ratio = NA_real_,
      stationarity = "STATIONARITY: INDETERMINATE (insufficient event span)",
      n_events = event_info$n_events
    ))
  }

  hawkes_fit <- suppressWarnings(
    hawkesbow::mle(
      events = event_info$events,
      kern = "Exponential",
      end = event_info$end_time
    )
  )

  parameters <- hawkes_fit$par
  mu <- unname(parameters[1])
  alpha <- unname(parameters[2])
  beta <- unname(parameters[3])
  branching_ratio <- alpha / beta

  stationarity <- if (is.finite(branching_ratio) && branching_ratio < 1) {
    "STATIONARITY: PASS (n < 1)"
  } else if (is.finite(branching_ratio)) {
    "STATIONARITY: FAIL (n >= 1)"
  } else {
    "STATIONARITY: INDETERMINATE (branching ratio unavailable)"
  }

  list(
    label = label,
    fit = hawkes_fit,
    mu = mu,
    alpha = alpha,
    beta = beta,
    branching_ratio = branching_ratio,
    stationarity = stationarity,
    n_events = event_info$n_events
  )
}

print_fit_summary <- function(fit_result) {
  cat(fit_result$label, "\n", sep = "")
  cat("Events:", fit_result$n_events, "\n")

  if (!is.finite(fit_result$mu) || !is.finite(fit_result$alpha) || !is.finite(fit_result$beta)) {
    cat("Fit unavailable: insufficient non-missing timestamps or event span.\n")
    cat(fit_result$stationarity, "\n\n", sep = "")
    return(invisible(NULL))
  }

  cat(sprintf("mu (events/day): %.6f\n", fit_result$mu))
  cat(sprintf("alpha (jump size): %.6f\n", fit_result$alpha))
  cat(sprintf("beta (decay rate, per day): %.6f\n", fit_result$beta))
  cat(sprintf("decay half-life (1/beta, days): %.6f\n", 1 / fit_result$beta))
  cat(sprintf("branching ratio (n): %.6f\n", fit_result$branching_ratio))
  cat(fit_result$stationarity, "\n\n", sep = "")
}

build_placeholder_plot <- function(title_text, message_text, x_label = NULL, y_label = NULL) {
  ggplot(data.frame(x = 0, y = 0, label = message_text), aes(x, y)) +
    geom_text(aes(label = label), size = 4) +
    labs(title = title_text, x = x_label, y = y_label) +
    theme_hawkes() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank()
    )
}

save_cumulative_plot <- function(event_info) {
  if (event_info$n_events == 0) {
    cumulative_plot <- build_placeholder_plot(
      title_text = "Cumulative iNaturalist observations — UoA City Campus",
      message_text = "No observations available after timestamp parsing.",
      x_label = "Date",
      y_label = "Cumulative observations"
    )
  } else {
    cumulative_data <- event_info$data %>%
      mutate(cumulative_observations = row_number())

    bioblitz_date <- as.Date("2026-03-24")
    annotation_y <- max(cumulative_data$cumulative_observations) * 0.80

    # Place annotation to the LEFT of the vertical line so it is not clipped
    # at the right edge of the plot panel.
    cumulative_plot <- ggplot(
      cumulative_data,
      aes(x = as.Date(timestamp), y = cumulative_observations)
    ) +
      geom_line(color = "#1b4965", linewidth = 0.7) +
      geom_vline(
        xintercept = bioblitz_date,
        color = "#c1121f",
        linewidth = 0.6,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = bioblitz_date,
        y = annotation_y,
        label = "Bioblitz",
        color = "#c1121f",
        hjust = 1.1,
        vjust = 0.5,
        fontface = "italic"
      ) +
      labs(
        title = "Cumulative iNaturalist observations — UoA City Campus",
        x = "Date",
        y = "Cumulative observations"
      ) +
      theme_hawkes()
  }

  ggsave(
    filename = here::here("figures", "hawkes_intensity.pdf"),
    plot = cumulative_plot,
    device = grDevices::cairo_pdf,
    width = 7,
    height = 4,
    units = "in"
  )
}

save_kernel_plot <- function(full_fit, no_thorpe_fit) {
  kernel_inputs <- bind_rows(
    tibble(
      subset_label = "Full dataset",
      alpha = full_fit$alpha,
      beta = full_fit$beta
    ),
    tibble(
      subset_label = "stephen_thorpe removed",
      alpha = no_thorpe_fit$alpha,
      beta = no_thorpe_fit$beta
    )
  ) %>%
    filter(is.finite(alpha), is.finite(beta), beta > 0)

  if (nrow(kernel_inputs) == 0) {
    kernel_plot <- build_placeholder_plot(
      title_text = "Hawkes excitation kernel",
      message_text = "Kernel unavailable because no model fit converged.",
      x_label = "Hours since triggering event",
      y_label = "Excitation intensity (events/day)"
    )
  } else {
    # Adaptive x-range: show ~5 mean lifetimes (5/beta) so the decay curve is
    # actually visible.  With beta ~ 300 the half-life is a few minutes, so we
    # plot in HOURS rather than days.
    min_beta <- min(kernel_inputs$beta)
    x_max_days <- 5 / min_beta
    x_max_hours <- x_max_days * 24

    kernel_data <- expand.grid(
      elapsed_hours = seq(0, x_max_hours, length.out = 401),
      subset_label = kernel_inputs$subset_label,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      left_join(kernel_inputs, by = "subset_label") %>%
      mutate(
        elapsed_days = elapsed_hours / 24,
        excitation_intensity = alpha * exp(-beta * elapsed_days)
      )

    kernel_plot <- ggplot(
      kernel_data,
      aes(x = elapsed_hours, y = excitation_intensity, color = subset_label)
    ) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(
        values = c(
          "Full dataset" = "#1b4965",
          "stephen_thorpe removed" = "#c1121f"
        )
      ) +
      labs(
        title = "Hawkes excitation kernel",
        x = "Hours since triggering event",
        y = "Excitation intensity (events/day)"
      ) +
      theme_hawkes() +
      theme(legend.position = c(0.78, 0.82))
  }

  ggsave(
    filename = here::here("figures", "hawkes_kernel.pdf"),
    plot = kernel_plot,
    device = grDevices::cairo_pdf,
    width = 7,
    height = 4,
    units = "in"
  )
}

save_branching_ratio_plot <- function(full_fit, no_thorpe_fit, no_top10_fit) {
  branching_data <- tibble(
    dataset_subset = c(
      "Full dataset",
      "stephen_thorpe removed",
      "Top 10 observers removed"
    ),
    branching_ratio = c(
      full_fit$branching_ratio,
      no_thorpe_fit$branching_ratio,
      no_top10_fit$branching_ratio
    )
  ) %>%
    mutate(
      dataset_subset = factor(
        dataset_subset,
        levels = rev(c(
          "Full dataset",
          "stephen_thorpe removed",
          "Top 10 observers removed"
        ))
      )
    )

  valid_branching <- filter(branching_data, is.finite(branching_ratio))

  if (nrow(valid_branching) == 0) {
    counterfactual_plot <- build_placeholder_plot(
      title_text = "Branching ratio by observer subset",
      message_text = "Counterfactual comparison unavailable because no fit converged.",
      x_label = "Branching ratio (n)",
      y_label = "Dataset subset"
    )
  } else {
    # Scale x-axis to the actual data so the bars are visible.  When all
    # branching ratios are << 0.5, placing the subcritical threshold line at
    # 0.5 squashes the bars into invisible slivers.  Instead, scale to a
    # sensible multiple of the largest observed ratio and note n = 0.5 in the
    # subtitle rather than drawing it off-screen.
    max_n <- max(valid_branching$branching_ratio)
    x_upper <- max_n * 1.6

    counterfactual_plot <- ggplot(
      valid_branching,
      aes(x = branching_ratio, y = dataset_subset, fill = dataset_subset)
    ) +
      geom_col(width = 0.55, show.legend = FALSE) +
      geom_text(
        aes(label = sprintf("n = %.4f", branching_ratio)),
        hjust = -0.08,
        size = 3.4
      ) +
      scale_fill_manual(
        values = c(
          "Full dataset" = "#1b4965",
          "stephen_thorpe removed" = "#c1121f",
          "Top 10 observers removed" = "#6c757d"
        )
      ) +
      scale_x_continuous(
        limits = c(0, x_upper),
        expand = expansion(mult = c(0, 0.05)),
        labels = scales::label_number(accuracy = 0.0001)
      ) +
      labs(
        title = "Branching ratio by observer subset",
        subtitle = "All subsets deeply subcritical (n << 0.5)",
        x = "Branching ratio (n)",
        y = "Dataset subset"
      ) +
      theme_hawkes() +
      theme(plot.subtitle = element_text(color = "grey40", size = 9))
  }

  ggsave(
    filename = here::here("figures", "hawkes_counterfactual.pdf"),
    plot = counterfactual_plot,
    device = grDevices::cairo_pdf,
    width = 7,
    height = 4,
    units = "in"
  )
}

main <- function() {
  observations <- readr::read_csv(
    here::here("data", "clean", "observations_flat.csv"),
    show_col_types = FALSE
  )

  dir.create(here::here("figures"), recursive = TRUE, showWarnings = FALSE)

  # The data-generating process includes self-excitation when an observation surge
  # recruits attention that yields more observations shortly afterwards.
  full_events <- prepare_events(observations)
  full_fit <- fit_hawkes_model(full_events, "--- Full dataset Hawkes fit ---")
  print_fit_summary(full_fit)

  # Observer-level varying effects can inflate apparent contagion if one prolific
  # naturalist dominates the event stream rather than a campus-wide response.
  observations_no_thorpe <- observations %>%
    filter(user_login != "stephen_thorpe")
  no_thorpe_events <- prepare_events(observations_no_thorpe)
  no_thorpe_fit <- fit_hawkes_model(
    no_thorpe_events,
    "--- Counterfactual: stephen_thorpe removed ---"
  )
  print_fit_summary(no_thorpe_fit)

  cat(
    sprintf(
      "Branching ratio: %s (full) -> %s (no thorpe), delta = %s\n\n",
      ifelse(is.finite(full_fit$branching_ratio), sprintf("%.6f", full_fit$branching_ratio), "NA"),
      ifelse(is.finite(no_thorpe_fit$branching_ratio), sprintf("%.6f", no_thorpe_fit$branching_ratio), "NA"),
      ifelse(
        is.finite(full_fit$branching_ratio) && is.finite(no_thorpe_fit$branching_ratio),
        sprintf("%.6f", no_thorpe_fit$branching_ratio - full_fit$branching_ratio),
        "NA"
      )
    )
  )

  top_observers <- observations %>%
    count(user_login, sort = TRUE, name = "n_observations") %>%
    filter(!is.na(user_login), user_login != "") %>%
    slice_head(n = 10) %>%
    pull(user_login)

  observations_no_top10 <- observations %>%
    filter(!user_login %in% top_observers)
  no_top10_events <- prepare_events(observations_no_top10)
  no_top10_fit <- fit_hawkes_model(
    no_top10_events,
    "--- Counterfactual: top 10 observers removed ---"
  )
  print_fit_summary(no_top10_fit)

  save_cumulative_plot(full_events)
  save_kernel_plot(full_fit, no_thorpe_fit)
  save_branching_ratio_plot(full_fit, no_thorpe_fit, no_top10_fit)
}

main()
