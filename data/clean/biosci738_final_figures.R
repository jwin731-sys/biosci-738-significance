## --- Section 0: packages + constants
# fig.width=12, fig.height=8
# Use cache=TRUE on the R Markdown chunk that contains the API calls below.

library(httr)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(lubridate)

bioblitz_date <- as.Date("2026-03-24")
thorpe_death <- as.Date("2024-08-01")
ua_palette <- c(
  "UoA Campus" = "#1A5276",
  "Dye Creek" = "#27AE60",
  "ArkState Querétaro" = "#8E44AD"
)
normal_colour <- "#7F8C8D"
user_agent_string <- "BIOSCI738-significance-article (R; educational use)"
api_base <- "https://api.inaturalist.org/v1"

## --- Section 1: iNat API helpers
# API helper for project-level observer rankings.
get_project_observers <- function(project_id, max_results = 500) {
  empty_df <- data.frame(
    rank = integer(),
    user_login = character(),
    obs_count = integer(),
    stringsAsFactors = FALSE
  )

  resp <- tryCatch(
    {
      GET(
        url = paste0(api_base, "/observations/observers"),
        query = list(project_id = project_id, per_page = max_results),
        timeout(30),
        user_agent(user_agent_string)
      )
    },
    error = function(e) {
      warning(sprintf("GET observers failed for project %s: %s", project_id, e$message))
      NULL
    }
  )
  Sys.sleep(0.5)

  if (is.null(resp)) {
    return(empty_df)
  }
  if (status_code(resp) != 200) {
    warning(sprintf("Observers endpoint returned status %s for project %s", status_code(resp), project_id))
    return(empty_df)
  }

  parsed <- tryCatch(
    {
      fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
    },
    error = function(e) {
      warning(sprintf("JSON parse failed for observers project %s: %s", project_id, e$message))
      NULL
    }
  )

  if (is.null(parsed) || is.null(parsed$results) || length(parsed$results) == 0) {
    return(empty_df)
  }

  rows <- lapply(parsed$results, function(x) {
    data.frame(
      user_login = if (!is.null(x$user$login)) as.character(x$user$login) else NA_character_,
      obs_count = suppressWarnings(as.integer(if (!is.null(x$observation_count)) x$observation_count else NA_integer_)),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(rows)
  if (nrow(out) == 0) {
    return(empty_df)
  }

  out <- out %>%
    filter(!is.na(user_login), !is.na(obs_count)) %>%
    arrange(desc(obs_count), user_login) %>%
    mutate(rank = row_number()) %>%
    select(rank, user_login, obs_count)

  if (nrow(out) == 0) {
    return(empty_df)
  }

  out
}

# API helper for annual counts of a focal user within a project.
get_user_annual_obs <- function(project_id, user_login) {
  empty_df <- data.frame(
    year = integer(),
    n_obs = integer(),
    stringsAsFactors = FALSE
  )

  per_page <- 200
  max_pages <- 50
  all_years <- vector("list", max_pages)
  stored <- 0

  for (p in seq_len(max_pages)) {
    resp <- tryCatch(
      {
        GET(
          url = paste0(api_base, "/observations"),
          query = list(project_id = project_id, user_id = user_login, per_page = per_page, page = p),
          timeout(30),
          user_agent(user_agent_string)
        )
      },
      error = function(e) {
        warning(sprintf("GET observations failed for project %s, user %s, page %s: %s", project_id, user_login, p, e$message))
        NULL
      }
    )
    Sys.sleep(0.5)

    if (is.null(resp)) {
      break
    }
    if (status_code(resp) != 200) {
      warning(sprintf("Observations endpoint returned status %s for project %s, user %s, page %s", status_code(resp), project_id, user_login, p))
      break
    }

    parsed <- tryCatch(
      {
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
      },
      error = function(e) {
        warning(sprintf("JSON parse failed for project %s, user %s, page %s: %s", project_id, user_login, p, e$message))
        NULL
      }
    )

    if (is.null(parsed) || is.null(parsed$results) || length(parsed$results) == 0) {
      break
    }

    this_years <- vapply(parsed$results, function(x) {
      observed_on_val <- if (!is.null(x$observed_on)) x$observed_on else NA_character_
      year_from_date <- suppressWarnings(as.integer(substr(observed_on_val, 1, 4)))
      if (!is.na(year_from_date)) {
        return(year_from_date)
      }
      year_from_details <- if (!is.null(x$observed_on_details$year)) suppressWarnings(as.integer(x$observed_on_details$year)) else NA_integer_
      year_from_details
    }, integer(1))

    stored <- stored + 1
    all_years[[stored]] <- data.frame(year = this_years, stringsAsFactors = FALSE)

    if (length(parsed$results) < per_page) {
      break
    }
  }

  if (stored == 0) {
    return(empty_df)
  }

  out <- bind_rows(all_years[seq_len(stored)])
  if (nrow(out) == 0) {
    return(empty_df)
  }

  out <- out %>%
    filter(!is.na(year)) %>%
    count(year, name = "n_obs") %>%
    arrange(year)

  if (nrow(out) == 0) {
    return(empty_df)
  }

  out
}

# Safe power-law fit on ranks 2..min(20, max_rank).
fit_power_law <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !all(c("rank", "obs_count") %in% names(df))) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }

  usable <- df %>%
    filter(rank >= 2, rank <= min(20L, max(rank, na.rm = TRUE)), obs_count > 0) %>%
    filter(is.finite(rank), is.finite(obs_count))

  if (nrow(usable) < 2) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }

  fit <- tryCatch(
    lm(log(obs_count) ~ log(rank), data = usable),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }

  coefs <- coef(fit)
  list(
    slope = unname(coefs[[2]]),
    intercept = unname(coefs[[1]])
  )
}

## --- Section 2: UoA data from CSV
# Load the flat observation export and derive all UoA-specific summaries from CSV.
obs_csv <- read.csv("data/clean/observations_flat.csv", stringsAsFactors = FALSE)
obs_csv$observed_on_date <- as.Date(obs_csv$observed_on, format = "%Y-%m-%d")

observer_ranking <- if (nrow(obs_csv) == 0) {
  data.frame(
    user_login = character(),
    obs_count = integer(),
    total_obs = numeric(),
    stringsAsFactors = FALSE
  )
} else {
  total_uoa_obs <- nrow(obs_csv)
  obs_csv %>%
    filter(!is.na(user_login), nzchar(user_login)) %>%
    count(user_login, name = "obs_count") %>%
    arrange(desc(obs_count), user_login) %>%
    mutate(total_obs = total_uoa_obs)
}

bioblitz_df <- if (nrow(obs_csv) == 0) {
  obs_csv[0, , drop = FALSE]
} else {
  obs_csv %>% filter(!is.na(observed_on_date) & observed_on_date == bioblitz_date)
}

thorpe_annual <- if (nrow(obs_csv) == 0) {
  data.frame(year = integer(), n_obs = integer(), stringsAsFactors = FALSE)
} else {
  obs_csv %>%
    filter(user_login == "stephen_thorpe") %>%
    mutate(year = suppressWarnings(as.integer(format(observed_on_date, "%Y")))) %>%
    filter(!is.na(year)) %>%
    count(year, name = "n_obs") %>%
    arrange(year)
}

first_taxon_dates <- if (nrow(obs_csv) == 0) {
  data.frame(taxon_name = character(), first_observed_on = as.Date(character()), first_on_bioblitz = logical(), stringsAsFactors = FALSE)
} else {
  obs_csv %>%
    filter(!is.na(taxon_name), nzchar(taxon_name), !is.na(observed_on_date)) %>%
    group_by(taxon_name) %>%
    summarise(first_observed_on = min(observed_on_date), .groups = "drop") %>%
    mutate(first_on_bioblitz = first_observed_on == bioblitz_date)
}
new_species_count <- if (nrow(first_taxon_dates) == 0) 0L else sum(first_taxon_dates$first_on_bioblitz, na.rm = TRUE)
message(sprintf("New species first recorded on bioblitz date: %s", comma(new_species_count)))

## --- Section 3: Cross-project pull
# Pull or derive cross-project observer rankings and focal annual counts.
projects <- list(
  "UoA Campus"         = list(id = 24826, r1 = "stephen_thorpe",    type = "thorpe"),
  "Dye Creek"          = list(id = 32783, r1 = "dyecreekscott",     type = "thorpe"),
  "ArkState Querétaro" = list(id = 27346, r1 = "shady_meerkat",     type = "thorpe"),
  "Georgia Tech"       = list(id = 25859, r1 = "featherenthusiast", type = "normal"),
  "West Valley"        = list(id = 17602, r1 = "merav",             type = "normal"),
  "Univ. Alabama"      = list(id = 6339,  r1 = "jdorris54",         type = "normal"),
  "Lane CC"            = list(id = 23714, r1 = "winterwren22",      type = "normal")
)

## --- BEGIN CACHEABLE API SECTION
all_observer_list <- vector("list", length(projects))
project_metrics_list <- vector("list", length(projects))
annual_list <- list()
project_names <- names(projects)

for (i in seq_along(projects)) {
  project_name <- project_names[[i]]
  project_info <- projects[[i]]

  if (project_name == "UoA Campus") {
    observers_df <- observer_ranking %>%
      mutate(rank = row_number()) %>%
      select(rank, user_login, obs_count, total_obs)
  } else {
    api_rank_df <- get_project_observers(project_info$id, max_results = 500)
    total_obs_project <- sum(api_rank_df$obs_count, na.rm = TRUE)
    observers_df <- if (nrow(api_rank_df) == 0) {
      data.frame(rank = integer(), user_login = character(), obs_count = integer(), total_obs = numeric(), stringsAsFactors = FALSE)
    } else {
      api_rank_df %>% mutate(total_obs = total_obs_project)
    }
  }

  observers_df <- observers_df %>%
    mutate(project = project_name, type = project_info$type) %>%
    select(project, type, rank, user_login, obs_count, total_obs)

  all_observer_list[[i]] <- observers_df

  total_obs_value <- if (nrow(observers_df) == 0) 0 else max(observers_df$total_obs, na.rm = TRUE)
  r1_row <- observers_df %>% filter(rank == 1)
  r2_row <- observers_df %>% filter(rank == 2)
  r1_n <- if (nrow(r1_row) == 0) NA_real_ else as.numeric(r1_row$obs_count[[1]])
  r2_n <- if (nrow(r2_row) == 0) NA_real_ else as.numeric(r2_row$obs_count[[1]])
  fit <- fit_power_law(observers_df)
  predicted_r1 <- if (is.na(fit$intercept)) NA_real_ else exp(fit$intercept)

  project_metrics_list[[i]] <- data.frame(
    project = project_name,
    type = project_info$type,
    r1_user = project_info$r1,
    r1_n = r1_n,
    r2_n = r2_n,
    total_obs = total_obs_value,
    r1r2_ratio = if (is.na(r1_n) || is.na(r2_n) || r2_n == 0) NA_real_ else r1_n / r2_n,
    pct_r1 = if (is.na(r1_n) || is.na(total_obs_value) || total_obs_value == 0) NA_real_ else r1_n / total_obs_value,
    power_law_slope = fit$slope,
    power_law_intercept = fit$intercept,
    predicted_r1 = predicted_r1,
    excess = if (is.na(r1_n) || is.na(predicted_r1) || predicted_r1 <= 0) NA_real_ else r1_n / predicted_r1,
    stringsAsFactors = FALSE
  )
}

annual_list[[1]] <- thorpe_annual %>% mutate(project = "UoA Campus", user_login = "stephen_thorpe")
annual_list[[2]] <- get_user_annual_obs(projects[["Dye Creek"]]$id, "dyecreekscott") %>% mutate(project = "Dye Creek", user_login = "dyecreekscott")
annual_list[[3]] <- get_user_annual_obs(projects[["ArkState Querétaro"]]$id, "shady_meerkat") %>% mutate(project = "ArkState Querétaro", user_login = "shady_meerkat")
## --- END CACHEABLE API SECTION

all_observers_df <- bind_rows(all_observer_list)
project_metrics <- bind_rows(project_metrics_list)
annual_df <- bind_rows(annual_list) %>%
  select(project, user_login, year, n_obs)

project_metrics <- project_metrics %>%
  mutate(
    project_short = c("UoA", "Dye", "ArkState", "GaTech", "West Valley", "Alabama", "Lane CC")[match(project, project_names)],
    plot_colour = if_else(type == "thorpe", ua_palette[project], normal_colour)
  )

all_observers_df <- all_observers_df %>%
  left_join(project_metrics %>% select(project, power_law_slope, power_law_intercept, predicted_r1, excess, pct_r1, r1r2_ratio, plot_colour), by = "project")

annual_df <- annual_df %>%
  mutate(
    project = factor(project, levels = c("UoA Campus", "Dye Creek", "ArkState Querétaro")),
    plot_colour = ua_palette[as.character(project)]
  )

## --- Section 4: blitz_5min
# Convert bioblitz timestamps to New Zealand time and summarise into 5-minute bins.
bioblitz_times <- bioblitz_df %>%
  mutate(time_utc = suppressWarnings(ymd_hms(time_observed_at, tz = "UTC", quiet = TRUE)))

dropped_time_rows <- if (nrow(bioblitz_times) == 0) 0L else sum(is.na(bioblitz_times$time_utc))
message(sprintf("Dropped %s bioblitz rows with unparseable time_observed_at.", comma(dropped_time_rows)))

blitz_5min <- if (nrow(bioblitz_times) == 0) {
  data.frame(bin = as.POSIXct(character(), tz = "Pacific/Auckland"), n = integer(), stringsAsFactors = FALSE)
} else {
  bioblitz_times %>%
    filter(!is.na(time_utc)) %>%
    mutate(
      time_nz = with_tz(time_utc, tzone = "Pacific/Auckland"),
      bin = floor_date(time_nz, unit = "5 minutes")
    ) %>%
    count(bin, name = "n") %>%
    arrange(bin)
}

## --- Section 5: Figure A
# fig.width=12, fig.height=8
# Multi-panel rank-abundance comparison across seven campus/nature projects.
figA_panels <- lapply(project_names, function(project_name) {
  df_proj <- all_observers_df %>% filter(project == project_name)
  metrics_proj <- project_metrics %>% filter(project == project_name)
  proj_type <- if (nrow(metrics_proj) == 0) "normal" else metrics_proj$type[[1]]
  panel_colour <- if (nrow(metrics_proj) == 0) normal_colour else metrics_proj$plot_colour[[1]]
  title_prefix <- if (identical(proj_type, "thorpe")) "★ " else ""
  title_text <- paste0(title_prefix, project_name)

  line_df <- if (nrow(df_proj) == 0 || nrow(metrics_proj) == 0 || is.na(metrics_proj$power_law_intercept[[1]]) || is.na(metrics_proj$power_law_slope[[1]])) {
    data.frame(rank = numeric(), obs_count = numeric())
  } else {
    data.frame(
      rank = seq(1, max(df_proj$rank, na.rm = TRUE)),
      obs_count = exp(metrics_proj$power_law_intercept[[1]]) * seq(1, max(df_proj$rank, na.rm = TRUE))^(metrics_proj$power_law_slope[[1]])
    )
  }

  pred_df <- if (nrow(metrics_proj) == 0 || is.na(metrics_proj$predicted_r1[[1]])) {
    data.frame(rank = numeric(), obs_count = numeric())
  } else {
    data.frame(rank = 1, obs_count = metrics_proj$predicted_r1[[1]])
  }

  annot_text <- if (nrow(metrics_proj) == 0) {
    "R1/R2: NA\nExcess: NA\nR1: NA"
  } else {
    sprintf(
      "R1/R2: %s\nExcess: %s\nR1: %s",
      ifelse(is.na(metrics_proj$r1r2_ratio[[1]]), "NA", number(metrics_proj$r1r2_ratio[[1]], accuracy = 0.1, suffix = "x")),
      ifelse(is.na(metrics_proj$excess[[1]]), "NA", number(metrics_proj$excess[[1]], accuracy = 0.1, suffix = "x")),
      ifelse(is.na(metrics_proj$pct_r1[[1]]), "NA", percent(metrics_proj$pct_r1[[1]], accuracy = 0.1))
    )
  }

  ggplot(df_proj, aes(x = rank, y = obs_count)) +
    geom_point(color = panel_colour, size = 1, alpha = 0.5, na.rm = TRUE) +
    geom_point(
      data = df_proj %>% filter(rank == 1),
      shape = 21,
      stroke = 0.7,
      fill = panel_colour,
      color = "black",
      size = 3,
      na.rm = TRUE
    ) +
    geom_line(
      data = line_df,
      aes(x = rank, y = obs_count),
      inherit.aes = FALSE,
      linewidth = 0.4,
      linetype = "dashed",
      color = "black",
      na.rm = TRUE
    ) +
    geom_point(
      data = pred_df,
      aes(x = rank, y = obs_count),
      inherit.aes = FALSE,
      shape = 24,
      size = 2.5,
      fill = "white",
      color = "black",
      stroke = 0.7,
      na.rm = TRUE
    ) +
    scale_x_log10(labels = label_number()) +
    scale_y_log10(labels = label_number()) +
    labs(title = title_text, x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = annot_text, hjust = 1.05, vjust = 1.1, size = 3) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      plot.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 8)
    )
})

figA <- wrap_plots(figA_panels, ncol = 3) +
  plot_annotation(
    title = "Observer rank-abundance patterns across seven iNaturalist projects",
    theme = theme_bw()
  ) &
  labs(x = "Observer rank", y = "Observations")

print(figA)

## --- Section 6: Figure B
# fig.width=12, fig.height=8
# Annual observation histories for the three thorpe-like top observers.
figure_b_titles <- c(
  "UoA Campus" = "stephen_thorpe | UoA Campus",
  "Dye Creek" = "dyecreekscott | Dye Creek",
  "ArkState Querétaro" = "shady_meerkat | ArkState Querétaro"
)
last_obs_years <- c(
  "stephen_thorpe" = if (nrow(thorpe_annual) == 0) NA_integer_ else max(thorpe_annual$year, na.rm = TRUE),
  "dyecreekscott" = NA_integer_,
  "shady_meerkat" = if (nrow(annual_df %>% filter(user_login == "shady_meerkat")) == 0) NA_integer_ else max((annual_df %>% filter(user_login == "shady_meerkat"))$year, na.rm = TRUE)
)

figB_panels <- lapply(c("UoA Campus", "Dye Creek", "ArkState Querétaro"), function(project_name) {
  df_proj <- annual_df %>% filter(as.character(project) == project_name)
  focal_login <- if (nrow(df_proj) == 0) projects[[project_name]]$r1 else df_proj$user_login[[1]]
  panel_colour <- ua_palette[[project_name]]
  total_label <- if (nrow(df_proj) == 0) "Total: 0" else sprintf("Total: %s", comma(sum(df_proj$n_obs, na.rm = TRUE)))
  last_year <- last_obs_years[[focal_login]]
  ymax <- if (nrow(df_proj) == 0) 1 else max(df_proj$n_obs, na.rm = TRUE)

  p <- ggplot(df_proj, aes(x = year, y = n_obs)) +
    geom_col(fill = panel_colour, width = 0.8, na.rm = TRUE) +
    labs(title = figure_b_titles[[project_name]], x = "Year", y = "Observations") +
    annotate("text", x = Inf, y = Inf, label = total_label, hjust = 1.05, vjust = 1.3, size = 3.2) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 11, face = "bold")
    )

  if (focal_login %in% c("stephen_thorpe", "shady_meerkat") && !is.na(last_year)) {
    p <- p +
      geom_vline(xintercept = last_year, color = "red3", linetype = "dashed", linewidth = 0.5) +
      annotate("text", x = last_year, y = ymax, label = "Last obs", color = "red3", vjust = -0.4, hjust = -0.05, size = 3)
  }

  p
})

figB <- wrap_plots(figB_panels, ncol = 3) +
  plot_annotation(
    title = "Annual contribution of devoted single observers within their bounded-site projects",
    theme = theme_bw()
  )

print(figB)

## --- Section 7: Figure C
# Comparative summary metrics across all seven projects.
metrics_long <- project_metrics %>%
  select(project, project_short, type, plot_colour, r1r2_ratio, excess, pct_r1) %>%
  pivot_longer(cols = c(r1r2_ratio, excess, pct_r1), names_to = "metric", values_to = "value")

metric_labels <- c(
  r1r2_ratio = "R1/R2 ratio",
  excess = "Observed vs predicted rank-1 excess",
  pct_r1 = "Share of project observations from rank-1 observer"
)

make_metric_plot <- function(metric_name) {
  df_metric <- metrics_long %>% filter(metric == metric_name)
  value_lab <- if (metric_name == "pct_r1") percent(df_metric$value, accuracy = 0.1) else number(df_metric$value, accuracy = 0.1)

  ggplot(df_metric, aes(x = factor(project_short, levels = project_metrics$project_short), y = value, fill = project)) +
    geom_col(width = 0.75, na.rm = TRUE) +
    geom_text(aes(label = value_lab), vjust = -0.3, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = setNames(project_metrics$plot_colour, project_metrics$project), guide = guide_legend(title = NULL)) +
    labs(title = metric_labels[[metric_name]], x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}

figC1 <- make_metric_plot("r1r2_ratio") + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4)
figC2 <- make_metric_plot("excess") + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4)
figC3 <- make_metric_plot("pct_r1") + scale_y_continuous(labels = percent_format(accuracy = 1))

figC <- (figC1 | figC2 | figC3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(figC)

## --- Section 8: Figure D
# fig.width=12, fig.height=8
# Species accumulation with and without Stephen Thorpe's observations.
# Each curve lives on its own observation-number axis (the two series have
# different lengths once Thorpe's 2,053 rows are removed). We therefore keep them
# in two separate long-format frames and overlay both lines directly. A row-wise
# join on observed_on_date would be wrong here: many observations share a date
# (up to ~790 on the bioblitz day), so joining on date fans out cartesian-style
# into hundreds of thousands of meaningless rows.
accum_all <- obs_csv %>%
  filter(!is.na(observed_on_date), !is.na(taxon_name), nzchar(taxon_name)) %>%
  arrange(observed_on_date, id) %>%
  mutate(
    observation_number = row_number(),
    cumulative_taxa = cumsum(!duplicated(taxon_name))
  ) %>%
  select(observation_number, observed_on_date, cumulative_taxa)

accum_no_thorpe <- obs_csv %>%
  filter(user_login != "stephen_thorpe", !is.na(observed_on_date), !is.na(taxon_name), nzchar(taxon_name)) %>%
  arrange(observed_on_date, id) %>%
  mutate(
    observation_number = row_number(),
    cumulative_taxa = cumsum(!duplicated(taxon_name))
  ) %>%
  select(observation_number, observed_on_date, cumulative_taxa)

# Stack the two curves into one long frame for a single ggplot call.
accum_long <- bind_rows(
  accum_all %>% mutate(series = "All observers"),
  accum_no_thorpe %>% mutate(series = "Excluding stephen_thorpe")
)

# Detection shadow = (final taxa with everyone) minus (final taxa without Thorpe).
final_taxa_all <- if (nrow(accum_all) == 0) 0 else max(accum_all$cumulative_taxa, na.rm = TRUE)
final_taxa_no_thorpe <- if (nrow(accum_no_thorpe) == 0) 0 else max(accum_no_thorpe$cumulative_taxa, na.rm = TRUE)
thorpe_exclusive_taxa <- final_taxa_all - final_taxa_no_thorpe
message(sprintf("Thorpe-exclusive taxa: %s", comma(thorpe_exclusive_taxa)))

# Reference x-positions (in observation-number space of the all-observers curve).
thorpe_death_x <- if (nrow(accum_all) == 0 || !any(accum_all$observed_on_date <= thorpe_death)) {
  NA_real_
} else {
  max(accum_all$observation_number[accum_all$observed_on_date <= thorpe_death], na.rm = TRUE)
}
bioblitz_x <- if (nrow(accum_all) == 0 || !any(accum_all$observed_on_date <= bioblitz_date)) {
  NA_real_
} else {
  max(accum_all$observation_number[accum_all$observed_on_date <= bioblitz_date], na.rm = TRUE)
}

# Shade the gap between the curves. The ribbon is defined on the all-observers
# x-axis: at each all-observers observation number, look up how many taxa the
# no-Thorpe curve had reached by the same calendar date (step-function lookup,
# carried forward). This is a faithful, non-fanned-out gap band.
no_thorpe_by_date <- accum_no_thorpe %>%
  group_by(observed_on_date) %>%
  summarise(taxa_no_thorpe = max(cumulative_taxa), .groups = "drop") %>%
  arrange(observed_on_date)

ribbon_df <- accum_all %>%
  rename(taxa_all = cumulative_taxa) %>%
  mutate(
    taxa_no_thorpe = if (nrow(no_thorpe_by_date) == 0) {
      0
    } else {
      # carry-forward step lookup: most recent no-Thorpe cumulative count at or
      # before each all-observers observation date.
      idx <- findInterval(observed_on_date, no_thorpe_by_date$observed_on_date)
      ifelse(idx == 0, 0, no_thorpe_by_date$taxa_no_thorpe[pmax(idx, 1)])
    }
  ) %>%
  select(observation_number, taxa_all, taxa_no_thorpe)

figD <- ggplot() +
  geom_ribbon(
    data = ribbon_df,
    aes(x = observation_number, ymin = taxa_no_thorpe, ymax = taxa_all),
    fill = ua_palette[["UoA Campus"]],
    alpha = 0.15,
    na.rm = TRUE
  ) +
  geom_line(
    data = accum_long,
    aes(x = observation_number, y = cumulative_taxa, color = series),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("All observers" = ua_palette[["UoA Campus"]], "Excluding stephen_thorpe" = normal_colour),
    name = NULL
  ) +
  labs(
    title = "UoA Campus species accumulation with and without Stephen Thorpe",
    x = "Observation number",
    y = "Cumulative taxa"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

if (!is.na(thorpe_death_x)) {
  figD <- figD + geom_vline(xintercept = thorpe_death_x, linetype = "dashed", color = "red3", linewidth = 0.5)
}
if (!is.na(bioblitz_x)) {
  figD <- figD + geom_vline(xintercept = bioblitz_x, linetype = "dashed", color = "dodgerblue4", linewidth = 0.5)
}

if (nrow(accum_all) > 0) {
  figD <- figD +
    annotate(
      "text",
      x = max(accum_all$observation_number, na.rm = TRUE),
      y = final_taxa_all,
      label = sprintf("Thorpe-exclusive taxa: %s", comma(thorpe_exclusive_taxa)),
      hjust = 1.02,
      vjust = -0.3,
      size = 3.2
    )
}

print(figD)
