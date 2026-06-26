suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Plot Sampling Pipeline (OSM-stratified, LCZ-rebalanced)
# ============================================================
# This script:
# 1) loads study area + LCZ/OSM rasters,
# 2) samples candidate points per city,
# 3) selects points stratified by OSM class,
# 4) rebalances LCZ globally by swapping points within same city + OSM,
# 5) exports points, 100x100 m squares, and statistics tables.

utils::globalVariables(c(
  "pct",
  "base_count",
  "id",
  "selected_n",
  "target",
  "deficit",
  "lcz_deficit",
  "osm_deficit",
  "label",
  "available",
  "scarcity"
))

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "R.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) return(NULL)
    current <- parent
  }
}

project_root <- find_project_root()
if (is.null(project_root)) stop("Could not locate project root (R.Rproj).")
setwd(project_root)

# -----------------------------
# User configuration
# -----------------------------
study_area_path <- "F:/KULeuven/BRANCH/Materials_and_methods/Sampling_design/Study_Area/BRANCH_StudyArea_31370.gpkg"
selected_cities <- c("BRUSSEL", "LUIK", "LEUVEN", "HASSELT")

n_plots_per_city <- 20
plot_size <- 100
min_distance <- 150
n_candidates <- 10000

# OSM stratification targets.
# If NULL, the script splits n_plots_per_city evenly over allowed_osm_labels.
n_plots_per_osm_class <- 5L

# Keep only these OSM classes (labels). Set to NULL to keep all non-NA OSM classes.
allowed_osm_labels <- c("High green", "Low green", "Mixed", "Mid green")

# Global LCZ rebalance settings (applied after OSM-balanced selection).
min_global_lcz_per_class <- 20L
max_lcz_rebalance_iter <- 200L

seed <- 42

export_folder <- "exports"
export_points <- file.path(export_folder, "veg_plot_points.gpkg")
export_squares <- file.path(export_folder, "veg_plot_squares.gpkg")
export_stats <- file.path(export_folder, "veg_plot_statistics.csv")
export_diagnostics <- file.path(export_folder, "veg_plot_diagnostics.csv")
export_validation_city <- file.path(export_folder, "validation_points_per_city.csv")
export_validation_city_osm <- file.path(export_folder, "validation_points_per_city_osm.csv")
export_validation_lcz <- file.path(export_folder, "validation_points_per_lcz.csv")

message("Project root: ", project_root)

# -----------------------------
# Data loading
# -----------------------------
Study_area <- vect(study_area_path)
LCZ <- rast("Output/LCZ_3classes.tif")
OSM <- rast("Output/OSM_raster_10m.tif")

names(LCZ) <- "LCZ_class"
names(OSM) <- "OSM_class"

# Normalize labels for robust class matching.
normalize_class <- function(x) {
  tolower(trimws(as.character(x)))
}

# Read raster level table so label-based filtering can also match numeric codes.
get_osm_level_map <- function(r) {
  lv <- levels(r)
  if (is.null(lv) || length(lv) == 0 || is.null(lv[[1]]) || nrow(lv[[1]]) == 0) return(NULL)
  tbl <- lv[[1]]

  num_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  chr_cols <- names(tbl)[vapply(tbl, function(col) is.character(col) || is.factor(col), logical(1))]
  if (length(num_cols) == 0 || length(chr_cols) == 0) return(NULL)

  out <- data.frame(
    code = as.character(tbl[[num_cols[1]]]),
    label = as.character(tbl[[chr_cols[1]]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$code) & !is.na(out$label), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out
}

osm_level_map <- get_osm_level_map(OSM)
allowed_osm_codes <- NULL
if (!is.null(allowed_osm_labels) && !is.null(osm_level_map)) {
  allowed_norm <- normalize_class(allowed_osm_labels)
  allowed_osm_codes <- unique(osm_level_map$code[normalize_class(osm_level_map[[2]]) %in% allowed_norm])
}

if (!identical(crs(OSM), crs(LCZ))) {
  OSM <- project(OSM, LCZ)
}
if (!identical(crs(Study_area), crs(LCZ))) {
  Study_area <- project(Study_area, LCZ)
}

# Pull a standardized city name from available municipality columns.
extract_city_label <- function(city_row) {
  city_df <- as.data.frame(city_row)
  candidates <- c(
    city_df$mun_name_1_2,
    city_df$mun_name_1,
    city_df$mun_name_n_2,
    city_df$mun_name_n,
    city_df$mun_name_f_2,
    city_df$mun_name_f
  )
  candidates <- as.character(candidates)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0) return("UNKNOWN_CITY")
  raw <- gsub("[\\[\\]']", "", candidates[1])
  city <- toupper(trimws(raw))
  if (city %in% c("WATERMAAL-BOSVOORDE", "WATERMAEL-BOITSFORT", "BRUSSEL", "BRUXELLES")) {
    city <- "BRUSSEL"
  }
  city
}

Study_area$city <- vapply(seq_len(nrow(Study_area)), function(i) extract_city_label(Study_area[i, ]), character(1))
Study_area_sub <- Study_area[Study_area$city %in% selected_cities, ]

if (nrow(Study_area_sub) == 0) {
  stop("No matching cities found for selected_cities: ", paste(selected_cities, collapse = ", "))
}

selected_labels <- unique(as.character(Study_area_sub$city))
message("Requested cities: ", paste(selected_cities, collapse = ", "))
message("Matched cities: ", paste(selected_labels, collapse = ", "))

# Build a 100x100 m square around a point center.
make_square <- function(x, y, size, crs_value) {
  half <- size / 2
  coords <- matrix(
    c(
      x - half, y - half,
      x + half, y - half,
      x + half, y + half,
      x - half, y + half,
      x - half, y - half
    ),
    ncol = 2,
    byrow = TRUE
  )
  vect(list(coords), type = "polygons", crs = crs_value)
}

# Compute percentage cover of classes inside a polygon for one raster.
cover_percentages <- function(r, poly, prefix) {
  vals <- values(mask(crop(r, poly), poly))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(data.frame())

  out <- prop.table(table(vals)) * 100
  out <- data.frame(class = names(out), pct = as.numeric(out))

  out |>
    mutate(var = paste0(prefix, "_", class)) |>
    dplyr::select(all_of(c("var", "pct"))) |>
    pivot_wider(
      names_from = var,
      values_from = all_of("pct"),
      values_fill = 0
    )
}

# Draw random candidate points inside city polygon (or inner buffer when possible).
sample_candidates_interior <- function(city_poly, n_candidates, half_size, seed_city) {
  set.seed(seed_city)
  city_inner <- try(buffer(city_poly, -half_size), silent = TRUE)
  if (inherits(city_inner, "try-error") || nrow(city_inner) == 0) city_inner <- city_poly

  pts <- try(
    spatSample(city_inner, size = n_candidates, method = "random"),
    silent = TRUE
  )

  if (inherits(pts, "try-error") || is.null(pts) || nrow(pts) == 0) return(NULL)
  pts
}

# Add points from a pool while keeping minimum distance to already selected points.
pick_additional_with_min_distance <- function(pool_df, existing_df, target_n, min_dist) {
  if (nrow(pool_df) == 0 || target_n <= 0) return(integer(0))

  if (is.null(existing_df)) {
    existing_df <- data.frame(x = numeric(0), y = numeric(0))
  }

  chosen <- integer(0)
  order_idx <- sample(seq_len(nrow(pool_df)))

  for (idx in order_idx) {
    if (length(chosen) >= target_n) break

    x0 <- pool_df$x[idx]
    y0 <- pool_df$y[idx]

    if (nrow(existing_df) > 0) {
      d_existing <- sqrt((existing_df$x - x0)^2 + (existing_df$y - y0)^2)
      if (!all(d_existing >= min_dist)) next
    }

    if (length(chosen) > 0) {
      chosen_df <- pool_df[pool_df$id %in% chosen, , drop = FALSE]
      d_chosen <- sqrt((chosen_df$x - x0)^2 + (chosen_df$y - y0)^2)
      if (!all(d_chosen >= min_dist)) next
    }

    chosen <- c(chosen, pool_df$id[idx])
  }

  chosen
}

# Last safety pass: drop spacing conflicts and refill from valid unused candidates.
sanitize_selected_spacing <- function(selected_df, candidates, target_n, min_dist) {
  if (is.null(selected_df) || nrow(selected_df) == 0) return(selected_df)

  keep_ids <- integer(0)
  order_idx <- sample(seq_len(nrow(selected_df)))

  for (idx in order_idx) {
    cand_row <- selected_df[idx, ]
    if (length(keep_ids) == 0) {
      keep_ids <- c(keep_ids, cand_row$id)
      next
    }

    kept <- selected_df[selected_df$id %in% keep_ids, , drop = FALSE]
    d <- sqrt((kept$x - cand_row$x)^2 + (kept$y - cand_row$y)^2)
    if (all(d >= min_dist)) keep_ids <- c(keep_ids, cand_row$id)
  }

  cleaned <- selected_df[selected_df$id %in% keep_ids, , drop = FALSE]

  if (nrow(cleaned) < target_n) {
    fill_pool <- candidates |>
      filter(!(.data$id %in% cleaned$id))

    extra_ids <- pick_additional_with_min_distance(
      pool_df = fill_pool,
      existing_df = cleaned,
      target_n = target_n - nrow(cleaned),
      min_dist = min_dist
    )

    if (length(extra_ids) > 0) {
      cleaned <- bind_rows(cleaned, candidates[candidates$id %in% extra_ids, , drop = FALSE])
    }
  }

  cleaned
}

sample_city_plots <- function(city_poly, city_id, seed_city) {
  half_size <- plot_size / 2

  pts <- sample_candidates_interior(city_poly, n_candidates, half_size, seed_city)
  if (is.null(pts)) {
    return(list(ok = FALSE, reason = "no_candidate_points"))
  }

  cand <- as.data.frame(pts, geom = "XY")
  if (nrow(cand) == 0) {
    return(list(ok = FALSE, reason = "empty_candidate_dataframe"))
  }

  cand$id <- seq_len(nrow(cand))

  ex_lcz <- terra::extract(LCZ, pts)
  ex_osm <- terra::extract(OSM, pts)
  cand$LCZ <- ex_lcz[, 2]
  cand$OSM <- ex_osm[, 2]

  n_raw <- nrow(cand)
  cand <- cand |>
    filter(!is.na(LCZ), !is.na(OSM))
  n_non_na <- nrow(cand)

  cand_lcz_counts <- paste(names(table(cand$LCZ)), as.integer(table(cand$LCZ)), sep = ":", collapse = "; ")
  cand_osm_counts <- paste(names(table(cand$OSM)), as.integer(table(cand$OSM)), sep = ":", collapse = "; ")

  if (!is.null(allowed_osm_labels)) {
    osm_chr <- as.character(cand$OSM)
    keep <- normalize_class(osm_chr) %in% normalize_class(allowed_osm_labels)

    if (!is.null(allowed_osm_codes) && length(allowed_osm_codes) > 0) {
      keep <- keep | (osm_chr %in% allowed_osm_codes)
    }

    cand <- cand[keep, , drop = FALSE]
  }
  n_allowed <- nrow(cand)

  if (nrow(cand) == 0) {
    return(list(
      ok = FALSE,
      reason = "all_candidates_filtered",
      diagnostic = data.frame(
        city_id = city_id,
        n_raw = n_raw,
        n_non_na = n_non_na,
        n_allowed = n_allowed,
        cand_lcz_counts = paste(cand_lcz_counts, collapse = "; "),
        cand_osm_counts = paste(cand_osm_counts, collapse = "; ")
      )
    ))
  }

  cand$OSM <- as.character(cand$OSM)

  # Determine per-OSM target count.
  if (!is.null(n_plots_per_osm_class)) {
    selected_per_class <- as.integer(n_plots_per_osm_class)
  } else {
    selected_per_class <- floor(n_plots_per_city / length(allowed_osm_labels))
  }

  selected_ids <- integer(0)

  for (osm_class in allowed_osm_labels) {
    pool <- cand |>
      filter(normalize_class(.data$OSM) == normalize_class(osm_class), !(.data$id %in% selected_ids))

    if (nrow(pool) == 0) next

    current_selected <- cand |>
      filter(.data$id %in% selected_ids)

    picked <- pick_additional_with_min_distance(
      pool_df = pool,
      existing_df = current_selected,
      target_n = selected_per_class,
      min_dist = min_distance
    )

    if (length(picked) < selected_per_class) {
      remaining_pool <- pool |>
        filter(!(.data$id %in% picked))

      if (nrow(remaining_pool) > 0) {
        already_for_class <- pool |>
          filter(.data$id %in% picked)

        extra_ids <- pick_additional_with_min_distance(
          pool_df = remaining_pool,
          existing_df = bind_rows(current_selected, already_for_class),
          target_n = selected_per_class - length(picked),
          min_dist = min_distance
        )

        picked <- c(picked, extra_ids)
      }
    }

    selected_ids <- c(selected_ids, picked)
  }

  if (length(selected_ids) < n_plots_per_city) {
    remaining <- cand |>
      filter(!(.data$id %in% selected_ids))

    if (nrow(remaining) > 0) {
      current_selected <- cand |>
        filter(.data$id %in% selected_ids)

      extra_ids <- pick_additional_with_min_distance(
        pool_df = remaining,
        existing_df = current_selected,
        target_n = n_plots_per_city - length(selected_ids),
        min_dist = min_distance
      )

      selected_ids <- c(selected_ids, extra_ids)
    }
  }

  selected_df <- cand |>
    filter(.data$id %in% selected_ids)

  selected_df <- sanitize_selected_spacing(
    selected_df = selected_df,
    candidates = cand,
    target_n = n_plots_per_city,
    min_dist = min_distance
  )

  if (nrow(selected_df) == 0) {
    return(list(ok = FALSE, reason = "no_selected_points"))
  }

  points <- vect(selected_df, geom = c("x", "y"), crs = crs(city_poly))
  points$city_id <- city_id

  sq_list <- lapply(seq_len(nrow(selected_df)), function(i) {
    sq <- make_square(selected_df$x[i], selected_df$y[i], plot_size, crs(city_poly))
    sq$plot_id <- paste0(city_id, "_", i)
    sq
  })
  squares <- do.call(rbind, sq_list)

  stats <- lapply(seq_len(nrow(squares)), function(i) {
    sq <- squares[i, ]

    lcz_cov <- cover_percentages(LCZ, sq, "LCZ")
    osm_cov <- cover_percentages(OSM, sq, "OSM")

    bind_cols(
      data.frame(
        plot_id = squares$plot_id[i],
        city_id = city_id,
        x = selected_df$x[i],
        y = selected_df$y[i],
        LCZ_center = selected_df$LCZ[i],
        OSM_center = selected_df$OSM[i]
      ),
      lcz_cov,
      osm_cov
    )
  })

  stats <- bind_rows(stats)

  list(
    ok = TRUE,
    points = points,
    squares = squares,
    stats = stats,
    selected_df = selected_df,
    candidates = cand,
    diagnostic = data.frame(
      city_id = city_id,
      n_raw = n_raw,
      n_non_na = n_non_na,
      n_allowed = n_allowed,
      n_selected = nrow(selected_df),
      cand_lcz_counts = paste(cand_lcz_counts, collapse = "; "),
      cand_osm_counts = paste(cand_osm_counts, collapse = "; ")
    )
  )
}

# Swap points to improve global LCZ balance while preserving city + OSM quotas.
rebalance_lcz_global <- function(city_results, min_per_class, max_iter) {
  for (iter in seq_len(max_iter)) {
    selected_all <- bind_rows(lapply(city_results, function(res) res$selected_df))
    if (nrow(selected_all) == 0) break

    lcz_counts <- sort(table(selected_all$LCZ))
    if (length(lcz_counts) == 0 || all(as.integer(lcz_counts) >= min_per_class)) break

    need_class <- names(which.min(lcz_counts))
    donor_class <- names(which.max(lcz_counts))
    if (as.integer(lcz_counts[[donor_class]]) <= min_per_class) break

    swapped <- FALSE
    for (city_id in names(city_results)) {
      selected_df <- city_results[[city_id]]$selected_df
      candidates <- city_results[[city_id]]$candidates
      if (is.null(selected_df) || is.null(candidates) || nrow(selected_df) == 0) next

      donor_rows <- selected_df |>
        filter(.data$LCZ == donor_class)
      if (nrow(donor_rows) == 0) next

      for (j in seq_len(nrow(donor_rows))) {
        donor_row <- donor_rows[j, ]
        other_selected <- selected_df |>
          filter(.data$id != donor_row$id)

        replacement_pool <- candidates |>
          filter(
            .data$OSM == donor_row$OSM,
            .data$LCZ == need_class,
            !(.data$id %in% selected_df$id),
            !((paste(.data$x, .data$y)) %in% (paste(other_selected$x, other_selected$y)))
          )

        if (nrow(replacement_pool) > 0 && nrow(other_selected) > 0) {
          keep_dist <- vapply(seq_len(nrow(replacement_pool)), function(k) {
            d <- sqrt((other_selected$x - replacement_pool$x[k])^2 + (other_selected$y - replacement_pool$y[k])^2)
            all(d >= min_distance)
          }, logical(1))
          replacement_pool <- replacement_pool[keep_dist, , drop = FALSE]
        }

        if (nrow(replacement_pool) == 0) next

        replacement_row <- replacement_pool[1, ]
        selected_df <- selected_df |>
          filter(.data$id != donor_row$id)
        selected_df <- bind_rows(selected_df, replacement_row)
        selected_df <- sanitize_selected_spacing(
          selected_df = selected_df,
          candidates = candidates,
          target_n = n_plots_per_city,
          min_dist = min_distance
        )
        city_results[[city_id]]$selected_df <- selected_df
        swapped <- TRUE
        break
      }

      if (swapped) break
    }

    if (!swapped) break
  }

  city_results
}

# Build final export layers/tables from selected points for one city.
build_city_outputs <- function(selected_df, city_id) {
  if (is.null(selected_df) || nrow(selected_df) == 0) return(NULL)

  selected_df$city_id <- city_id
  selected_df$plot_id <- paste0(city_id, "_", seq_len(nrow(selected_df)))

  points <- vect(selected_df, geom = c("x", "y"), crs = crs(Study_area))
  points$city_id <- city_id

  sq_list <- lapply(seq_len(nrow(selected_df)), function(i) {
    sq <- make_square(selected_df$x[i], selected_df$y[i], plot_size, crs(Study_area))
    sq$plot_id <- selected_df$plot_id[i]
    sq
  })
  squares <- do.call(rbind, sq_list)

  stats <- lapply(seq_len(nrow(selected_df)), function(i) {
    sq <- squares[i, ]
    lcz_cov <- cover_percentages(LCZ, sq, "LCZ")
    osm_cov <- cover_percentages(OSM, sq, "OSM")
    bind_cols(
      data.frame(
        plot_id = selected_df$plot_id[i],
        city_id = city_id,
        x = selected_df$x[i],
        y = selected_df$y[i],
        LCZ_center = selected_df$LCZ[i],
        OSM_center = selected_df$OSM[i]
      ),
      lcz_cov,
      osm_cov
    )
  })

  list(
    points = points,
    squares = squares,
    stats = bind_rows(stats)
  )
}

city_results <- list()
all_diag <- list()

for (i in seq_len(nrow(Study_area_sub))) {
  city <- Study_area_sub[i, ]
  city_id <- extract_city_label(city)

  res <- sample_city_plots(city_poly = city, city_id = city_id, seed_city = seed + i)

  if (!is.null(res$diagnostic)) all_diag[[city_id]] <- res$diagnostic

  if (!isTRUE(res$ok)) {
    message("City ", city_id, " skipped: ", res$reason)
    next
  }

  city_results[[city_id]] <- res

  message("City ", city_id, ": selected ", nrow(res$stats), " plots")
}

city_results <- rebalance_lcz_global(
  city_results,
  min_per_class = min_global_lcz_per_class,
  max_iter = max_lcz_rebalance_iter
)

all_points <- NULL
all_squares <- NULL
all_stats <- list()

for (city_id in names(city_results)) {
  rebuilt <- build_city_outputs(city_results[[city_id]]$selected_df, city_id)
  if (is.null(rebuilt)) next
  all_points <- if (is.null(all_points)) rebuilt$points else rbind(all_points, rebuilt$points)
  all_squares <- if (is.null(all_squares)) rebuilt$squares else rbind(all_squares, rebuilt$squares)
  all_stats[[city_id]] <- rebuilt$stats
}

all_stats <- bind_rows(all_stats)
all_diag <- bind_rows(all_diag)

if (!dir.exists(export_folder)) dir.create(export_folder)

if (!is.null(all_points)) {
  writeVector(all_points, export_points, overwrite = TRUE)
  message("Exported points: ", export_points)
} else {
  warning("No point plots to export.")
}

if (!is.null(all_squares)) {
  writeVector(all_squares, export_squares, overwrite = TRUE)
  message("Exported squares: ", export_squares)
} else {
  warning("No square plots to export.")
}

if (nrow(all_stats) > 0) {
  write.csv(all_stats, export_stats, row.names = FALSE)
  message("Exported statistics: ", export_stats)
} else {
  warning("No plot statistics to export.")
}

if (nrow(all_diag) > 0) {
  write.csv(all_diag, export_diagnostics, row.names = FALSE)
  message("Exported diagnostics: ", export_diagnostics)
}

if (nrow(all_stats) > 0) {
  validation_city <- all_stats |>
    dplyr::count(.data$city_id, name = "n_points") |>
    arrange(.data$city_id)

  validation_city_osm <- all_stats |>
    dplyr::count(.data$city_id, .data$OSM_center, name = "n_points") |>
    arrange(.data$city_id, .data$OSM_center)

  validation_lcz <- all_stats |>
    dplyr::count(.data$LCZ_center, name = "n_points") |>
    arrange(.data$LCZ_center)

  write.csv(validation_city, export_validation_city, row.names = FALSE)
  write.csv(validation_city_osm, export_validation_city_osm, row.names = FALSE)
  write.csv(validation_lcz, export_validation_lcz, row.names = FALSE)

  message("Validation - points per city:")
  print(validation_city)
  message("Validation - points per city and OSM:")
  print(validation_city_osm)
  message("Validation - points per LCZ:")
  print(validation_lcz)
}

if (!is.null(all_points)) {
  message("Total plots generated: ", nrow(all_points))
} else {
  message("Total plots generated: 0")
}
