# Cartogram deformation: pluggable methods behind deform_regions().
#
# "dcn" / "gn"  -> cartogramR (C implementation; Suggests)
# "damped"      -> built-in pure-R Dougenik with per-iteration mass clamping
# "none"        -> pass-through (geometry already deformed by the caller)

# ── Internal helpers (damped Dougenik) ────────────────────────────────────────

.region_stats <- function(gdf, values) {
  areas       <- pmax(as.numeric(sf::st_area(gdf)), 1e-9)
  total_area  <- sum(areas)

  cents <- suppressWarnings(sf::st_centroid(sf::st_geometry(gdf)))
  cx    <- vapply(cents, `[[`, numeric(1), 1)
  cy    <- vapply(cents, `[[`, numeric(1), 2)
  radius <- sqrt(areas / pi)

  desired_area   <- values / sum(values) * total_area
  desired_radius <- sqrt(desired_area / pi)
  mass           <- desired_radius - radius

  size_error      <- pmax(areas, desired_area) / pmin(areas, desired_area)
  force_reduction <- 1.0 / (1.0 + mean(size_error))

  list(cx = cx, cy = cy, radius = radius, mass = mass,
       force_reduction = force_reduction)
}

# Returns a closure that displaces an N×2 coordinate matrix.
.make_displacer <- function(cx, cy, radius, mass, frf) {
  force(cx); force(cy); force(radius); force(mass); force(frf)

  function(coords) {
    N  <- nrow(coords)
    px <- coords[, 1]; py <- coords[, 2]

    dx        <- outer(px, cx, "-")   # N × J
    dy        <- outer(py, cy, "-")
    dist      <- sqrt(dx^2 + dy^2)
    dist_safe <- pmax(dist, 1e-9)

    R_mat <- matrix(radius, N, length(radius), byrow = TRUE)
    m_mat <- matrix(mass,   N, length(mass),   byrow = TRUE)

    outside <- dist > R_mat
    force   <- ifelse(outside,
                      m_mat * R_mat / dist_safe,
                      m_mat * dist_safe / pmax(R_mat, 1e-9))

    new_x <- px + frf * rowSums(force * (dx / dist_safe))
    new_y <- py + frf * rowSums(force * (dy / dist_safe))

    cbind(new_x, new_y)
  }
}

# Apply a coordinate-displacement function to a single sfg POLYGON or
# MULTIPOLYGON. Returns a new sfg of the same type.
.transform_geom <- function(geom, fun) {
  # class(sfg)[1] is "XY", not the geometry type — use inherits()
  if (!inherits(geom, "POLYGON") && !inherits(geom, "MULTIPOLYGON")) return(geom)

  coords <- sf::st_coordinates(geom)   # matrix: X Y L1 [L2 [L3]]
  new_xy <- fun(coords[, c("X", "Y"), drop = FALSE])

  close_ring <- function(m) {
    if (!isTRUE(all.equal(m[1, ], m[nrow(m), ])))
      m <- rbind(m, m[1, , drop = FALSE])
    m
  }

  if (inherits(geom, "POLYGON")) {
    ring_ids <- coords[, "L1"]
    rings <- lapply(sort(unique(ring_ids)), function(l) {
      idx <- ring_ids == l
      close_ring(cbind(new_xy[idx, 1], new_xy[idx, 2]))
    })
    sf::st_polygon(rings)

  } else {  # MULTIPOLYGON
    part_ids <- coords[, "L2"]
    ring_ids <- coords[, "L1"]
    parts <- lapply(sort(unique(part_ids)), function(l2) {
      pidx   <- part_ids == l2
      p_xy   <- new_xy[pidx, , drop = FALSE]
      p_ring <- ring_ids[pidx]
      lapply(sort(unique(p_ring)), function(l1) {
        ridx <- p_ring == l1
        close_ring(cbind(p_xy[ridx, 1], p_xy[ridx, 2]))
      })
    })
    sf::st_multipolygon(parts)
  }
}

# ── Public functions ──────────────────────────────────────────────────────────

#' Damped Dougenik cartogram (pure R)
#'
#' Dougenik–Chrisman–Niemeyer (1985) with two stability guards: per-iteration
#' mass clamping (`max_step_fraction`) and a vertex budget that reverts to the
#' last good geometry if `sf::st_make_valid()` inflates the vertex count.
#' Deliberately trades convergence speed for robustness — downstream quota
#' fitting corrects residual area error exactly.
#'
#' @param gdf sf data frame (projected, equal-area) with a `region_name`
#'   column and polygon/multipolygon geometries.
#' @param values Named list/vector of target weights, or a numeric vector
#'   aligned with `gdf` rows.
#' @param value_key Column used to look up `values` when a named vector is passed.
#' @param iterations Number of rubber-sheet iterations (15–25 for ~20 regions).
#' @param max_step_fraction Per-iteration mass cap as fraction of mean radius.
#' @param verbose Print per-iteration area error.
#' @return Copy of `gdf` with deformed geometries.
#' @export
dougenik_cartogram <- function(gdf,
                               values,
                               value_key         = "region_name",
                               iterations        = 20L,
                               max_step_fraction = 0.35,
                               verbose           = TRUE) {
  out <- gdf

  if (!is.null(names(values))) {
    val <- pmax(vapply(as.character(out[[value_key]]), function(n) {
      v <- values[[n]]; if (is.null(v)) 0 else as.numeric(v)
    }, numeric(1)), 1e-9)
  } else {
    val <- pmax(as.numeric(values), 1e-9)
  }

  start_verts <- tryCatch(nrow(sf::st_coordinates(sf::st_geometry(out))), error = function(e) 0L)
  vert_budget <- max(2000L, start_verts * 4L)
  prev_geoms  <- sf::st_geometry(out)

  for (it in seq_len(iterations)) {
    stats        <- .region_stats(out, val)
    step_cap     <- max_step_fraction * mean(stats$radius)
    mass_clamped <- pmax(pmin(stats$mass, step_cap), -step_cap)

    displace <- .make_displacer(
      stats$cx, stats$cy, stats$radius, mass_clamped, stats$force_reduction
    )

    new_geoms <- sf::st_sfc(
      lapply(sf::st_geometry(out), .transform_geom, fun = displace),
      crs = sf::st_crs(out)
    )
    new_geoms <- sf::st_make_valid(new_geoms)

    nverts <- tryCatch(nrow(sf::st_coordinates(new_geoms)), error = function(e) 0L)
    if (nverts > vert_budget) {
      if (verbose) message(sprintf(
        "cartogram stop at iter %d/%d - vertices %d exceed budget %d",
        it, iterations, nverts, vert_budget
      ))
      sf::st_geometry(out) <- prev_geoms
      break
    }

    prev_geoms <- new_geoms
    sf::st_geometry(out) <- new_geoms

    if (verbose) {
      areas   <- as.numeric(sf::st_area(out))
      desired <- val / sum(val) * sum(areas)
      err     <- pmax(areas, desired) / pmin(pmax(areas, 1e-9), desired)
      message(sprintf(
        "cartogram iter %2d/%d - mean area error %.3f, max %.3f",
        it, iterations, mean(err), max(err)
      ))
    }
  }

  out
}

#' Deform region polygons so area becomes proportional to values
#'
#' Dispatcher over the available cartogram implementations. `"dcn"` and
#' `"gn"` require the suggested \pkg{cartogramR} package; `"damped"` is the
#' built-in pure-R fallback; `"none"` returns the input unchanged. `"auto"`
#' picks `"dcn"` when \pkg{cartogramR} is installed, else `"damped"`.
#'
#' @param gdf Projected sf data frame with polygon/multipolygon geometries.
#' @param values Numeric vector of weights, one per row of `gdf`.
#' @param method One of `"auto"`, `"dcn"`, `"gn"`, `"damped"`, `"none"`.
#' @param iterations,max_step_fraction Passed to [dougenik_cartogram()]
#'   (`"damped"` method only).
#' @param options Named list passed to `cartogramR::cartogramR()`'s `options`
#'   argument (`"dcn"`/`"gn"` only).
#' @param verbose Print progress messages.
#' @return Copy of `gdf` with deformed geometries.
#' @export
deform_regions <- function(gdf,
                           values,
                           method            = c("auto", "dcn", "gn", "damped", "none"),
                           iterations        = 20L,
                           max_step_fraction = 0.35,
                           options           = NULL,
                           verbose           = TRUE) {
  method <- match.arg(method)
  stopifnot(length(values) == nrow(gdf))

  if (method == "auto") {
    method <- if (requireNamespace("cartogramR", quietly = TRUE)) "dcn" else "damped"
    if (verbose) message("deform method: ", method)
  }

  if (method == "none") return(gdf)

  if (method == "damped") {
    return(dougenik_cartogram(gdf, as.numeric(values),
                              iterations        = iterations,
                              max_step_fraction = max_step_fraction,
                              verbose           = verbose))
  }

  # dcn / gn via cartogramR
  if (!requireNamespace("cartogramR", quietly = TRUE)) {
    stop("method = \"", method, "\" requires the 'cartogramR' package. ",
         "Install it or use method = \"damped\".")
  }

  dat <- gdf
  # cartogramR requires strictly positive counts and POLYGON/MULTIPOLYGON.
  dat$.popgridR_count <- pmax(as.numeric(values), max(as.numeric(values)) * 1e-6)
  geom_types <- as.character(sf::st_geometry_type(dat))
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    dat <- sf::st_collection_extract(sf::st_make_valid(dat), "POLYGON")
  }
  dat <- sf::st_cast(dat, "MULTIPOLYGON")

  res <- cartogramR::cartogramR(dat, count = ".popgridR_count", method = method,
                                options = options)

  deformed <- if (inherits(res, "cartogramR")) res$cartogram else res
  out <- gdf
  sf::st_geometry(out) <- sf::st_make_valid(sf::st_geometry(sf::st_as_sf(deformed)))
  out
}

#' Per-region (target_share, actual_share, ratio) for validation
#'
#' @param gdf Projected sf data frame.
#' @param values Named list/vector of weights keyed by `value_key`.
#' @param value_key Column used to look up `values`.
#' @return data.frame with region, target_share, actual_share, ratio.
#' @export
area_error <- function(gdf, values, value_key = "region_name") {
  areas  <- as.numeric(sf::st_area(gdf))
  val    <- vapply(as.character(gdf[[value_key]]), function(n) {
    v <- values[[n]]; if (is.null(v)) 0 else as.numeric(v)
  }, numeric(1))
  tshare <- val / sum(val)
  ashare <- areas / sum(areas)
  data.frame(
    region       = gdf[[value_key]],
    target_share = tshare,
    actual_share = ashare,
    ratio        = ifelse(ashare > 0, tshare / ashare, NA_real_),
    stringsAsFactors = FALSE
  )
}

#' Scale all geometries uniformly so total area equals target_area
#'
#' @param gdf Projected sf data frame.
#' @param target_area Target total area in squared projected units.
#' @return Copy of `gdf` with scaled geometries.
#' @export
scale_to_area <- function(gdf, target_area) {
  cur <- sum(as.numeric(sf::st_area(gdf)))
  if (cur <= 0 || target_area <= 0) return(gdf)

  fac      <- sqrt(target_area / cur)
  centroid <- suppressWarnings(sf::st_centroid(sf::st_union(sf::st_geometry(gdf))))
  ox <- centroid[[1]][1]; oy <- centroid[[1]][2]

  out <- gdf
  sf::st_geometry(out) <- sf::st_sfc(
    lapply(sf::st_geometry(gdf), function(g) {
      .transform_geom(g, function(coords) {
        cbind(ox + (coords[, 1] - ox) * fac,
              oy + (coords[, 2] - oy) * fac)
      })
    }),
    crs = sf::st_crs(gdf)
  )
  out
}
