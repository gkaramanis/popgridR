# Main entry point: sf regions + weights -> mosaic cartogram cells.

#' Build a mosaic (block) cartogram
#'
#' Takes any sf polygon data frame and a numeric weight per region, and
#' returns an sf data frame of equal-sized square cells where each region's
#' cell count is proportional to its weight (fixed exactly by Hamilton
#' largest-remainder allocation).
#'
#' Pipeline: project (equal-area) -> allocate quotas -> detect landmasses ->
#' per landmass: scale to quota share, deform ([deform_regions()]),
#' rasterise ([generate_cells()]) -> fit exact quotas ([fit_quota()]) ->
#' repair adjacency ([enforce_topology()]) -> fill interior gaps
#' ([fill_gaps()]).
#'
#' Rows sharing the same value in `region_col` are treated as one region:
#' their weights pool (see `values_agg`) and their cells count against one
#' quota. Sub-unit polygons therefore need no dissolving beforehand.
#'
#' @param x sf data frame with POLYGON/MULTIPOLYGON geometries in any CRS.
#'   Geographic (lon/lat) input is auto-projected to a data-centred Lambert
#'   Azimuthal Equal-Area CRS unless `crs` is given.
#' @param values Weights controlling each region's cell share. Either the
#'   name of a numeric column in `x`, a named numeric vector keyed by
#'   `region_col` values, or an unnamed numeric vector aligned with `x` rows.
#' @param region_col Name of the column identifying regions. Default
#'   `"region_name"`.
#' @param n Total number of cells to allocate. Ignored when
#'   `value_per_cell` is given.
#' @param value_per_cell Alternative to `n`: weight units represented by one
#'   cell (e.g. `50000` for "1 cell = 50 000 trips").
#' @param values_agg How to pool `values` when several rows share a region.
#'   `"auto"` (default) sums per-row quantities but takes a value that is
#'   constant within every multi-row region only once — the case where a
#'   share or rate has been joined onto sub-unit polygons, where summing
#'   would inflate regions in proportion to their row count. `"sum"` and
#'   `"first"` force the behaviour. Ignored when `values` is a named vector,
#'   which is already per-region.
#' @param method Deformation method, see [deform_regions()]. Default
#'   `"auto"`: \pkg{cartogramR}'s `"dcn"` when installed, else `"damped"`.
#' @param crs Optional CRS (anything `sf::st_crs()` accepts) to project into
#'   instead of the automatic LAEA. Must be projected, equal-area strongly
#'   recommended.
#' @param cluster_distance_m Polygons whose centroids are within this
#'   distance (metres) cluster into one landmass. Default 200 km.
#' @param oversize Pre-rasterisation area multiplier compensating perimeter
#'   cell losses. Default 1.18.
#' @param iterations,max_step_fraction Tuning for `method = "damped"`, see
#'   [dougenik_cartogram()].
#' @param deform_options Named list forwarded to
#'   `cartogramR::cartogramR(options = )` for `method = "dcn"`/`"gn"`.
#' @param topology_passes Maximum repair passes in [enforce_topology()].
#' @param fill Run [fill_gaps()] after topology repair. Default TRUE.
#' @param verbose Print progress messages.
#'
#' @return sf data frame of square cells with columns `cx`, `cy`,
#'   `region_id`, `region_name`, `landmass_id` and, when `region_col` is not
#'   `"region_name"`, a copy of the region key under that name. Attributes:
#'   `quotas` (named integer vector), `cell_size` (metres),
#'   `value_per_cell` (weight units per cell), `n` (requested total).
#'
#' @examples
#' # Four synthetic regions in a projected CRS
#' library(sf)
#' r <- function(x0) st_polygon(list(cbind(c(x0, x0 + 1, x0 + 1, x0, x0) * 1e5,
#'                                          c(0, 0, 1e5, 1e5, 0))))
#' regions <- st_sf(name = c("A", "B", "C", "D"),
#'                  geometry = st_sfc(lapply(0:3, r), crs = 3857))
#' cells <- popgrid(regions, values = c(A = 40, B = 30, C = 20, D = 10),
#'                  region_col = "name", n = 100, method = "damped",
#'                  verbose = FALSE)
#' table(cells$region_name)
#' @export
popgrid <- function(x,
                    values,
                    region_col         = "region_name",
                    n                  = 1000L,
                    value_per_cell     = NULL,
                    values_agg         = c("auto", "sum", "first"),
                    method             = c("auto", "dcn", "gn", "damped", "none"),
                    crs                = NULL,
                    cluster_distance_m = 200000,
                    oversize           = 1.18,
                    iterations         = 20L,
                    max_step_fraction  = 0.35,
                    deform_options     = NULL,
                    topology_passes    = 30L,
                    fill               = TRUE,
                    verbose            = TRUE) {
  method     <- match.arg(method)
  values_agg <- match.arg(values_agg)

  # ── Validate input ──────────────────────────────────────────────────────
  if (!inherits(x, "sf")) {
    stop("`x` must be an sf data frame of region polygons, not ",
         class(x)[1], ".\n",
         "  Pass the sf object itself; name its columns with `values` and ",
         "`region_col`:\n",
         "    popgrid(my_sf, values = \"weight_col\", region_col = \"name_col\")")
  }
  regions <- x
  if (!(region_col %in% names(regions)))
    stop("`region_col` = '", region_col, "' is not a column in `x`. ",
         "Available: ", paste(setdiff(names(regions), attr(regions, "sf_column")),
                              collapse = ", "))
  if (is.na(sf::st_crs(regions)))
    stop("`x` has no CRS. Assign one with sf::st_set_crs() first.")

  region_name <- as.character(regions[[region_col]])
  if (anyNA(region_name)) stop("`", region_col, "` contains NA values.")

  # ── Drop rows with no geometry ──────────────────────────────────────────
  # Statistical boundary files routinely carry placeholder rows (offshore,
  # migratory, "no usual address") with empty geometry. They cannot be
  # tiled, and they poison the projection origin, so remove them before
  # anything else reads the geometry.
  nonempty <- !sf::st_is_empty(sf::st_geometry(regions))
  if (!all(nonempty)) {
    if (verbose)
      message(sum(!nonempty), " row(s) with empty geometry dropped.")
    regions     <- regions[nonempty, ]
    region_name <- region_name[nonempty]
    if (nrow(regions) == 0) stop("No non-empty geometries in `x`.")
  }

  # ── Resolve values to one numeric per row ───────────────────────────────
  if (is.character(values) && length(values) == 1) {
    if (!(values %in% names(regions)))
      stop("`values` = '", values, "' is not a column in `x`. ",
           "Available: ", paste(setdiff(names(regions), attr(regions, "sf_column")),
                                collapse = ", "))
    vals <- as.numeric(regions[[values]])
  } else if (!is.null(names(values))) {
    missing_nms <- setdiff(unique(region_name), names(values))
    if (length(missing_nms) > 0)
      stop("No value supplied for: ", paste(missing_nms, collapse = ", "))
    vals <- as.numeric(values[region_name])
  } else if (is.numeric(values) && length(values) == nrow(regions)) {
    vals <- as.numeric(values)
  } else {
    stop("`values` must be one of:\n",
         "  - a column name in `x`, e.g. values = \"population\"\n",
         "  - a named numeric vector keyed by `region_col` values\n",
         "  - an unnamed numeric vector with one value per row of `x` (",
         nrow(regions), " needed, got ", length(values), ")")
  }
  if (anyNA(vals)) stop("`values` contains NA for region(s): ",
                        paste(unique(region_name[is.na(vals)]), collapse = ", "))
  if (any(vals < 0)) stop("`values` must be non-negative.")

  # ── Pool weights per region ─────────────────────────────────────────────
  # Rows sharing a region name pool into one weight. A named vector is
  # already per-region. For a column or per-row vector, `values_agg` decides:
  # "auto" sums per-row quantities (populations of sub-units) but takes a
  # value that is constant within every multi-row region only once (a share
  # or rate repeated across sub-units), since summing those would inflate
  # regions in proportion to their row count.
  if (!is.null(names(values)) && !is.character(values)) {
    weights <- setNames(as.numeric(values[unique(region_name)]),
                        unique(region_name))
  } else {
    grp       <- split(vals, region_name)
    multi_row <- any(lengths(grp) > 1L)
    constant  <- all(vapply(grp, function(v) length(unique(v)) == 1L, logical(1)))

    agg <- values_agg
    if (agg == "auto") agg <- if (multi_row && constant) "first" else "sum"

    w <- vapply(grp, if (agg == "first") function(v) v[[1]] else sum, numeric(1))
    weights <- setNames(as.numeric(w), names(w))

    if (verbose && multi_row) {
      message(sprintf(
        "%d rows pooled into %d regions by `%s`; values %s.",
        nrow(regions), length(weights), region_col,
        if (agg == "first")
          "constant within each region, taken once (values_agg = \"sum\" to add them)"
        else "summed per region"
      ))
    }
  }

  # ── Cell count ──────────────────────────────────────────────────────────
  if (!is.null(value_per_cell)) {
    if (value_per_cell <= 0) stop("`value_per_cell` must be positive.")
    n <- max(1L, as.integer(round(sum(weights) / value_per_cell)))
  }
  n <- as.integer(n)
  if (n < 1) stop("`n` must be at least 1.")

  # ── Project ─────────────────────────────────────────────────────────────
  if (!is.null(crs)) {
    regions_proj <- sf::st_transform(regions, sf::st_crs(crs))
    if (isTRUE(sf::st_is_longlat(regions_proj)))
      stop("`crs` must be a projected CRS (units metres), not geographic.")
  } else if (isTRUE(sf::st_is_longlat(regions))) {
    if (verbose) message("Projecting to data-centred LAEA.")
    regions_proj <- project_regions(regions)$regions
  } else {
    regions_proj <- regions
  }

  regions_proj <- sf::st_make_valid(regions_proj)
  keep         <- !sf::st_is_empty(sf::st_geometry(regions_proj))
  regions_proj <- regions_proj[keep, ]
  region_name  <- region_name[keep]
  vals         <- vals[keep]

  regions_proj$region_name <- region_name
  regions_proj$region_id   <- match(region_name, unique(region_name)) - 1L

  # ── Quotas ──────────────────────────────────────────────────────────────
  quotas <- hamilton_quotas(weights, n)
  if (verbose) {
    zero <- names(quotas)[quotas == 0]
    if (length(zero) > 0)
      message("Regions with 0 cells at n = ", n, ": ",
              paste(zero, collapse = ", "))
  }

  # ── Landmasses & adjacency ──────────────────────────────────────────────
  landmasses <- detect_landmasses(regions_proj, cluster_distance_m)
  if (verbose) message(length(landmasses), " landmass(es) detected.")

  lm_names_list <- lapply(landmasses, function(lm)
    unique(regions_proj$region_name[lm$region_idx]))
  dup_names <- unlist(lm_names_list)[duplicated(unlist(lm_names_list))]
  if (length(dup_names) > 0)
    warning("Region(s) split across landmasses (quota applied in each): ",
            paste(unique(dup_names), collapse = ", "),
            ". Consider a larger `cluster_distance_m`.")

  real_adj <- real_adjacency(regions_proj)

  total_area <- sum(as.numeric(sf::st_area(regions_proj)))
  cell_size  <- cell_size_for_n(total_area, n)

  # ── Per landmass: scale -> deform -> rasterise ──────────────────────────
  cell_frames <- lapply(landmasses, function(lm) {
    lm_regions <- regions_proj[lm$region_idx, ]
    lm_vals    <- vals[lm$region_idx]

    lm_quota_total <- sum(quotas[unique(lm_regions$region_name)], na.rm = TRUE)
    if (lm_quota_total == 0) return(NULL)  # too small for any cell at this n

    # Target area = quota-weighted share of total (x oversize for perimeter loss)
    target_area <- lm_quota_total / n * total_area * oversize
    lm_scaled   <- scale_to_area(lm_regions, target_area)

    lm_deformed <- if (nrow(lm_regions) > 1) {
      deform_regions(lm_scaled, pmax(lm_vals, 1e-9),
                     method            = method,
                     iterations        = iterations,
                     max_step_fraction = max_step_fraction,
                     options           = deform_options,
                     verbose           = verbose)
    } else {
      lm_scaled  # single region: scaling is enough
    }

    lm_union <- sf::st_union(sf::st_geometry(lm_deformed))[[1]]
    generate_cells(lm_union, lm_deformed, cell_size, lm$id)
  })

  cells <- do.call(rbind, Filter(Negate(is.null), cell_frames))
  if (is.null(cells) || nrow(cells) == 0)
    stop("No cells generated - check `n` and input geometries.")
  if (verbose) message(nrow(cells), " raw cells generated.")

  # ── Post-processing ─────────────────────────────────────────────────────
  cells <- fit_quota(cells, as.list(quotas), cell_size)
  cells <- enforce_topology(cells, real_adj, as.list(quotas), cell_size,
                            passes = topology_passes)
  if (fill) cells <- fill_gaps(cells, cell_size)
  if (verbose) message(nrow(cells), " cells after post-processing.")

  if (region_col != "region_name") cells[[region_col]] <- cells$region_name

  attr(cells, "quotas")         <- quotas
  attr(cells, "cell_size")      <- cell_size
  attr(cells, "value_per_cell") <- sum(weights) / n
  attr(cells, "n")              <- n
  cells
}
