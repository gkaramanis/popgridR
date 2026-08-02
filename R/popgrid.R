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
#' their weights pool and their cells count against one quota.
#'
#' @param regions sf data frame with POLYGON/MULTIPOLYGON geometries in any
#'   CRS. Geographic (lon/lat) input is auto-projected to a data-centred
#'   Lambert Azimuthal Equal-Area CRS unless `crs` is given.
#' @param values Weights controlling each region's cell share. Either the
#'   name of a numeric column in `regions`, a named numeric vector keyed by
#'   `region_col` values, or an unnamed numeric vector aligned with
#'   `regions` rows.
#' @param region_col Name of the column identifying regions. Default
#'   `"region_name"`.
#' @param n Total number of cells to allocate. Ignored when
#'   `value_per_cell` is given.
#' @param value_per_cell Alternative to `n`: weight units represented by one
#'   cell (e.g. `50000` for "1 cell = 50 000 trips").
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
popgrid <- function(regions,
                    values,
                    region_col         = "region_name",
                    n                  = 1000L,
                    value_per_cell     = NULL,
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
  method <- match.arg(method)

  # ── Validate input ──────────────────────────────────────────────────────
  if (!inherits(regions, "sf")) stop("`regions` must be an sf data frame.")
  if (!(region_col %in% names(regions)))
    stop("Column '", region_col, "' not found in `regions`.")
  if (is.na(sf::st_crs(regions)))
    stop("`regions` has no CRS. Assign one with sf::st_set_crs() first.")

  region_name <- as.character(regions[[region_col]])
  if (anyNA(region_name)) stop("`", region_col, "` contains NA values.")

  # ── Resolve values to one numeric per row ───────────────────────────────
  if (is.character(values) && length(values) == 1) {
    if (!(values %in% names(regions)))
      stop("Values column '", values, "' not found in `regions`.")
    vals <- as.numeric(regions[[values]])
  } else if (!is.null(names(values))) {
    missing_nms <- setdiff(unique(region_name), names(values))
    if (length(missing_nms) > 0)
      stop("No value supplied for: ", paste(missing_nms, collapse = ", "))
    vals <- as.numeric(values[region_name])
  } else if (is.numeric(values) && length(values) == nrow(regions)) {
    vals <- as.numeric(values)
  } else {
    stop("`values` must be a column name, a named numeric vector, ",
         "or an unnamed numeric vector with one value per row.")
  }
  if (anyNA(vals)) stop("`values` contains NA for region(s): ",
                        paste(unique(region_name[is.na(vals)]), collapse = ", "))
  if (any(vals < 0)) stop("`values` must be non-negative.")

  # Weights pooled per region name (handles multiple rows per region).
  # A named vector is a per-region total, so take it once per name;
  # column/positional values are per-row and sum.
  if (!is.null(names(values)) && !is.character(values)) {
    weights <- values[unique(region_name)]
    weights <- setNames(as.numeric(weights), unique(region_name))
  } else {
    w       <- tapply(vals, region_name, sum)
    weights <- setNames(as.numeric(w), names(w))
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
