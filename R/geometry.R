# Equal-area projection and landmass detection.

# ── CRS helpers ───────────────────────────────────────────────────────────────

#' Lambert Azimuthal Equal-Area CRS centred on (lon_0, lat_0), units metres
#'
#' @param lon_0,lat_0 Projection origin in degrees.
#' @return An `sf::st_crs()` object.
#' @export
laea_crs <- function(lon_0, lat_0) {
  sf::st_crs(sprintf(
    "+proj=laea +lat_0=%.6f +lon_0=%.6f +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
    lat_0, lon_0
  ))
}

#' Reproject a geographic sf data frame to a data-centred LAEA CRS
#'
#' Uses the median of each polygon's bbox centre as the projection origin,
#' so a single large outlying territory (e.g. Alaska) doesn't pull the centre
#' away from the main inhabited area.
#'
#' @param regions_wgs84 sf data frame in a geographic (lon/lat) CRS.
#' @return Named list: `regions` (projected sf) and `crs` (the LAEA crs object).
#' @export
project_regions <- function(regions_wgs84) {
  # Per-polygon bbox centre
  poly_bounds <- do.call(rbind, lapply(sf::st_geometry(regions_wgs84), function(g) {
    bb <- sf::st_bbox(g)
    c((bb["xmin"] + bb["xmax"]) / 2,
      (bb["ymin"] + bb["ymax"]) / 2)
  }))

  lon_0 <- stats::median(poly_bounds[, 1], na.rm = TRUE)
  lat_0 <- stats::median(poly_bounds[, 2], na.rm = TRUE)
  if (is.na(lon_0) || is.na(lat_0))
    stop("Cannot derive a projection origin - all geometries are empty.")

  crs       <- laea_crs(lon_0, lat_0)
  projected <- sf::st_transform(regions_wgs84, crs)

  list(regions = projected, crs = crs)
}

# ── Landmass detection ────────────────────────────────────────────────────────

# Flatten any sfg geometry to a list of POLYGON sfg objects.
# Note: class(sfg)[1] is "XY", not the geometry type — use inherits().
.explode_to_polygons <- function(geom) {
  if (inherits(geom, "POLYGON")) return(list(geom))
  if (inherits(geom, "MULTIPOLYGON")) {
    return(lapply(seq_along(geom), function(i) sf::st_polygon(geom[[i]])))
  }
  if (inherits(geom, "GEOMETRYCOLLECTION")) {
    unlist(lapply(seq_along(geom), function(i) .explode_to_polygons(geom[[i]])),
           recursive = FALSE)
  } else {
    list()
  }
}

# Connected components of the "centroids within max_dist_m" graph.
# Returns a list of integer-index vectors, one per cluster.
#
# The neighbour search is delegated to sf::st_is_within_distance, which uses
# a spatial index. A pairwise distance loop is O(n^2) and untenable here:
# a national coastline unions into thousands of polygons (Australia gives
# ~6,600, i.e. ~22 million pairs).
.union_find_cluster <- function(polygons, max_dist_m) {
  n <- length(polygons)
  if (n == 1) return(list(1L))

  cents <- suppressWarnings(sf::st_centroid(do.call(sf::st_sfc, polygons)))
  nb    <- sf::st_is_within_distance(cents, dist = max_dist_m)

  # Union-Find with path halving.
  parent <- seq_len(n)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  for (i in seq_len(n)) {
    for (j in nb[[i]]) {
      if (j <= i) next
      ra <- find_root(i); rb <- find_root(j)
      if (ra != rb) parent[ra] <- rb
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  unname(split(seq_len(n), roots))
}

#' Identify separate landmass clusters from a projected sf data frame
#'
#' @param regions_proj Projected (equal-area) sf data frame.
#' @param cluster_distance_m Polygon centroids within this distance (metres) are
#'   merged into the same landmass. Default 200 km works for most countries.
#'
#' @return List of landmass objects (each a list with `id`, `geom`, `region_idx`),
#'   sorted by area descending so index 1 is always the main landmass.
#' @export
detect_landmasses <- function(regions_proj, cluster_distance_m = 200000) {
  full_union <- sf::st_union(sf::st_geometry(regions_proj))

  if (is.null(full_union) || sf::st_is_empty(full_union)) {
    stop("Geometry union is empty - check input data.")
  }

  raw_polygons <- .explode_to_polygons(full_union[[1]])
  n_poly <- length(raw_polygons)

  if (n_poly == 1) {
    return(list(list(
      id         = 0L,
      geom       = raw_polygons[[1]],
      region_idx = seq_len(nrow(regions_proj))
    )))
  }

  poly_clusters <- .union_find_cluster(raw_polygons, cluster_distance_m)

  landmasses <- lapply(seq_along(poly_clusters), function(ci) {
    idx          <- poly_clusters[[ci]]
    cluster_sfc  <- do.call(sf::st_sfc, c(raw_polygons[idx],
                                          list(crs = sf::st_crs(regions_proj))))
    cluster_geom <- sf::st_union(cluster_sfc)[[1]]

    # Which regions intersect this cluster?
    cluster_sf  <- sf::st_sf(geometry = sf::st_sfc(cluster_geom, crs = sf::st_crs(regions_proj)))
    region_idx  <- which(
      !sf::st_is_empty(sf::st_geometry(regions_proj)) &
      lengths(sf::st_intersects(regions_proj, cluster_sf)) > 0
    )

    list(id = as.integer(ci - 1L), geom = cluster_geom, region_idx = region_idx)
  })

  # Sort largest-first: index 1 is always the main landmass
  areas <- vapply(landmasses, function(lm) {
    as.numeric(sf::st_area(sf::st_sfc(lm$geom, crs = sf::st_crs(regions_proj))))
  }, numeric(1))
  landmasses <- landmasses[order(-areas)]

  # Deduplicate: a region spanning multiple clusters (e.g. a coastal province
  # with a tiny offshore fragment that clusters separately) must belong to
  # exactly one landmass — the one it overlaps most by area.
  .assign_regions_to_landmasses(landmasses, regions_proj)
}

# Each region is assigned to the single landmass whose geometry it overlaps most.
.assign_regions_to_landmasses <- function(landmasses, regions_proj) {
  crs       <- sf::st_crs(regions_proj)
  n_regions <- nrow(regions_proj)

  # Build sfc for each landmass geometry once
  lm_sfcs <- lapply(landmasses, function(lm) {
    sf::st_sfc(lm$geom, crs = crs)
  })

  # For each region find which landmass has maximum intersection area
  best_lm <- vapply(seq_len(n_regions), function(ri) {
    reg_sfc <- sf::st_geometry(regions_proj)[ri]
    overlaps <- vapply(lm_sfcs, function(lm_sfc) {
      inter <- suppressWarnings(sf::st_intersection(reg_sfc, lm_sfc))
      if (length(inter) == 0 || all(sf::st_is_empty(inter))) return(0)
      as.numeric(sf::st_area(inter))
    }, numeric(1))
    which.max(overlaps)
  }, integer(1))

  # Rebuild region_idx lists
  for (i in seq_along(landmasses)) {
    landmasses[[i]]$region_idx <- which(best_lm == i)
  }

  # Drop landmasses that ended up with no regions
  Filter(function(lm) length(lm$region_idx) > 0, landmasses)
}
