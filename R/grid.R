# Rasterise a landmass into N equal-area square cells.

# ── Cell-size helper ──────────────────────────────────────────────────────────

#' Cell side length so a perfect rasterisation yields exactly n cells
#'
#' `cell_size = sqrt(total_area / n)`
#'
#' @param total_area_m2 Total projected area in square metres.
#' @param n Target cell count (positive integer).
#' @return Cell side length in projected units (metres).
#' @export
cell_size_for_n <- function(total_area_m2, n) {
  if (n <= 0)             stop("n must be positive, got ", n)
  if (total_area_m2 <= 0) stop("total_area_m2 must be positive, got ", total_area_m2)
  sqrt(total_area_m2 / n)
}

# ── Forced single-cell fallback ───────────────────────────────────────────────

.single_centroid_cell <- function(landmass_geom, regions_proj, cell_size, landmass_id) {
  crs     <- sf::st_crs(regions_proj)
  lm_sfc  <- sf::st_sfc(landmass_geom, crs = crs)
  cent    <- sf::st_centroid(lm_sfc)
  cx      <- sf::st_coordinates(cent)[1, "X"]
  cy      <- sf::st_coordinates(cent)[1, "Y"]
  half    <- cell_size / 2.0

  pt_sf       <- sf::st_sf(data.frame(cx = cx, cy = cy, row_id = 1L),
                            geometry = cent, crs = crs)
  nearest_idx <- sf::st_nearest_feature(pt_sf, regions_proj)

  box <- sf::st_polygon(list(matrix(
    c(cx - half, cy - half,
      cx + half, cy - half,
      cx + half, cy + half,
      cx - half, cy + half,
      cx - half, cy - half),
    ncol = 2, byrow = TRUE
  )))

  sf::st_sf(
    cx          = cx,
    cy          = cy,
    region_id   = regions_proj$region_id[nearest_idx],
    region_name = regions_proj$region_name[nearest_idx],
    landmass_id = landmass_id,
    geometry    = sf::st_sfc(box, crs = crs)
  )
}

# ── Core grid generator ───────────────────────────────────────────────────────

#' Rasterise a single landmass into grid cells and assign each to a region
#'
#' Algorithm:
#' 1. Tile the landmass bounding box with candidate cell centres.
#' 2. Filter to centres strictly inside the landmass.
#' 3. `sf::st_join` (within) assigns each centre to its region.
#'    `sf::st_nearest_feature` handles boundary cells that fall through.
#' 4. Replace point geometry with square box polygons.
#'
#' @param landmass_geom sfg POLYGON/MULTIPOLYGON for the landmass (projected).
#' @param regions_proj  sf data frame for regions belonging to this landmass;
#'   must contain `region_id` (int) and `region_name` (chr) columns.
#' @param cell_size     Square cell side length in projected units (metres).
#' @param landmass_id   Integer tag written into every output row.
#'
#' @return sf data frame with columns cx, cy, region_id, region_name,
#'   landmass_id, geometry (square box per cell).
#' @export
generate_cells <- function(landmass_geom, regions_proj, cell_size, landmass_id) {
  if (nrow(regions_proj) == 0) {
    warning("Landmass ", landmass_id, ": no regions supplied - skipping.")
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(regions_proj))))
  }

  crs  <- sf::st_crs(regions_proj)
  bbox <- sf::st_bbox(sf::st_sfc(landmass_geom, crs = crs))
  half <- cell_size / 2.0

  # Guard: seq() errors when from > to with positive step (territory < cell_size)
  xs <- if (bbox["xmin"] + half <= bbox["xmax"])
    seq(bbox["xmin"] + half, bbox["xmax"], by = cell_size) else numeric(0)
  ys <- if (bbox["ymin"] + half <= bbox["ymax"])
    seq(bbox["ymin"] + half, bbox["ymax"], by = cell_size) else numeric(0)

  # Territory smaller than one cell
  if (length(xs) == 0 || length(ys) == 0) {
    return(.single_centroid_cell(landmass_geom, regions_proj, cell_size, landmass_id))
  }

  # ── Step 1: candidate centres ─────────────────────────────────────────────
  grid   <- expand.grid(cx = xs, cy = ys)
  pts_sf <- sf::st_as_sf(
    data.frame(cx = grid$cx, cy = grid$cy, row_id = seq_len(nrow(grid))),
    coords = c("cx", "cy"),
    remove = FALSE,
    crs    = crs
  )

  # ── Step 2: point-in-landmass filter ──────────────────────────────────────
  lm_sf  <- sf::st_sf(geometry = sf::st_sfc(landmass_geom, crs = crs))
  inside <- lengths(sf::st_within(pts_sf, lm_sf)) > 0

  if (!any(inside)) {
    warning("Landmass ", landmass_id,
            ": no cell centres fell inside landmass (thin/irregular) - forcing centroid cell.")
    return(.single_centroid_cell(landmass_geom, regions_proj, cell_size, landmass_id))
  }

  pts_in <- pts_sf[inside, ]

  # ── Step 3: spatial join — assign centres to regions ──────────────────────
  reg_sub <- regions_proj[c("region_id", "region_name")]
  joined  <- sf::st_join(pts_in, reg_sub, join = sf::st_within, left = TRUE)

  # Drop duplicates from points that touch shared boundaries
  joined <- joined[!duplicated(joined$row_id), ]

  # Nearest-region fallback for any still-unassigned boundary points
  miss <- is.na(joined$region_id)
  if (any(miss)) {
    ni <- sf::st_nearest_feature(joined[miss, ], regions_proj)
    joined$region_id[miss]   <- regions_proj$region_id[ni]
    joined$region_name[miss] <- regions_proj$region_name[ni]
  }

  # ── Step 4: replace point geometry with square boxes ──────────────────────
  cx_f <- joined$cx
  cy_f <- joined$cy

  box_geoms <- sf::st_sfc(
    mapply(function(x, y) {
      sf::st_polygon(list(matrix(
        c(x - half, y - half,
          x + half, y - half,
          x + half, y + half,
          x - half, y + half,
          x - half, y - half),
        ncol = 2, byrow = TRUE
      )))
    }, cx_f, cy_f, SIMPLIFY = FALSE),
    crs = crs
  )

  sf::st_geometry(joined) <- box_geoms
  joined$landmass_id      <- landmass_id
  joined$row_id           <- NULL

  joined
}
