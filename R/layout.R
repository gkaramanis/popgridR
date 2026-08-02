# Landmass layout. Ported from popgrid/core.py (_fill_outline_holes,
# _reanchor_secondary, and the global lattice snap in AreaGrid.build).
#
# Without these, each landmass is rasterised at its own geographic position
# after being scaled independently, so a small landmass whose weight share
# far exceeds its area share grows into the mainland and the two cell blocks
# overlap.

# Drop interior rings from a polygon/multipolygon outline.
#
# Deformed region polygons leave thin gaps between one another, so their
# union has holes and cells falling in those holes are lost during
# rasterisation. Filling them lets every interior cell be rasterised and
# assigned to its nearest region.
.fill_outline_holes <- function(geom) {
  if (is.null(geom)) return(geom)

  if (inherits(geom, "POLYGON")) {
    return(sf::st_polygon(geom[1]))          # exterior ring only
  }
  if (inherits(geom, "MULTIPOLYGON")) {
    parts <- lapply(geom, function(p) sf::st_polygon(p[1]))
    if (length(parts) == 0) return(geom)
    merged <- sf::st_union(do.call(sf::st_sfc, parts))
    if (length(merged) == 0) return(geom)
    return(merged[[1]])
  }
  geom
}

# Rebuild square cell geometry from cx/cy.
.cells_to_boxes <- function(cells, cell_size) {
  half <- cell_size / 2
  sf::st_geometry(cells) <- sf::st_sfc(
    mapply(function(x, y) {
      sf::st_polygon(list(matrix(
        c(x - half, y - half, x + half, y - half,
          x + half, y + half, x - half, y + half,
          x - half, y - half),
        ncol = 2, byrow = TRUE
      )))
    }, cells$cx, cells$cy, SIMPLIFY = FALSE),
    crs = sf::st_crs(cells)
  )
  cells
}

# Reposition non-mainland cells so they hold the same position relative to the
# *deformed* mainland that they held relative to the *geographic* mainland.
#
# Each landmass moves as a rigid block: its centroid is remapped into the
# deformed frame and the whole block translated by that offset. Cells are not
# rescaled individually, which would push them off the lattice.
.reanchor_secondary <- function(cells, mainland_id, geo_bbox, def_bbox, cell_size) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)
  others <- setdiff(unique(cells$landmass_id), mainland_id)
  if (length(others) == 0) return(cells)

  gw <- geo_bbox[3] - geo_bbox[1]; if (gw == 0) gw <- 1
  gh <- geo_bbox[4] - geo_bbox[2]; if (gh == 0) gh <- 1

  out <- cells
  for (lm_id in others) {
    m   <- out$landmass_id == lm_id
    cx0 <- mean(out$cx[m]); cy0 <- mean(out$cy[m])
    rx  <- def_bbox[1] + (cx0 - geo_bbox[1]) / gw * (def_bbox[3] - def_bbox[1])
    ry  <- def_bbox[2] + (cy0 - geo_bbox[2]) / gh * (def_bbox[4] - def_bbox[2])
    out$cx[m] <- out$cx[m] + (rx - cx0)
    out$cy[m] <- out$cy[m] + (ry - cy0)
  }
  .cells_to_boxes(out, cell_size)
}

# Snap every cell onto one global lattice (the mainland's grid phase) so the
# landmass blocks line up instead of sitting a fraction of a cell off. Each
# landmass is already internally cell-spaced, so this is a uniform per-landmass
# shift and cannot create collisions.
.snap_to_lattice <- function(cells, mainland_id, cell_size) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)

  ref <- cells[cells$landmass_id == mainland_id, ]
  if (nrow(ref) == 0) ref <- cells
  ox <- min(ref$cx); oy <- min(ref$cy)

  cells$cx <- ox + round((cells$cx - ox) / cell_size) * cell_size
  cells$cy <- oy + round((cells$cy - oy) / cell_size) * cell_size
  .cells_to_boxes(cells, cell_size)
}
