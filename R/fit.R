# Post-processing on the rasterised cell grid: exact quota fitting,
# topology repair against real adjacency, interior gap filling.
#
# The sparse integer-lattice grid is a named character vector:
# names are "i_j" keys, values are region names.

# ── Grid key helpers ──────────────────────────────────────────────────────────

.key_enc <- function(i, j) paste(i, j, sep = "_")

.key_dec <- function(k) {
  p <- strsplit(k, "_", fixed = TRUE)[[1]]
  c(as.integer(p[1]), as.integer(p[2]))
}

.nbrs <- function(k) {
  p <- strsplit(k, "_", fixed = TRUE)[[1]]
  i <- as.integer(p[1]); j <- as.integer(p[2])
  c(.key_enc(i + 1L, j), .key_enc(i - 1L, j),
    .key_enc(i, j + 1L), .key_enc(i, j - 1L))
}

# Reconstruct an sf data frame of square cells from a grid (named char vector).
.grid_to_sf <- function(grid, region_id_of, lm_id, minx, miny, cell_size, crs) {
  if (length(grid) == 0) return(NULL)

  ks  <- names(grid)
  ij  <- do.call(rbind, lapply(ks, .key_dec))
  cx  <- minx + ij[, 1] * cell_size
  cy  <- miny + ij[, 2] * cell_size

  half  <- cell_size / 2.0
  geoms <- sf::st_sfc(
    mapply(function(x, y) {
      sf::st_polygon(list(matrix(
        c(x - half, y - half,
          x + half, y - half,
          x + half, y + half,
          x - half, y + half,
          x - half, y - half),
        ncol = 2, byrow = TRUE
      )))
    }, cx, cy, SIMPLIFY = FALSE),
    crs = crs
  )

  sf::st_sf(
    cx          = cx,
    cy          = cy,
    region_id   = region_id_of[as.character(grid)],
    region_name = as.character(grid),
    landmass_id = lm_id,
    geometry    = geoms,
    stringsAsFactors = FALSE
  )
}

#' Adjust each region's cell count to exactly its quota
#'
#' Trimming removes least-connected boundary cells; growth adds most-connected
#' adjacent empty cells. Operates per landmass. Regions are keyed by
#' `region_name`, so multiple input polygons sharing a name pool into one
#' quota.
#'
#' @param cells sf data frame with cx, cy, region_id, region_name, landmass_id.
#' @param quotas Named list/integer vector: target cell count per region name.
#' @param cell_size Side length of each square cell.
#' @return sf data frame of cells with counts fitted to quotas.
#' @export
fit_quota <- function(cells, quotas, cell_size) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)

  out_frames <- list()

  for (lm_id in unique(cells$landmass_id)) {
    lm_cells <- cells[cells$landmass_id == lm_id, ]
    minx <- min(lm_cells$cx); miny <- min(lm_cells$cy)

    keys_vec     <- .key_enc(round((lm_cells$cx - minx) / cell_size),
                             round((lm_cells$cy - miny) / cell_size))
    occupied     <- setNames(lm_cells$region_name, keys_vec)
    region_id_of <- setNames(lm_cells$region_id,   lm_cells$region_name)
    reg_keys     <- tapply(keys_vec, lm_cells$region_name, unique, simplify = FALSE)

    # ── Drop disconnected specks ──────────────────────────────────────────
    if (length(reg_keys) > 1) {
      for (nm in names(reg_keys)) {
        kset <- reg_keys[[nm]]
        if (length(kset) <= 1) next

        seen <- character(0); components <- list()
        for (start in kset) {
          if (start %in% seen) next
          comp  <- start; stack <- start; seen <- c(seen, start)
          while (length(stack) > 0) {
            cur   <- stack[length(stack)]; stack <- stack[-length(stack)]
            for (nb in .nbrs(cur)) {
              if (nb %in% kset && !(nb %in% seen)) {
                seen <- c(seen, nb); comp <- c(comp, nb); stack <- c(stack, nb)
              }
            }
          }
          components <- c(components, list(comp))
        }

        if (length(components) > 1) {
          largest <- components[[which.max(lengths(components))]]
          drop    <- unlist(components)[!(unlist(components) %in% largest)]
          occupied <- occupied[!(names(occupied) %in% drop)]
          reg_keys[[nm]] <- largest
        }
      }
    }

    # ── Trim over-quota regions ───────────────────────────────────────────
    for (nm in names(reg_keys)) {
      quota <- as.integer(quotas[[nm]] %||% length(reg_keys[[nm]]))

      while (length(reg_keys[[nm]]) > max(quota, 1L)) {
        cand <- NULL; cand_score <- NULL

        for (k in reg_keys[[nm]]) {
          nb_k    <- .nbrs(k)
          empty   <- sum(!(nb_k %in% names(occupied)))
          if (empty == 0) next
          same    <- sum(occupied[nb_k[nb_k %in% names(occupied)]] == nm, na.rm = TRUE)
          score   <- c(empty, -same)

          if (is.null(cand_score) ||
              score[1] > cand_score[1] ||
              (score[1] == cand_score[1] && score[2] > cand_score[2])) {
            cand_score <- score; cand <- k
          }
        }

        if (is.null(cand)) break
        reg_keys[[nm]] <- setdiff(reg_keys[[nm]], cand)
        occupied       <- occupied[names(occupied) != cand]
      }
    }

    # ── Grow under-quota regions ──────────────────────────────────────────
    for (nm in names(reg_keys)) {
      quota <- as.integer(quotas[[nm]] %||% length(reg_keys[[nm]]))
      guard <- 0L

      while (length(reg_keys[[nm]]) < quota && guard < 50000L) {
        guard <- guard + 1L

        # Collect empty neighbours -> best-connected wins
        cands <- setNames(integer(0), character(0))
        for (k in reg_keys[[nm]]) {
          for (nb in .nbrs(k)) {
            if (!(nb %in% names(occupied))) {
              sc <- as.integer(sum(.nbrs(nb) %in% reg_keys[[nm]]))
              if (is.na(cands[nb])) cands[nb] <- sc
              else cands[nb] <- max(cands[nb], sc)
            }
          }
        }

        if (length(cands) == 0) break
        best           <- names(which.max(cands))
        reg_keys[[nm]] <- c(reg_keys[[nm]], best)
        occupied[best] <- nm
      }
    }

    # ── Rebuild ───────────────────────────────────────────────────────────
    recs <- do.call(rbind, lapply(names(reg_keys), function(nm) {
      ks <- reg_keys[[nm]]
      if (length(ks) == 0) return(NULL)
      ij <- do.call(rbind, lapply(ks, .key_dec))
      data.frame(cx = minx + ij[, 1] * cell_size,
                 cy = miny + ij[, 2] * cell_size,
                 region_id   = region_id_of[[nm]],
                 region_name = nm,
                 landmass_id = lm_id,
                 stringsAsFactors = FALSE)
    }))
    if (is.null(recs) || nrow(recs) == 0) next

    frame <- .grid_to_sf(
      setNames(recs$region_name, .key_enc(
        round((recs$cx - minx) / cell_size),
        round((recs$cy - miny) / cell_size)
      )),
      region_id_of, lm_id, minx, miny, cell_size, sf::st_crs(cells)
    )
    if (!is.null(frame)) out_frames[[length(out_frames) + 1]] <- frame
  }

  if (length(out_frames) == 0) return(cells[0, ])
  do.call(rbind, out_frames)
}

#' Ground-truth border graph
#'
#' @param regions_proj Projected sf data frame with a `region_name` column.
#' @return Named list: for each region name, a character vector of neighbours.
#' @export
real_adjacency <- function(regions_proj) {
  nms    <- regions_proj$region_name
  scale  <- sqrt(mean(as.numeric(sf::st_area(regions_proj)))) * 0.02
  buffed <- sf::st_buffer(sf::st_geometry(regions_proj), scale)

  ints <- sf::st_intersects(buffed, buffed)
  adj  <- setNames(vector("list", length(nms)), nms)
  for (nm in nms) adj[[nm]] <- character(0)

  for (i in seq_along(nms)) {
    for (j in ints[[i]]) {
      if (j == i) next
      adj[[nms[i]]] <- union(adj[[nms[i]]], nms[j])
      adj[[nms[j]]] <- union(adj[[nms[j]]], nms[i])
    }
  }
  adj
}

#' Repair grid adjacencies so regions touch only if they touch geographically
#'
#' @param cells sf data frame with cx, cy, region_id, region_name, landmass_id.
#' @param real_adj Named list from [real_adjacency()].
#' @param quotas Named list/integer vector of target cell counts.
#' @param cell_size Side length of each square cell.
#' @param passes Maximum number of repair passes.
#' @return sf data frame of repaired cells.
#' @export
enforce_topology <- function(cells, real_adj, quotas, cell_size, passes = 30L) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)

  out_frames <- list()

  for (lm_id in unique(cells$landmass_id)) {
    lm_cells <- cells[cells$landmass_id == lm_id, ]
    present  <- unique(lm_cells$region_name)

    if (length(present) < 2) {
      out_frames[[length(out_frames) + 1]] <- lm_cells
      next
    }

    minx <- min(lm_cells$cx); miny <- min(lm_cells$cy)
    keys_vec     <- .key_enc(round((lm_cells$cx - minx) / cell_size),
                             round((lm_cells$cy - miny) / cell_size))
    grid         <- setNames(lm_cells$region_name, keys_vec)
    region_id_of <- setNames(lm_cells$region_id,   lm_cells$region_name)

    reg_keys_of <- function(region) names(grid)[grid == region]

    connected_without <- function(region, drop) {
      keys <- setdiff(reg_keys_of(region), drop)
      if (length(keys) == 0) return(FALSE)
      seen <- keys[1]; stack <- keys[1]
      while (length(stack) > 0) {
        cur <- stack[length(stack)]; stack <- stack[-length(stack)]
        for (nb in .nbrs(cur)) {
          if (nb %in% keys && !(nb %in% seen)) {
            seen <- c(seen, nb); stack <- c(stack, nb)
          }
        }
      }
      length(seen) == length(keys)
    }

    legal_at <- function(region, k) {
      adj_R <- real_adj[[region]] %||% character(0)
      for (nb in .nbrs(k)) {
        o <- grid[nb][[1L]]   # single-bracket: NA for absent keys, no crash
        if (!is.null(o) && !is.na(o) && o != region && !(o %in% adj_R))
          return(FALSE)
      }
      TRUE
    }

    for (pass in seq_len(passes)) {
      changed <- 0L
      counts  <- as.list(table(grid))   # list for NULL-safe access

      # ── 1. Repair illegal adjacencies ──────────────────────────────────
      for (k in names(grid)) {
        A     <- grid[[k]]; if (is.null(A) || is.na(A)) next
        adj_A <- real_adj[[A]] %||% character(0)

        for (nb_k in .nbrs(k)) {
          B <- grid[nb_k][[1L]]   # single-bracket: NA for absent keys, no crash
          if (is.null(B) || is.na(B) || B == A || B %in% adj_A) next

          conn <- intersect(intersect(adj_A, real_adj[[B]] %||% character(0)), present)
          if (length(conn) == 0) next

          C      <- conn[which.min(vapply(conn, function(c) as.integer(counts[[c]] %||% 0L), integer(1)))]
          over_A <- as.integer(counts[[A]] %||% 0L) - as.integer(quotas[[A]] %||% 0L)
          over_B <- as.integer(counts[[B]] %||% 0L) - as.integer(quotas[[B]] %||% 0L)

          if (over_A >= over_B) { target_k <- k;    target_r <- A
          } else                { target_k <- nb_k; target_r <- B }

          if (legal_at(C, target_k) && connected_without(target_r, target_k)) {
            grid[[target_k]] <- C
            counts[[target_r]] <- as.integer(counts[[target_r]] %||% 0L) - 1L
            counts[[C]]        <- as.integer(counts[[C]]        %||% 0L) + 1L
            changed <- changed + 1L
            break
          }
        }
      }

      # ── 2. Topology-safe gradient quota flow ───────────────────────────
      surplus_order <- sort(vapply(present, function(r) {
        -(as.integer(counts[[r]] %||% 0L) - as.integer(quotas[[r]] %||% 0L))
      }, numeric(1)))

      for (region in names(surplus_order)) {
        guard <- 0L
        err_R <- as.integer(counts[[region]] %||% 0L) - as.integer(quotas[[region]] %||% 0L)

        while (err_R > 0 && guard < 4000L) {
          guard <- guard + 1L; moved <- FALSE

          for (k in reg_keys_of(region)) {
            for (nb_k in .nbrs(k)) {
              S <- grid[nb_k][[1L]]   # single-bracket: NA for absent keys
              if (is.null(S) || is.na(S) || !(S %in% present) || S == region) next
              err_S <- as.integer(counts[[S]] %||% 0L) - as.integer(quotas[[S]] %||% 0L)

              if (err_R - err_S >= 2 && legal_at(S, k) && connected_without(region, k)) {
                grid[[k]] <- S
                counts[[region]] <- as.integer(counts[[region]] %||% 0L) - 1L
                counts[[S]]      <- as.integer(counts[[S]]      %||% 0L) + 1L
                err_R <- err_R - 1L; moved <- TRUE; break
              }
            }
            if (moved) break
          }
          if (!moved) break
        }
      }

      if (changed == 0L) break
    }

    frame <- .grid_to_sf(grid, region_id_of, lm_id, minx, miny, cell_size,
                          sf::st_crs(cells))
    if (!is.null(frame)) out_frames[[length(out_frames) + 1]] <- frame
  }

  if (length(out_frames) == 0) return(cells[0, ])
  do.call(rbind, out_frames)
}

#' Close empty cells inside the filled area (seam filling)
#'
#' @param cells sf data frame with cx, cy, region_id, region_name, landmass_id.
#' @param cell_size Side length of each square cell.
#' @param quotas Optional named integer vector (unused in filling, kept for API parity).
#' @param min_neighbours Minimum occupied neighbours to fill an exterior cell.
#' @param passes Maximum flood-fill passes.
#' @return sf data frame of cells with interior pockets filled.
#' @export
fill_gaps <- function(cells, cell_size, quotas = NULL,
                      min_neighbours = 3L, passes = 12L) {
  if (is.null(cells) || nrow(cells) == 0) return(cells)

  out_frames <- list()

  for (lm_id in unique(cells$landmass_id)) {
    lm_cells <- cells[cells$landmass_id == lm_id, ]
    minx <- min(lm_cells$cx); miny <- min(lm_cells$cy)

    keys_vec     <- .key_enc(round((lm_cells$cx - minx) / cell_size),
                             round((lm_cells$cy - miny) / cell_size))
    occ          <- setNames(lm_cells$region_name, keys_vec)
    region_id_of <- setNames(lm_cells$region_id,   lm_cells$region_name)
    orig_counts  <- as.list(table(occ))

    for (pass in seq_len(passes)) {
      all_ij <- do.call(rbind, lapply(names(occ), .key_dec))
      imin <- min(all_ij[, 1]) - 1L; imax <- max(all_ij[, 1]) + 1L
      jmin <- min(all_ij[, 2]) - 1L; jmax <- max(all_ij[, 2]) + 1L

      # ── Flood-fill exterior from bounding-box border ──────────────────
      border <- unique(c(
        unlist(lapply(seq(imin, imax), function(i)
          c(.key_enc(i, jmin), .key_enc(i, jmax)))),
        unlist(lapply(seq(jmin, jmax), function(j)
          c(.key_enc(imin, j), .key_enc(imax, j))))
      ))

      exterior <- character(0)
      stack    <- border[!(border %in% names(occ))]

      while (length(stack) > 0) {
        ck <- stack[length(stack)]; stack <- stack[-length(stack)]
        if (ck %in% exterior || ck %in% names(occ)) next
        exterior <- c(exterior, ck)
        ij <- .key_dec(ck); ci <- ij[1]; cj <- ij[2]
        for (nb in list(c(ci+1L,cj), c(ci-1L,cj), c(ci,cj+1L), c(ci,cj-1L))) {
          ni <- nb[1]; nj <- nb[2]
          if (ni >= imin && ni <= imax && nj >= jmin && nj <= jmax) {
            nb_k <- .key_enc(ni, nj)
            if (!(nb_k %in% names(occ))) stack <- c(stack, nb_k)
          }
        }
      }

      # ── Collect enclosed empty cells ──────────────────────────────────
      MAX_POCKET <- 50L
      enclosed <- character(0)
      for (i in seq(imin, imax)) {
        for (j in seq(jmin, jmax)) {
          k <- .key_enc(i, j)
          if (!(k %in% names(occ)) && !(k %in% exterior))
            enclosed <- c(enclosed, k)
        }
      }

      if (length(enclosed) == 0) break

      enclosed_set <- enclosed
      seen_p  <- character(0)
      to_fill <- setNames(character(0), character(0))

      for (start in enclosed) {
        if (start %in% seen_p) next

        pocket <- start; seen_p <- c(seen_p, start); qi <- 1L
        while (qi <= length(pocket)) {
          ij <- .key_dec(pocket[qi]); ci <- ij[1]; cj <- ij[2]; qi <- qi + 1L
          for (nb in list(c(ci+1L,cj), c(ci-1L,cj), c(ci,cj+1L), c(ci,cj-1L))) {
            nb_k <- .key_enc(nb[1], nb[2])
            if (nb_k %in% enclosed_set && !(nb_k %in% seen_p)) {
              seen_p <- c(seen_p, nb_k); pocket <- c(pocket, nb_k)
            }
          }
        }

        if (length(pocket) > MAX_POCKET) next

        for (k in pocket) {
          hit <- .nbrs(k); hit <- hit[hit %in% names(occ)]
          if (length(hit) > 0) {
            winner <- names(sort(table(occ[hit]), decreasing = TRUE))[1]
            to_fill[k] <- winner
          }
        }
      }

      if (length(to_fill) == 0) break

      occ[names(to_fill)] <- to_fill

      # ── Count-neutral: shed outer silhouette cells to compensate ──────
      now_counts <- as.list(table(occ))

      for (region in names(now_counts)) {
        surplus <- as.integer(now_counts[[region]] %||% 0L) -
                   as.integer(orig_counts[[region]] %||% 0L)

        for (s in seq_len(max(surplus, 0L))) {
          cand_k <- NULL; cand_e <- 0L
          for (k in names(occ)[occ == region]) {
            e <- sum(!(.nbrs(k) %in% names(occ)))
            if (e > cand_e) { cand_e <- e; cand_k <- k }
          }
          if (is.null(cand_k) || cand_e == 0L) break
          occ <- occ[names(occ) != cand_k]
        }
      }
    }

    frame <- .grid_to_sf(occ, region_id_of, lm_id, minx, miny, cell_size,
                          sf::st_crs(cells))
    if (!is.null(frame)) out_frames[[length(out_frames) + 1]] <- frame
  }

  if (length(out_frames) == 0) return(cells[0, ])
  do.call(rbind, out_frames)
}
