# Optional convenience: admin-1 boundaries via rnaturalearth (Suggests).

#' Load admin-1 regions for a country as an sf data frame (WGS84)
#'
#' Convenience helper around \pkg{rnaturalearth}. Not required by
#' [popgrid()] — any sf polygon data frame works as input.
#'
#' @param country_iso3 ISO 3166-1 alpha-3 code, e.g. "ESP", "DEU".
#' @param dissolve_by  Optional column to dissolve sub-units into parent regions.
#'   Typical value: "region_group" or "region" (Natural Earth column names).
#'   NULL (default) returns raw admin-1 rows.
#'
#' @return sf data frame (WGS84) with columns region_id, region_name, geometry.
#' @export
load_country_regions <- function(country_iso3, dissolve_by = NULL) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE))
    stop("load_country_regions() requires the 'rnaturalearth' package.")
  if (!requireNamespace("rnaturalearthdata", quietly = TRUE))
    stop("load_country_regions() requires the 'rnaturalearthdata' package.")

  iso <- toupper(trimws(country_iso3))

  # ne_states() doesn't accept ISO3. Use ne_countries() — which is keyed by
  # adm0_a3 — to resolve the country name, then pass that to ne_states().
  countries   <- rnaturalearth::ne_countries(returnclass = "sf")
  country_row <- countries[toupper(trimws(countries$adm0_a3)) == iso, ]
  if (nrow(country_row) == 0) stop("Unknown ISO3 code '", iso, "'")
  country_name <- country_row$name[1]

  regions <- rnaturalearth::ne_states(country = country_name, returnclass = "sf")

  if (nrow(regions) == 0) stop("No regions found for ISO3 code '", iso, "'")

  regions <- regions[!sf::st_is_empty(regions$geometry), ]
  regions$geometry <- sf::st_make_valid(regions$geometry)

  if (!is.null(dissolve_by)) {
    if (!(dissolve_by %in% names(regions))) {
      warning("Column '", dissolve_by, "' not found - returning raw admin-1 rows.")
    } else {
      filled <- !is.na(regions[[dissolve_by]]) & trimws(regions[[dissolve_by]]) != ""
      regions <- regions[filled, ]
      regions[[dissolve_by]] <- trimws(regions[[dissolve_by]])

      nms <- sort(unique(trimws(regions[[dissolve_by]])))
      dissolved_list <- lapply(nms, function(nm) {
        sub  <- regions[trimws(regions[[dissolve_by]]) == nm, ]
        geom <- sf::st_make_valid(sf::st_union(sf::st_geometry(sub)))
        sf::st_sf(region_name = nm, geometry = geom)
      })
      dissolved <- do.call(rbind, dissolved_list)
      dissolved$region_id <- seq_len(nrow(dissolved)) - 1L
      return(sf::st_sf(
        region_id   = dissolved$region_id,
        region_name = dissolved$region_name,
        geometry    = sf::st_geometry(dissolved)
      ))
    }
  }

  regions$region_id <- seq_len(nrow(regions)) - 1L

  # Prefer the GeoNames English label (gn_name) for readable region names,
  # falling back to the native `name` where gn_name is absent/blank.
  region_name <- regions$name
  if ("gn_name" %in% names(regions)) {
    en       <- trimws(regions$gn_name)
    use_en   <- !is.na(en) & en != ""
    region_name[use_en] <- en[use_en]
  }
  regions$region_name <- region_name

  sf::st_sf(
    region_id   = regions$region_id,
    region_name = regions$region_name,
    geometry    = regions$geometry
  )
}
