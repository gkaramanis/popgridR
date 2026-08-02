# popgridR

Mosaic (block) cartograms in R: regions are deformed so that area is
proportional to a numeric weight, then rasterised into equal-area square
cells. Each region receives an exact cell count set by Hamilton
(largest-remainder) allocation, so one cell always represents the same
number of weight units (people, trips, euros).

popgridR is an R port of the Python package
[popgrid](https://github.com/databites-tech/popgrid) by
[Josep Ferrer](https://github.com/iamjosepferrer),
generalised to accept any `sf` polygon input and any numeric weight instead
of country codes with bundled population data. The port was written largely
with Claude Code.

Development version, not on CRAN.

## Installation

```r
devtools::install_github("gkaramanis/popgridR")

# optional, for the faster C deformation backends
install.packages("cartogramR")
```

## Usage

```r
library(popgridR)
library(sf)
library(ggplot2)

# any sf polygon data frame; any numeric weight
regions <- read_sf("regions.gpkg")   # user-supplied boundaries

cells <- popgrid(
  regions,
  values     = "weight",   # column name, named vector, or per-row vector
  region_col = "region",
  n          = 500         # or value_per_cell = 50000
)

ggplot(cells) +
  geom_sf(aes(fill = region), color = NA) +
  theme_void()

# what one cell represents
attr(cells, "value_per_cell")
# exact cell counts per region
attr(cells, "quotas")
```

## Example: Greece population grid

Admin-1 boundaries fetched with the bundled helper (requires
rnaturalearth); population figures from the 2021 ELSTAT census. Also in
`inst/examples/greece.R`.

```r
library(popgridR)
library(ggplot2)

# Natural Earth's admin-1 `name` column already is the 13 Greek peripheries.
regions <- load_country_regions("GRC")

# Mount Athos is counted separately by ELSTAT and too small to earn a cell
# at n = 1000 - drop it for a clean 13-region map.
regions <- regions[regions$region_name != "Mount Athos", ]

population <- c(
  "Attica"                    = 3814064,
  "Central Macedonia"         = 1795669,
  "Thessaly"                  =  688255,
  "West Greece"               =  643449,
  "Crete"                     =  624408,
  "East Macedonia and Thrace" =  562201,
  "Peloponnese"               =  539535,
  "Central Greece"            =  508254,
  "Epirus"                    =  319991,
  "South Aegean"              =  309015,
  "West Macedonia"            =  254595,
  "Ionian Islands"            =  204532,
  "North Aegean"              =  194943
)

cells <- popgrid(regions, values = population, n = 1000)

ggplot(cells) +
  geom_sf(aes(fill = region_name), color = NA) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank()) +
  labs(
    title    = "Greece - population grid",
    subtitle = sprintf("1 cell ~ %.0f people | %d cells total",
                       attr(cells, "value_per_cell"), nrow(cells))
  )
```

![Mosaic cartogram of Greece: 1000 square cells, one per ~10,500 people, coloured by region](man/figures/greece.png)

## Deformation methods

The `method` argument selects the cartogram backend:

| method     | backend                              | notes                                        |
|------------|--------------------------------------|----------------------------------------------|
| `"auto"`   | `"dcn"` if cartogramR is installed, else `"damped"` | default                     |
| `"dcn"`    | cartogramR (Dougenik et al. 1985, C) | fast, converges fully                        |
| `"gn"`     | cartogramR (Gastner–Newman 2004, C)  | better behaved on extreme area/weight ratios |
| `"damped"` | built-in pure R Dougenik             | no compiled dependency; mass clamped per iteration for stability |
| `"none"`   | pass-through                         | input geometry is used as-is                 |

Exact per-region cell counts are enforced after rasterisation regardless of
method, so residual area error from the deformation step does not affect
counts, only shape.

## Pipeline

`popgrid()` runs: equal-area projection (data-centred LAEA unless `crs` is
given) → Hamilton quotas → landmass clustering → per landmass: scale to
quota share, deform, rasterise → exact quota fitting → adjacency repair
against true neighbours → interior gap filling.

Every step is exported (`project_regions()`, `detect_landmasses()`,
`deform_regions()`, `generate_cells()`, `fit_quota()`,
`enforce_topology()`, `fill_gaps()`) for use outside the wrapper.

## References

- Dougenik, J. A., Chrisman, N. R., & Niemeyer, D. R. (1985). An algorithm
  to construct continuous area cartograms. *The Professional Geographer*,
  37(1), 75–81.
- Gastner, M. T., & Newman, M. E. J. (2004). Diffusion-based method for
  producing density-equalizing maps. *PNAS*, 101(20), 7499–7504.
- Original Python implementation:
  [databites-tech/popgrid](https://github.com/databites-tech/popgrid) by
  [Josep Ferrer](https://github.com/iamjosepferrer).
