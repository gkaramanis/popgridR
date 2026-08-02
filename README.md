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

It covers the geometry pipeline, not the rendering: `popgrid()` returns cells
as `sf` for plotting with ggplot2, where the Python package returns a
matplotlib figure. The styled output, colour palettes, per-country
recommended settings and detached-territory inset panels are not ported.

Development version, not on CRAN. See Limitations below.

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
given) → Hamilton quotas → drop regions whose quota rounds to zero →
landmass clustering → per landmass: scale to quota share, deform, rasterise
→ exact quota fitting → adjacency repair against true neighbours → interior
gap filling → re-anchor detached landmasses to the deformed mainland → snap
all blocks to one lattice.

Every step is exported (`project_regions()`, `detect_landmasses()`,
`deform_regions()`, `generate_cells()`, `fit_quota()`,
`enforce_topology()`, `fill_gaps()`) for use outside the wrapper.

## Input data

Any `sf` polygon data frame works, but two properties of the input matter a
lot in practice.

**Resolution.** Generalised boundaries (Natural Earth) union into tens of
polygons; official statistical boundaries at survey resolution union into
thousands, nearly all of them sub-square-kilometre islets that cannot
receive a cell. Dropping them first is worth it:

```r
regions <- regions |>
  dplyr::filter(!sf::st_is_empty(geometry)) |>
  sf::st_cast("MULTIPOLYGON") |>   # normalise before casting to POLYGON
  sf::st_cast("POLYGON") |>
  dplyr::filter(as.numeric(sf::st_area(geometry)) > 1e6) |>
  dplyr::group_by(region) |>
  dplyr::summarise(.groups = "drop")
```

The `st_cast("MULTIPOLYGON")` step is required: casting a column that mixes
POLYGON and MULTIPOLYGON straight to POLYGON silently keeps only each
feature's first part.

**Sub-units sharing a region name.** Rows sharing a value in `region_col`
pool into one region, so sub-unit polygons need no dissolving. A weight that
is constant within a region (a share or rate joined onto sub-units) is taken
once rather than summed; see `values_agg`.

**Remote territories** are placed at their real position relative to the
mainland, which can leave a single cell thousands of kilometres out and
stretch the plot. Filter them out before building if you don't want them.

## Limitations

- **Speed.** `fit_quota()` and `enforce_topology()` index the cell lattice
  with character keys and linear scans, where the Python original uses hash
  lookups. Cost grows sharply with `n` and with the number of regions.
  Full-resolution national boundaries at `n = 1000` are impractical; reduce
  `n`, or simplify the input as above.
- **Nested enclaves.** Regions fully enclosed by another (a capital city
  inside its surrounding region) constrain topology repair heavily, and it
  may exhaust `topology_passes` without converging.
- **`cartogramR`'s `gsm` method is not exposed**, only `dcn` and `gn`.

## References

- Dougenik, J. A., Chrisman, N. R., & Niemeyer, D. R. (1985). An algorithm
  to construct continuous area cartograms. *The Professional Geographer*,
  37(1), 75–81.
- Gastner, M. T., & Newman, M. E. J. (2004). Diffusion-based method for
  producing density-equalizing maps. *PNAS*, 101(20), 7499–7504.
- Original Python implementation:
  [databites-tech/popgrid](https://github.com/databites-tech/popgrid) by
  [Josep Ferrer](https://github.com/iamjosepferrer).
