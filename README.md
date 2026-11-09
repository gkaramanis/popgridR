# popgridR

popgridR draws mosaic cartograms: maps made of equal squares, where each
region gets as many squares as its weight calls for. A region with twice the
population gets twice the squares, and every square stands for the same
amount (people, trips, euros).

It works in two steps. First the regions are deformed so their areas match
their weights, then they are filled with square cells. The number of cells
per region is fixed in advance with Hamilton (largest-remainder) allocation,
so the counts are exact even when the shapes are not.

It is an R port of [popgrid](https://github.com/databites-tech/popgrid), a
Python package by [Josep Ferrer](https://github.com/iamjosepferrer). The
Python version takes a country code and uses its own population data. This
one takes any `sf` polygons and any numeric column. I wrote the port largely
with Claude Code.

popgridR stops at the geometry. `popgrid()` returns the cells as an `sf`
data frame for you to plot, with ggplot2 or anything else. The Python
package's styled figures, palettes, per-country settings and inset panels
are not ported.

It is a development version and not on CRAN.

## Installation

```r
devtools::install_github("gkaramanis/popgridR")

# Optional. Adds the faster deformation methods written in C.
install.packages("cartogramR")
```

## Usage

```r
library(popgridR)
library(sf)
library(ggplot2)

regions <- read_sf("regions.gpkg")

cells <- popgrid(
  regions,
  values     = "weight",   # a column name, a named vector or one value per row
  region_col = "region",
  n          = 500         # total cells; or value_per_cell = 50000
)

ggplot(cells) +
  geom_sf(aes(fill = region), color = NA) +
  theme_void()

# What one cell stands for
attr(cells, "value_per_cell")

# Cells per region
attr(cells, "quotas")
```

## Example: Greece

The boundaries come from Natural Earth through `load_country_regions()`,
which needs rnaturalearth. The population figures are from the 2021 ELSTAT
census. The same script is in `inst/examples/greece.R`.

```r
library(popgridR)
library(ggplot2)

# Natural Earth's admin-1 regions for Greece are the 13 peripheries
regions <- load_country_regions("GRC")

# Mount Athos is too small to get a cell at n = 1000
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

`method` picks how the regions are deformed before they are filled with
cells.

| method     | backend                                 | notes                                         |
|------------|-----------------------------------------|-----------------------------------------------|
| `"auto"`   | `"dcn"` with cartogramR, else `"damped"` | the default                                   |
| `"dcn"`    | cartogramR, Dougenik et al. (1985), C   | fast and converges fully                      |
| `"gn"`     | cartogramR, Gastner and Newman (2004), C | copes better with very uneven weights         |
| `"gsm"`    | cartogramR, Gastner, Seguy and More (2018), C | flow-based successor to `"gn"`, faster |
| `"damped"` | Dougenik in plain R                     | no compiled dependency, slower but stable     |
| `"none"`   | no deformation                          | the input shapes are used as they are         |

The method changes the shapes only. Cell counts are corrected after the
cells are drawn, so every method gives the same counts.

## What `popgrid()` does

1. Projects the regions to an equal-area CRS (Lambert azimuthal centred on
   the data, unless you pass `crs`).
2. Allocates cells to regions and drops those that get none.
3. Groups the regions into landmasses.
4. For each landmass, scales and deforms the regions and fills them with
   cells.
5. Adds or removes cells until every region has its exact count.
6. Moves cells so that regions border the same neighbours as on the real
   map, and fills gaps inside the mosaic.
7. Places detached landmasses next to the mainland and puts all cells on
   one grid.

Each step is also exported on its own: `project_regions()`,
`detect_landmasses()`, `deform_regions()`, `generate_cells()`,
`fit_quota()`, `enforce_topology()` and `fill_gaps()`.

## Input data

Any `sf` polygons work, but three things about them change the result.

Detail costs time. Runtime grows with the number of polygon parts, not the
number of regions. Small islands and coastline fragments that are too small
to get a cell can be dropped first without changing the map.

Rows that share a value in `region_col` are treated as one region, so you
don't need to dissolve sub-units first. If the weight is repeated on every
row of a region, it is counted once rather than summed. See `values_agg`.

Detached territories are placed next to the mainland, not in an inset. A
far-off one can end up well outside the main map, so filter it out first if
you don't want it.

## Limitations

- Runtime depends on the number of vertices, not on `n`. Generalised
  boundaries are quick (Greece at `n = 1000` takes about three seconds), but
  full-resolution official boundaries can take minutes. `simplify = TRUE`
  helps with those and isn't worth it otherwise.
- A region enclosed by another makes the neighbour repair hard. It may use
  up `topology_passes` without settling.
- The Python original is much faster on the same input.

## References

- Dougenik, J. A., Chrisman, N. R., & Niemeyer, D. R. (1985). An algorithm
  to construct continuous area cartograms. *The Professional Geographer*,
  37(1), 75–81.
- Gastner, M. T., & Newman, M. E. J. (2004). Diffusion-based method for
  producing density-equalizing maps. *PNAS*, 101(20), 7499–7504.
- Gastner, M. T., Seguy, V., & More, P. (2018). Fast flow-based algorithm
  for creating density-equalizing map projections. *PNAS*, 115(10),
  E2156–E2164.
- Original Python implementation:
  [databites-tech/popgrid](https://github.com/databites-tech/popgrid) by
  [Josep Ferrer](https://github.com/iamjosepferrer).
