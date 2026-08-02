# Population grid cartogram for Greece.
# Requires: rnaturalearth, rnaturalearthdata, ggplot2.

library(popgridR)
library(ggplot2)

# Natural Earth's admin-1 `name` column already IS the 13 Greek peripheries.
regions <- load_country_regions("GRC")

# Mount Athos is counted separately by ELSTAT and far too small to earn a
# cell at n = 1000 — drop it for a clean 13-region map.
regions <- regions[regions$region_name != "Mount Athos", ]

# 2021 ELSTAT census
population <- c(
  "Attica"                    = 3814064,
  "Central Macedonia"         = 1795669,
  "West Greece"               =  643449,
  "Thessaly"                  =  688255,
  "Crete"                     =  624408,
  "East Macedonia and Thrace" =  562201,
  "Peloponnese"               =  539535,
  "Central Greece"            =  508254,
  "South Aegean"              =  309015,
  "Epirus"                    =  319991,
  "West Macedonia"            =  254595,
  "North Aegean"              =  194943,
  "Ionian Islands"            =  204532
)

cells <- popgrid(regions, values = population, n = 1000)

ggplot(cells) +
  geom_sf(aes(fill = region_name), color = NA) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank()) +
  labs(
    title    = "Greece - population grid",
    subtitle = sprintf("1 cell ~ %.0f people  |  %d cells total",
                       attr(cells, "value_per_cell"), nrow(cells))
  )
