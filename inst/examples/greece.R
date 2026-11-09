# Population grid cartogram for Greece.
# Requires: rnaturalearth, rnaturalearthdata, ggplot2, ggrepel, RColorBrewer.

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

# Five ColorBrewer Set2 colours, assigned by hand so that no two
# neighbouring regions share one. Labels name the regions.
set2 <- RColorBrewer::brewer.pal(5, "Set2")

region_colour <- c(
  "Thessaly"                  = set2[1],
  "East Macedonia and Thrace" = set2[1],
  "Peloponnese"               = set2[1],
  "Central Greece"            = set2[2],
  "Epirus"                    = set2[2],
  "South Aegean"              = set2[2],
  "West Greece"               = set2[3],
  "West Macedonia"            = set2[3],
  "Crete"                     = set2[3],
  "Central Macedonia"         = set2[4],
  "Attica"                    = set2[4],
  "Ionian Islands"            = set2[5],
  "North Aegean"              = set2[5]
)

blocks <- aggregate(cells["region_name"], by = list(region = cells$region_name), FUN = \(x) x[1])

ggplot(cells) +
  geom_sf(aes(fill = region_name), color = "white", linewidth = 0.1) +
  ggrepel::geom_text_repel(data = blocks, aes(label = region, geometry = geometry), stat = "sf_coordinates", fun.geometry = sf::st_point_on_surface, size = 3, bg.color = "white", bg.r = 0.15, min.segment.length = 0, seed = 1) +
  scale_fill_manual(values = region_colour, guide = "none") +
  theme_void() +
  labs(
    title    = "Greece - population grid",
    subtitle = sprintf("1 cell ~ %.0f people  |  %d cells total",
                       attr(cells, "value_per_cell"), nrow(cells))
  )
