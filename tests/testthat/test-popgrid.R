# Synthetic end-to-end tests. No network, no cartogramR (method = "damped").

make_regions <- function() {
  # Four adjacent 100 km squares in a row, projected CRS
  r <- function(x0) sf::st_polygon(list(cbind(
    c(x0, x0 + 1, x0 + 1, x0, x0) * 1e5,
    c(0, 0, 1e5, 1e5, 0)
  )))
  sf::st_sf(
    name     = c("A", "B", "C", "D"),
    w        = c(40, 30, 20, 10),
    geometry = sf::st_sfc(lapply(0:3, r), crs = 3857)
  )
}

test_that("popgrid returns exact per-region quotas", {
  regions <- make_regions()
  cells <- popgrid(regions, values = "w", region_col = "name",
                   n = 100L, method = "damped", verbose = FALSE)

  expect_s3_class(cells, "sf")
  counts <- table(cells$region_name)
  quotas <- attr(cells, "quotas")
  expect_equal(sum(counts), 100L)
  expect_equal(as.integer(counts[names(quotas)]), unname(quotas))
  # user's region column is present in the output
  expect_true("name" %in% names(cells))
})

test_that("values as named vector matches values as column", {
  regions <- make_regions()
  c1 <- popgrid(regions, values = "w", region_col = "name",
                n = 60L, method = "damped", verbose = FALSE)
  c2 <- popgrid(regions, values = c(A = 40, B = 30, C = 20, D = 10),
                region_col = "name", n = 60L, method = "damped",
                verbose = FALSE)
  expect_equal(attr(c1, "quotas"), attr(c2, "quotas"))
})

test_that("value_per_cell overrides n", {
  regions <- make_regions()
  cells <- popgrid(regions, values = "w", region_col = "name",
                   value_per_cell = 2, method = "damped", verbose = FALSE)
  # sum(w) = 100, value_per_cell = 2 -> n = 50
  expect_equal(attr(cells, "n"), 50L)
  expect_equal(sum(table(cells$region_name)), 50L)
})

test_that("rows sharing a region name pool into one quota", {
  regions <- make_regions()
  regions$name <- c("A", "A", "B", "B")   # two rows per region
  cells <- popgrid(regions, values = "w", region_col = "name",
                   n = 50L, method = "damped", verbose = FALSE)
  quotas <- attr(cells, "quotas")
  expect_equal(sort(names(quotas)), c("A", "B"))
  expect_equal(sum(quotas), 50L)
  # A pools 40 + 30 = 70, B pools 20 + 10 = 30
  expect_equal(unname(quotas["A"]), 35L)
})

test_that("missing values fail loudly", {
  regions <- make_regions()
  expect_error(
    popgrid(regions, values = c(A = 1, B = 2, C = 3), region_col = "name",
            n = 10L, method = "damped", verbose = FALSE),
    "No value supplied for: D"
  )
  expect_error(
    popgrid(regions, values = "nope", region_col = "name", n = 10L),
    "not a column in"
  )
  expect_error(popgrid(regions$name, values = 1:4), "must be an sf data frame")
})

test_that("a value constant within multi-row regions is not summed", {
  regions <- make_regions()
  regions$name  <- c("A", "A", "B", "B")
  regions$share <- c(0.7, 0.7, 0.3, 0.3)   # per-region share repeated per row

  auto <- popgrid(regions, values = "share", region_col = "name",
                  n = 100L, method = "damped", verbose = FALSE)
  expect_equal(unname(attr(auto, "quotas")[["A"]]), 70L)

  # forcing sum double-counts both regions; here they scale equally
  summed <- popgrid(regions, values = "share", region_col = "name",
                    n = 100L, values_agg = "sum", method = "damped",
                    verbose = FALSE)
  expect_equal(unname(attr(summed, "quotas")[["A"]]), 70L)

  # per-row quantities that differ within a region are still summed
  regions$pop <- c(40, 30, 20, 10)
  q <- attr(popgrid(regions, values = "pop", region_col = "name", n = 100L,
                    method = "damped", verbose = FALSE), "quotas")
  expect_equal(unname(q[["A"]]), 70L)
})

test_that("empty geometries are dropped", {
  regions <- make_regions()
  empty   <- regions[1, ]
  sf::st_geometry(empty) <- sf::st_sfc(sf::st_polygon(), crs = 3857)
  empty$name <- "E"; empty$w <- 5
  regions <- rbind(regions, empty)

  cells <- popgrid(regions, values = "w", region_col = "name",
                   n = 50L, method = "damped", verbose = FALSE)
  expect_false("E" %in% cells$region_name)
})

test_that("method = none tiles undeformed geometry", {
  regions <- make_regions()
  cells <- popgrid(regions, values = "w", region_col = "name",
                   n = 40L, method = "none", verbose = FALSE)
  expect_equal(sum(table(cells$region_name)), 40L)
})
