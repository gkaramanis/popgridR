test_that("cell_size_for_n is sqrt(area / n)", {
  expect_equal(cell_size_for_n(100, 4), 5)
  expect_error(cell_size_for_n(100, 0))
  expect_error(cell_size_for_n(0, 4))
})

test_that("generate_cells tiles a square exactly", {
  sq <- sf::st_polygon(list(cbind(c(0, 10, 10, 0, 0), c(0, 0, 10, 10, 0)) * 1e4))
  regions <- sf::st_sf(
    region_id   = 0L,
    region_name = "A",
    geometry    = sf::st_sfc(sq, crs = 3857)
  )
  cells <- generate_cells(sq, regions, cell_size = 1e4, landmass_id = 0L)
  expect_s3_class(cells, "sf")
  expect_equal(nrow(cells), 100L)
  expect_true(all(cells$region_name == "A"))
})
