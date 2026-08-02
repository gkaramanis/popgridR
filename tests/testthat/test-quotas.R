test_that("hamilton_quotas sums to n and keeps names", {
  w <- c(a = 10, b = 20, c = 70)
  q <- hamilton_quotas(w, 10L)
  expect_equal(sum(q), 10L)
  expect_named(q, c("a", "b", "c"))
  expect_equal(q[["c"]], 7L)
})

test_that("hamilton_quotas handles zero weights", {
  q <- hamilton_quotas(c(a = 0, b = 0), 5L)
  expect_equal(sum(q), 0L)
})

test_that("hamilton_quotas resolves remainders deterministically", {
  q <- hamilton_quotas(c(a = 1, b = 1, c = 1), 10L)
  expect_equal(sum(q), 10L)
  expect_true(all(q >= 3L))
})
