# Integer seat allocation.

#' Hamilton (largest-remainder) quota allocation
#'
#' Allocates `n` integer seats from floating-point weights, guaranteeing
#' `sum(result) == n`.
#'
#' @param weights Named numeric vector of weights.
#' @param n       Total integer seats to allocate.
#' @return Named integer vector.
#' @export
hamilton_quotas <- function(weights, n) {
  nms     <- names(weights)           # pmax() drops names - capture first
  weights <- pmax(as.numeric(weights), 0)
  total_w <- sum(weights)
  if (total_w == 0) return(setNames(rep(0L, length(weights)), nms))

  raw     <- weights / total_w * n
  floors  <- floor(raw)
  deficit <- n - sum(floors)
  remaind <- raw - floors

  order_r <- order(-remaind)
  floors[order_r[seq_len(as.integer(deficit))]] <-
    floors[order_r[seq_len(as.integer(deficit))]] + 1L

  setNames(as.integer(floors), nms)
}
