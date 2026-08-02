#' @keywords internal
#' @importFrom stats median setNames
"_PACKAGE"

`%||%` <- function(x, y) if (!is.null(x)) x else y
