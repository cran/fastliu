#' Print Method for Liu Regression Statistics
#'
#' \code{statliu} computes the statistics related to the Liu regression.
#'
#' @param x A \code{statliu} object
#' @param digits Number of decimal places in the data frame of Liu regression statistics.
#' @param ... Other parameters related to \code{print}.
#'
#' @return The return object is the statistics relatec to the Liu regression.
#' @author Murat Genç
#' @method print statliu
#' @export
#'
#' @seealso [liureg()], [summary()], [pressliu()], [residuals()]
#'
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' stats <- statliu(liu.mod)
#' print(stats)
print.statliu <- function(x, digits = 5, ...) {
  cat("\nLiu Regression Statistics:\n\n")
  res <- do.call(cbind.data.frame, x)
  row.names(res) <- attr(x, "row.names")
  print(round(res, digits), ...)
}
