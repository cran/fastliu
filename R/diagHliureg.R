#' Diagonal Elements of the Hat Matrix
#'
#' For each value of the regularization parameter lambda,
#'    \code{diagHliu} returns the diagonal elements of the hat matrix.
#'    Unlike the \code{hatliu} function, only the diagonal
#'    elements of the hat matrix are calculated, thus the
#'    computation of diagonal elements is faster than \code{hatliu}.
#'
#' @param obj A \code{liureg} object
#'
#' @return The returned object is a matrix whose columns are the
#'    diagonal elements of the hat matrix for each value of the
#'    lambda regularization parameter.
#' @author Murat Genç
#' @export
#'
#' @seealso [liureg()], [summary()], [pressliu()], [residuals()]
#'
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' diagHliu(liu.mod)
diagHliu <- function(obj){
  diagH <- diagHcpp(obj)
  rownames(diagH) <- obj$rnames
  colnames(diagH) <- obj$lnames
  diagH
}



