#' Hat matrix of Liu Regression
#'
#' For each value of the regularization parameter lambda,
#'    \code{hatliu} returns the hat matrix of Liu regression.
#'    The hat matrix for Liu regression is computed using the formula
#'    \eqn{\mathbf{H}=\mathbf{X}\left(\mathbf{X}^{T}\mathbf{X}+\mathbf{I}_{p}\right)^{-1}
#'    \left(\mathbf{X}^{T}\mathbf{X}+\lambda\mathbf{I}_{p}\right)\left(\mathbf{X}^{T}
#'    \mathbf{X}\right)^{-1}\mathbf{X}^{T}.}

#'
#' @param obj A \code{liureg} object.
#'
#' @return The returned object is a list of matrices whose elements are
#'    the hat matrices for the values of the
#'    \code{lambda} regularization parameter.
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
#' # Hat matrix list
#' hatlist <- hatliu(liu.mod)
#' # Hat matrix for third regularization parameter
#' hatlist[[3]]
hatliu <- function(obj){
  H <- hatcpp(obj)
  names(H) <- obj$lnames
  for(i in 1:length(H)){
    colnames(H[[i]]) <- rownames(H[[i]]) <- obj$rnames
  }
  H
}
