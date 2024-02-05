#' Covariance matrix based on a fitted liureg object.
#'
#' For a scalar or vector tuning parameter lambda,
#'    the \code{covliureg} computes the covariance matrix
#'    for the estimates of a Liu regression model.
#'
#' @param obj A \code{liureg} object.
#' @return The returned object is a list of the matrix of estimated covariances.
#' @author Murat Genç and Ömer Özbilen
#' @export
#'
#' @seealso [liureg()], [coef()], [predict()], [summary()], [pressliu()], [residuals()]
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.01)
#' liu.mod <- liureg(X, y, lam)
#' # List of covariance matrices for 101 lambda values
#' cov.mat <- covliu(liu.mod)
#' print(cov.mat$lam1)
covliu <- function(obj){
  covliu <- covliucpp(obj)
  names(covliu) <- obj$lnames
  for(i in 1:length(covliu)){
    rownames(covliu[[i]]) <- colnames(covliu[[i]]) <- obj$cnames
  }
  covliu
}




















































