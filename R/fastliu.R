#' Fit a Liu Regression Coefficients Path
#'
#' \code{liureg} fits coefficients paths for Liu regression models
#'    over a grid of values for the regularization (biasing)
#'    parameter \code{lambda}. The returned object is of class \code{liureg}.
#'
#' @param X The design matrix of features. \code{fastliu}
#' standardizes the data and includes an intercept term by default.
#' @param y The response vector.
#' @param lambda User-specified values of \code{lambda}. The default
#'    value is 1, which corresponds to the least squares estimator.
#'    A \code{lambda} sequence can be entered to generate multiple models.
#' @param scale Scaling type of the design matrix. \code{"ulength"} corresponds
#'    to unit-length scaling. In this scaling the scaled design matrix is
#'    in the form of a correlation matrix. \code{"unormal"} scales the features
#'    to have unit variance (using \eqn{1/n} rather than \eqn{1/(n-1)} formula). \code{"none"}
#'    does not make scaling and computations are done on centered features.
#' @param ... Not used in this implementation.
#'
#' @details
#' The sequence of Liu regression models indexed by the tuning parameter.
#' \eqn{\lambda} are obtained by
#' \deqn{\hat{\boldsymbol{\beta}}^{liu}\left(\lambda\right)=
#' \left(\mathbf{X}^{T}\mathbf{X}+\mathbf{I}_{p}\right)^{-1}
#' \left(\mathbf{X}^{T}\mathbf{y}+\lambda\hat{\boldsymbol{\beta}}^{ls}\right),}
#' where \eqn{\hat{\boldsymbol{\beta}}^{ls}} is the ordinary least squares
#' estimator.To obtain the models, the singular value decomposition (SVD)
#' of the matrix \eqn{\mathbf{X}} is used. This SVD is done once and
#' is used to generate all models.
#'
#'
#' Explanatory variables in the design matrix are always centered
#' before fitting a model in the \code{fastliu} package.
#' For scaling, two options are possible:
#' unit-length and unit-normal scaling. In unit-length scaling,
#' the matrix of explanatory variables has correlation form.
#' In unit-normal scaling, the explanatory variables have zero
#' mean and unit variance.
#' Both Coefficient estimates based on the scaled data and
#' in original scale are presented.
#' The intercept of the model is not penalized and computed by
#' \eqn{\bar{y}-\bar{X}\boldsymbol{\hat{\beta}}_1}, where \eqn{\bar{X}}
#' is the row vector of the explanatory variables and \eqn{\boldsymbol{\hat{\beta}}_1}
#' is computed based on centered design matrix.
#'
#' The returned \code{liureg} object
#' is used for statistical testing of Liu coefficients,
#' plotting method and computing the Liu regression related statistics.
#'
#' @return Fitted Liu regression object with the class of \code{liureg}
#' @author Murat Genç and Ömer Özbilen
#' @export
#'
#' @seealso [coef()], [predict()], [summary()], [pressliu()], [residuals()]
#' @examples
#' data("Hitters")
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' lam <- seq(0, 1, 0.05)
#' liu.mod <- liureg(X, y, lam)
liureg <- function(X, y, lambda = 1, scale = c("ulength","unormal","none"), ...){
  scale <- match.arg(scale)
  # credit to the package ncvreg
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if (!is.null(ncol(y)) && ncol(y) > 1) stop("y should be a vector of responses, not a matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }

  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to liureg function.", call.=FALSE)

  if(!(inherits(lambda,"numeric") | inherits(lambda,"matrix"))) {
    stop("lambda must be a sequence or a matrix with one column.")
  }

  if(inherits(lambda,"matrix")){
    if(ncol(lambda) > 1){
      stop("lambda must be a sequence or a matrix with one column.")
    }
  }

  lambda <- drop(lambda)
  l <- length(lambda)

  n <- nrow(X)
  p <- ncol(X)

  if(is.null(colnames(X))) colnames(X) <- paste0("V", 1:p)
  if(is.null(rownames(X))) rownames(X) <- 1:n
  if(is.null(names(lambda))) names(lambda) <- paste0("lam", 1:l)

  fit <- liuregcpp(X, y, lambda, scale)
  fit$cnames <- colnames(X)
  fit$rnames <- rownames(X)
  fit$lnames <- names(lambda)
  res <- structure(fit, class = "liureg")
  res
}
