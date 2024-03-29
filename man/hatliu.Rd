% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hatliureg.R
\name{hatliu}
\alias{hatliu}
\title{Hat matrix of Liu Regression}
\usage{
hatliu(obj)
}
\arguments{
\item{obj}{A \code{liureg} object.}
}
\value{
The returned object is a list of matrices whose elements are
the hat matrices for the values of the
\code{lambda} regularization parameter.
}
\description{
For each value of the regularization parameter lambda,
\code{hatliu} returns the hat matrix of Liu regression.
The hat matrix for Liu regression is computed using the formula
\eqn{\mathbf{H}=\mathbf{X}\left(\mathbf{X}^{T}\mathbf{X}+\mathbf{I}_{p}\right)^{-1}
   \left(\mathbf{X}^{T}\mathbf{X}+\lambda\mathbf{I}_{p}\right)\left(\mathbf{X}^{T}
   \mathbf{X}\right)^{-1}\mathbf{X}^{T}.}
}
\examples{
data("Hitters")
Hitters <- na.omit(Hitters)
X <- model.matrix(Salary ~ ., Hitters)[, -1]
y <- Hitters$Salary
lam <- seq(0, 1, 0.01)
liu.mod <- liureg(X, y, lam)
# Hat matrix list
hatlist <- hatliu(liu.mod)
# Hat matrix for third regularization parameter
hatlist[[3]]
}
\seealso{
\code{\link[=liureg]{liureg()}}, \code{\link[=summary]{summary()}}, \code{\link[=pressliu]{pressliu()}}, \code{\link[=residuals]{residuals()}}
}
\author{
Murat Genç
}
