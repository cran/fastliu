% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagHliureg.R
\name{diagHliu}
\alias{diagHliu}
\title{Diagonal Elements of the Hat Matrix}
\usage{
diagHliu(obj)
}
\arguments{
\item{obj}{A \code{liureg} object}
}
\value{
The returned object is a matrix whose columns are the
diagonal elements of the hat matrix for each value of the
lambda regularization parameter.
}
\description{
For each value of the regularization parameter lambda,
\code{diagHliu} returns the diagonal elements of the hat matrix.
Unlike the \code{hatliu} function, only the diagonal
elements of the hat matrix are calculated, thus the
computation of diagonal elements is faster than \code{hatliu}.
}
\examples{
data("Hitters")
Hitters <- na.omit(Hitters)
X <- model.matrix(Salary ~ ., Hitters)[, -1]
y <- Hitters$Salary
lam <- seq(0, 1, 0.01)
liu.mod <- liureg(X, y, lam)
diagHliu(liu.mod)
}
\seealso{
\code{\link[=liureg]{liureg()}}, \code{\link[=summary]{summary()}}, \code{\link[=pressliu]{pressliu()}}, \code{\link[=residuals]{residuals()}}
}
\author{
Murat Genç
}
