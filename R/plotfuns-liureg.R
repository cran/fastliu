plt.liucoefpath <- function(obj, abline=TRUE, legend.pos="topright", ...) {
  lambda <- obj$lambda
  coefs <- obj$coefliu
  p <- nrow(coefs)
  l <- ncol(coefs)
  cols <- c("#000000", "#0082c8", "#e6194b", "#3cb44b", "#ffe119", "#f58231",
            "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#0000aa",
            "#aa6e28", "#800000", "#aaffc3", "#808000", "#ffd8b1",
            "#000080", "#808080")
  cols <- rep(cols, p %/% length(cols) + 1)
  lends <- rownames(coefs)
  if(is.null(lends)) lends <- paste0("V",1:p)

  LIUSTATS <- statliu(obj)
  MSE <- LIUSTATS$MSE
  msemin <- min(MSE)
  lambda_min <- lambda[which.min(MSE)]

  xlab <- expression(lambda)
  ylab <- expression(hat(beta))
  main <- "Liu Coefficient Path"


  if (length(lambda) == 1) {
    plot(x = rep(lambda,length(coefs)), y = drop(coefs),
         xlab = xlab,ylab = "", main = main,
         col = "#4682b4", pch = 19, ...)
    mtext(ylab, 2, 3, las=1)
    legend("topright",
           legend = bquote(lambda*" = "*.(lambda)),
           cex = .85,
           col = "#4682b4", pch = 19,
           x.intersp=0.5,
           bty="n",
           bg="transparent")
  } else {
    matplot(x = lambda, y = t(coefs),
            xlab = xlab, ylab = "", main = main,
            col = cols[1:p],
            type = "l", lty = 1:p, lwd = 1, ...)
    mtext(ylab, 2, 3, las=1)

    legend(legend.pos,
           legend = lends,
           lwd = 1,
           cex = .75,
           bty = "n",
           bg="transparent",
           col = cols[1:p],
           x.intersp=0.5, y.intersp=0.9,
           merge=TRUE)
  }

  if (abline) {
    abline(v = lambda_min,
           lty = 2,
           col = "#4682b4")
    text(lambda_min,
         max(coefs),
         bquote("min(MSE) = "*.(round(msemin,3))),
         col="#4682b4",
         cex=0.85,
         pos=4)
    text(lambda_min,
         max(coefs),
         bquote("at "*lambda*" = "*.(lambda_min)),
         col="#4682b4",
         cex=0.85,
         adj=c(-0.22,1.75))
  }
}



plt.biasliureg <- function(obj, abline=TRUE, legend.pos = "topright", ...){
  lambda <- obj$lambda
  LIUSTATS <- statliu(obj)
  BIAS2 <- LIUSTATS$BIAS2
  VAR <- LIUSTATS$VAR
  MSE <- LIUSTATS$MSE
  MSEMIN <- min(MSE)
  lambda_min <- lambda[which.min(MSE)]
  cols <- cbind("#b4464b", "#46b478", "#4682b4")
  yvals <- cbind(VAR, BIAS2, MSE)

  main <- "Bias-Variance Trade-off"
  xlab <- expression(lambda)

  if(length(lambda)==1){
    plot(x=rep(lambda, length(yvals)), y=yvals,
         main=main, xlab=xlab, ylab=" ",
         pch=19, col=cols)
    legend("topright",
           legend=c("Var", expression("Bias"^"2"), "MSE"),
           cex=0.75,
           pch=19, col=cols,
           x.intersp=0.5, y.intersp=0.9,
           bg="transparent",
           bty = "n")
  } else{
    matplot(x=lambda, y=yvals,
            main=main, xlab=xlab, ylab=" ",
            col=cols, lwd = 1, lty=c(2,4,1), type="l")

    legend(legend.pos,
           legend=c("Var", expression("Bias"^"2"), "MSE"),
           col=cols,
           cex=0.75,
           lwd = 1,lty=c(2,4,1),
           x.intersp=0.5, y.intersp=0.9,
           bg = "transparent",
           bty = "n")
  }

  if(abline){
    abline(v=lambda_min, lty=2)
    text(lambda_min,
         max(LIUSTATS$MSE),
         bquote("min(MSE) = "*.(round(MSEMIN,3))),
         col="#4682b4",
         cex=0.85,
         pos=4)
    text(lambda_min,
         max(LIUSTATS$MSE),
         bquote("at "*lambda*" = "*.(lambda_min)),
         col="#4682b4",
         cex=0.85,
         adj=c(-0.2,1.75))
  }
}


plt.infoliureg <- function(obj, abline=TRUE, legend.pos="topright", ...){

  DF <- colSums(diagHliu(obj))
  INFO <- infoliu(obj)
  LIUSTATS <- statliu(obj)
  lambda <- obj$lambda

  aic <- INFO$AIC
  bic <- INFO$BIC
  cols <- c("#46b478", "#4682b4")
  mse <- LIUSTATS$MSE
  msemin <- min(mse)
  df_min <- DF[which.min(mse)]
  mselect <- cbind(aic, bic)

  main <- "Model Selection Criteria"
  xlab <- "DF"
  ylab <- ""
  legendlab <- c("AIC", "BIC")

  if(length(lambda)==1){
    plot(x=rep(DF, length(mselect)), y=mselect,
         main=main, xlab=xlab, ylab=ylab,
         col=cols, pch=19)
    legend("topright",
           legend=legendlab,
           cex=0.75,
           pch=19, col=cols,
           x.intersp=0.5, y.intersp=0.9,
           bg="transparent",
           bty = "n")
  }
  else{
    matplot(DF, mselect,
            main=main, xlab=xlab, ylab=ylab,
            type="l", col=cols, lwd = 1, lty=c(1,2))
    legend(legend.pos,
           legend=legendlab,
           col=cols,
           cex=0.75,
           lwd = 1,lty=c(2,4,1),
           x.intersp=0.5, y.intersp=0.9,
           bg = "transparent",
           bty = "n")
  }


  if(abline){

    abline(v=df_min, lty=2)
    text(df_min,
         max(mselect),
         bquote("min(MSE) = "*.(round(msemin,3))),
         col="#4682b4",
         cex=0.85,
         pos=4)
    text(df_min,
         max(mselect),
         bquote("at df = "*.(round(df_min, 3))),
         col="#4682b4",
         cex=0.85,
         adj=c(-0.13,1.75))
  }

}


#' Plot method for liureg objects
#'
#' @param x A \code{liureg} object.
#' @param type What to plot on the vertical axis.  \code{coefpath} plots the
#' coefficient path of the Liu regression; \code{biasvar} generates a bias-variaance plot,
#' \code{info} plots the information criteria corresponding the regularization parameter values.
#' @param ... Other graphical parameters to \code{plot}.
#'
#' @method plot liureg
#' @export
#'
#' @seealso [liureg()], [predict()], [summary()]
#' @return No return value.
#' @author Murat Genç
#' @examples
#' Hitters <- na.omit(Hitters)
#' X <- model.matrix(Salary ~ ., Hitters)[, -1]
#' y <- Hitters$Salary
#' liu.mod <- liureg(X, y, seq(0, 1, 0.01))
#'
#' # Liu coefficient paths
#' plot(liu.mod)
#'
#' # Bias-variance trade-off
#' plot(liu.mod, type="biasvar")
plot.liureg <- function(x, type=c("coefpath","biasvar","info"), ...){

  type <- match.arg(type)

  if(type == "coefpath"){
    plt.liucoefpath(x, ...)
  } else if(type == "biasvar"){
    plt.biasliureg(x, ...)
  } else {
    plt.infoliureg(x, ...)
  }
}

