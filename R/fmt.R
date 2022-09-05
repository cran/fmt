#' Variance estimation of FMT method (Fully Moderated T-statistic)
#'
#' This function computes posterior residual variances to be used
#' in the denominator of a moderated t-statistic from a linear model
#' analysis of gene expression data. It is an extension of the moderated
#' t-statistic originally proposed by Smyth (Statistical Applications in
#' Genetics and Molecular Biology, 2004;3:Article3). LOESS local regression
#' and empirical Bayesian method are used to estimate gene specific prior
#' degrees of freedom and prior variance based on average gene intensity levels.
#' The posterior residual variance in the denominator is a weighted average of
#' prior and residual variance and the weights are prior degrees of freedom
#' and residual variance degrees of freedom. The degrees of freedom of the
#' moderated t-statistic is simply the sum of prior and residual variance
#' degrees of freedom.
#'
#' @param Amean     average log intensity levels of all genes
#' @param sigmasq   residual variances of all genes
#' @param df        degrees of freedom for sigmasq
#' @param span1     span parameter in LOESS smoothing function, default is 0.5
#' @param span2     span parameter in LOESS smoothing function, default is 0.95
#' @param iter1     iteration number in LOESS smoothing function, default is 4
#' @param iter2     iteration number in LOESS smoothing function, default is 4
#' @param b         number of genes on either side of moving average
#'   window when calculating variance of log residual variances,
#'   default is 20
#' @return A data frame with the following components:
#'
#'  'df.prior' the estimated prior degrees of freedom.
#'
#'  'df.post' the estimated posterior degrees of freedom.
#'
#'  's2.prior' the estimated prior variance.
#'
#'  's2.post' the estimated posterior variance.
#'
#'  'Ameansort' intermediate result.
#'
#'  'eg' intermediate result.
#'
#'  'egpred' intermediate result.
#'
#'  'MAvar' intermediate result.
#'
#'  'tri.d0' intermediate result.
#' @examples
#' ## Simulate gene expression data for 1000 genes and 10 samples in two groups.
#' exp <- rnorm(1000,8,2)
#' sd <- 0.5*sqrt(4/rchisq(1000, df=7))
#' y <- matrix(rnorm(1000*10, exp, sd),1000,10)
#' rownames(y) <- paste("Gene",1:1000)
#' design <- cbind(Grp1=1, Grp2vs1=c(0,0,0,0,0,1,1,1,1,1))
#'
#' ## limma fit
#' fit <- lmFit(y,design)
#'
#' ## fmt fit
#' fmt.fit <- fmt(fit$Amean, fit$sigma, fit$df.residual)
#'
#' @importFrom limma loessFit trigammaInverse
#' @export

fmt <- function(Amean, sigmasq, df, span1=0.5, span2=0.95, iter1=4, iter2=4, b=20) {
	eg <- log(sigmasq) - digamma(df/2) + log(df/2)
	egpred <- loessFit(eg, Amean, iterations=iter1, span=span1)$fitted
  N <- length(Amean)
	mat <- cbind(Amean,(eg - egpred)^2)
	order <- sort(Amean,index.return=TRUE)$ix
	matsort <- mat[order,]
	MAvar <- NULL
	for (i in 1:b) {
 	      MAvar[i] <- mean(matsort[1:(i+b),2])
	}
	for (i in (b+1):(N-b)) {
	      MAvar[i] <- mean(matsort[(i-b):(i+b),2])
	}
	for (i in (N-b+1):N) {
	      MAvar[i] <- mean(matsort[(i-b):N,2])
	}
	tri.d0 <- loessFit(MAvar, matsort[,1], iterations=iter2, span=span2)$fitted - trigamma(df/2)
	tri.d0[tri.d0<=0] <- 0.04
	df.prior.sort <- 2*trigammaInverse(tri.d0)
  df.prior <- df.prior.sort[sort(order,index.return=TRUE)$ix]
	df.post <- df.prior + df
  s2.prior <- exp(egpred + digamma((df.prior)/2) - log(df.prior/2))
	s2.post <- (df.prior*s2.prior + df*sigmasq) / df.post
	output <- data.frame(df = df,
	                     df.prior = df.prior,
	                     df.post = df.post,
	                     s2.prior = s2.prior,
	                     s2.post = s2.post,
	                     Amean = Amean,
	                     Ameansort = matsort[,1],
	                     eg = eg,
	                     egpred = egpred,
	                     MAvar = MAvar,
	                     tri.d0 = tri.d0)
	return(output)
}






#' Plot fitting results from fmt function
#'
#' This function provides fitting plots of fmt function.
#'
#' @param x     return list from fmt function
#' @param type  type of plots
#' @param ...   extra parameters for plot function
#' @return No return value.
#'
#' @importFrom graphics par points
#' @export

plotFMT <- function(x, type, ...) {
    if (type=="all") {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mfrow = c(2,2))
      plot(x$Amean, x$df.prior, cex=0.5, col="black",
           xlab="Average Log Intensity", ylab="Prior Degrees of Freedom",
           type="p", ...)
      plot(x$Ameansort, x$MAvar, cex=0.5,
           xlab="Average Log Intensity", ylab="Variance of Log Variance",
           pch=1, ...)
      points(x$Ameansort, x$tri.d0+trigamma(x$df/2), cex=0.5, col="red")
      plot(x$Amean, x$s2.prior, cex=0.5, col="black",
           xlab="Average Log Intensity", ylab="Prior Variance",
           type="p", ...)
      plot(x$Amean, x$eg, cex=0.5,
           xlab="Average Log Intensity", ylab="Log Variance",
           pch=1, ...)
      points(x$Amean, x$egpred, col="red", cex=0.5)
      par(mfrow = c(1,1))
    }
    if (type=="priordf") {
      plot(x$Amean, x$df.prior, cex=0.5, col="black",
           xlab="Average Log Intensity", ylab="Prior Degrees of Freedom",
           type="p", ...)
    }
    if (type=="varoflogvar") {
      plot(x$Ameansort, x$MAvar, cex=0.5,
           xlab="Average Log Intensity", ylab="Variance of Log Variance",
           pch=1, ...)
      points(x$Ameansort, x$tri.d0+trigamma(x$df/2), cex=0.5, col="red")
    }
  if (type=="priorvar") {
     plot(x$Amean, x$s2.prior, cex=0.5, col="black",
          xlab="Average Log Intensity", ylab="Prior Variance",
          type="p", ...)
  }
  if (type=="logvar") {
    plot(x$Amean, x$eg, cex=0.5,
         xlab="Average Log Intensity", ylab="Log Variance",
         pch=1, ...)
    points(x$Amean, x$egpred, col="red", cex=0.5)
  }
}



