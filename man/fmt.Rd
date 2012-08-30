\name{fmt}
\alias{fmt}
\alias{fmt.default}
\alias{plot.fmt}
\title{Variance estimation of FMT method (Fully Moderated t-statistic)}
\description{
  This function computes posterior residual variances to be used in the denominator of a moderated t-statistic from a linear model analysis of microarray data.  It is an extension of the moderated t-statistic original proposed by Smyth (2004). LOESS local regression and empirical Bayesian method are used to estimate gene specific prior degrees of freedom and prior variance based on average gene intensity level. The posterior residual variance in the denominator is a weighted average of prior and residual variance and the weights are prior degrees of freedom and residual variance degrees of freedom. The degrees of freedom of the moderated t-statistic is simply the sum of prior and residual variance degrees of freedom.
}
\usage{
fmt(Amean, sigmasq, df, ...)

\method{fmt}{default}(Amean, sigmasq, df, span1 = 0.5, span2 = 0.95, iter1 = 4, iter2 = 4, b = 20, ...)
\method{plot}{fmt}(x, type, ...)
}
\arguments{
  \item{Amean}{ vector of average log intensity levels of all genes. }
  \item{sigmasq}{ vector of residual variances of all genes. }
  \item{df}{ degrees of freedom for sigmasq. }
  \item{span1}{ span parameter in LOESS smoothing function. }
  \item{span2}{ span parameter in LOESS smoothing function. }
  \item{iter1}{ iteration number in LOESS smoothing function. }
  \item{iter2}{ iteration number in LOESS smoothing function. }
  \item{b}{ number of genes on either side of moving average window when calculating variance of log residual variances. }
  \item{x}{ object of class "fmt". }
  \item{type}{ type of plot. }
  \item{\dots}{ not used. }
}
\value{
  \item{df.prior}{ estimated prior degrees of freedom. }
  \item{df.post}{ estimated posterior degrees of freedom. }
  \item{s2.prior}{ estimated prior variance. }
  \item{s2.post}{ estimated posterior variance. }
  \item{Ameansort}{ intermediate result for plotting. }
  \item{eg}{ intermediate result for plotting. }
  \item{egpred}{ intermediate result for plotting. }
  \item{MAvar}{ intermediate result for plotting. }
  \item{tri.d0}{ intermediate result for plotting. }
}
\references{Yu L, Gulati P, Fernandez S, Pennell M, Kirschner L, Jarjoura D. Fully moderated t-statistic for small sample size gene expression arrays. Statistical Applications in Genetics and Molecular Biology, 2011; 10:article 42.}
\author{Lianbo Yu}
\keyword{variance smoothing}
\keyword{fully moderated t-statistic}
