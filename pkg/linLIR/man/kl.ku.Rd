\name{kl.ku}
\alias{kl.ku}
\title{
Finding \code{k.l} and \code{k.u}
}
\description{
Internal function to determine the quantities \code{k.l} and \code{k.u} on the basis of \code{n}, \code{p}, \code{bet}. This function is used within the (internal) functions \code{gen.lms}, \code{undom.a} and \code{undom.para}.
}
\usage{
kl.ku(n, p = 0.5, bet, epsilon = 0)
}
\arguments{
  \item{n}{Number of observations.
}
  \item{p}{
Quantile of the abolute residuals' distribution to be used as loss function in the LIR analysis. (0.5 corresponds to the median.)
}
  \item{bet}{
Cutoff-point for the normalized profile likelihood function.
}
  \item{epsilon}{
Fraction of coarsening errors considered.
}
}
\value{
A vector with 2 elements, \code{k.l} and \code{k.u}.
}
\references{
M. Cattaneo, A. Wiencierz (2012c). On the implementation of LIR: the case of simple linear regression with interval data. Technical Report No. 127. Department of Statistics. LMU Munich.

A. Wiencierz, M. Cattaneo (2012b). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. In: R. Kruse et al. (Eds.). Advances in Intelligent Systems and Computing. Vol. 190. Springer. pp. 293-301.

M. Cattaneo, A. Wiencierz (2012a). Likelihood-based Imprecise Regression. International Journal of Approximate Reasoning. Vol. 53. pp. 1137-1154.
}
\seealso{
\code{\link{s.linlir}},
\code{\link{gen.lms}},
\code{\link{undom.para}}
}
