\name{gen.lms}
\alias{gen.lms}
\title{
Finding LRM regression line(s)
}
\description{
Function within the \code{s.linlir}-function that determines the Likelihood-based Region Minimax (LRM) line(s). 
}
\usage{
gen.lms(dat, p = 0.5, bet, epsilon = 0, k.u = 0)
}
\arguments{
  \item{dat}{
An \code{n}x4 \code{data.frame} containing the imprecise data of the analyzed variables. Columns 1 and 2 correspond to the interval-valued observations of the regressor variable, columns 3 and 4 to those of the dependent variable.
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
  \item{k.u}{
As default \code{k.u} is calculated on the basis of \code{p}, \code{bet} and \code{epsilon}.
}
}
\details{
The \code{gen.lms}-function implements the first part of the exact algorithm for the simple linear LIR analysis with interval data developed in M. Cattaneo, A. Wiencierz (2012c). This first part of the algorithm can be seen as a generalization of the basic algorithm for Least Median of Squares Regression (see, e.g., Steele / Steiger (1986) and Rousseeuw / Leroy (1987)). 
}
\value{
A list of two components.
\item{lrm}{A vector (or a matrix) of the intercept and slope parameter values of the LRM line(s) and the value of the p-quantile of the absolute residuals associated with the LRM line(s).}
\item{b.val}{A vector of all slope values considered for the LRM line(s).}
\item{dat.unique}{A data frame of all unique observations in the analyzed data set.}
\item{ind.dat}{A matrix of the indices of all combinations of 2 unique data points.}
}
\references{
M. Cattaneo, A. Wiencierz (2012c). On the implementation of LIR: the case of simple linear regression with interval data. Technical Report No. 127. Department of Statistics. LMU Munich.

A. Wiencierz, M. Cattaneo (2012b). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. In: R. Kruse et al. (Eds.). Advances in Intelligent Systems and Computing. Vol. 190. Springer. pp. 293-301.

M. Cattaneo, A. Wiencierz (2012a). Likelihood-based Imprecise Regression. International Journal of Approximate Reasoning. Vol. 53. pp. 1137-1154.

P. Rousseeuw, A. Leroy (1987). Robust Regression and Outlier Detection. Wiley.

J. Steele, W. Steiger (1986). Algorithms and complexity for least median of squares regression. Discret Appl Math 14. 93-100.
}
\seealso{
\code{\link{s.linlir}},
\code{\link{kl.ku}},
\code{\link{undom.para}}
}
