\name{gen.lms}
\alias{gen.lms}
\title{
Finding Likelihood-based Region Minimax (LRM) line(s)
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
Fraction of errors considered.
}
  \item{k.u}{
As default \code{k.u} is calculated on the basis of \code{p}, \code{bet} and \code{epsilon}.
}
}
\details{
The exact algorithm implemented by the function gem.lms can be seen as a generalization of the basic algorithm for Least Median of Squares Regression (see, e.g., Steele / Steiger (1986) and Rousseeuw / Leroy (1987)). 
}
\value{
A list of two components.
\item{lrm}{A vector (or a matrix) of the intercept and slope parameter values of the LRM line(s) and the value of the p-quantile of the absolute residuals associated with the LRM line(s).}
\item{max.b}{The maximal absolute value of the considered slopes for the LRM line(s).}
}
\references{
M. Cattaneo, A. Wiencierz (2012b). On the implementation of LIR: the case of simple linear regression with interval data. Technical Report No. 127. Department of Statistics. LMU Munich.

A. Wiencierz, M. Cattaneo (2012). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. (Accepted for the 6th International Conference on Soft Methods in Probability and Statistics (SMPS 2012). Publication in the series Advances in Intelligent and Soft Computing. Springer-Verlag.)

M. Cattaneo, A. Wiencierz (2012a). Likelihood-based Imprecise Regression. International Journal of Approximate Reasoning. In Press. (Uncorrected Proof: \url{http://www.sciencedirect.com/science/article/pii/S0888613X12000862})

P. Rousseeuw, A. Leroy (1987). Robust Regression and Outlier Detection. Wiley

J. Steele, W. Steiger (1986). Algorithms and complexity for least median of squares regression. Discret Appl Math 14. 93-100.
}
\seealso{
\code{\link{s.linlir}},
\code{\link{kl.ku}},
\code{\link{undom.para}}
}
