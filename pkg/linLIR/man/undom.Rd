\name{undom}
\alias{undom.a}
\alias{undom.para}
\title{
Finding undominated parameters
}
\description{
Functions within the \code{s.linlir}-function that determine the parameter combinations corresponding to undominated regression lines. The \code{undom.a}-function finds the set of undominated intercept values associated with a given slope and the \code{undom.para}-function finds the set of undominated intercept values associated with all slope values in a given range.
}
\usage{
undom.a(dat, b, q.lrm, p = 0.5, bet, epsilon = 0)

undom.para(dat, b.range, b.extra = 0, b.grid = 1000, q.lrm, p = 0.5, bet, epsilon = 0)
}
\arguments{
  \item{dat}{
An \code{n}x4 \code{data.frame} containing the imprecise data of the analyzed variables. Columns 1 and 2 correspond to the interval-valued observations of the regressor variable, columns 3 and 4 to those of the dependent variable.
}
  \item{b}{
A given value for the slope of a regression line.
}
  \item{q.lrm}{
Value of the p-quantile of the absolute residuals associated with the LRM line(s).
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
  \item{b.range}{
Considered interval of slope values.
}
  \item{b.extra}{
If it is likely that the set of undominated regression parameters is unbounded, one can consider very small or large values in addition to the \code{b.range}.
}
  \item{b.grid}{
1/\code{b.grid} is the width between two points of the considered grid of slope values within \code{b.range}.
}
}
\value{
The \code{undom.a}-function returns a list of 2 components:
\item{result1}{A 2-column matrix of possibly degenerate intervals for the undominated intercept values associated with the given slope b.}
\item{result2}{The information of \code{result1} reduced to the fewest intervals possible.}

The \code{undom.para}-function returns a list of 3 components:
\item{a.undom}{Range of intercept values of the undominated regression lines.}
\item{b.undom}{Range of slope values of the undominated regression lines.}
\item{undom.para}{A matrix of undominated parameter combinations approximating the entire set of parameters corresponding to the set of undominated regression lines.}
}
\references{
M. Cattaneo, A. Wiencierz (2012b). On the implementation of LIR: the case of simple linear regression with interval data. Technical Report No. 127. Department of Statistics. LMU Munich.

A. Wiencierz, M. Cattaneo (2012). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. (Accepted for the 6th International Conference on Soft Methods in Probability and Statistics (SMPS 2012). Publication in the series Advances in Intelligent and Soft Computing. Springer-Verlag.)

M. Cattaneo, A. Wiencierz (2012a). Likelihood-based Imprecise Regression. International Journal of Approximate Reasoning. In Press. (Uncorrected Proof: \url{http://www.sciencedirect.com/science/article/pii/S0888613X12000862})
}
\seealso{
\code{\link{s.linlir}},
\code{\link{gen.lms}},
\code{\link{kl.ku}}
}
