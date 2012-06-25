\name{s.linlir}
\alias{s.linlir}
\alias{print.s.linlir}
\alias{summary.s.linlir}
\alias{plot.s.linlir}
\title{
Simple linear Likelihood-based Imprecise Regression
}
\description{
Conducts a LIR analysis for 2 variables with interval-valued observations whose relation is assumed to be linear. 
}
\usage{
s.linlir(dat.idf, var = NULL, p = 0.5, bet, epsilon = 0, b.grid = 1000)

\method{print}{s.linlir}(x, ...)
\method{summary}{s.linlir}(object, ...)
\method{plot}{s.linlir}(x, y=NULL, ..., typ, para.typ = "polygon", b.range = c(-1e-05, 1e-05), b.grid = 1000, nb.func = 1000, seed.func = NULL, pl.lrm = TRUE, pl.band = FALSE, pl.dat = FALSE, k.x = 1, k.y = 1, p.cex = 1, x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, y.padj = 0, x.lim = c(0, 0), y.lim = c(0, 0), x.lab = " ", y.lab = " ")
}
\arguments{
  \item{dat.idf}{
The \code{idf}-object to be analyzed.
}
  \item{var}{
Names of the two variables out of the \code{idf}-object to be analyzed.
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
  \item{b.grid}{
Particular parameter of the internal function undom.para determining the undominated parameter combinations.
}

  \item{x}{
Argument of the generic functions \code{plot} and \code{print}. Here \code{x} is the \code{s.linlir}-object to be plotted or printed.
}
  \item{...}{
Argument of the generic functions \code{plot}, \code{print} and \code{summary}: Other parameters.
}

\item{object}{
The \code{s.linlir}-object to be summarized.
}

  \item{y}{
Argument of the generic function \code{plot}. Here \code{y=NULL}.
}


  \item{typ}{Type of the plot. Possible values are \code{"para"}: plot undominated parameter set, \code{"lrm"}: plot \code{f.lrm}, \code{"func"}: plot undominated regression functions.
}
  \item{para.typ}{Options for plot of \code{typ="para"} are \code{"polygon"} (default) or \code{"points"} (approximation).
}
  \item{b.range}{Considered interval of slopes for plots of \code{typ="para"}. If default is chosen, then \code{x.lim} ist taken as \code{b.range}.
}
  \item{nb.func}{Number of (randomly chosen) plotted undominated lines for plots of \code{typ="func"}.}
  \item{seed.func}{Set seed for the random selection of plotted regression lines for plots of \code{typ="func"}. (Optional)}
  \item{pl.lrm}{Logical for plots of \code{typ=c("para","func")}. If \code{pl.lrm=TRUE} (default), \code{f.lrm} is highlighted in the plot.}
  \item{pl.band}{Logical for plots of \code{typ="func"}. If \code{pl.band=TRUE}, the band around \code{f.lrm} is added to the plot.}
  \item{pl.dat}{Logical for plots of \code{typ=c("lrm","func")}. If \code{pl.dat=TRUE}, the data are plotted in the background of the plot.}
  \item{k.x}{
Particular plot function parameter. 1/\code{k.x} is the step width along the abscissa.
}
  \item{k.y}{
Particular plot function parameter. 1/\code{k.y} is the step width along the ordinate.
}

  \item{p.cex}{
Particular plot function parameter. \code{p.cex} is the point size to fill the rectangles with grey color.
}
  \item{x.adj}{
Horizontal position of the text for the abscissa. 
}
  \item{x.padj}{
Vertical position of the text for the abscissa. 
}
  \item{y.las}{
Orientation of the text for the ordinate. \code{y.las=1} will turn the axis labels and the text in reading direction.
}
  \item{y.adj}{
\code{y.adj} regulates the position of the text for the ordinate in reading direction, i.e. if \code{y.las=0} it sets the vertical position and if \code{y.las=1} the horizontal position.
}
  \item{y.padj}{
\code{y.padj} regulates the position of the text for the ordinate orthogonal to the reading direction, i.e. if \code{y.las=0} it sets the horizontal position and if \code{y.las=1} the vertical position.
}

  \item{x.lim}{
The limits for the abscissa of the plot.
}
  \item{y.lim}{
The limits for the ordinate of the plot.
}
  \item{x.lab}{
Title of the abscissa.
}
  \item{y.lab}{
Title of the ordinate.
}
}
\value{
\item{f.lrm}{Intercept and slope of the Likelihood-based Region Minimax (LRM) line(s).}
\item{q.lrm}{Value of the p-quantile of the absolute residuals associated with the LRM line(s).}
\item{a.undom}{Range of intercept values of the undominated regression lines.}
\item{b.undom}{Range of slope values of the undominated regression lines.}
\item{undom.para}{A matrix of undominated parameter combinations approximating the entire set of parameters corresponding to the set of undominated regression lines.}
\item{config}{A list containing information about the settings of the LIR analysis.}
\item{dat}{An \code{n}x4 \code{data.frame} containing the imprecise data of the analyzed variables. Columns 1 and 2 correspond to the interval-valued observations of the regressor variable, columns 3 and 4 to those of the dependent variable.}
\item{n}{Number of observations.}
\item{call}{Call of the function \code{s.linlir}.}

}
\references{
A. Wiencierz, M. Cattaneo (2012). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. (Accepted for the 6th International Conference on Soft Methods in Probability and Statistics (SMPS 2012). Publication in the series Advances in Intelligent and Soft Computing. Springer-Verlag.)

M. Cattaneo, A. Wiencierz (2012). Likelihood-based Imprecise Regression. (Accepted for publication in the International Journal of Approximate Reasoning. A preliminary version of the paper is available as a research report at: \url{http://epub.ub.uni-muenchen.de/12450/}.)
}
\seealso{
\code{\link{idf.create}},
\code{\link{gen.lms}},
\code{\link{kl.ku}},
\code{\link{undom.para}}
}
\examples{
data('toy.smps')
toy.idf <- idf.create(toy.smps, var.labels=c("x","y"))

test <- s.linlir(toy.idf, bet=0.5)
test

summary(test)

plot(test, typ="para", x.adj=0.7, y.las=1, y.adj=6, y.padj=-3)
plot(test, typ="func", pl.lrm=FALSE, x.adj=0.7, y.adj=0.7, y.padj=-3)
plot(test, typ="lrm", pl.dat=TRUE, pl.band=TRUE, k.x=10, k.y=10, y.las=1, y.adj=6)
}