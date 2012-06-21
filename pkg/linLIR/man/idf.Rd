\name{idf}
\alias{idf.create}
\alias{summary.idf}
\alias{plot.idf}
\title{
Interval data frame
}
\description{
Create an interval data frame (\code{idf}-object), summarize its content and visualize subsets of two variables.
}
\usage{
idf.create(dat, var.labels = NULL)

\method{summary}{idf}(dat.idf)
\method{plot}{idf}(dat.idf, var = NULL, k.x = 1, k.y = 1, x.lim = c(0,0), y.lim = c(0,0), x.lab = "X", y.lab = "Y")
}
\arguments{
  \item{dat}{
A \code{data.frame} containing 2 neighboring columns for each variable, the first column for the left endpoints of the interval observations, the second for the right endpoints.
}
  \item{var.labels}{
Names of the variables corresponding to the interval-valued observations in the \code{data.frame}.
}

  \item{dat.idf}{
The \code{idf}-object to be summarized or plotted.
}

  \item{var}{
Names of the two variables out of the \code{idf}-object to be plotted. (Optional)
}
  \item{k.x}{
Particular plot function parameter. 1/\code{k.x} is the step width along the abscissa.
}
  \item{k.y}{
Particular plot function parameter. 1/\code{k.y} is the step width along the ordinate.
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
\details{
Within the LIR framework all types of interval data are possible, including the particular cases of actually precise data (i.e. lower endpoint = upper endpoint) or missing data (i.e., in case of a real valued variable, lower endpoint = -Inf and upper endpoint = Inf). For the LIR analysis it makes no difference if the intervals are closed or not, therefore, the created \code{idf}-object does not contain this information.
}
\value{
An \code{idf}-object of \code{m} variables, which is a list of \code{m+1} entries.
\item{Var1 ... varm}{\code{m} different \code{data.frame}s with 2 columns each, one for each of the 1st to \code{m}th variables.}
\item{n}{Number of observations.}
}
\references{
A. Wiencierz, M. Cattaneo (2012). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. (Accepted for the 6th International Conference on Soft Methods in Probability and Statistics (SMPS 2012). Publication in the series Advances in Intelligent and Soft Computing. Springer-Verlag.)

M. Cattaneo, A. Wiencierz (2012). Likelihood-based Imprecise Regression. (Accepted for publication in the International Journal of Approximate Reasoning. A preliminary version of the paper is available as a research report at: \url{http://epub.ub.uni-muenchen.de/12450/}.)
}
\examples{
data('toy.smps')
toy.idf <- idf.create(toy.smps, var.labels=c("x","y"))

summary(toy.idf)

plot(toy.idf, k.x=10, k.y=10)
}