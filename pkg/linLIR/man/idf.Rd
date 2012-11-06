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

\method{summary}{idf}(object, ...)
\method{plot}{idf}(x, y=NULL, ..., var = NULL, typ="hist", k.x = 1, k.y = 1, inf.margin=10, p.cex=1, col.lev=15, plot.grid=FALSE, x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, y.padj = 0, x.lim = c(0, 0), y.lim = c(0, 0), x.lab = "X", y.lab = "Y")
}
\arguments{
  \item{dat}{
A \code{data.frame} containing 2 neighboring columns for each variable, the first column for the left endpoints of the interval observations, the second for the right endpoints.
}
  \item{var.labels}{
Names of the variables corresponding to the interval-valued observations in the \code{data.frame}.
}

  \item{object}{
The \code{idf}-object to be summarized.
}
  \item{...}{
Argument of the generic functions \code{plot} and \code{summary}: Other parameters.
}

  \item{x}{
Argument of the generic function \code{plot}. Here \code{x} is the \code{idf}-object to be plotted.
}
  \item{y}{
Argument of the generic function \code{plot}. Here \code{y=NULL}.
}

  \item{var}{
Names of the two variables out of the \code{idf}-object to be plotted. (Optional)
}
  \item{typ}{Type of the plot. Possible values are \code{"hist"}: plot 2-dim. histogram (default) and \code{"draft"}.
 
}
  \item{k.x}{
Particular plot function parameter. 1/\code{k.x} is the step width along the abscissa.
}
  \item{k.y}{
Particular plot function parameter. 1/\code{k.y} is the step width along the ordinate.
}
  \item{inf.margin}{
Particular parameter for plot type \code{"draft"}. \code{inf.margin} is the number of steps that the infinite observations are drawn beyond the limits of the plot.
}
  \item{p.cex}{
Particular parameter for plot type \code{"draft"}. \code{p.cex} is the point size to fill the rectangles with grey color.
}
  \item{col.lev}{
Particular parameter for plot type \code{"hist"} indicating the number of different grey levels in the 2-dim. histogram. 
}
  \item{plot.grid}{
Logical for plot type \code{"hist"}. If \code{plot.grid=TRUE} dashed lines are added to the plot to indicate the location of the interval endpoints. This is particularly useful for categorized data.
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
\details{
Within the LIR framework all types of interval data are possible, including the particular cases of actually precise data (i.e., lower endpoint = upper endpoint) or missing data (i.e., in case of a real valued variable, lower endpoint = -Inf and upper endpoint = Inf). For the LIR analysis it makes practically no difference if the intervals are closed or not, therefore, the created \code{idf}-object does not contain this information.
}
\value{
An \code{idf}-object of \code{m} variables, which is a list of \code{m+1} entries.
\item{Var1 ... varm}{\code{m} different \code{data.frame}s with 2 columns each, one for each of the 1st to \code{m}th variables.}
\item{n}{Number of observations.}
}
\references{
M. Cattaneo, A. Wiencierz (2012c). On the implementation of LIR: the case of simple linear regression with interval data. Technical Report No. 127. Department of Statistics. LMU Munich.

A. Wiencierz, M. Cattaneo (2012b). An exact algorithm for Likelihood-based Imprecise Regression in the case of simple linear regression with interval data. In: R. Kruse et al. (Eds.). Advances in Intelligent Systems and Computing. Vol. 190. Springer. pp. 293-301.

M. Cattaneo, A. Wiencierz (2012a). Likelihood-based Imprecise Regression. International Journal of Approximate Reasoning. Vol. 53. pp. 1137-1154.
}
\examples{
data('toy.smps')
toy.idf <- idf.create(toy.smps, var.labels=c("x","y"))

summary(toy.idf)

plot(toy.idf, typ="draft", k.x=10, k.y=10, p.cex=1.5, y.las=1, y.adj=6) 
plot(toy.idf, typ="draft", k.x=10, k.y=10, x.adj=0.7, y.las=1, y.adj=6, y.padj=-3)
plot(toy.idf, k.x=10, k.y=10, x.adj=0.7, x.padj=4, y.adj=0.7, y.padj=-4)

data('pm10')
pm.idf <- idf.create(pm10)

summary(pm.idf)

plot(pm.idf, typ="draft", k.x=10, k.y=20, p.cex=0.35, x.adj=0.5, x.padj=4, y.las=0, y.adj=0.5, y.padj=-4, x.lab="temperature", y.lab="particulate matter concentration")
}