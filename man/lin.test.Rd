\name{lin.test}

\Rdversion{1.1}

\alias{lin.test}

\title{
  Lack-of-fit test for the mean structure based on cumulated residuals
}

\description{
  The function provides a lack-of-fit test for the mean structure based on cumulated residuals from the model fit.
}

\usage{
  lin.test(object, noksSim = 20, seed = 20070325, plotit = TRUE, 
  log = "", bp = 0.01, xlab, ylab, ylim, ...)
}

\arguments{
  \item{object}{
  object of class 'drc'.
}
  \item{noksSim}{
  numeric specifying the number of simulations used to obtain the p-value.
}
  \item{seed}{
  numeric specifying the seed value for the random number generator.
}
  \item{plotit}{
  logical indicating whether or not the observed cumulated residual process should be plotted. Default is to
  plot the process.
}
  \item{log}{
  character string which should contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be 
  logarithmic and '"xy"' or '"yx"' if both axes are to be logarithmic. The default is "x". 
  The empty string "" yields the original axes.
}
  \item{bp}{
  numeric value specifying the break point below which the dose is zero (the amount of stretching on 
  the dose axis above zero in order to create the visual illusion of a logarithmic scale \emph{including} 0).
}
  \item{xlab}{
  string character specifying an optional label for the x axis.
}
  \item{ylab}{
  character string specifying an optional label for the y axis.
}
  \item{ylim}{
  numeric vector of length two, containing the lower and upper limit for the y axis.
}
  \item{\dots}{
  additional arguments to be passed further to the basic \code{\link{plot}} method.
}
}

\details{
  The function provides a graphical model checking of the mean structure in a dose-response model. The graphical
  display is supplemented by a p-value based on a supremum-type test.
  
  The test is applicable even in cases where data are non-normal or exhibit variance heterogeneity.
}

\value{
  A p-value for test of the null hypothesis that the mean structure is appropriate. 
  Ritz and Martinussen (2009) provide the details.
}
\references{
  Ritz, C and Martinussen, T. (2009) 
  Lack-of-fit tests for assessing mean structures for continuous dose-response data, 
  \emph{Submitted manuscript}
}
\author{
 Christian Ritz
}

%\note{
%}

\seealso{
 Other available lack-of-fit tests are the Neill test (\code{\link{neill.test}}) and 
 ANOVA-based test (\code{\link{modelFit}}).
}

\examples{

## Fitting a log-logistic model to the dataset 'etmotc'
etmotc.m1<-drm(rgr1~dose1, data=etmotc[1:15,], fct=LL.4())

## Test based on umulated residuals
lin.test(etmotc.m1, 1000)
#lin.test(etmotc.m1, 10000, plotit = FALSE)  # more precise

## Fitting an exponential model to the dataset 'O.mykiss'
O.mykiss.m1<-drm(weight~conc, data=O.mykiss, fct=EXD.2(), na.action=na.omit)

## ANOVA-based test
modelFit(O.mykiss.m1)

## Test based on umulated residuals
lin.test(O.mykiss.m1, log = "", cl = 0.2, xlab = "Dose (mg/l)", main = "B", ylim = c(-0.6, 0.6))
#lin.test(O.mykiss.m1, noksSim = 10000, plotit = FALSE)  # more precise

}

\keyword{models}
\keyword{nonlinear}
