\name{isobole}
\alias{isobole}

\title{Creating isobolograms}

\description{
  'isobole' displays isobole based on EC/ED50 estimates from a log-logistic model.
  Additionally isoboles determined by the concentration addition model, Hewlett's
  model and Voelund's model can be added to the plot.
}

\usage{
isobole(object1, object2, exchange = 1, cifactor = 2, ename = "e", 
xaxis = "100", xlab, ylab, xlim, ylim, ...)
}

\arguments{
  \item{object1}{object of class 'drc' where EC/ED50 parameters vary freely.}
  \item{object2}{object of class 'drc' where EC/ED50 parameters vary according to Hewlett's model.}
  \item{ename}{character string. The name of the EC/ED50 variable.}
  \item{xaxis}{character string. Is the mixture "0:100" or "100:0" on the x axis?}
  \item{exchange}{numeric. The exchange rate between the two substances.}
  \item{cifactor}{numeric. The factor to be used in the confidence intervals. 
  Default is 2, but 1 has been used in publications.}  
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{\dots}{Additional graphical parameters.}
}

\details{
  The model fits to be supplied as first and optionally second argument are obtained
  using \code{\link{mixture}} and \code{\link{drm}}.
}

\value{
  No value is returned. Only used for the side effect: the isobologram shown.
}

\references{
  Ritz, C. and Streibig, J. C. (2014)
  From additivity to synergism - A modelling perspective
  \emph{Synergy}, \bold{1}, 22--29.
}

\author{Christian Ritz}

%\note{}

\seealso{
  The examples in \code{\link{acidiq}}, \code{\link{glymet}} and \code{\link{mecter}}.
}

%\examples{}

\keyword{models}
\keyword{nonlinear}
