\name{EXD}

\alias{EXD.2}
\alias{EXD.3}

\title{Exponential decay model}

\description{
  Exponential decay model with or without a nonzero lower limit.
}

\usage{
  EXD.2(fixed = c(NA, NA), names = c("d", "e"), ...)

  EXD.3(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The default parameter names are: init, plateau, k.}
  \item{...}{additional arguments to be passed from the convenience functions.}   
}

\details{
  The exponential decay model is a three-parameter model with mean function:
  
  \deqn{f(x) = c + (d-c)(\exp(-x/e))}
  
  The parameter init is the upper limit (attained at \eqn{x=0}), the parameter plateau is the lower limit 
  reached for x going to infinity and the parameter \eqn{e>0} is determining the steepness of the 
  decay. The curve is monotonously decreasing in \eqn{x}.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

\references{
  Organisation for Economic Co-operation and Development (OECD) (2006)
  \emph{Current approaches in the statistical analysis of ecotoxicity data: A guidance to application - annexes},
  Paris: OECD (p. 80).
}

\author{Christian Ritz}

%\note{}

\seealso{
  Similar models giving exponential increasing curves are \code{\link{AR.2}} and \code{\link{AR.3}}.
}

\examples{

## Fitting an exponential decay model
ryegrass.m1<-drm(rootl~conc, data=ryegrass, fct=EXD.3())

plot(ryegrass.m1)

summary(ryegrass.m1)

}
\keyword{models}
\keyword{nonlinear}
