\name{confint.drc}

\alias{confint.drc}

\title{Confidence Intervals for model parameters}

\description{
  Computes confidence intervals for one or more parameters in a model of class 'drc'.
}

\usage{

  \method{confint}{drc}(object, parm, level = 0.95, pool = TRUE, ...)

}

\arguments{
  \item{object}{a model object of class 'drc'.}
  \item{parm}{a specification of which parameters are to be given
  confidence intervals, either a vector of numbers or a vector
  of names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
%  \item{type}{the type of confidence interval: based on the standard normal distribution or 
%  based on a t-distribution (default).}  
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.}  
  \item{\dots}{additional argument(s) for methods. Not used.}
}

\details{
  For binomial and Poisson data the confidence intervals are based on the normal distribution, whereas \emph{t} distributions
  are used of for continuous/quantitative data.
}

\value{
  A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as
  (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Confidence intervals for all parameters
confint(ryegrass.m1)

## Confidence interval for a single parameter
confint(ryegrass.m1, "e")

}
\keyword{models}
\keyword{nonlinear}
