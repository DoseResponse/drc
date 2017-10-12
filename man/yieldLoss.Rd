\name{yieldLoss}

\alias{yieldLoss}

\title{Calculating yield loss parameters}

\description{
  Calculation of parameters in the re-parameterization of the Michaelis-Menten model that is commonly
  used to assess yield loss (the rectangular hyperbola model)
}

\usage{
  yieldLoss(object, interval = c("none", "as"), level = 0.95, display = TRUE)
}

\arguments{
  \item{object}{object of class 'drc}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "as" for asymptotically-based confidence intervals.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
}

\details{
  The rectangular hyperbola model is a reparameterization of the Michaelis-Menten in terms of parameters
  \eqn{A} and \eqn{I}
  
  \deqn{ Y_L = \frac{Id}{1+Id/A}}
  
  where \eqn{d} denotes the weed density and \eqn{Y_L} the resulting yield loss.
}

\value{
  For each of the two parameters, a matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
}

\references{
  Cousens, R. (1985). A simple model relating yield loss to weed density,
  \emph{Ann. Appl. Biol.}, \bold{107}, 239--252.  
}

\author{Christian Ritz}

\note{
  This function is only for use with model fits based on Michaelis-Menten models.
}

%\seealso{}

\examples{

## Fitting Michaelis-Menten model
met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.3(), 
pmodels = list(~1, ~factor(product), ~factor(product)))

## Yield loss parameters with standard errrors
yieldLoss(met.mm.m1)

## Also showing confidence intervals
yieldLoss(met.mm.m1, "as")

}

\keyword{models}
\keyword{nonlinear}

\concept{rectangular hyperbola model}
