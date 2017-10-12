\name{residuals.drc}

\alias{residuals.drc}

\title{Extracting residuals from the fitted dose-response model}

\description{
  'residuals' extracts different types of residuals from an object of class 'drc'.
}

\usage{

  \method{residuals}{drc}(object, typeRes = c("working", "standardised", "studentised"), 
  trScale = TRUE, ...)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{typeRes}{character string specifying the type of residual to be returned: raw/working residuals, 
  residuals standardised using the estimated residual standard error,
  or studentised residuals based on the H matrix of partial derivatives of the model function.}
  \item{trScale}{logical value indicating whether or not to return residuals on the transformed scale (in case a Box-Cox transformation was applied).}
  \item{...}{additional arguments.}
}

\value{
  The raw (also called working) residuals or some kind of scaled residuals extracted from 'object'. 
}

\details{
  Standardised residuals are the raw residuals divided by a scale estimate (if available). 
  
  Studentised residuals are obtained by dividing by a scale estimate and in
  addition a correction factor (square root of 1 minus h with h is a diagonal element in the hat matrix). 
}  

\note{
  The 'standardised' residuals are available for least squares estimation
  with or without Box-Cox transformation or variance as a power of the
  mean.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())

## Displaying the residual plot (raw residuals)
plot(fitted(ryegrass.m1), residuals(ryegrass.m1))

## Using the standardised residuals
plot(fitted(ryegrass.m1), residuals(ryegrass.m1, typeRes = "standard"))

## Overlayering the studentised residuals ... not much of a difference
points(fitted(ryegrass.m1), residuals(ryegrass.m1, typeRes = "student"), col = 2)

}
\keyword{models}
\keyword{nonlinear}
