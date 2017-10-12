\name{compParm}

\alias{compParm}

\title{Comparison of parameters}

\description{
  Compare parameters from different assays, either by means of ratios or differences.
}

\usage{
  compParm(object, strVal, operator = "/", vcov. = vcov, od = FALSE, 
  pool = TRUE, display = TRUE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{strVal}{a name of parameter to compare.}
  \item{operator}{a character. If equal to "/" (default) parameter ratios are compared. If equal to "-" parameter differences are compared.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default, 
  but \code{sandwich} is also an option (for obtaining robust standard errors).}  
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}  
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, values of t-statistics and p-values for the null hypothesis that 
  the ratio equals 1 or that the difference equals 0 (depending on the \code{operator} argument). 
}

\details{
  The function compares actual parameter estimates, and therefore the results depend on the parameterisation used. Probably it is most useful
  in combination with the argument \code{collapse} in \code{\link{drm}} for specifying parameter constraints in models, either through 
  data frames or lists with formulas without intercept (-1). 
}

\author{Christian Ritz}

\examples{

# Fitting a model with names assigned to the parameters!
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach,
fct = LL.4(names = c("b", "lower", "upper", "ed50")))

## Calculating ratios of parameter estimates for the parameter named "ed50"
compParm(spinach.m1, "ed50")      

## Calculating differences between parameter estimates for the parameter named "ed50"
compParm(spinach.m1, "ed50", "-")  

}
\keyword{models}
\keyword{nonlinear}
