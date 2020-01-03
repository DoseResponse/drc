\name{ED.drc}

\alias{ED}
\alias{ED.drc}

\title{Estimating effective doses}

\description{

  \code{ED} estimates effective concentration or doses for one or more specified absolute or relative response levels.
}

\usage{
  \method{ED}{drc}(object, respLev, interval = c("none", "delta", "fls", "tfls", "inv"), 
  clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"), type = c("relative", "absolute"), lref, uref,
  bound = TRUE, vcov. = vcov, display = TRUE, logBase = NULL, 
  multcomp = FALSE, intType = "confidence", ...)  
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  See Details below for more explanation.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates
  are shown for all curves.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{lref}{numeric value specifying the lower limit to serve as reference.}
  \item{uref}{numeric value specifying the upper limit to serve as reference (e.g., 100\%).}  
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default, 
  but \code{sandwich} is also an option (for obtaining robust standard errors).}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with independently fitted curves as specified in \code{\link{drm}}.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}  
  \item{multcomp}{logical to switch on output for use with the package multcomp (which needs to be activated first). Default is FALSE (corresponding to the original output).}
  \item{intType}{string specifying the type of interval to use with the predict method in case the type of confidence
  interval chosen with the argument "type" is "inverse regression."}
  \item{...}{see the details section below.}
}

\details{
  There are several options for calculating confidence intervals through the argument \code{interval}. The option "delta" results in asymptotical Wald-type confidence intervals (using the delta method and the normal or t-distribution depending on the type of response). The option "fls" produces (possibly skewed) confidence intervals through back-transformation from the logarithm scale (only meaningful in case the parameter in the model is log(ED50) as for the \code{\link{llogistic2}}) models. The option "tfls" is for transforming back and forth from log scale (experimental). The option "inv" results in confidence intervals obtained through inverse regression.

  For hormesis models (\code{\link{braincousens}} and \code{\link{cedergreen}}), the additional 
  arguments \code{lower} and \code{upper} may be supplied. These arguments specify the lower and upper limits 
  of the bisection method used to find the ED values.  The lower and upper limits need to be smaller/larger 
  than the EDx level to be calculated. The default limits are 0.001 and 1000 for \code{braincousens} and 
  0.0001 and 10000 for \code{cedergreen} and \code{\link{ucedergreen}}, but this may need to be modified 
  (for \code{\link{cedergreen}} the upper limit may need to be increased and for \code{\link{ucedergreen}} 
  the lower limit may need to be increased). Note that the lower limit should not be set to 0 (use instead 
  something like 1e-3, 1e-6, ...).
}

\value{
  An invisible matrix containing the shown matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
  Or, alternatively, a list with elements that may be plugged directly into \code{parm} 
  in the package \emph{multcomp} (in case the argument \code{multcomp} is TRUE).
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\seealso{
  \code{\link{backfit}}, \code{\link{isobole}}, and \code{\link{maED}} use \code{\link{ED}} for specific calculations involving estimated ED values.

  The related function \code{\link{EDcomp}} may be used for estimating differences and ratios of ED values, 
  whereas \code{\link{compParm}} may be used to compare other model parameters. 
}

\examples{

## Fitting 4-parameter log-logistic model
ryegrass.m1 <- drm(ryegrass, fct = LL.4())

## Calculating EC/ED values
ED(ryegrass.m1, c(10, 50, 90)) 
## first column: the estimates of ED10, ED50 and ED90
## second column: the corresponding estimated standard errors 

### How to use the argument 'ci'

## Also displaying 95% confidence intervals
ED(ryegrass.m1, c(10, 50, 90), interval = "delta")

## Comparing delta method and back-transformed 
##  confidence intervals for ED values

## Fitting 4-parameter log-logistic 
##  in different parameterisation (using LL2.4)
ryegrass.m2 <- drm(ryegrass, fct = LL2.4())  

ED(ryegrass.m1, c(10, 50, 90), interval = "fls")
ED(ryegrass.m2, c(10, 50, 90), interval = "delta")


### How to use the argument 'bound'

## Fitting the Brain-Cousens model
lettuce.m1 <- drm(weight ~ conc, 
data = lettuce, fct = BC.4())

### Calculating ED[-10]

# This does not work
#ED(lettuce.m1, -10)  

## Now it does work
ED(lettuce.m1, -10, bound = FALSE)  # works
ED(lettuce.m1, -20, bound = FALSE)  # works

## The following does not work for another reason: ED[-30] does not exist 
#ED(lettuce.m1, -30, bound = FALSE)  

}

\keyword{models}
\keyword{nonlinear}
