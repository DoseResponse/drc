\name{EDcomp}

\alias{EDcomp}
\alias{relpot}

\title{Comparison of relative potencies between dose-response curves}

\description{
  Relative potencies (also called selectivity indices) for arbitrary doses are compared between
  fitted dose-response curves.
}

\usage{
  EDcomp(object, percVec, percMat = NULL, compMatch = NULL, od = FALSE, vcov. = vcov, 
  reverse = FALSE, 
  interval = c("none", "delta", "fieller", "fls"), 
  level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"),  
  type = c("relative", "absolute"),
  display = TRUE, pool = TRUE, logBase = NULL, 
  multcomp = FALSE, ...)
  
  relpot(object, plotit = TRUE, compMatch = NULL, percVec = NULL, interval = "none", 
  type = c("relative", "absolute"), 
  scale = c("original", "percent", "unconstrained"), ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{percVec}{a numeric vector of dosage values.}
  \item{percMat}{a matrix with 2 columns providing the pairs of indices \code{percVec} to be compared. 
  By default all pairs are compared.}
  \item{compMatch}{an optional character vector of names of assays to be compared. If not specified all comparisons are supplied.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used. This argument only makes a difference for
  binomial data.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default, 
  but \code{sandwich} is also an option (for obtaining robust standard errors).}
  \item{reverse}{logical. If TRUE the order of comparison of two curves is reversed.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fieller" for confidence intervals based on Fieller's theorem (with help from the delta method).
  Use "fls" for confidence interval back-transformed from logarithm scale (in case the parameter in the model fit is 
  log(ED50) as is the case for the \code{\link{logistic}} or \code{\link{llogistic2}} models); currently
  the argument \code{logBase} then also needs to be specified.}
  \item{level}{numeric. The level for the confidence intervals. Default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}  
  \item{type}{character string specifying whether absolute or relative response levels are supplied.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.}
  \item{multcomp}{logical to switch on output for use with the package multcomp (which needs to be activated first). Default is FALSE (corresponding to the original output).}  
  \item{...}{In \code{SI}: additional arguments to the function doing the calculations. 
  For instance the upper limit for the bisection method 
  needs to be larger than the ED values used in the required relative pontency.
  In \code{relpot}: additional graphical parameters.}  
  \item{plotit}{logical. If TRUE the relative potencies are plotted as a function of the response level.}    
  \item{scale}{character string indicating the scale to be used on the x axis: original or percent response level
  (only having an effect for type="relative").}    
}

\details{
  The function \code{relpot} is a convenience function, which is useful for assessing how the relative potency
  changes as a function of the response level (e.g., for plotting as outlined by Ritz \emph{et al} (2006)).
  
  Fieller's theorem is incorporated using the formulas provided by Kotz and Johnson (1983) and Finney (1978).

  For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper' argument
  or the 'interval' argument. The 'upper' argument specifies the upper limit of the bisection method. 
  The upper limits needs to be larger than the EDx level to be calculated. The default limit is 1000.
  The 'interval' argument should specify a rough interval in which the dose
  yielding the maximum hormetical response lies. The default interval is 'c(0.001, 1000)'. 
  Notice that the lower limit should not be set to 0 (use something like 1e-3, 1e-6, ...).
}

\value{  
  An invisible matrix containing the shown matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
  Or, alternatively, a list with elements that may be plugged directly into \code{parm} 
  in the package \emph{multcomp} (in case the argument \code{multcomp} is TRUE).
}

\references{
  Finney, D. J. (1978) \emph{Statistical method in Biological Assay}, London: Charles Griffin House, 
  3rd edition (pp. 80--82).

  Kotz, S. and Johnson, N. L. (1983) \emph{Encyclopedia of Statistical Sciences Volume 3}, 
  New York: Wiley \& Sons (pp. 86--87).
  
  Ritz, C. and Cedergreen, N. and Jensen, J. E. and Streibig, J. C. (2006) 
  Relative potency in nonsimilar dose-response curves, \emph{Weed Science}, \bold{54}, 407--412.
}

\author{Christian Ritz}

\note{
  This function only works for the following built-in functions available in the package \emph{drc}: 
  \code{\link{braincousens}}, \code{\link{cedergreen}}, \code{\link{ucedergreen}}, \code{\link{llogistic}},
  and \code{\link{weibull1}}.
}

\seealso{
  A related function is \code{\link{ED.drc}} (used for calculating effective doses).
}

\examples{

spinach.LL.4 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

EDcomp(spinach.LL.4, c(50,50))
EDcomp(spinach.LL.4, c(10,50))
EDcomp(spinach.LL.4, c(10,50), reverse = TRUE)

## Using the package multcomp
#sires <- SI(spinach.LL.4, c(25, 50, 75))
#library(multcomp)
#summary(glht(parm(sires[[2]][[1]], sires[[2]][[2]]), rhs = 1))

## Comparing specific ratios: 25/25, 50/50, 75/75
#sires2 <- SI(spinach.LL.4, c(25, 50, 75), matrix(c(1, 1, 2, 2, 3, 3), 3, 2, byrow = TRUE))
#library(multcomp)
#summary(glht(parm(sires2[[2]][[1]], sires2[[2]][[2]]), rhs = 1))


## Relative potency of two herbicides
m2 <- drm(DryMatter~Dose, Herbicide, 
data = S.alba, fct = LL.3())

EDcomp(m2, c(50, 50))  
EDcomp(m2, c(50, 50), interval = "delta")
EDcomp(m2, c(50, 50), interval = "fieller")

## Comparison based on an absolute
##  response level

m3 <- drm(SLOPE~DOSE, CURVE,
data = spinach, fct = LL.4())

EDcomp(m3, c(0.5,0.5), compMatch = c(2,4), type = "absolute", interval = "fieller")

EDcomp(m3, c(55,80), compMatch = c(2,4))  
# same comparison using a relative response level


## Relative potency transformed from log scale
m4 <- drm(drymatter~log(dose), treatment, data=G.aparine[-c(1:40), ], 
pmodels = data.frame(treatment,treatment,1,treatment), fct = LL2.4())

EDcomp(m4, c(50,50), interval = "fls", logBase = exp(1))

}

\keyword{models}
\keyword{nonlinear}
