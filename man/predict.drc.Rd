\name{predict.drc}

\alias{predict.drc}
%\alias{predict.mrdrc}

\title{Prediction}

\description{
  Predicted values for models of class 'drc'.
}

\usage{

  \method{predict}{drc}(object, newdata, se.fit = FALSE, 
  interval = c("none", "confidence", "prediction", "ssd"), 
  level = 0.95, na.action = na.pass, od = FALSE, vcov. = vcov, 
  ssdSEfct = NULL, constrain = TRUE, checkND = TRUE, ...)

%  \method{predict}{mrdrc}(object, newdata, se.fit = FALSE, 
%  interval = c("none", "confidence", "prediction"), 
%  level = 0.95, pava = FALSE, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{newdata}{An optional data frame in which to look for variables with which to predict.
  If omitted, the fitted values are used.}
  \item{se.fit}{logical. If TRUE standard errors are required.}
  \item{interval}{character string. Type of interval calculation: "none", "confidence", "prediction", or "ssd".}
  \item{level}{Tolerance/confidence level.}
  \item{na.action}{function determining what should be done with missing values in 'newdata'. 
  The default is to predict 'NA'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default, 
  but \code{sandwich} is also an option (for obtaining robust standard errors).}  
  \item{ssdSEfct}{specifies the function for interpolating standard errors between observed standard errors.
  The default is linear interpolation on log-log scale (back-transformed). See Details for more explanation.}
  \item{constrain}{logical. If TRUE (default) predicted values are truncated within meaningful limits, i.e., 
  0 and, possibly, 1.}
%  \item{pava}{logical. If TRUE the fit is monotoniosed using pool adjacent violators algorithm.}
  \item{checkND}{logical indicating whether or not names in "newdata" data frame match 
  the names in the original data frame (used for fitting the model). Default is TRUE.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  For the built-in log-logistic, log-normal, and Weibull-type models standard errors and confidence/prediction
  intervals can be calculated. For other built-in models it may not yet be implemented (drop us an e-mail if 
  you need them). 
  
  The function for interpolating standard errors of estimates, which may be used when fitting an SSD, should 
  have 3 arguments: observed estimates and corresponding standard errors and future estimates and should return
  interpolated standard errors corresponding to the future estimates provided.
}

\value{
  A matrix with as many rows as there are dose values provided in 'newdata' or in the original dataset 
  (in case 'newdata' is not specified) and, at most, 4 columns containing fitted, standard errors, lower and 
  upper limits of confidence/prediction intervals.
}

%\references{}

\author{Christian Ritz}

%\note{}

\seealso{For details are found in the help page for \code{\link{predict.lm}}.}

\examples{

## Fitting a model
spinach.model1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

## Predicting values a dose=2 (with standard errors)
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")), se.fit = TRUE)

## Getting confidence intervals
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")), 
interval = "confidence")

## Getting prediction intervals
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")), 
interval = "prediction")

}

\keyword{models}
\keyword{nonlinear}

\concept{prediction}