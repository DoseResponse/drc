\name{maED}

\Rdversion{1.1}

\alias{maED}

\title{Estimation of ED values using model-averaging}

\description{
 Estimates and confidence intervals for ED values are estimated using model-averaging.
}

\usage{
maED(object, fctList = NULL, respLev, interval = c("none", "buckland", "kang"), 
linreg = FALSE, clevel = NULL, level = 0.95, type = c("relative", "absolute"), 
display = TRUE, na.rm = FALSE, extended = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{fctList}{a list of non-linear functions to be compared.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  The choices "buckland" and "kang" are explained in the Details section.}
  \item{linreg}{logical indicating whether or not additionally a simple linear regression model 
  should be fitted.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates
  are shown for all curves.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{na.rm}{logical indicating whether or not NA occurring during model fitting should be left out of
  subsequent calculations.}
  \item{extended}{logical specifying whether or not an extended output (including fit summaries) should be
  returned.}
}
\details{
  Model-averaging of individual estimates is carried out as described by Buckland \emph{et al.} (1997) and 
  Kang \emph{et al.} (2000) using AIC-based weights. The two approaches differ w.r.t. the calculation of confidence
  intervals: Buckland \emph{et al.} (1997) provide an approximate variance formula under the assumption of
  perfectly correlated estimates (so, confidence intervals will tend to be too wide). 
  Kang \emph{et al.} (2000) use the model weights to calculate confidence limits as weighted means of 
  the confidence limits for the individual fits; this procedure corresponds to using the standard error in Equation (3) 
  given by Buckland \emph{et al.} (1997) (assuming symmetric confidence intervals based on the same percentile).
}

\value{
  A matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
}

\references{ 
  Buckland, S. T. and Burnham, K. P. and Augustin, N. H. (1997)
  Model Selection: An Integral Part of Inference,
  \emph{Biometrics} \bold{53}, 603--618.
  
  Kang, Seung-Ho and Kodell, Ralph L. and Chen, James J. (2000)
  Incorporating Model Uncertainties along with Data Uncertainties in Microbial Risk Assessment,
  \emph{Regulatory Toxicology and Pharmacology} \bold{32}, 68--72.
}

\author{Christian Ritz}

%\note{}

\seealso{
  The function \code{\link{mselect}} provides a summary of fit statistics for several models fitted to the same data.
}

\examples{

## Fitting an example dose-response model
ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
  
## Comparing models (showing the AIC values)  
mselect(ryegrass.m1, 
list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)))
  
## Doing the actual model-averaging  
maED(ryegrass.m1, 
list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)), 
c(10, 50, 90))  
  
## With confidence intervals according to Buckland et al. (1997)  
maED(ryegrass.m1, 
list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)), 
c(10, 50, 90), "buckland")

## With confidence intervals according to Kang et al. (2000)
maED(ryegrass.m1, 
list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)), 
c(10, 50, 90), "kang")

## Comparing to model-averaged ED values with simple linear regression included
maED(ryegrass.m1, 
list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)), 
c(10, 50, 90), interval = "buckland", linreg = TRUE)



## Example with a model fit involving two compounds/curves
S.alba.m1 <- drm(DryMatter~Dose, Herbicide, data=S.alba, fct = LL.4(),
pmodels=data.frame(Herbicide,1,1,Herbicide)) 

## Model-averaged ED50 for both compounds
maED(S.alba.m1, list(LL.3(), LN.4()), 50)

## Model-averaged ED50 only for one compound (glyphosate)
maED(S.alba.m1, list(LL.3(), LN.4()), 50, clevel="Glyphosate")

## With confidence intervals
maED(S.alba.m1, list(LL.3(), LN.4()), 50, interval="buckland")

## For comparison model-specific confidence intervals
ED(S.alba.m1, 50, interval="delta")  # wider!

}

\keyword{models}
\keyword{nonlinear}