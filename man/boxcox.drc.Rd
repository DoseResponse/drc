\name{boxcox.drc}

\alias{boxcox.drc}

\title{Transform-both-sides Box-Cox transformation}

\description{
  Finds the optimal Box-Cox transformation for non-linear regression.
}

\usage{
\method{boxcox}{drc}(object, lambda = seq(-2, 2, by = 0.25), plotit = TRUE, bcAdd = 0, 
method = c("ml", "anova"),  level = 0.95, eps = 1/50, 
xlab = expression(lambda), ylab = "log-Likelihood", ...)
}

\arguments{
  \item{object}{object of class \code{drc}.}
  \item{lambda}{numeric vector of lambda values; the default is (-2, 2) in steps of 0.25.}
  \item{plotit}{logical which controls whether the result should be plotted.}
  \item{bcAdd}{numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}
  \item{method}{character string specifying the estimation method for lambda: maximum likelihood or ANOVA-based 
  (optimal lambda inherited from more general ANOVA model fit.}
  \item{eps}{numeric value: the tolerance for lambda = 0; defaults to 0.02.}
  \item{level}{numeric value: the confidence level required.}
  \item{xlab}{character string: the label on the x axis, defaults to "lambda".}
  \item{ylab}{character string: the label on the y axis, defaults to "log-likelihood".}
  \item{\dots}{additional graphical parameters.}
}

\details{
  The optimal lambda value is determined using a profile likelihood approach: 
  For each lambda value the dose-response regression model is fitted and the lambda value (and corresponding model fit) resulting in the largest
  value of the log likelihood function is chosen.
}

\value{
  An object of class "drc" (returned invisibly). 
  If plotit = TRUE a plot of loglik vs lambda is shown indicating a confidence interval (by default 95%) about 
  the optimal lambda value.
}

\references{
    Carroll, R. J. and Ruppert, D. (1988) \emph{Transformation and Weighting in Regression}, 
    New York: Chapman and Hall (Chapter 4).
}

\author{Christian Ritz}

%\note{}

\seealso{
  For linear regression the analogue is \code{\link[MASS]{boxcox}}. 
}

\examples{

## Fitting log-logistic model without transformation
ryegrass.m1 <- drm(ryegrass, fct = LL.4())
summary(ryegrass.m1)

## Fitting the same model with the optimal Box-Cox transformation
ryegrass.m2 <- boxcox(ryegrass.m1)
summary(ryegrass.m2)

}

\keyword{models}
\keyword{nonlinear}
