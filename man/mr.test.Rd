\name{mr.test}

\Rdversion{1.1}

\alias{mr.test}

\title{
  Mizon-Richard test for dose-response models 
}

\description{
  The function provides a lack-of-fit test for the mean structure based on the Mizon-Richard test as compared to a
  specific alternative model.
}

\usage{
mr.test(object1, object2, object, x, var.equal = TRUE, component = 1)
}

\arguments{
  \item{object1}{
  object of class 'drc' (null model).
}
  \item{object2}{
  object of class 'drc' (alternative model).
}
  \item{object}{
  object of class 'drc' (fitted model under alternative).
}
  \item{x}{
  numeric vector of dose values.
}
  \item{var.equal}{
  logical indicating whether or not equal variances can be assumed across doses.
}
  \item{component}{
  numeric vector specifying the component(s) in the parameter vector to use in the test.
}
}

\details{
  The function provides a p-value indicating whether or not the mean structure is appropriate.
  
  The test is applicable even in cases where data are non-normal or exhibit variance heterogeneity.
}

\value{
  A p-value for test of the null hypothesis that the chosen mean structure is appropriate as compared
  to the alternative mean structure provided (see Ritz and Martinussen (2011) for a detailed explanation).
}

\references{
  
  Ritz, C and Martinussen, T. (2011)   
  Lack-of-fit tests for assessing mean structures for continuous dose-response data,
 \emph{Environmental and Ecological Statistics}, \bold{18}, 349--366
}

\author{
  Christian Ritz
}

\note{
  This functionality is still experimental: Currently, the null and alternative models are hardcoded!
  In the future the function will be working for null and alternative models specified by the user.
}

\seealso{
  See also \code{\link{modelFit}} for details on the related lack-of-fit test against an ANOVA model.
}
\examples{

## Fitting log-logistic and Weibull models 
## The Weibull model is the alternative
etmotc.m1<-drm(rgr1~dose1, data=etmotc[1:15,], fct=LL.4())
etmotc.m2 <- update(etmotc.m1, fct=W1.4())

## Fitting the fitted model (using the alternative model)
etmotc.m3 <- drm(fitted(etmotc.m1)~dose1, data=etmotc[1:15,], fct=W1.4())

## Handling missing values
xVec <- etmotc[1:15,]$dose1
xVec[1:8] <- 1e-10  # avoiding 0's

## Obtaining the Mizon-Richard test
mr.test(etmotc.m1, etmotc.m2, etmotc.m3, xVec, var.equal = FALSE)

}

\keyword{models}
\keyword{nonlinear}
