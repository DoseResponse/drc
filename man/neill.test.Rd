\name{neill.test}

\alias{neill.test}

\title{
  Neill's lack-of-fit test for dose-response models
}

\description{
  'neill.test' provides a lack-of-fit test for non-linear regression models. It is applicable both in cases 
  where there are replicates (in which case it reduces to the standard lack-of-fit test against an ANOVA
  model) and in cases where there are no replicates, though then a grouping has to be provided.
}

\usage{
  neill.test(object, grouping, method = c("c-finest", "finest", "percentiles"), 
  breakp = NULL, display = TRUE)
}

\arguments{
  \item{object}{
  object of class 'drc' or 'nls'.
}
  \item{grouping}{
  character or numeric vector that provides the grouping of the dose values.
}
  \item{method}{
  character string specifying the method to be used to generate a grouping of the dose values.
}
  \item{breakp}{
  numeric vector of break points for generating dose intervals that form a grouping.
}
  \item{display}{
  logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).
}
}

\details{
  The functions used the methods \code{\link{df.residual}} and \code{\link{residuals}} and the 'data'
  component of \code{object} (only for determining the number of observations).
}

\value{
  The function returns an object of class anova which is displayed using \code{print.anova}.
}

\references{
  Neill, J. W. (1988) Testing for lack of fit in nonlinear regression, 
  \emph{Ann. Statist.}, \bold{16}, 733--740  
}

\author{Christian Ritz}

\note{
  A clustering technique could be employed to determine the grouping to be used in cases where there are
  no replicates. There should at most be ceiling(n/2) clusters as otherwise some observations will not be used
  in the test. At the other end there need to be more clusters than parameters in the model.
}

\seealso{
  See also \code{\link{modelFit}} for details on the lack-of-fit test against an ANOVA model.}

\examples{

### Example with 'drc' object

## Lack-of-fit test against ANOVA
ryegrass.m1 <-drm(rootl~conc, data = ryegrass, fct = LL.4())
modelFit(ryegrass.m1)

## The same test using 'neill.test'
neill.test(ryegrass.m1, ryegrass$conc)

## Generating a grouping
neill.test(ryegrass.m1, method="c-finest")
neill.test(ryegrass.m1, method="finest")
neill.test(ryegrass.m1, method="perc")


}

\keyword{models}
\keyword{nonlinear}
