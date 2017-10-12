\name{MAX}

\alias{MAX}

\title{Maximum mean response}

\description{
  \code{MAX} estimates the maximum mean response and the dose at which it occurs.
}

\usage{
MAX(object, lower = 1e-3, upper = 1000, pool = TRUE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{lower}{numeric. Lower limit for bisection method. Need to be smaller than EDx level to be calculated.}  
  \item{upper}{numeric. Upper limit for bisection method. Need to be larger than EDx level to be calculated.}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.} 
}

\details{
  This function is only implemented for the built-in functions of class \code{\link{braincousens}} and 
  \code{\link{cedergreen}}.
  
  This function was used for obtaining the results on hormesis effect size reported in Cedergreen et al. (2005).
}

\value{
  A matrix with one row per curve in the data set and two columns: 
  one containing the dose at which the maximum occurs 
  and one containing the corresponding maximum response.
}

\references{
  Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved empirical models describing hormesis, 
  \emph{Environmental Toxicology and Chemistry} \bold{24}, 3166--3172. 
}

\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens' and 'mlogistic'.}

%\seealso{The related function \code{\link{SI}}.}

\examples{

## Fitting a Cedergreen-Ritz-Streibig model
lettuce.m1 <- drm(weight~conc, data = lettuce, fct = CRS.4c())

## Finding maximum average response and the corrresponding dose
MAX(lettuce.m1)

}
\keyword{models}
\keyword{nonlinear}
