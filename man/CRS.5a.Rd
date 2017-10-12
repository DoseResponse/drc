\name{CRS.5a}

\alias{CRS.5a}
\alias{CRS.5b}
\alias{CRS.5c}

\alias{ml4a}
\alias{ml4b}
\alias{ml4c}

\alias{UCRS.5a}
\alias{UCRS.5b}
\alias{UCRS.5c}

\alias{uml4a}
\alias{uml4b}
\alias{uml4c}


\title{Cedergreen-Ritz-Streibig dose-reponse model for describing hormesis}

\description{
  'CRS.5a', 'CRS.5b' and 'CRS.5c' provide the Cedergreen-Ritz-Streibig modified log-logistic model for describing
  (inverse u-shaped or j-shaped) hormesis.
  
  'UCRS.5a', 'UCRS.5b' and 'UCRS.5c' provide the Cedergreen-Ritz-Streibig modified log-logistic model for
  describing u-shaped hormesis. 
}

\usage{
  CRS.5a(names = c("b", "c", "d", "e", "f"), ...)

  UCRS.5a(names = c("b", "c", "d", "e", "f"), ...)
}

\arguments{
  \item{names}{a vector of character strings giving the names of the parameters.}
  \item{...}{additional arguments to be passed from the convenience functions.}  
}

\details{
  The model function for inverse u-shaped hormetic patterns is
  
  \deqn{ f(x) = c + \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}},

  which is a five-parameter model. It is a modification of the four-parameter log-logistic curve 
  to take hormesis into account.
  
  The parameters have the following interpretations
  \itemize{
    \item \eqn{b}: Not direct interpretation
    \item \eqn{c}: Lower horizontal asymptote
    \item \eqn{d}: Upper horizontal asymptote
    \item \eqn{e}: Not direct interpretation
    \item \eqn{f}: Size of the hormesis effect: the larger the value the larger is the hormesis effect. \eqn{f=0}
    corresponds to no hormesis effect and the resulting model is the four-parameter log-logistic model. 
    This parameter should be positive in order for the model to make sense.
  }  
  
  The model function for u-shaped hormetic patterns is
  \deqn{ f(x) = c + d - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
  This model also simplifies to the four-parameter log-logistic model in case \eqn{f=0} (in a slightly
  different parameterization as compared to the one used in \code{\link{LL.4}}).
  
  The models denoted a,b,c are obtained by fixing the alpha parameter at 1, 0.5 and 0.25, respectively.  
}

\value{
  See \code{\link{cedergreen}}.
}

\references{
  See the reference under \code{\link{cedergreen}}.
}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.  
}

\seealso{
  Similar functions are \code{\link{CRS.4a}} and \code{\link{UCRS.4a}}, but with the
  lower limit (the parameter \eqn{c}) fixed at 0 (one parameter less to be estimated).
}

\examples{

## Modified logistic model
lettuce.m1 <- drm(lettuce[,c(2,1)], fct=CRS.5a())
summary(lettuce.m1)
ED(lettuce.m1, c(50))

lettuce.m2 <- drm(lettuce[,c(2,1)], fct=CRS.5b())
summary(lettuce.m2)
ED(lettuce.m2, c(50))

lettuce.m3 <- drm(lettuce[,c(2,1)], fct=CRS.5c())
summary(lettuce.m3)
ED(lettuce.m3, c(50))

}
\keyword{models}
\keyword{nonlinear}
