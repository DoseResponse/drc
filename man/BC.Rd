\name{BC.5}

\alias{BC.5}
\alias{BC.4}

\title{The Brain-Cousens hormesis models}

\description{
  'BC.4' and 'BC.5' provide the Brain-Cousens modified log-logistic models for describing u-shaped hormesis.
}

\usage{
  BC.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)

  BC.4(fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"), ...)
}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at which values they are fixed. 
  NAs designate parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters.}
  \item{...}{additional arguments to be passed on to convenience functions.} 
}

\details{
  The model function for the Brain-Cousens model is
  
  \deqn{ f(x, b,c,d,e,f) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}.
  
  It is a five-parameter model, obtained by extending the four-parameter log-logistic model (\code{\link{LL.4}}
  to take into account inverse u-shaped hormesis effects (Brain and Cousens, 1989).
  
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
  
  Fixing the lower limit at 0 results in the four-parameter model: 
  
  \deqn{ f(x) = 0 + \frac{d-0+fx}{1+\exp(b(\log(x)-\log(e)))}}
  
  used by van Ewijk and Hoekstra (1993).
}

\value{
  See \code{\link{braincousens}}.
}

\references{
  Brain, P. and Cousens, R. (1989) 
  An equation to describe dose responses where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.

  van Ewijk, P. H. and Hoekstra, J. A. (1993) 
  Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
  \emph{Ecotoxicology and Environmental Safety}, \bold{25}, 25--32.  
}

\author{Christian Ritz}

\note{This model function is for use with the function \code{\link{drm}}.}

\seealso{
  More details are found for the general model function \code{\link{braincousens}}.
}

\examples{

## Fitting a five-parameter BC model to data from van Ewijk and Hoekstra (1993)
lettuce.BC.m1 <- drm(weight ~ conc, data = lettuce, fct = BC.5())
modelFit(lettuce.BC.m1)
plot(lettuce.BC.m1)

lettuce.BC.m2 <- drm(weight ~ conc, data = lettuce, fct = BC.4())
summary(lettuce.BC.m2)
ED(lettuce.BC.m2, c(10, 20, 50))  
# compare the parameter estimate and its estimated standard error 
#   to the values in the paper by van Ewijk and Hoekstra (1993)


## Fitting a five-parameter BC model to the ryegrass data
ryegrass.BC.m1 <- drm(rootl ~ conc, data = ryegrass, fct = BC.5())
summary(ryegrass.BC.m1)
# No hormesis effect: the estimate of f is equal 0.357 (p=0.61)
ED(ryegrass.BC.m1, c(10, 20, 50, 80, 90))

## Fitting a model with the constraint b>3
ryegrass.BC.m2 <- drm(rootl ~ conc, data = ryegrass, fct = BC.5(), 
lower = c(3, -Inf, -Inf, -Inf, -Inf), control = drmc(constr = TRUE))

summary(ryegrass.BC.m2)

## Fitting a model with the constraint f>0 
## (no effect as the estimate of f is positive anyway)
ryegrass.BC.m3 <- drm(rootl ~conc, data = ryegrass, fct = BC.5(), 
lower = c(-Inf, -Inf, -Inf, -Inf, 0), control = drmc(constr = TRUE))

summary(ryegrass.BC.m3)

}

\keyword{models}
\keyword{nonlinear}
