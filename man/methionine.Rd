%\encoding{latin1}

\name{methionine}

\alias{methionine}

\docType{data}

\title{Weight gain for different methionine sources}

\description{
  Data consist of average body weight gain of chickens being treated 
  with one of the two methionine sources DLM and HMTBA.
}

\usage{data(methionine)}

\format{
  A data frame with 9 observations on the following 3 variables:
  \describe{
    \item{\code{product}}{a factor with levels \code{control}, \code{DLM} and \code{MHA} denoting the treatments}
    \item{\code{dose}}{a numeric vector of methionine dose}
    \item{\code{gain}}{a numeric vector of average body weight gain}
  }
}

\details{
  The dataset contains a common control measurement for the two treatments. More examples using this
  dataset are found under \code{\link{AR.2}} and \code{\link{MM.2}}.
}

\source{
  Kratzer. D. D. and Littell, R. C. (2006) Appropriate Statistical Methods to Compare 
  Dose Responses of Methionine Sources, \emph{Poultry Science}, \bold{85}, 947--954.
}

%\references{
%  Schutte, J. B. and de Jong, J. (1996) Biological efficacy of DL-methinonine hydroxy analog free acid compared
%  to DL-methionine in broiler chicks as determined by performance and breast meat yield,
%  \emph{Agribiol. Res.}, \bold{49}, 74--82.
%}

\examples{

## Fitting model with constraint on one parameter 
met.ar.m1 <- drm(gain~dose, product, data = methionine, 
fct = AR.3(), pmodels = list(~1, ~factor(product), ~factor(product)), 
upperl = c(Inf, Inf, 1700, Inf, Inf)) 

plot(met.ar.m1, xlim=c(0,0.3), ylim=c(1450, 1800))
abline(h=1700, lty=1)

summary(met.ar.m1)

}
\keyword{datasets}
