\name{W1.4}

\alias{W1.4}
\alias{w4}
\alias{W2.4}
\alias{W2x.4}

\title{The four-parameter Weibull functions}

\description{
  'W1.4' and 'W2.4' provide the four-parameter Weibull functions, self starter function and 
  names of the parameters.   
}

\usage{
  W1.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
  
  W2.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{additional arguments to be passed from the convenience functions.}   
}

\details{
  The equations for the mean functions are given at \code{\link{weibull1}}.
}

\value{
  See \code{\link{weibull1}}.
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).
  
  Ritz, C (2009)
  Towards a unified approach to dose-response modeling in ecotoxicology
  \emph{To appear in Environ Toxicol Chem}.
}

\author{Christian Ritz}

\note{This function is for use with the model fitting function \code{\link{drm}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{W1.3}}. A more flexible function, allowing 
fixing or constraining parameters, is available through \code{\link{weibull1}}.}

\examples{

## Fitting a four-parameter Weibull (type 1) model
terbuthylazin.m1 <- drm(rgr~dose, data = terbuthylazin, fct = W1.4())
summary(terbuthylazin.m1)

## Fitting a first-order multistage model
## to data from BMDS by EPA
## (Figure 3 in Ritz (2009))
bmds.ex1 <- data.frame(ad.dose=c(0,50,100), dose=c(0, 2.83, 5.67), 
num=c(6,10,19), total=c(50,49,50))

bmds.ex1.m1<-drm(num/total~dose, weights=total, data=bmds.ex1, 
fct=W2.4(fixed=c(1,NA,1,NA)), type="binomial")

modelFit(bmds.ex1.m1)  # same as in BMDS

summary(bmds.ex1.m1)  # same background estimate as in BMDS

logLik(bmds.ex1.m1)

## BMD estimate identical to BMDS result
## BMDL estimate differs from BMDS result (different method)
ED(bmds.ex1.m1, 10, ci="delta")

## Better fit

bmds.ex1.m2<-drm(num/total~dose, weights=total, data=bmds.ex1, 
fct=W1.4(fixed=c(-1,NA,1,NA)), type="binomial")
modelFit(bmds.ex1.m2)
summary(bmds.ex1.m2)

ED(bmds.ex1.m2, 50, ci = "delta")

## Creating Figure 3 in Ritz (2009)
bmds.ex1.m3 <- drm(num/total~dose, weights=total, data=bmds.ex1, 
fct=LL.4(fixed=c(-1,NA,1,NA)), type="binomial")

plot(bmds.ex1.m1, ylim = c(0.05, 0.4), log = "", lty = 3, lwd = 2, 
xlab = "Dose (mg/kg/day)", ylab = "",
cex=1.2, cex.axis=1.2, cex.lab=1.2)

mtext("Tumor incidence", 2, line=4, cex=1.2)  # tailored y axis label

plot(bmds.ex1.m2, ylim = c(0.05, 0.4), log = "", add = TRUE, lty = 2, lwd = 2)

plot(bmds.ex1.m3, ylim = c(0.05, 0.4), log = "", add = TRUE, lty = 1, lwd = 2)

arrows(2.6 , 0.14, 2, 0.14, 0.15, lwd=2)
text(2.5, 0.14, "Weibull-1", pos=4, cex=1.2)

}

\keyword{models}
\keyword{nonlinear}
