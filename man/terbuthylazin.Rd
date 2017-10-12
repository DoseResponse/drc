\name{terbuthylazin}

\alias{terbuthylazin}

\docType{data}

\title{The effect of terbuthylazin on growth rate}

\description{
  Test on the effect of terbuthylazin on  \emph{Lemna minor}, performed on an aseptic
  culture according to the OECD-guidelines.  
}

\usage{data(terbuthylazin)}

\format{
  A data frame with 30 observations on the following 2 variables.
  \describe{
    \item{dose}{a numeric vector of dose values.}
    \item{rgr}{a numeric vector of relative growth rates.}
  }
}

\details{
  Dose is \deqn{\mu l^{-1}} and rgr is the relative growth rate of \emph{Lemna}.  
}

\source{
  Cedergreen N. (2004). Unpublished bioassay data.
}

\examples{

## displaying first 6 rows of the data set
head(terbuthylazin)

## Fitting log-logistic model
terbuthylazin.m1 <- drm(rgr~dose, data = terbuthylazin, fct = LL.4())
summary(terbuthylazin.m1)

## Fitting log-logistic model
##  with Box-Cox transformation
terbuthylazin.m2 <- boxcox(terbuthylazin.m1, method = "anova")
summary(terbuthylazin.m2)

}
\keyword{datasets}
