\name{drm}

\alias{drm}

\title{Fitting dose-response models}

\description{
  A general model fitting function for analysis of various types of dose-response data.
}

\usage{
  drm(formula, curveid, pmodels, weights, data = NULL, subset, fct, 
  type = c("continuous", "binomial", "Poisson", "negbin1", "negbin2", "event", "ssd"), 
  bcVal = NULL, bcAdd = 0, 
  start, na.action = na.omit, robust = "mean", logDose = NULL, 
  control = drmc(), lowerl = NULL, upperl = NULL, separate = FALSE, 
  pshifts = NULL, varcov = NULL)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. Either of the form 'response \eqn{~} dose' 
  or as a data frame with response values in first column and dose values in second column.}
  
  \item{curveid}{a numeric vector or factor containing the grouping of the data.}
  
  \item{pmodels}{a data frame with a many columns as there are parameters in the non-linear function. 
  Or a list containing a formula for each parameter in the nonlinear function.}
  
  \item{weights}{
    a numeric vector containing weights. For continuous/quantitative responses weights are multiplied inside the squared errors (see the details below).
    For binomial reponses weights provide information about the total number of binary observations used to obtain the response (which is a proportion): 
    1/2 and 10/20 lead to different analyses due to the different totals (2 vs. 20) even though the proportion in both cases is 0.5.
  }
  
  \item{data}{
    an optional data frame containing the variables in the model.}
  
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  
  \item{fct}{a list with three or more elements specifying the non-linear
  function, the accompanying self starter function, the names of the parameter in the non-linear function and, 
  optionally, the first and second derivatives as well as information used for calculation of ED values.
  Currently available functions include, among others, the four- and five-parameter log-logistic models 
  \code{\link{LL.4}}, \code{\link{LL.5}} and the Weibull model \code{\link{W1.4}}. Use
  \code{\link{getMeanFunctions}} for a full list.}
     
  \item{type}{a character string specifying the distribution of the data (parameter estimation will depend on the assumed distribution as different log likelihood functions will be used). 
  The default is "continuous", corresponding to assuming a normal distribution. The choices "binary" and "event" imply a binomial and multinomial distribution, respectively. The choice "ssd" is for fitting a species sensitivity
  distribution.}
  
  \item{bcVal}{a numeric value specifying the lambda parameter to be used in the Box-Cox transformation.}
  
  \item{bcAdd}{a numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}  
  
  \item{start}{an optional numeric vector containing starting values for all mean parameters in the model. 
  Overrules any self starter function.}
   
  \item{na.action}{a function for treating mising values ('NA's). Default is \code{\link{na.omit}}.}
   
  \item{robust}{a character string specifying the rho function for robust estimation. Default is non-robust 
  least squares estimation ("mean"). Available robust methods are: median estimation ("median"), 
  least median of squares ("lms"), least trimmed squares ("lts"), metric trimming ("trimmed"), 
  metric winsorizing ("winsor") and Tukey's biweight ("tukey").}
  
  \item{logDose}{a numeric value or NULL. If log doses value are provided the base of the logarithm should be specified (exp(1) for the natural logarithm and 10 for the base 10 logarithm).}

  \item{control}{a list of arguments controlling constrained optimisation (zero as boundary), 
  maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}
  
  \item{lowerl}{a numeric vector of lower limits for all parameters in the model 
  (the default corresponds to minus infinity for all parameters).}  
  
  \item{upperl}{a numeric vector of upper limits for all parameters in the model 
  (the default corresponds to plus infinity for all parameters).}
  
  \item{separate}{logical value indicating whether curves should be fit separately (independent of each other).}
  
  \item{pshifts}{a matrix of constants to be added to the matrix of parameters. Default is no shift for all parameters.}
  
  \item{varcov}{an optional user-defined known variance-covariance matrix for the responses. Default is the identity matrix (NULL), corresponding to independent response values with a common standard deviation, which will be estimated from the data.}  
} 

\details{
  Response values are assumed to be mutually independent unless the argument \code{varcov} is specified (in which the correlation structure is assumed to be completely known). 

  This function relies on the general optimiser function \code{\link{optim}} for the minimisation of negative log likelihood function. For a continuous response this reduces to least squares estimation, which is carried out by minimising the following sums of squares
  
  \deqn{\sum_{i=1}^N [w_i (y_i-f_i)]^2}
  
  where \eqn{y_i}, \eqn{f_i}, and \eqn{w_i} correspond to the observed value, expected value, and the weight respectively, for the ith observation (from 1 to \eqn{N}).

  For robust estimation MAD (median abslolute deviance) is used to estimate the residual variance.  
  
  Setting \code{lowerl} and/or \code{upperl} automatically invokes constrained optimisation.
  
  The columns of a data frame argument to \code{pmodels} are automatically converted into factors. 
  This does not happen if a list is specified.   

  The control arguments are specified using the function \code{\link{drmc}}.
  
}

\value{
  An object of (S3) class 'drc'. 
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz and Jens C. Streibig}

%\note{}

%\examples{}

\keyword{models}
\keyword{nonlinear}
