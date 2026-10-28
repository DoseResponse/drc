"llogistic" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL, 
fctName, fctText)
{
    ## Matching 'adjust' argument
    method <- match.arg(method)
    
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

#    ## Defining the basic non-linear function
#    bfct <- function(x, parm)
#    {
#        parm[2] + (parm[3]-parm[2])/((1+(x/parm[4])^parm[1]))^parm[5]
#    }

    ## Defining the model function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
#        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose)-log(parmMat[, 4]))))^parmMat[, 5])
        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose/parmMat[, 4]))))^parmMat[, 5])
    }
       
    ## Defining the model function adjusted for scaling
    retFct <- function(doseScaling, respScaling)
    {          
        fct <- function(dose, parm) 
        {
            parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), 
                              nrow(parm), numParm, byrow = TRUE)
            parmMat[, notFixed] <- parm
        
            cParm <- parmMat[, 2]
            cParm + (parmMat[, 3] - cParm)/((1 + exp(parmMat[, 1]*(log(dose / parmMat[, 4]))))^parmMat[, 5])
        }
        fct        
    }

    ## Defining the derivative in x adjusted for scaling
    retFctDx <- function(doseScaling, respScaling)
    {          
      fct <- function(dose, parm) 
      {
#        print(dose[1:10])
        parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), 
                          nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
#        bNeg <- parmMat[, 1] < 0
#        parmMat[bNeg, 1] <- -parmMat[bNeg, 1]
        
        temp1 <- dose/parmMat[, 4]          
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec <- (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
#        retVec[bNeg] <- -retVec[bNeg]
        retVec
      }
      fct        
    }
    
        
if (FALSE) {  ## will work once plotFct does not depend on drcFct
    ## Defining the model function adjusted for scaling
    retFct <- function(doseScaling, respScaling, numObs)
    {   
        parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), numObs, numParm, byrow = TRUE)
        
        fct <- function(dose, parm) 
        {        
            parmMat[, notFixed] <- parm        
            cParm <- parmMat[, 2]
            cParm + (parmMat[, 3] - cParm)/((1 + exp(parmMat[, 1]*(log(dose / parmMat[, 4]))))^parmMat[, 5])
        }
        fct        
    }
}
        
    ## Defining the scale function
    scaleFct <- function(doseScaling, respScaling)
    {        
        c(1, respScaling, respScaling, doseScaling, 1)[notFixed]
    }    

    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- llogistic.ssf(method, fixed)
    }
   

    ## Defining names
    names <- names[notFixed]
    
    ## Defining the first derivatives (in the parameters)        
    deriv1 <- function(dose, parm) {

        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
#        t3 <- (1 + t2)^(2*parmMat[, 5])
#        t4 <- parmMat[, 5]*((1 + t2)^(parmMat[, 5] - 1))
##        t3 <- parmMat[, 5]*((1 + t2)^(-parmMat[, 5] - 1))
        t5 <- (1 + t2)^parmMat[, 5]                  

        cbind( -t1 * xlogx(dose/parmMat[, 4], parmMat[, 1], parmMat[, 5] + 1) * parmMat[, 5],  # *t4/t3, 
               1 - 1/t5, 
               1/t5, 
#               t1*t2*t4*parmMat[, 1]/parmMat[, 4]/t3, 
               t1 * parmMat[, 5] * divAtInf(t2, (1 + t2)^(parmMat[, 5] + 1)) * parmMat[, 1] / parmMat[, 4], 
               -t1 * divAtInf(log(1+t2), t5) )[, notFixed]
    }
       
       
    deriv2 <- NULL


    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
#        bNeg <- parmMat[, 1] < 0
#        parmMat[bNeg, 1] <- -parmMat[bNeg, 1]
                  
        temp1 <- x/parmMat[, 4]          
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec <- (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
#        retVec[bNeg] <- -retVec[bNeg]
        retVec
    }

  
    ## The three definitions below are not needed in future ('drm')

#    ## Defining parameter to be scaled
#    if (is.na(fixed[4]))  #  (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }
#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos=conPos, fct=confct2))
#    }
#    ## Defining flag to indicate if more general ANOVA model
##    anovaYes <- list(bin = !any(is.na(fixed[c(2,3,5)])) , cont = TRUE)
#    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}
#    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm

        ## Converting absolute to relative if needed
        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- log(1 - p/100)

#        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])
#
#        EDder <- EDp * c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
#                         0, 0, 1 / parmVec[4], 
#        exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))
#
#        EDderxFct0 <- deriv(~e*(exp(-log((d - p)/(d - c))/f)-1)^(1/b), "p", function(p,b,c,d,e,f){})
#        EDderxFct <- function(p, parm) {
#            EDderxTemp <- as.vector(attr(EDderxFct0(p, parm[1], parm[2], parm[3], parm[4], parm[5]), "gradient"))
#            EDderxTemp * c(0, (p - parm[3]) / ((parm[3] - parm[2])^2), (parm[2] - p) / ((parm[3] - parm[2])^2), 0, 0)
#            # based on the chain rule
#        }
#        if (type == "absolute") {EDder <- EDder + EDderxFct(respl, parmVec)}

        # Using: deriv(~e * (exp(-log((y-c)/(d-c))/f)-1)^(1/b), 
        #        c("b", "c", "d", "e", "f"), function(b,c,d,e,f,y){})
        derED <- function (b, c, d, e, f) {
    
            #.expr1 <- y - c
            .expr2 <- d - c
            
            .expr4 <- tempVal  # log(.expr3)
            .expr3 <- exp(.expr4)  # .expr1/.expr2
            .expr7 <- exp(-.expr4/f)
            .expr8 <- .expr7 - 1
            .expr9 <- 1/b
            .expr10 <- .expr8^.expr9
            .expr20 <- .expr8^(.expr9 - 1)
            .expr23 <- .expr3 / .expr2  # .expr1/.expr2^2
            
            .value <- e * .expr10
            .grad <- array(0, c(length(.value), 5L), 
                           list(NULL, c("b", "c", "d", "e", "f")))
    
            .grad[, "b"] <- -(e * (.expr10 * (log(.expr8) * (1/b^2))))
            .grad[, "c"] <- e * (.expr20 * (.expr9 * (.expr7 * ((1/.expr2 - .expr23)/.expr3/f))))
            .grad[, "d"] <- e * (.expr20 * (.expr9 * (.expr7 * (.expr23/.expr3/f))))
            .grad[, "e"] <- .expr10
            .grad[, "f"] <- e * (.expr20 * (.expr9 * (.expr7 * (.expr4/f^2))))
    
            attr(.value, "gradient") <- .grad
            .value
        }

        EDp <- derED(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDder <- attr(EDp, "gradient")
        if (!identical(type, "absolute")) {EDder[c(2, 3)] <- 0}

        return(list(as.vector(EDp), as.vector(EDder[notFixed])))
    }
    

    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))
    } 
    

    ## Defining functions returning lower and upper limit and monotonicity
    lowerAs <- pickParm(parmVec, notFixed, 2)
    upperAs <- pickParm(parmVec, notFixed, 3)
    monoton <- monoParm(parmVec, notFixed, 1, -1)


    ## Returning the function with self starter and names
    returnList <- list(
        fct = fct, ssfct = ssfct, names = names, 
        deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
        edfct = edfct, inversion = invfct, scaleFct = scaleFct,
#    scaleInd = scaleInd, confct=confct, anovaYes=anovaYes, lowerc=lowerLimits, upperc=upperLimits,
        name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text = ifelse(missing(fctText), "Log-logistic (ED50 as parameter)", fctText), 
        noParm = sum(is.na(fixed)), lowerAs = lowerAs, upperAs = upperAs, monoton = monoton,
        retFct = retFct, fixed = fixed, retFctDx = retFctDx)
    # the 4th last line is not needed in the future ('drm')
    #     , sxInd = sxInd, syInd = syInd,    
    
    class(returnList) <- "llogistic"
    invisible(returnList)
}


"LL.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...) {
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, upper, fixed[2], 1), 
    names = c(names[1], "c", "d", names[2], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Log-logistic (ED50 as parameter)", upper),
    ...) )
}


"LL.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...) {
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, fixed[2:3], 1), 
    names = c(names[1], "c", names[2:3], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Log-logistic (ED50 as parameter)"), 
    ...) )
}


"LL.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...) {
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1:2], upper, fixed[3], 1), 
    names = c(names[1:2], "d", names[3], "f"),
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Log-logistic (ED50 as parameter)", upper), 
    ...) )
}


"LL.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...) {
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed, 1), names = c(names, "f"),
    fctName = as.character(match.call()[[1]]), ...) )
}


"LL.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...) {
    return( llogistic(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]),
    fctText = "Generalized log-logistic (ED50 as parameter)", ...) )
}


"MM.2" <-
function(fixed = c(NA, NA), names = c("d", "e"), ...) {
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(-1, 0, fixed[1:2], 1), names = c("b", "c", names[1:2], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Michaelis-Menten", 
    ...) )
}


"MM.3" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...) {
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(-1, fixed[1:3], 1), names = c("b", names[1:3], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Shifted Michaelis-Menten", 
    ...) )
}
