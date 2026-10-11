"ucedergreen" <- function(
    fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
    method = c("1", "2", "3", "4"), ssfct = NULL,
    alpha, fctName, fctText) {

    numParm <- 5
    if (!is.character(names) || !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}
    
    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec

    ## Defining the function
    fct <- function(dose, parm) {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        numTerm <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*exp(-1/dose^alpha)
        denTerm <- 1 + exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        parmMat[, 3] - numTerm/denTerm
    }

    ## Defining self starter function
    if (!is.null(ssfct)){
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe) {
            initval <- llogistic(method = method)$ssfct(dframe)   
            initval[1] <- -initval[1]
            initval[5] <- 0  # better solution?
    
            return(initval[notFixed])
        }
    }   


    ## Specifying the derivatives    
    deriv1OLD <- function(dose, parm)
    # not used anymore
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t0 <- exp(-1/(dose^alpha))
        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*t0
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1/t3, 
               1 - 1/t3, 
               -t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               -t0/t3 )[, notFixed]
    }


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }


    ## Defining names
    names <- names[notFixed]


    ## Defining the first derivative (in the parameters)
    deriv1 <- function(dose, parm){
        
        deriv1FctTemp <- function(dose, pmRow) {
            deriv1Fct <- deriv(~d-(d-c+f*exp(-1/x^alpha))/(1+(x/e)^b), c("b", "c", "d", "e", "f"), 
                               function(x,b,c,d,e,f){})
            dVal <- attr(deriv1Fct(dose, pmRow[1], pmRow[2], pmRow[3], pmRow[4], pmRow[5]), "gradient")
            dVal[is.na(dVal)] <- 0
            # NaN's for b and e are due to a power-log term not being well-defined; it should return 0
            dVal
        }

        #notFixed <- rep(TRUE, 5)
        #parmVec <- c(NA, NA, NA, NA, NA)
        #numParm <- 5

        nrpar <- nrow(parm)
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        derivMat <- matrix(NA, nrpar, numParm)
        for (i in 1:nrpar) {
            derivMat[i, ] <- deriv1FctTemp(dose[i], parmMat[i, ])
        }
        derivMat[, notFixed]
    }
        
    deriv2 <- NULL


    ## Defining the first derivative (in the dose)
    derivx <- function(dose, parm) {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t2 <- 1 + t1
        t3 <- exp(-1/(dose^alpha)) * alpha * (dose^(-alpha - 1))

        -(parmMat[, 5] * t3 / t2 - (parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*exp(-1/(dose^alpha))) * t1 * parmMat[, 1] / (t2^2 * dose) + parmMat[, 5] * t3 / t2)
    }
    # Identical to derivx in cedergreen.R, but with a - in front of the last term


    ## Defining the ED function    
    edfct <- function(parm, respl, reference, type, lower = NULL, upper = NULL, ...) {

        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type, TRUE)       
        tempVal <- p / 100
        
        ## Defining the numerator of derivative of the BC model function wrt. dose
        derx <- function(dose) {
            expVal1 <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            expVal2 <- exp(-1/(dose^alpha)) * alpha * (dose^(-alpha - 1))

            -(parmVec[5]*(1+expVal1)*expVal2 - (parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))*expVal1*parmVec[1]/dose)
        
        }
        tempVal2 <- parmVec[3] - parmVec[2]
        nearzepVal1 <- parmVec[4] * ((tempVal2) / parmVec[5])^(1 / parmVec[1])
        nearzepval2 <- ifelse((tempVal2)/2 - parmVec[5]*exp(-(parmVec[4]/2)^(-alpha)) > 0, 
        parmVec[4]/2, ifelse((tempVal2)/2 - parmVec[5]*exp(-parmVec[4]^(-alpha)) > 0, 
        parmVec[4], NA))
        print(c(nearzepVal1, nearzepval2))
        nearzepVal <- min(nearzepVal1, nearzepval2, na.rm = TRUE)
        print(nearzepVal)

        if (is.null(lower)) {
            xvalues <- 10^(runif(100, -10, log10(nearzepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec>0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        }
      
        # Checking if the curve is effectively decreasing
        if (length(lower) == 0) {maxAt <- 1e-20} else {
            maxAt <- uniroot(derx, lower = lower, upper = nearzepVal)$root 
        } 
        if (is.null(upper)) {
            upper <- parmVec[4] * 10000
        }
        #print(maxAt)

        eqn <- function(dose) {
            (1 - tempVal)*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*exp(-1/(dose^alpha))/(parmVec[3]-parmVec[2]))}
        EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root

        tempVal1 <- exp(parmVec[1] * (log(EDp) - log(parmVec[4])))
        #tempVal2 <- parmVec[3] - parmVec[2]

        derNum <- c(-(1 - tempVal) * tempVal1 * (log(EDp) - log(parmVec[4])), 
                     -1 / tempVal2, 
                     -(tempVal * tempVal2 + parmVec[5] * exp(-1/(EDp^alpha))) / ((1 - tempVal) * tempVal2^2),
                     (1 - tempVal) * tempVal1 * parmVec[1] / parmVec[4],
                     exp(-1/(EDp^alpha)) / tempVal2)
        derDenom <- (1 - tempVal) * tempVal1 *parmVec[1] / EDp - (exp(-1/(EDp^alpha)) * alpha * (EDp^(-alpha - 1))) * parmVec[5] / tempVal2
        EDder <- derNum / derDenom
        
        return(list(EDp, EDder[notFixed]))
    }


    ## Finding the maximal hormesis
    maxfct <- function(parm, lower = NULL, noRan = 100) {
        # Note that the code is mostly copied from the edfct function above, 
        #   but with some changes to find the maximum hormesis

        parmVec[notFixed] <- parm
        # No parameter check yet, but should be added
        
         ## Defining the numerator of derivative of the CRS model function wrt. dose
        derx <- function(dose) {
            expVal1 <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            expVal2 <- exp(-1/(dose^alpha)) * alpha * (dose^(-alpha - 1))

            -(parmVec[5]*(1+expVal1)*expVal2 - (parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))*expVal1*parmVec[1]/dose)
        }
        nearzepVal <- parmVec[4] * ((parmVec[3]-parmVec[2]) / parmVec[5])^(1 / parmVec[1])

        if (is.null(lower)) {
            xvalues <- 10^(runif(noRan, -10, log10(nearzepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec<0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        } 
        maxAt <- uniroot(derx, lower = lower, upper = nearzepVal)$root

        return(c(maxAt, fct(maxAt, matrix(parm, 1, length(names)))))
    }
    # Identical to maxfct in cedergreen.R, but with a - in front of the derivative function derx


    returnList <- list(
        fct = fct, ssfct = ssfct, names = names, 
        deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
        edfct = edfct, maxfct = maxfct,
        name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text = ifelse(missing(fctText), "U-shaped Cedergreen-Ritz-Streibig", fctText),     
        noParm = sum(is.na(fixed)), fixed = fixed)
    
    class(returnList) <- "UCRS"
    invisible(returnList)
}


"UCRS.4a" <-
function(names = c("b", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 1, ...))
}

"UCRS.4b" <-
function(names = c("b", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.5, ...))
}

"UCRS.4c" <-
function(names = c("b", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.25, ...))
}

"UCRS.5a" <-
function(names = c("b", "c", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 1, ...))
}

"UCRS.5b" <-
function(names = c("b", "c", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.5, ...))
}

"UCRS.5c" <-
function(names = c("b", "c", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.25, ...))
}
