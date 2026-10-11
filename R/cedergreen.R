"cedergreen" <- function(
    fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
    method = c("1", "2", "3", "4"), ssfct = NULL, 
    alpha, fctName, fctText) {

    ## Checking arguments
    numParm <- 5
    if (!is.character(names) || !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}

    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm) {

        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2] + (parmMat[,3] - parmMat[,2] + parmMat[,5]*exp(-1/(dose^alpha)))/(1 + exp(parmMat[,1]*(log(dose) - log(parmMat[,4]))))
    }


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
#
#
#    ## Defining flag to indicate if more general ANOVA model
##    anovaYes <- TRUE
#    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}
#    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
#        ssfct <- cedergreen.ssf(method, fixed, alpha)       
        ssfct <- function(dframe)
        {
            initval <- llogistic(method = method)$ssfct(dframe)   
            initval[5] <- (2*(median(dframe[, 2])-initval[2])-(initval[3]-initval[2]))*exp(1/(initval[4]^alpha))
    
            return(initval[notFixed])
        }        
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

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1 - 1/t3, 
               1/t3, 
               t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               t0/t3 )[, notFixed]
    }

    deriv1 <- function(dose, parm) {
    
        # Obtaining the derivatives
        deriv1FctTemp <- function(dose, pmRow)
        {
            deriv1Fct <- deriv(~c+(d-c+f*exp(-1/x^alpha))/(1+(x/e)^b), c("b", "c", "d", "e", "f"), 
                                function(x, b, c, d, e, f){})
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
        for (i in 1:nrpar)
        {
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

        parmMat[, 5] * t3 / t2 - (parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*exp(-1/(dose^alpha))) * t1 * parmMat[, 1] / (t2^2 * dose) + parmMat[, 5] * t3 / t2
    }


    ## Defining the ED function    
    edfct <- function(parm, respl, reference, type, lower = NULL, upper = NULL, ...) {

        # lower = 1e-4, upper = 10000
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type, TRUE)  # FALSE)  Changed 2010-06-02 after e-mail from Claire       
        tempVal <- 1 - p / 100
        
        ## Defining the numerator of derivative of the BC model function wrt. dose
        derx <- function(dose) {
            expVal1 <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            expVal2 <- exp(-1/(dose^alpha)) * alpha * (dose^(-alpha - 1))

            parmVec[5]*(1+expVal1)*expVal2 - (parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))*expVal1*parmVec[1]/dose
        }
        nearzepVal <- parmVec[4] * (parmVec[5] / (parmVec[3]-parmVec[2]))^(1 / parmVec[1])

        if (is.null(lower)) {
            xvalues <- 10^(runif(100, -10, log10(nearzepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec<0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        }
      
        # Checking if the curve is effectively decreasing
        if (length(lower) == 0) {maxAt <- 1e-20} else {
            maxAt <- uniroot(derx, lower = lower, upper = nearzepVal)$root 
        } 
        if (is.null(upper)) {
            upper <- parmVec[4] * 10000
        }

        eqn <- function(dose) {
            tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*exp(-1/(dose^alpha))/(parmVec[3]-parmVec[2]))}
        EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root

        tempVal1 <- exp(parmVec[1] * (log(EDp) - log(parmVec[4])))
        tempVal2 <- parmVec[3] - parmVec[2]

        derNum <- c(-tempVal * tempVal1 * (log(EDp) - log(parmVec[4])), 
                     ((1-tempVal) * tempVal2 + parmVec[5] * exp(-1/(EDp^alpha))) / (tempVal * tempVal2^2),
                     1 / tempVal2, 
                     tempVal * tempVal1 * parmVec[1] / parmVec[4],
                     exp(-1/(EDp^alpha)) / tempVal2)
        derDenom <- tempVal * tempVal1 *parmVec[1] / EDp - (exp(-1/(EDp^alpha)) * alpha * (EDp^(-alpha - 1))) * parmVec[5] / tempVal2
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

            parmVec[5]*(1+expVal1)*expVal2 - (parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))*expVal1*parmVec[1]/dose
        }
        nearzepVal <- parmVec[4] * (parmVec[5] / (parmVec[3]-parmVec[2]))^(1 / parmVec[1])

        if (is.null(lower)) {
            xvalues <- 10^(runif(noRan, -10, log10(nearzepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec<0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        }  
        maxAt <- uniroot(derx, lower = lower, upper = nearzepVal)$root

        return(c(maxAt, fct(maxAt, matrix(parm, 1, length(names)))))
    }

    
    returnList <- list(
        fct = fct, ssfct = ssfct, names = names, 
        deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
        edfct = edfct, maxfct = maxfct, 
        name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text = ifelse(missing(fctText), "Cedergreen-Ritz-Streibig", fctText),     
        noParm = sum(is.na(fixed)), fixed = fixed)

    class(returnList) <- "mllogistic"
    invisible(returnList)
}

"CRS.4a" <-
function(names = c("b", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 1, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1)",
    ...))
}

"CRS.4b" <-
function(names = c("b", "d", "e", "f"), ...) {
    
    ## Checking arguments
    if (!is.character(names) || !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 0.5, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=.5)",
    ...))
}

"CRS.4c" <-
function(names = c("b", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 0.25, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=.25)",
    ...))
}

"CRS.5a" <-
function(names = c("b", "c", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names)==5)) {stop("Not correct 'names' argument")}
 
    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 1, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=1)",
    ...))
}

"CRS.5b" <-
function(names = c("b", "c", "d", "e", "f"), ...) {
    
    ## Checking arguments
    if (!is.character(names) || !(length(names)==5)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 0.5, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=.5)",
    ...))
}

"CRS.5c" <-
function(names = c("b", "c", "d", "e", "f"), ...) {

    ## Checking arguments
    if (!is.character(names) || !(length(names)==5)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 0.25, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=.25)",
    ...))
}


