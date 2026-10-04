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
    
#    ## Setting the names of the parameters
#    names <- names[notFixed]


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }


    ## Defining derivatives
    
#    ## Constructing a helper function
#    xlogx <- function(x, p)
#    {
#        lv <- (x < 1e-12)
#        nlv <- !lv
#        
#        rv <- rep(0, length(x))
#        
#        xlv <- x[lv] 
#        rv[lv] <- log(xlv^(xlv^p[lv]))
#        
#        xnlv <- x[nlv]
#        rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
#    
#        rv
#    }
    
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

    ## Setting limits
#    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, lower = 1e-4, upper = 10000, ...)
    {    
        cedergreen(fixed =  fixed, names = names, alpha = alpha)$edfct(parm, 100 - respl, reference, type, lower, upper, ...) 
    }

#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair, upper = 10000, interval = c(1e-4, 10000))
#    {
#        cedergreen(alpha = alpha)$sifct(parm1, parm2, 100-pair, upper, interval)
#    }    

    ## Finding the maximal hormesis
    maxfct <- function(parm, upper, interval) {
       retVal <- cedergreen(fixed =  fixed, names = names, alpha = alpha)$maxfct(parm, upper, interval)
#       retVal[2] <- (parm[2] + parm[3]) - (retVal[2] - parm[2])
       retVal[2] <- (parm[2] + parm[3]) - retVal[2]
              
       return(retVal)
    }

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2,  # lowerc=lowerLimits, upperc=upperLimits, 
    edfct = edfct, maxfct = maxfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "U-shaped Cedergreen-Ritz-Streibig", fctText),     
    noParm = sum(is.na(fixed)))
    #list(fct = fct, ssfct = ssfct, names = names[notFixed], edfct = edfct, maxfct = maxfct,
    #name = "ucedergreen",
    #text = "U-shaped Cedergreen-Ritz-Streibig", 
    #noParm = sum(is.na(fixed)))
    
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
