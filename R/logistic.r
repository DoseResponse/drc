"logistic" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText) {

    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
           
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2]+(parmMat[,3]-parmMat[,2])/((1+exp(parmMat[,1]*(dose-parmMat[,4])))^parmMat[,5])
    }

    ## Defining self starter function
    if (!is.null(ssfct)) {
        ssfct <- ssfct
    } else {
        ssfct <- logistic.ssf(method, fixed)     
    }

    ## Defining names
    names <- names[notFixed]

    ## Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(dose - parmMat[, 4]))
        t3 <- (1 + t2)^(2*parmMat[, 5])
        t4 <- parmMat[, 5]*((1 + t2)^(-parmMat[, 5] - 1))
        t5 <- (1 + t2)^(parmMat[, 5])                  

        cbind( -t1*t2*t4*(dose - parmMat[ , 4]), 
               1 - 1/t5, 
               1/t5, 
               t1*t2*t4*parmMat[, 1], 
               -t1*log(1+t2)/t5 )[, notFixed]
    }

    deriv2 <- NULL
    
    ## Defining the first derivatives (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        temp1 <- exp(parmMat[, 1]*(x - parmMat[, 4]))       
        
        (-parmMat[, 5]*(parmMat[, 3] - parmMat[, 2])*temp1*parmMat[, 1])/((1 + temp1)^(parmMat[, 5] + 1))
    }

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...) {

        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- 1 - p/100
            
        # deriv(~e + log(((d-c)/(y-c))^(1/f) - 1) / b, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
        EDderFct <- function (b, c, d, e, f) {

            .expr1 <- d - c
            .expr2 <- tempVal * .expr1
            .expr3 <- 1 / tempVal  # .expr1/.expr2
            .expr4 <- 1/f
            .expr5 <- .expr3^.expr4
            .expr6 <- .expr5 - 1
            .expr7 <- log(.expr6)
            .expr14 <- .expr3^(.expr4 - 1)
            .expr15 <- 1/.expr2
            .value <- e + .expr7/b
            .grad <- array(0, c(length(.value), 5L), 
            list(NULL, c("b", "c", "d", "e", "f")))
            .grad[, "b"] <- -(.expr7/b^2)
            .grad[, "c"] <- -(.expr14 * (.expr4 * (.expr15 - .expr1/.expr2^2))/.expr6/b)
            .grad[, "d"] <- .expr14 * (.expr4 * .expr15)/.expr6/b
            .grad[, "e"] <- 1
            .grad[, "f"] <- -(.expr5 * (log(.expr3) * (1/f^2))/.expr6/b)
            attr(.value, "gradient") <- .grad
            .value
        } 

        EDcalc <- EDderFct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDp <- as.numeric(EDcalc)
        EDder <- attr(EDcalc, "gradient")
        if (!identical(type, "absolute")) {EDder[c(2, 3)] <- 0}

        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + parmVec[4]
    }

    ## Defining return list
    returnList <- list(
        fct = fct, ssfct = ssfct, names = names, 
        deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, 
        edfct = edfct, inversion = invfct,
        name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text = ifelse(missing(fctText), "Logistic (ED50 as parameter)", fctText),    
        noParm = sum(is.na(fixed)), fixed = fixed)

    class(returnList) <- "logistic"  # "Boltzmann"
    invisible(returnList)
}

#"boltzmann" <- logistic

"L.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Logistic (ED50 as parameter) with lower limit fixed at 0", ...))
}


"L.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed, 1), names = c(names, "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Logistic (ED50 as parameter)", ...))
}


"L.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return(logistic(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), 
    fctText = "Generalised logistic (ED50 as parameter)", ...))
}

