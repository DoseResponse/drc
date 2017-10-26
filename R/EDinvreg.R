"EDinvreg" <- function(object, respLev, catLev = NA, intType = "confidence", level, type, extFactor = 10)
{
  if (!is.na(catLev))
  {  
    EDval <- ED(object, respLev, clevel = catLev, type = type, display = FALSE)
  } else {
    EDval <- ED(object, respLev, type = type, display = FALSE)
  }
  EDval1.1 <- EDval[1, 1]
  newData0 <- data.frame(EDval[, 1], catLev)
  
  objDL <- object[["dataList"]][["names"]]
  colnames(newData0) <- c(objDL[["dName"]], objDL[["cName"]])
  yval <- predict(object, newData0)
  #    print(yval)
  
  rootFct1 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
#    print(c(x, predict(object, newData, interval = intType, level = level)[2] - yval))
    predict(object, newData, interval = intType, level = level)[2] - yval
  }
  
  rootFct2 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
#    print(c(x, predict(object, newData, interval = intType, level = level)[3] - yval))
    predict(object, newData, interval = intType, level = level)[3] - yval
  }
  
  maxdose <- extFactor * max(object[["dataList"]][["dose"]])
  uroot1 <- try(uniroot(rootFct1, c(EDval1.1, maxdose)), silent = TRUE)
  if (inherits(uroot1, "try-error"))  # an error happens in case of a decreasing curve
  {
    #      print(c(0, EDval1.1)) 
    uroot2 <- try(uniroot(rootFct1, c(0, EDval1.1)), silent = TRUE)
    #if (inherits(uroot2, "try-error")) {lowlim <- 0} else {lowlim <- uroot2[["root"]]}
    uroot1 <- try(uniroot(rootFct2, c(EDval1.1, maxdose)), silent = TRUE)
    #if (inherits(uroot1, "try-error")) {uplim <- Inf} else {uplim <- uroot2[["root"]]}
  } else {
    uroot2 <- try(uniroot(rootFct2, c(0, EDval1.1)), silent = TRUE)
  }
  if (inherits(uroot1, "try-error")) {uplim <- Inf} else {uplim <- uroot1[["root"]]}
  if (inherits(uroot2, "try-error")) {lowlim <- 0} else {lowlim <- uroot2[["root"]]}
  
  
  #return(c(uroot2[["root"]], uroot1[["root"]]))
  return(c(lowlim, uplim))
}


EDinvreg1 <- Vectorize(EDinvreg, "respLev")
