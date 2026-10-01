absToRel <- function(parmVec, respl, typeCalc) {
    
    #stop("Ooops")
    # Convert absolute to relative
    if (typeCalc == "absolute") {
        100 * (parmVec[3] - respl) / (parmVec[3] - parmVec[2])
    } else {
        respl
    }
}