"resPrint" <- function(resMat, headerText, interval, intervalLabel, display)
{
# Note: arguments "interval", "intervalLabel" no longer used
    if (display)
    {
        cat("\n")
        cat(paste(headerText, "\n", sep = ""))
        # if (!identical(interval, "none"))
        # {
        #     intervalText <- paste("(", intervalLabel, "-based confidence interval(s))\n", sep = "")
        #     cat(intervalText)
        # }
        cat("\n")
        printCoefmat(resMat, cs.ind = 1:ncol(resMat), tst.ind = NULL, has.Pvalue = FALSE)
    }
#    invisible(resMat)
}