formatF <- function(x, w, d) {
    
    # Format x as a Fixed Point number with fixed-width w and d places to the right of the decimal pt. If the number is too wide to fit, fill the
    # field with '*'s.  If d=0, do not print the decimal point.  Pad left with blanks. Pad right with zeros.  x can be a vector.  All elements of
    # the vector will be formatted the same. x cannot be a matrix or data frame.
    wholePart <- as.integer(x)
    wholeStrings <- as.character(wholePart)
    wholeStrings <- ifelse((wholeStrings == 0) & (x < 0), paste("-", wholeStrings, sep = ""), wholeStrings)
    wholeWidths <- ifelse(d > 0, w - d - 1, w)
    leftPad <- wholeWidths - nchar(wholeStrings)
    decimalPart <- round(x - wholePart, digits = d)
    decimalStrings <- as.character(decimalPart)
    decimalStrings <- substring(decimalStrings, regexpr("\\.", decimalStrings) + 1)
    rightPad <- ifelse(rep(d, length(x)) > 0, d - nchar(decimalStrings), 0)
    for (i in seq(along = wholeStrings)) {
        if (leftPad[i] >= 0) {
            wholeStrings[i] <- paste(paste(rep(" ", leftPad[i]), collapse = ""), wholeStrings[i], sep = "")
        } else {
            wholeStrings[i] <- paste(rep("*", wholeWidths[i]), collapse = "")
        }
        if (rightPad[i] > 0) {
            decimalStrings[i] <- paste(decimalStrings[i], paste(rep("0", rightPad[i]), collapse = ""), sep = "")
        }
    }
    decimalPoint <- ifelse(d > 0, ".", "")
    if (d > 0) 
        paste(wholeStrings, decimalPoint, decimalStrings, sep = "") else paste(wholeStrings)
}

printF <- function(x, w, d, matrix = F) {
    if (matrix) 
        x <- as.matrix(x)
    if (!is.numeric(x)) 
        stop("x must be numeric in printF")
    if (length(x) == 1) {
        cat(formatF(x, w = w, d = d))
    } else if (!matrix) {
        cat(formatF(as.vector(x), w = w, d = d))
    } else {
        apply(as.matrix(x), 1, FUN = function(y, w, d) {
            cat(formatF(y, w = w, d = d))
        }, w = w, d = d)
    }
    invisible()
}
