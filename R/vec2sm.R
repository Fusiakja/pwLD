`vec2sm` <-
function (vec, diag = FALSE, order = NULL) 
{
    n <- (sqrt(1 + 8 * length(vec)) + 1)/2
    if (diag == TRUE) 
        n <- n - 1
    if (ceiling(n) != floor(n)) 
        stop("Length of vector incompatible with symmetric matrix")
    m <- matrix(NA, nrow = n, ncol = n)
    lo <- lower.tri(m, diag)
    if (is.null(order)) {
        m[lo] <- vec
    }
    else {
        vec.in.order <- rep(NA, length(order))
        vec.in.order[order] <- vec
        m[lo] <- vec.in.order
    }
    for (i in 1:(n - 1)) for (j in (i + 1):n) m[i, j] <- m[j, 
        i]
    return(m)
}

