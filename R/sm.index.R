`sm.index` <-
function (m, diag = FALSE) 
{
    m.dim <- length(diag(m))
    if (diag == TRUE) 
        num.entries <- m.dim * (m.dim + 1)/2
    else num.entries <- m.dim * (m.dim - 1)/2
    index1 <- rep(NA, num.entries)
    index2 <- rep(NA, num.entries)
    if (diag == TRUE) 
        delta <- 0
    else delta <- 1
    z <- 1
    for (i in 1:(m.dim - delta)) for (j in (i + delta):m.dim) {
        index1[z] <- i
        index2[z] <- j
        z <- z + 1
    }
    return(cbind(index1, index2))
}

