`sm2vec` <-
function (m, diag = FALSE) 
{
    return(as.vector(m[lower.tri(m, diag)]))
}

