perc <- function(arr, pc.wert) {
  neu <- arr
  neu <- sort(neu)
  n <- length(neu)
  exakt.perz <- (n+1)*pc.wert/100
  hilfswert <- floor(exakt.perz);
  untergr <- neu[hilfswert];
  obergr <- neu[hilfswert + 1];
  if ((hilfswert < 1) || ((hilfswert+1) > n))
  {
    print("Percentile is out of range");
    return(invisible())
  }
  perzentil <- untergr + (obergr - untergr)*(exakt.perz - floor(exakt.perz))
  return(perzentil)
}