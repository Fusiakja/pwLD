\name{cardano.complex}
\alias{cardano.complex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analytic solution of a cubic polynomial}
\description{Analytic solution of a cubic polynomial by Cardano's formula.}
\usage{
cardano.complex(coeff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coeff}{numeric vector of length four indicating the coefficients of the polynomial in descending order. See \_Details\_ below.}
}
\details{The function calculates the roots of the polynomial within the set of complex numbers and hence returns complex roots. The coefficients are determined by parameter \code{coeff} and correspond to the standard form of cubic polynomial



	y = a \eqn{x^3} + b \eqn{x^2} + c \eqn{x} + d



in the following order: coeff = c(a,b,c,d).
}
\value{
  \item{D}{discriminant of Cardano's formula}
  \item{x}{complex vector containing the roots}
}
% \references{ ~put references to the literature/web site here ~ }
\author{Karsten Krug <\email{karsten\_krug@gmx.net}>}
% \note{ ~~further notes~~ 
% 
%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
% \seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
cardano.complex(c(4,6,3,1))
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
