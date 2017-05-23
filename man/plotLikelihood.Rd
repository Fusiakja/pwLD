\name{plotLikelihood}
\alias{plotLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of the log-likelihood function of pairwise haplotype frequencies}
\description{
The function visualises the log-likelihood function of haplotype frequencies of a given SNP pair.
}
\usage{
plotLikelihood(snps, code = c(0, 1, 2, 3), 
	paradigm = c("freq", "bayes"), dirich = rep(1,9), main = "", 
	ylim = NULL, legend = F, tol = .Machine$double.eps^0.6, 
	plot.tol=.Machine$double.eps^.25, digits=12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{\eqn{2} x \eqn{N} matrix  representing two loci with \emph{N} observed genotypes. The entries are items of parameter \code{code}.}
  \item{code}{vector of length four that indicates the code used in \code{snps} to mark allele combinations. The ordering is NOT arbitrary and is stated as follows: allele 1 homozygote, heterozygote, allele 2 homozygote, missing data. For instance, if there are two alleles, A and B, and missing values are denoted as NN, the parameter has to be \code{code=c("AA", "AB", "BB", "NN")}.}
  \item{paradigm}{character indicating the statistical paradigm to use in order to estimate genotypic frequencies. If \code{freq}, genotype frequencies are estimated by dividing absolute counts by \code{N}. If \code{bayes}, a Dirichlet prior is assumed and the mean of posterior distribution serves as estimate.}
  \item{dirich}{numeric vector of length 9 indicating the shape parameters of Dirichlet prior.}
  \item{main}{character, main title for the plot.}
  \item{ylim}{numeric vector of length two indicating the y limits of the plot.}
  \item{legend}{logical, if \code{TRUE} a legend will be plotted.}
  \item{plot.tol}{numeric, tolerance at the margins in order to avoid NaN values. (can be ignored by the user)}
  \item{tol}{numeric, accuracy value for numerical comparison of floating point numbers. See the help pages of \code{LD.cardano} for further details.}
  \item{digits}{numeric, specifies the number of digits for rounding of the cubic polynomial's coefficients. See the help pages of \code{LD.cardano} for further details.}
  }
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
 \references{Further informations and documentation can be found on \href{https://github.com/Fusiakja/pwLD}{github} }
\author{Karsten Krug <\email{karsten\_krug@gmx.net}>}
% \note{ ~~further notes~~ 
% }
\seealso{\code{LD.ml}, \code{LD.cardano}}
\examples{
## load sample HapMap data
data(HapMap)

## plot the log likelihood and add a legend
##plotLikelihood(HapMap_geno[c(45,46),], legend=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
