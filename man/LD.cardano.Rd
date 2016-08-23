\name{LD.cardano}
\alias{LD.cardano}
\alias{LD.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of pairwise linkage disequilibrium measures}
\description{
Estimation of pairwise haplotype frequencies and linkage disequilibrium measures from genotyping SNP data.
\code{LD.cardano} determines the fix-point of the EM algorithm analytically in order to infer double heterozygote haplotype frequencies.
\code{LD.all} estimates LD measures of all pairwise SNP combinations contained in \code{data}.
}
\usage{
LD.cardano(snps, code = c(0, 1, 2, 3), 
	LD = c("Dprime", "Q", "r2"), CI = T, 
	alpha = .05, n.sim = 5000, returnLDdist=F, 
	paradigm = c("freq", "bayes"),
	dirich = rep(1, 9), all.solutions = F, 
	tol=.Machine$double.eps^.6, digits=12)

LD.all(data, code = c(0, 1, 2, 3),  
	LD = c("D", "Dprime", "Q", "r2", "OR", "MI", "chi2"), MAF=0, 
	paradigm = c("freq", "bayes"), dirich = rep(1,9), verbose=T,
	tol=.Machine$double.eps^.6, digits=12 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{\eqn{2} x \eqn{N} matrix  representing two loci with \code{N} observed genotypes. The entries are items of the parameter \code{code}.}
  \item{data}{\eqn{P} x \eqn{N} matrix representing \code{P} loci with \code{N} observed genotypes. The entries are items of the parameter \code{code}.}
  \item{code}{vector of length four that indicates the code used in \code{snps} to mark allele combinations. The ordering is NOT arbitrary and is stated as follows: allele 1 homozygote, heterozygote, allele 2 homozygote, missing data. For instance, if there are two alleles, A and B, and missing values are denoted as NN, the parameter has to be \code{code=c("AA", "AB", "BB", "NN")}.}
  \item{LD}{character vector indicating the LD measures to estimate. See \_Details\_ below.}
  \item{MAF}{numeric specifying the minor allele frequency. All loci with allele frequency \eqn{<=} MAF are excluded from \code{data}. If \code{MAF = 0} and \code{paradigm = "bayes"} no loci are excluded.} 
  \item{CI}{logical, if \code{TRUE} either a frequentistic confidence interval or a Bayesian credible interval is estimated, depending on parameter \code{paradigm}.}
  \item{alpha}{numerical, confidence level of confidence or credible interval, respectively.}
  \item{n.sim}{numerical indicating the number of bootstrap or Dirichlet replicates to use in order to determine confidence or credible intervals, respectively. Default value is set to 5000.}
  \item{paradigm}{character indicating the statistical paradigm to use in order to estimate genotypic frequencies. If \code{freq}, genotype frequencies are estimated by dividing absolute counts by \code{N}. If \code{bayes}, a Dirichlet prior is assumed and the mean of posterior distribution serves as an estimate.}
  \item{dirich}{numeric vector of length 9 indicating the shape parameters of the Dirichlet prior.}
  \item{returnLDdist}{logical, if \code{TRUE} the distribution of each LD measure over \code{n.sim} simulations is returned.}
  \item{verbose}{logical, if \code{TRUE} some informations about the progress of function \code{LD.all} will be displayed.}
  \item{all.solutions}{logical, if \code{TRUE} all solutions of the cubic polynomial will be returned.}
  \item{tol}{numeric, accuracy value for numerical comparison of floating point numbers. See \_Details\_ below. }
  \item{digits}{numeric, specifies the number of digits for rounding of the cubic polynomial's coefficients. See \_Details\_ below.}
}
\details{
% In order to infer haplotype frequencies caused by double heterozygotes an analytical solution of the EM algorithm is used. The conditional probability of the contribution of double heterozygotes results in a cubic polynome that is solved by the use of Cardano's Formula and the solution maximizing the likelihood of the data is chosen.
% 
Pairwise haplotype frequencies are represented by a \eqn{2 \times 2}{2 x 2} contingency table with following notation:
 \tabular{lccrr}{
 \tab alleles	\tab 0		\tab 1 		\tab Sum  \cr
 \tab 0		\tab \eqn{p_{00}}{P00}	\tab \eqn{p_{01}}{P01}	\tab \eqn{p_{0.}}{p0.} \cr
 \tab 1		\tab \eqn{p_{10}}{P10}	\tab \eqn{p_{11}}{P11}	\tab \eqn{p_{1.}}{p1.} \cr
 \tab Sum 		\tab \eqn{p_{.0}}{p.0}	\tab \eqn{p_{.1}}{p.1}	\tab 1 \cr
 }
The rows correspond to the first, the columns to the second SNP. The marginal distributions of this table are the allele frequencies that are directly estimated from genotypic frequencies. Haplotype frequencies are estimated via maximum-likelihood by determining the fix-point of the EM algorithm analytically, i.e. the resulting cubic polynomial is solved by the use of Cardano's formulae and the solution that is maximum-likelihood serves as an estimate.
The measures of linkage disequilibrium specified by parameter \code{LD} are:
 \tabular{ll}{
}
The parameters \code{tol} and \code{digits} are used to avoid numerical problems during the computation. Parameter \code{tol} specifies an epsilon for a numerically robust comparison of floating point numbers, parameter \code{digits} guarantees a robust case discrimination within the Cardano formulae. The chosen default values should be unaffected by the user. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Karsten Krug <\email{karsten\_krug@gmx.net}>}
\seealso{\code{LD.all}, \code{LD.cardano}}
\examples{
## load sample HapMap data
##data(HapMap)
## estimate all implemented LD measures between the fourth and fifth SNP

##res <- LD.cardano(HapMap_geno[c(4,5),  ])

## estimate Q based on Bayesian genotypic frequencies and return 
## the respective credible interval 

##res.ci <- LD.cardano(HapMap_geno[c(4,5), ], LD="Q", paradigm="bayes", CI=TRUE)

## estimate LD between all SNP pairs

##res.all <- LD.all(HapMap_geno)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}
