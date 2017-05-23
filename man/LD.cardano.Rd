\name{LD.cardano}
\alias{LD.cardano}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of pairwise linkage disequilibrium measures}
\description{
Estimation of pairwise haplotype frequencies and linkage disequilibrium measures from genotyping SNP data.
\code{LD.cardano} determines the fix-point of the EM algorithm analytically in order to infer double heterozygote haplotype frequencies.}
\usage{
LD.cardano(snps, code = c(0, 1, 2, 3), LD = c("Dprime", "Q", "r", "Y", "HS"), 
HSweight = 4, CI = T, strategy = c("jackknife", "bootstrap", "zapata"), 
alpha = 0.05, n.sim = 5000, returnLDdist = F, paradigm = c("freq", "bayes", "fullbayes"),
dirich = rep(1, 9), all.solutions = F, tol = .Machine$double.eps^0.6, 
digits = 12, seed = F, mc = 1000, intervall = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{
\eqn{2} x \eqn{N} matrix  representing two loci with \code{N} observed genotypes. The entries are items of the parameter \code{code}.}
  \item{code}{
vector of length four that indicates the code used in \code{snps} to mark allele combinations. The ordering is NOT arbitrary and is stated as follows: allele 1 homozygote, heterozygote, allele 2 homozygote, missing data. For instance, if there are two alleles, A and B, and missing values are denoted as NN, the parameter has to be \code{code=c("AA", "AB", "BB", "NN")}.}
  \item{LD}{
character vector indicating the LD measures to estimate. See \_Details\_ below.}
  \item{HSweight}{
numeric, parameter n of the HS measure for up-and downweighting}
  \item{CI}{
logical, if \code{TRUE} either a frequentistic confidence interval or a Bayesian credible interval is estimated, depending on parameter \code{paradigm}.}
  \item{strategy}{
haracter, if \code{bootstrap} a bootstrap confidence interval is constructed, if \code{jackknife} a jackknife confidence interval is constructed}
  \item{alpha}{
numerical, confidence level of confidence or credible interval, respectively.}
  \item{n.sim}{
numerical indicating the number of bootstrap or Dirichlet replicates to use in order to determine confidence or credible intervals, respectively. Default value is set to 5000.}
  \item{returnLDdist}{
logical, if \code{TRUE} the distribution of each LD measure over \code{n.sim} simulations is returned.}
  \item{paradigm}{
character indicating the statistical paradigm to use in order to estimate genotypic frequencies. If \code{freq}, genotype frequencies are estimated by dividing absolute counts by \code{N}. If \code{bayes}, a Dirichlet prior is assumed and the mean of posterior distribution serves as an estimate.}
  \item{dirich}{
numeric vector of length 9 indicating the shape parameters of the Dirichlet prior.}
  \item{all.solutions}{
logical, if \code{TRUE} all solutions of the cubic polynomial will be returned.}
  \item{tol}{
numeric, accuracy value for numerical comparison of floating point numbers. See \_Details\_ below}
  \item{digits}{
numeric, specifies the number of digits for rounding of the cubic polynomial's coefficients. See \_Details\_ below.}
  \item{seed}{
logical, if \code{TRUE} random seed of the CI construction was set}
  \item{mc}{
numberic, numbers of iterations of the monte carlo simulation for the full bayesian haplotype frequency estimator}
  \item{intervall}{
integer, if \code{strategy=="bootstrap"&&interval=="0"} bootstrap quantile method is used, if \code{strategy=="bootstrap"&&interval=="1"} bootstrap standard error method is used, if \code{strategy=="jackknife"&&interval=="0"} jackknife leave one out method is used, if \code{strategy=="jackknife"&&interval=="1"} jackknife pseudo value mehtod is used}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
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

The output of \code{LD.cardano} is a list of haplotype frequencies, chi-.
If the paramter:
If \code{CI} ist \code{FALSE} a list containing haplotype frequencies, chi square value, number of genotype counts, number of double heterozygotes, names of the SNPs LD is measured between and the defined LD measrues
If \code{CI} a list containing haplotype frequencies, chi square value, number of genotype counts, number of double heterozygotes, names of the SNPs LD is measured between and the defined LD measrues, CIs for the LD measures with the wanted CI method
If \code{all.solutions} is \code{TRUE} the solutions of maximizing the likelihood of the data is printed out too.
If \code{LDdist} a list of tables of used by the Bootstrap/Jackknife confidence interval estimation and the LD measrues for these tables are printed out too.
}
\references{
Further informations and documentation can be found on \href{https://github.com/Fusiakja/pwLD}{github} }
\author{
Karsten Krug <\email{karsten\_krug@gmx.net}>, Jakub Fusiak <\email{fusiakjakub@gmx.net}>}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LD.all}}
}
\examples{
## load sample HapMap data
data(HapMap)

##Calculation of pairwise LD measures with default values, no conifidence interval
##res.ci <- LD.cardano(snps=HapMap_geno[c(4,5), ], LD=c("Dprime","Q","r","Y","HS"), 
##paradigm="bayes", CI=FALSE,strategy = "bootstrap", 
##intervall = 1,alpha = 0.05,HSweight = 4,returnLDdist = FALSE,
##mc = 1000,seed = FALSE,all.solutions = FALSE ,n.sim = 1000)

##Calculation of pairwise LD measures D', Q, r, Y and HS with the bayesian plug in estimator
##res.ci <- LD.cardano(snps=HapMap_geno[c(4,5), ], LD=c("Dprime","Q","r","Y","HS"), 
##paradigm="bayes", CI=TRUE,strategy = "bootstrap", 
##intervall = 1,alpha = 0.05,HSweight = 4,returnLDdist = FALSE,
##mc = 1000,seed = FALSE,all.solutions = FALSE ,n.sim = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
