\name{LD.all}
\alias{LD.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of pairwise linkage disequilibrium measures}
\description{
Estimation of pairwise haplotype frequencies and linkage disequilibrium measures from genotyping SNP data.
\code{LD.all}estimates LD measures of all pairwise SNP combinations contained in \code{data}.}
\usage{
LD.all(data, code = c(0, 1, 2, 3), 
LD = c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), MAF = 0, 
paradigm = c("freq", "bayes", "fullbayes"), 
strategy = c("bootstrap", "jackknife"), dirich = rep(1, 4), 
verbose = T, tol = .Machine$double.eps^0.6, digits = 12, 
CI = F, HSweight = 4, alpha = 0.2, nSim = 1000, seed = F, intervall = c(0, 1), mc=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
\eqn{P} x \eqn{N} matrix representing \code{P} loci with \code{N} observed genotypes. The entries are items of the parameter \code{code}.}
  \item{code}{
vector of length four that indicates the code used in \code{snps} to mark allele combinations. The ordering is NOT arbitrary and is stated as follows: allele 1 homozygote, heterozygote, allele 2 homozygote, missing data. For instance, if there are two alleles, A and B, and missing values are denoted as NN, the parameter has to be \code{code=c("AA", "AB", "BB", "NN")}.}
  \item{LD}{
character vector indicating the LD measures to estimate. See \_Details\_ below.}
  \item{MAF}{
numeric specifying the minor allele frequency. All loci with allele frequency \eqn{<=} MAF are excluded from \code{data}. If \code{MAF = 0} and \code{paradigm = "bayes"} no loci are excluded.}
  \item{paradigm}{
character indicating the statistical paradigm to use in order to estimate haplotype frequencies. If \code{freq}, haplotype frequencies are estimated by dividing absolute counts by \code{N}. If \code{bayes}, a Dirichlet prior is assumed and the mean of posterior distribution serves as an estimate.If \code{fullbayes}, a Dirichlet prior is assumed and the mean of posterior distribution serves as an estimate using a monte carlo simulation.}
  \item{strategy}{
character, if \code{bootstrap} a bootstrap confidence interval is constructed, if \code{jackknife} a jackknife confidence interval is constructed}
  \item{dirich}{
numeric vector of length 9 indicating the shape parameters of the Dirichlet prior.}
  \item{verbose}{
logical, if \code{TRUE} some informations about the progress of function \code{LD.all} will be displayed.}
  \item{tol}{
numeric, accuracy value for numerical comparison of floating point numbers. See \_Details\_ below. }
  \item{digits}{
numeric, specifies the number of digits for rounding of the cubic polynomial's coefficients. See \_Details\_ below.}
  \item{CI}{
logical, if \code{TRUE} either a frequentistic confidence interval or a Bayesian credible interval is estimated, depending on parameter \code{paradigm}.}
  \item{HSweight}{
numeric, parameter n of the HS measure for up-and downweighting}
  \item{alpha}{
numerical, confidence level of confidence or credible interval, respectively.}
  \item{nSim}{
numerical, indicating the number of bootstrap or Dirichlet replicates to use in order to determine confidence or credible intervals, respectively. Default value is set to 1000.}
  \item{seed}{
logical, if \code{TRUE} random seed of the CI construction was set}
  \item{intervall}{
integer, if \code{strategy=="bootstrap"&&intervall=="0"} bootstrap quantile method is used, if \code{strategy=="bootstrap"&&intervall=="1"} bootstrap standard error method is used, if \code{strategy=="jackknife"&&intervall=="1"} jackknife leave one out method is used, if \code{strategy=="jackknife"&&intervall=="0"} jackknife pseudo value mehtod is used.}
  \item{mc}{
  integer, number of Monte Carlo iterations for full bayesian estimator.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The output of \code{LD.all} is a list of correlation matrizes.
If the paramter:
\code{CI}is \code{FALSE}, a list of one entry containing the correlation matrix of all defined LD measures between all genotyped SNPs.
\code{CI}is \code{TRUE}, a list of three entries containing the correlation matrizies of all defined LD measures between all genotyped SNPs, the lower bounds of confidence interval for all defined LD measures between all genotyped SNPs, the upper bounds of confidence interval for all defined LD measures between all genotyped SNPs.
}
\references{
Further informations and documentation can be found on \href{https://github.com/Fusiakja/pwLD}{github} 
}
\author{
Karsten Krug <\email{karsten\_krug@gmx.net}>, Jakub Fusiak <\email{fusiakjakub@gmx.net}>}
%%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LD.cardano}}
}
\examples{
## load sample HapMap data
data(HapMap)

## estimate LD measures Dprime, Q, r, OR, Y and HS between all SNP pairs:
#res.all <- LD.all(data = HapMap_geno, code = c(0, 1, 2, 3), 
#LD = c("Dprime", "Q", "r", "OR", "Y", "HS"), MAF = 0.1, 
#paradigm = c("freq"), strategy = c("jackknife"), verbose = TRUE, 
#tol = .Machine$double.eps^0.6, digits = 12, CI = TRUE, HSweight = 4, 
#alpha = 0.05, intervall = 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SNPs }% use one of  RShowDoc("KEYWORDS")
\keyword{ all }% __ONLY ONE__ keyword per line
