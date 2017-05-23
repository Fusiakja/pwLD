\name{LD.corrplot}
\alias{LD.corrplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of pairwise linkage disequilibrium measures
}
\description{
Graphical output of a correlation matrix consists of the LD measures between all genotyped SNPs. This data is the result of the \code{LD.all} command.}
\usage{
LD.corrplot(valuematrix, values = c("values", "upperCI", "lowerCI"), 
LD = c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), 
snp.zoom = FALSE, snp.target, snp.region)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{valuematrix}{output of the \code{LD.all} command - 
symmatric matrix of the LD measure between all genotyped SNPs.}
  \item{values}{
Parameter to define wich values shoul be plotted.
If \code{vales} is \code{"values"} the LD values are taken for the plotting.
If \code{values} is \code{"upperCI"} the upper bounds of the confidence interval estimation are taken for the plotting.
If \code{values} is \code{"lowerCI"} the lower bounds of the confidence interval estimation are taken for the plotting.}
  \item{LD}{
Defines which LD measure will be plotted}
  \item{snp.zoom}{
logical, if \code{FALSE} all data will be plotted
if \code{TRUE} only part of the data will be plotted}
  \item{snp.target}{
Definition of the SNP, that is wanted to be investigated, e.g. "rs16984366"}
  \item{snp.region}{
Definition of the region of the SNP, that is wanted to be investigated, e.g. 5 SNPs before and after the "rs16984366" SNP}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Output is a correlation matrix plot
}
\references{
Further informations and documentation can be found on \href{https://github.com/Fusiakja/pwLD}{github} }
\author{
Karsten Krug <\email{karsten\_krug@gmx.net}>, Jakub Fusiak <\email{fusiakjakub@gmx.net}>
}
\note{
This command uses the package \code{corrplot} to give out the graphical output.
}
\seealso{
\code{\link{LD.snpplot}}
}
\examples{
## load sample HapMap data
data(HapMap)

##Creating correlation matrix
res.all <- LD.all(data = HapMap_geno, code = c(0, 1, 2, 3), 
LD = c("Dprime", "Q", "r", "OR", "Y", "HS"), MAF = 0.1, 
paradigm = c("freq"), strategy = c("jackknife"), verbose = TRUE, 
tol = .Machine$double.eps^0.6, digits = 12, CI = TRUE, HSweight = 4, 
alpha = 0.05, inervall = 0)

##Plotting correlation matrix of D' without zoom
LD.corrplot(res.all, values = c("values"), LD = c("Dprime"), snp.zoom = FALSE)

##Plotting correlation matrix of D' with zoom on rs16984366 around a region of 5 SNPs
LD.corrplot(res.all, values = c("values"), LD = c("Dprime"), 
snp.zoom = TRUE, snp.target="rs16984366", snp.region=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
