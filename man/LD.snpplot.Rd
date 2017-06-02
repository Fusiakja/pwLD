\name{LD.snpplot}
\alias{LD.snpplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of pairwise linkage disequilibrium measures
}
\description{
Graphical output of a scatterplot consists of the LD measures between all genotyped SNPs. This data is the result of the \code{LD.all} command.}
\usage{
LD.snpplot(valuematrix, annotation, values = c("values", "upperCI", "lowerCI"),
LD = c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), 
snp.zoom = FALSE, snp.target, snp.region, unit = c("Mb", "Kb", "bp"), legends=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{valuematrix}{
object, output of the \code{LD.all} command - 
symmatric matrix of the LD measure between all genotyped SNPs.}
  \item{annotation}{
object, numeric vector containing the genomic positions of SNPs.}
  \item{values}{
character, Parameter to define wich values shoul be plotted.
If \code{values} is \code{"values"} the LD values are taken for the plotting.
If \code{values} is \code{"upperCI"} the upper bounds of the confidence interval estimation are taken for the plotting.
If \code{values} is \code{"lowerCI"} the lower bounds of the confidence interval estimation are taken for the plotting.}
  \item{LD}{
character, Defines which LD measure will be plotted. One of the following values is possible: \code{Dprime}, \code{r}, \code{Q}, \code{Y}, \code{HS} }
  \item{snp.zoom}{
logical, if \code{FALSE} all data will be plotted
if \code{TRUE} only part of the data will be plotted}
  \item{snp.target}{
character, definition of the SNP, that is wanted to be investigated, e.g. "rs16984366"}
  \item{snp.region}{
integer, definition of the region of the SNP, that is wanted to be investigated, e.g. 5 SNPs before and after the "rs16984366" SNP}
  \item{unit}{
character, indicating the scale used for the plot. One of the following values is possible: \code{Mb}, \code{Kb}, \code{bp}.}
  \item{legends}{
  logical, show legend? \code{TRUE} Yes. \code{FALSE} No.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Output is a scatterplot of LD values of a region of SNPs against the position of the Chromosome in base pairs
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
\code{\link{LD.corrplot}}
}
\examples{
## load sample HapMap data
data(HapMap)

## Creating correlation matrix
res.all <- LD.all(data = HapMap_geno, code = c(0, 1, 2, 3), 
LD = c("Dprime", "Q", "r", "OR", "Y", "HS"), MAF = 0.1, paradigm = c("freq"), 
strategy = c("jackknife"), verbose = TRUE, tol = .Machine$double.eps^0.6, 
digits = 12, CI = TRUE, HSweight = 4, alpha = 0.05, intervall = 0)

## Plotting scatterplot for D' without zoom
LD.snpplot(valuematrix =  res.all, annotation = HapMap_anno, 
values = "values",LD ="Dprime",snp.zoom = FALSE, unit = "Mb", 
snp.target = "rs16984366", snp.region = 5)

## Plotting scatterplot for D' with zoom on rs16984366 around a region of 5 SNPs
LD.snpplot(valuematrix =  res.all, annotation = HapMap_anno, 
values = "values",LD ="Dprime",snp.zoom = TRUE, snp.target = "rs16984366", 
snp.region = 5,unit = "Mb")
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line