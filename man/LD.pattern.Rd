\name{LD.pattern}
\alias{LD.pattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualisation of LD patterns along chromosomes or chromosomal regions.}
\description{The function visualises patterns of LD along chromosomes or chromosomal regions. 
}
\usage{
LD.pattern(LD.mat, map = NULL, window = 1.7e6, 
	minSNPs = 10, scale = c("Mb", "Kb", "bp"), plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LD.mat}{object returned by function \code{LD.all}.}
  \item{map}{numeric vector containing the genomic positions of SNPs.}
  \item{window}{numeric indicating the size of the sliding window in bp. See \_Details\_ below. }
  \item{minSNPs}{minimum number of SNPs that have to be within a window.}
  \item{scale}{character indicating the scale used for the plot. One of the following values is possible: \code{Mb}, \code{Kb}, \code{bp}.}
  \item{plot}{logical, if \code{TRUE} the plot will be produced.}
}
\details{
The absolute LD values of SNPs within a sliding window of size \code{window} are averaged and plotted against the middle position of the corresponding interval. Parameter \code{minSNPs} indicates the minimal number of SNPs that have to be within the window in order to average the LD values.
}
\value{
  \item{LD.mean}{vector containing the averaged LD values of each sliding window, i.e. the y-axis of the plot.}
  \item{BP.mean}{vector containing the middle position of each sliding window, i.e. the x-axis of the plot.}
  \item{#SNPs}{vector containing the number of SNPs within each sliding window.}
  \item{info}{some general informations.}
}
% \references{ ~put references to the literature/web site here ~ }
\author{Karsten Krug <\email{karsten\_krug@gmx.net}>}
% \note{ ~~further notes~~ 
% 
%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{LD.all}
\examples{
## load sample HapMap data
##data(HapMap)

## estimate pairwise LD between all SNPs with MAF > 0
##res.LD <- LD.all(HapMap_geno, LD=c("Dprime","Q"), MAF=0) 

## make the plot
##res.pattern <- LD.pattern(res.LD, window=3e5, 
##map=HapMap_anno[rownames(res.LD[[1]]), "pos"], scale="Kb" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
