\name{agg.gs.tseries}
\alias{agg.gs.tseries}
\title{Information from genetic sampling events for use in a boundary setting algorithm (BSA)}
\description{Genetic sample information is aggregated across years for use use by a boundary-
setting algorithm (BSA).  This function will likely be used only by those developing a BSA.}
\usage{agg.gs.tseries(nlocal=sys.parent())}
\value{
\item{agg.gs}{Genetic samples (alleles at all loci for all individuals sampled) aggregated across sampling polygons for each genetic sampling year (list of length == number of sampling polygons, each list element containing an array of dimension (gs.years*n.samples) X n.loci X 2.}
\item{agg.gs.gtypes}{Identical to \code{agg.gs} except the format for \code{agg.gs.gtypes} is slightly different, 
with each list element comprised of a matrix of dimension (gs.years*n.samples) X (n.haploid.loci+2*n.diploid.loci).}
\item{n.loci}{The number of loci used in the current simulation (numeric vector of length == 1).}
\item{n.areas}{The number of sampling polygons in the current simulation (numeric vector of length == 1).}
\item{n.alleles}{the number of alleles that have occurred over the entire simulation at each locus in each sampling polygon (vector of length equal to number of loci).}
\item{coords}{Individual ID, birth year and spatial coordinates for all individuals sampled in a \code{sample.poly}}
}
\seealso{\code{\link{hyptest.network.BSA}},\code{\link{run.tossm}}}
\keyword{internal}