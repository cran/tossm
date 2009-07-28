\name{def.genetic.sampler}
\alias{def.genetic.sampler}
\title{Default method for drawing samples from sampling polygons}
\description{\code{def.genetic.sampler} draws a random sample of individuals from each sampling polygon and
assigns coordinates for each sample consistent with a random spatial distribution.}
\usage{def.genetic.sampler(rland,sample.polys,bp.polys,n.samples,n.in.ss,ID.tracker=NULL)}
\arguments{
\item{rland}{an \pkg{rmetasim} landscape object}
\item{sample.polys}{polygon(s) defining where genetic samples are to be collected}
\item{bp.polys}{polygons defining the number and geographic extent of breeding populations}
\item{n.samples}{number of genetic samples taken from each \code{sample.poly} in each genetic sampling year of the simulation}
\item{n.in.ss}{number of individuals present in a sampling polygon}
\item{ID.tracker}{internal placekeeping argument}
}
\details{
\code{def.genetic.sampler} generates random coordinates within each \code{sample.polys}, and draws genetic samples
from each of these points. \cr
Other sampling strategies can be developed through revision of the \code{def.genetic.sampler} code.
}
\seealso{\code{\link{tossm-package}} \code{\link{run.tossm}}}
\keyword{internal}