\name{MU.poly.generator}
\alias{MU.poly.generator}
\title{Divide the landscape into management units (MUs) based on the MU membership of the sampling polygons}
\description{

A function to assist BSAs that work at the sampling polygon level. Once the BSA designates management unit
membership to each of the \code{sampling.polys}, \code{MU.poly.generator} uses this designation to split the entire study
area into management units.  
}
\usage{

MU.poly.generator(MU.from.BSA,sample.polys,landscape.poly)

}
\arguments{
\item{ MU.from.BSA}{a list of length equal to the number of MUs, the elements of which are vectors indicating which 
\code{sample.polys} are contained in each MU.}
\item{ sample.poly}{a list of length equal to the number of sampling polygons, each element being a sampling 
polygon of class \code{gpc.poly}}
\item{ landscape.poly}{the bounding box of the breeding population polygons}
}
\value{\code{munits}--A list with each element consisting of a management unit, represented by a polygon of class 
\code{gpc.poly}. }

\details{

A BSA that works at the sampling polygon level can be assisted by \code{MU.poly.generator}, which divides 
the entire space of the study area (=\code{landscape.poly}) into management units based on, 1) the nearest distance to 
each of the \code{sample.polys} and, 2) the management unit under which each of the \code{sample.polys} should be placed. 
The function \code{MU.poly.generator} should be included within (near the end) of the BSA.

\code{MU.poly.generator} generates a bounding box around the breeding population polygons, and within this box a 60X60 
grid of cells. Each cell is assigned the management unit of the \code{sample.poly} which it is closest to.  Grid cells 
that share a management unit are then joined together.  The resulting management units meet the requirements of a BSA, as they
are non-overlapping and cover the entire simulated study area.

}
\seealso{ \code{\link{tossm-package}},\code{\link{run.tossm}},\code{\link{hyptest.network.BSA}}}
\examples{
# In this example, run.tossm is used with the BSA set to hyptest.network.BSA, 
#which uses MU.poly.generator to return management unit polygons. 
#For details on the two functions, see their code.

#set up sampling and management schedule using def.make.schedule
schedule=def.make.schedule(n.pre.RMP=5,n.RMP=10,n.post.RMP=5,abund.gap=2)

#specify historic removals from each historic harvest area for 5 pre-modern
#management years
removals.h1<-rep(75,5)
removals.h2<-rep(0,5)
historic.removals<-cbind(removals.h1,removals.h2)

example<-run.tossm(rland = rland.end, bp.polys=bp.polys,
schedule = schedule,n.samples = 25, sample.polys=sample.polys, 
historic.removals=historic.removals,historic.polys = historic.polys, 
BSA = hyptest.network.BSA, BSA.args = list(sig.level=.05),
harvest.interval=10)
}
\keyword{internal}