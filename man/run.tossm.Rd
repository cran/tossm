\name{run.tossm}
\alias{run.tossm}
\title{Simulate harvest of a spatially structured population}
\description{
A function for simulating the population dynamics and genetics of a harvested, spatially structured group of populations. The initial population input into the simulation should be a valid \code{\link[<rmetasim>:<rmetasim>-package]{<rmetasim>}} landscape object.  The other required arguments are introduced in help for the \pkg{tossm} package (\code{\link{tossm-package}}).  If the user is new to the \pkg{tossm} package, the introduction provided by \code{tossm-package} may be valuable.
}
\usage{
run.tossm(rland=NULL, bp.polys, schedule=NULL, n.samples=NULL, sample.polys, 
	initial.depletion=.30, historic.removals=NULL, historic.polys=NULL, 
	BSA=fixed.MU.BSA, BSA.args=list(n.mus=1), 
	stop.years=schedule$stop.years, gs.years=schedule$gs.years, 
	abund.est.years=schedule$abund.est.years, 
	pre.RMP.years=schedule$pre.RMP.years, CLA.years=schedule$CLA.years, 
	post.RMP.years=schedule$post.RMP.years, BSA.years=schedule$BSA.years, 
	harvest.interval, genetic.sampler=def.genetic.sampler,
	abund.for.10pc.CV=70000, quota.calc=PBR, 
	quota.args=list(r.max=0.04,F.r=1,multiplier=1), CLA.prog=NULL, 
	CLA.dir=NULL, plot.polys=FALSE, seed=-1) 
}
\arguments{
\item{rland}{A landscape from \code{rmetasim}}
\item{bp.polys}{A list of polygons (class \code{gpc.poly}) representing the spatial extent of breeding populations.}
\item{n.samples}{A vector indicating the number of genetic samples to collect within each sample polygon in each of \code{gs.years}.   If a single number is given, then that number of samples is taken from every sample polygon.}
\item{sample.polys}{A list of polygons (class \code{gpc.poly}) representing the geographic extent of sampling areas.
These polygons must not overlap, and must be fully contained within one or more breeding population polygons represented by \code{bp.polys}.}
\item{initial.depletion}{A vector indicating the level to which populations are depleted at the beginning of the modern management phase.  The length of the vector should equal the number of populations.  If not, the end of the vector will be padded with 0.99.  Thus, the default is that population 1 is depleted to 0.30 of carrying capacity, while all other populations are at 0.99 of K.  The initial depletion argument is only used if \code{historic.removals} is null, otherwise the \code{initial.depletion} argument will be ignored.}
\item{historic.removals}{This argument is used to specify in detail any historic harvest before the period in which modern management takes place.  It consists of a matrix with each row representing a single \code{pre.RMP.year} and each column representing an area of historic harvest.  The entries are the number of animals killed in each historic harvest area in each \code{pre.RMP.year}.  Note that values of 0 in this matrix will be automatically changed to 1 by \code{run.tossm} due to idiosyncracies of the CLA program, MANAGE-D.exe.  If \code{historic.removals} is set, the argument \code{initial.depletion} will be ignored.}
\item{historic.polys}{A list of polygons (class \code{gpc.poly}) representing the geographic extent of historic harvest areas.  As with the sampling polygons, these historic harvest polygons must not overlap each other, and must be fully contained within the \code{bp.polys}.}
\item{BSA}{A function used to set management unit boundaries.  Leaving the \code{BSA} argument of \code{run.tossm} as the default (\code{BSA=fixed.MU.BSA} and \code{BSA.args= list(n.mus=1)}) results in all breeding population polygons being combined into a single management unit. \cr See \code{\link{hyptest.network.BSA}} for further details.}
\item{BSA.args}{Extra arguments passed to the boundary-setting algorithm (see \code{\link{hyptest.network.BSA}}).}
\item{schedule}{A main prerequisite is setting up a schedule of events for the simulation. This schedule defines a simulation timeline.  The simulation schedule includes some terminology specific to International Whaling Commission (IWC) issues. The schedule must include the following components: \cr 
\code{stop.years}: The total number of years in the simulation \cr
\code{gs.years}: The years in which genetic samples are taken \cr
\code{abund.est.years}: The years in which abundance estimates are obtained \cr
\code{pre.RMP.years}: A period that predates modern management (may include historic harvest levels) \cr
\code{CLA.years}: Years in which quotas are calculated via the catch-limit algorithm \cr
\code{post.RMP.years}: A recovery period after modern management during which no harvest occurs(optional) \cr
\code{BSA.years}: Years in which the boundary-setting algorithm is called to detect spatial structure and define management units.

These different events are contingent on each other, in the following way:


\bold{1)} Genetic samples (\code{gs.years}) must be obtained before boundaries are set using the BSA (\code{BSA.years}) (collecting them in the same simulation year is fine). 

\bold{2)} Boundaries (\code{BSA.years}) must be set and abundance estimates (\code{abund.est.years}) must be obtained before an initial catch limit is set (\code{CLA.years}).  In other words, the first element in \code{BSA.years} and \code{abund.est.years} must be less than or equal to the first element in \code{CLA.years}.  Subsequent calls to the CLA, however, can be made without additional calls to the BSA or abundance estimates. 

\bold{3)} All \code{CLA.years} take place during the RMP phase. 

Though it is possible to set up a simulation schedule manually by creating vectors for each of the schedule components 
and passing them as arguments to \code{run.tossm}, the function \code{\link{def.make.schedule}}
 is a convenient way of setting up a schedule of all simulation events.  It is recommended to use 
\code{def.make.schedule} initially to run simulations; the user can experiment with customizing this schedule once more 
comfortable with the package. \cr
}
 
\item{stop.years, gs.years, abund.est.years, pre.RMP.years,\cr CLA.years, post.RMP.years, BSA.years}{integer vectors dictating the timing of events in \code{schedule}.}

\item{harvest.interval}{A vector of length one, which specifies the width of the harvest intervals into which the study area 
is divided.  This argument controls the degree to which harvest is spatially concentrated within MU.  Specifying a harvest 
interval width that is very small relative to the x-range of the breeding populations will result maximize the degree to 
which harvest is spatially concentrated.  Setting the harvest interval width equal to the x-range of the entire simulated 
landscape (i.e., all breeding populations combined) will result in the harvest being taken uniformly across each management 
unit.  See \code{\link{tossm-package}} for more details on harvest intervals.}
\item{ genetic.sampler}{The genetic sampler dictates the sampling design used for obtaining genetic samples. The default (\code{def.genetic.sampler}) draws a sample that is random with respect to age, sex, and all other demographic characteristics.  The distribution of samples within sampling sites is proportional to density.  If density is uniform across a sampling sites, samples will be uniformly distributed.  Otherwise, the density of samples will reflect the density of animals across the sampling site.  The default sampler should be fine for most applications, but the user can supply an alternate sampler if desired.}
\item{ abund.for.10pc.CV}{The coefficient of variation (CV) of the abundance estimates is assumed to be proportional to 1/sqrt( abund).  The constant of proportionality is chosen so that the CV is 0.1 at the value specified by this argument.  The default value for this argument (70,000) results in a CV of 0.30 when abundance equals 7,500.}
\item{ quota.calc}{The algorithm to use for calculating catch quotas.  There are two quota calculating algorithms included with the \code{TOSSM} package - \code{CLA} (IWC, 1994), which is the algorithm used by the International Whaling Commission to set quotas, and \code{PBR} (Taylor et al. 2000), which is the algorithm for calculating quotas under the U.S. Marine Mammal Protection Act.}
\item{ quota.args}{Extra arguments passed to the quota calculating algorithm. \code{multiplier} is used to scale the quota-- default is 1X.  \code{multiplier} is used by both \code{CLA} and \code{PBR}.  \code{r.max} and \code{F.r} are only used by \code{PBR}.  \code{F.r} is the recovery factor.  \code{r.max} is the maximum population growth rate for the species in question.  The default value for \code{r.max} (0.04) is the default value used in PBR calculations for any cetacean species for which a species-specific estimate is not available.}
\item{ CLA.prog}{The path to the Fortran program that implements the CLA. The default should be correct unless MANAGE-D.exe or MANAGE-D (in the case of unix-alikes) has been moved subsequent to installation of the \pkg{tossm} package.}
\item{ CLA.dir}{The location of some temporary files created by the CLA program}
\item{ plot.polys}{A flag indicating whether or not the breeding polygons, sampling polygons, and historic removal polygons should be plotted at the beginning of each simulation.  Defaults to FALSE, however, if a problem with the polygons is detected by \code{tossm.diagnostics}, the polygons will be plotted even if this argument is set to false.}
\item{ seed}{Integer. If \code{seed>0}, it ensures a reproducible sequence of datasets.}
}
\value{
  A list of class \code{tossm.obj}. These components can be accessed through the tossm object resulting from a simulation (eg. x<-run.tossm(...), x component).
  
  \item{abund.b }{a matrix tracking the abundance of the breeding populations through the simulation.}
  \item{catches }{a matrix tracking the catch that takes place in each breeding population polygon through the simulation.}
  \item{effort }{a vector tracking the average distance (in the x-dimension only) that hunters had to travel in order to catch the animals killed during the simulation (assuming that hunters are based at the left-hand edge of the landscape).}
  \item{mu.hist }{a list consisting of the Management Unit/historic polygons, and the catches taken from each MU/historic polygon, in each year in which animals are harvested}
  \item{est.abund.mu }{the abundance estimates for each Managment Unit throughout the simulation.}
  \item{var.abund.mu }{the variance of abundance estimates for each Management Unit throughout the simulation.}
  \item{gs }{a list containing all of the data from genetic samples taken during the simulation by year and by sampling area polygon.  Each sample is also attributed with IDs of the individuals sampled as well as coordinate information. 
This information can be accessed using \code{unclass}.}
  \item{agg.gs, agg.gtypes }{the same as \code{gs} but with genetic samples from each sampling polygon aggregated by year.  See \code{\link{agg.gs.tseries}} for description of formats.  Each aggregated sample is also attributed with IDs of the individuals sampled,the birth year of each individual, as well as coordinate information. 
This information can be accessed using \code{unclass}.}
  \item{agg.gfreq }{a single 3-D array with dimension (number of sampling polygons, number of loci, max number of alleles at a locus) containing allele frequencies for each locus within each sampling polygon.}
  \item{call }{the original call to run.tossm}
  \item{seed }{the seed used for the random number generator.}  
}
\seealso{
\code{\link{tossm-package}}, \code{\link{hyptest.network.BSA}}}
}
\references{

IWC. 1994  The Revised Management Procedure (RMP) for Baleen Whales. \emph{Rep. Int. Whal. Commn.} \bold{44}, 145-167. 

Taylor, B.L., P.R. Wade, D.P. DeMaster, and J. Barlow. 2000  Incorporating uncertainty into management models for marine 
mammals. \emph{Conservation Biology} \bold{14}, 1243-1252. 
}
\examples{

#set up sampling and management schedule using def.make.schedule
schedule=def.make.schedule(n.pre.RMP=5,n.RMP=10,n.post.RMP=5,abund.gap=2)

#Example 1: Use initial.depletion to specify pre-modern management depletion
example1<-run.tossm(rland = rland.end, bp.polys=bp.polys,
schedule = schedule,n.samples = 25, sample.polys=sample.polys, 
initial.depletion=.30,historic.removals=NULL, 
BSA = fixed.MU.BSA, BSA.args = list(n.mus=5),harvest.interval=10,plot.polys=TRUE)

\dontrun{
#Example 2: Use historic.removals to specify pre-modern management depletion

#specify historic removals from each historic harvest area for 10 pre-modern
#management years
removals.h1<-rep(75,5)
removals.h2<-rep(0,5)
historic.removals<-cbind(removals.h1,removals.h2)

example2<-run.tossm(rland = rland.end, bp.polys=bp.polys,
schedule = schedule,n.samples = 25, sample.polys=sample.polys, 
historic.removals=historic.removals,historic.polys = historic.polys, 
BSA = fixed.MU.BSA, BSA.args = list(n.mus=5),harvest.interval=10,plot.polys=T)

}
}
\keyword{internal}