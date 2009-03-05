\name{def.make.schedule}
\alias{def.make.schedule}
\title{Schedule constructor for run.tossm}
\description{Sets up a schedule of events that will take place during the years that are simulated by \code{\link{run.tossm}}.}
\usage{def.make.schedule(n.pre.RMP, n.RMP,n.post.RMP, abund.gap=5)}
\arguments{
\item{ n.pre.RMP, n.RMP, n.post.RMP}{number of years for historic harvest phase, managed harvest phase, and post-harvest recovery phase}
\item{ abund.gap}{the interval (in years) at which abundance estimates will be generated during the managed harvest phase}
}
\value{
A list with these components:
\item{stop.years}{total number of years in the simulation}
\item{gs.years}{year(s) in which genetic samples will be taken}
\item{abund.est.years}{year(s) in which abundance estimates will be taken}
\item{pre.RMP.years}{number of years simulation will be run prior to RMP phase}
\item{CLA.years}{year(s) in which catch-limit algorithm will be applied to breeding populations}
\item{post.RMP.years}{number of years of population recovery after the RMP phase}
\item{ BSA.years}{year(s) in which BSA should be implemented}
}
\details{

Though it is possible to set up a simulation schedule manually by creating vectors for each of the schedule components and passing them as arguments to \code{run.tossm} individually, the \code{def.make.schedule} function is a convenient way of setting up a schedule of all simulation events. 

The different events in the schedule are contingent on each other.  For instance, a harvest quota cannot be calculated until after an abundance estimate is generated.  Using \code{def.make.schedule} to set up a schedule of events ensures these contingencies are met. 

\code{def.make.schedule} generates a schedule in which genetic samples are taken (\code{gs.years}) and boundaries are set (\code{BSA.years}) in the last pre-RMP year only.  Abundance estimates are obtained (\code{abund.est.years}) and a quota is calculated (\code{CLA.years}) in the first year of managed harvest and every \code{abund.gap} years thereafter until the end of the managed harvest phase.  Thus, for the following inputs to \code{def.make.schedule}:

\code{my.schedule <- def.make.schedule(n.pre.RMP=10,n.RMP=20,n.post.RMP=10,abund.gap=5)}

the function returns a list with the following components: \cr

\code{>my.schedule} \cr
\code{$stop.years} \cr
\code{[1] 40} \cr
\code{$gs.years} \cr
\code{[1] 10} \cr
\code{$abund.est.years} \cr
\code{[1] 11 16 21 26} \cr
\code{$pre.RMP.years} \cr
\code{[1]  1  2  3  4  5  6  7  8  9 10} \cr
\code{$CLA.years} \cr
\code{[1] 11 16 21 26} \cr
\code{$post.RMP.years} \cr
\code{[1] 31 32 33 34 35 36 37 38 39 40} \cr
\code{$BSA.years} \cr
\code{[1] 10} \cr


By default, there is only one collection of genetic data, at the end of the pre-RMP phase: so \code{gs.years} will be the same as \code{n.pre.RMP}.
If \code{n.RMP==0}, genetic data will be collected in the final year, but the boundary setter will not be called.

If the user wants to set up a different schedule (e.g. multi-year genetic samples, and/or multiple boundary setting opportunities), some or all of the above components to the schedule can be passed directly to \code{\link{run.tossm}}.  
When manually inputting a schedule of events, rather than using \code{def.make.schedule}, the user should be aware of the following dependencies:

1) Genetic samples (\code{gs.years}) must be obtained before boundaries are set using the BSA (\code{BSA.years}) (collecting them in the same simulation year is fine).

2) Boundaries (\code{BSA.years}) must be set and abundance estimates (\code{abund.est.years}) must be obtained before an initial catch limit is set (\code{CLA.years}).  In other words, the first element in \code{BSA.years} and \code{abund.est.years} must be less than or equal to the first element in \code{CLA.years}.  Subsequent calls to the CLA, however, can be made without additional calls to the BSA or new abundance estimates.

3) All CLA.years take place during the RMP phase.



}
\examples{

#Example 1: just using def.make.schedule to generate a schedule

#set up sampling and management schedule using def.make.schedule
schedule=def.make.schedule(n.pre.RMP=5,n.RMP=10,n.post.RMP=5,abund.gap=2)

example1<-run.tossm(rland = example.landscape, bp.polys=example.bp.polys,
schedule = schedule,n.samples = 25, sample.polys=example.sample.polys, 
initial.depletion=0.30,historic.removals=NULL,
BSA = fixed.MU.BSA, BSA.args = list(n.mus=5),harvest.interval=10)

\dontrun{
#Example 2: specify extra genetic sampling years beyond the def.make.schedule default
#using gs.years
example2<-run.tossm(rland = example.landscape, gs.years=c(10,12,14),bp.polys=example.bp.polys,
schedule = schedule,n.samples = 25, sample.polys=example.sample.polys, 
initial.depletion=0.30,historic.removals=NULL, 
BSA = fixed.MU.BSA, BSA.args = list(n.mus=5),harvest.interval=10)
}
}
\seealso{
\code{\link{run.tossm}},\code{\link{tossm-package}}
}
\keyword{internal}