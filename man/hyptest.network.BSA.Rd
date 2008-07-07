\name{hyptest.network.BSA}
\alias{hyptest.network.BSA}
\alias{fixed.MU.BSA}
\title{Boundary setting algorithms}
\description{

\code{hyptest.network.BSA} and \code{fixed.MU.BSA} are examples of boundary-setting algorithms (BSAs) used to set
 management boundaries in a \pkg{tossm} simulation.  
}
\usage{

hyptest.network.BSA(gs,abund,var,C,landscape.poly,sample.polys,sig.level) 
fixed.MU.BSA(gs,abund,var,C,landscape.poly,sample.polys,n.mus)
}
\arguments{
\item{gs}{genetic samples; a list of arrays of the alleles at each locus for each individual animal from which genetic
 samples were taken, arranged by genetic sampling year and sampling polygon. Each list is a list of length = number of sampling polygons
 (\code{sample.poly}), and each element of these lists is an array of dimensions (number of genetic samples per
 \code{sample.poly}) X (number haploid loci + number of diploid loci) X (2 haplotypes)}
\item{abund}{matrix of abundance estimates for each management unit sampling year}
\item{var}{variance of these abundance estimates}
\item{C}{catches by simulation year and management unit}
\item{landscape.poly}{a landscape polygon object defining the study area extent (the bounding box of the 
\code{bp.polys})}
\item{sample.polys}{a list of sampling polygons to be lumped or split into separate management units by the BSA}
\item{sig.level}{significance value threshold for setting a management boundary}
\item{n.mus}{number of MUs to create}
}
\value{\code{munits}--A list with each element consisting of a management unit, represented by a polygon of class 
\code{gpc.poly}.}

\details{


\code{hyptest.network.BSA} uses the genetic data in \code{gs} to create networks of sampling polygons which are joined 
together by non-significant statistical tests (G-statistic, Goudet et al., 1996).  This method was proposed by Waples and 
Gagiotti (2006), and is illustrated further there (p.1425).
\cr
\code{hyptest.network.BSA} checks for pairwise significance between all sampling polygons.  Sampling polygons are joined together in a 
network if they do not show significant genetic differentiation (at the alpha level chosen by the user).  Network membership
extends to all sampling polygons connected by at least one non-significant test; ie. a sampling polygon might not be
directly connected to all sampling polygons in its assigned network. Networks are separated from each other by exclusive
significant tests with all other sampling polygon networks.  Each resulting network constitutes a management unit.
\cr
\cr
\code{fixed.MU.BSA} simply splits the study region into equally-sized management units, and disregards all genetic
 information.  The number of management units is specified by the user.  When initially running TOSSM simulations, the 
user may find it useful to use \code{fixed.MU.BSA} for testing purposes.
\cr

\strong{\emph{How to write a BSA that will run within the \code{run.tossm} code}}\cr


The BSA in essence needs to do two things: 1) analyze genetic and/or abundance data from \code{run.tossm} and,  2) set 
management units. The management units output from the BSA must completely cover the extent of the breeding population 
polygons, yet not overlap each other.  They must be returned from the BSA in the form of a list, each element consisting
of one management unit polygon.  These polygons should be of class \code{gpc.poly}, implemented via \pkg{gpclib}.

The analysis is likely to be done exclusively by an external program not written in R.  If your BSA is going to
call outside genetic software, you will probably need to write the genetic data to disk in some specified format, then
call the software, then read the results back in to R and reorganize them for return to \code{run.tossm}.  
Regardless, there needs to be at least some R code written, even if just a wrapper that allows the analytic method 
and \code{run.tossm} to interface.

There are a number of examples of mixing and migration models that have been proposed for use as analytic methods 
(Pritchard et al., 2000; Dupanloup et al, 2002; Guillot et al., 2005).  It is foreseen that some of these methods can 
stand alone or be combined to supply the information desired to analyze \pkg{tossm} simulated data and define management 
units. 

\bold{Note--} there is an important distinction between two types of potential BSAs: 

1) Those BSAs that simply assign each \code{sample.poly} to a management unit.  In this type of BSA, statistical analysis 
of genetic samples will be done at the level of the sampling polygon.  For such BSAs, the helper function 
\code{MU.poly.generator} described below may be of use to the BSA developer.  

2) Those BSAs that work at the level of individual genetic samples.  This type of BSA might disregard which sampling 
polygon samples come from.  In this case, management unit boundaries could quite possibly bisect sampling polygons and/or, 
conversely, assign individuals from different sampling polygons to the same management unit.  \code{MU.poly.generator} 
will not be of use to those developing such BSAs, and the developer is tasked with ensuring that a list of 
non-overlapping management unit polygons that cover the entire study area (=\code{landscape.poly}) is returned from the BSA.

Whatever the method used, the R code written for the BSA takes the following arguments:


\code{My.BSA<- function(genetic samples, abundance estimates, variances,  catches,}\cr
\code{    landscape.poly, sample.polys, optional param1, optional param2)}

The first six of these arguments are passed automatically to the BSA function by \code{run.tossm}, so the function 
must have these six arguments whether it uses the values contained in them or not.  They correspond to the first six 
arguments to \code{hyptest.network.BSA}, and are described in the arguments section above.

The last two arguments can optionally be used to pass any additional information that might be required to run the
 BSA function.  For example, \code{hyptest.network.BSA} uses only the first optional parameter, an alpha level 
for genetic comparisons it performs, and \code{fixed.MU.BSA} also uses only one optional argument, \code{n.mus}, the number
of management units the user wishes to define.  Optional arguments for the BSA are specified within the \code{run.tossm} 
argument  \code{BSA.args}.  Each of the arguments (supplied by \code{run.tossm} or user-defined) can be named as the user
 wishes.  The optional arguments can be of any type. Thus, if a BSA requires more than two additional parameters, 
the optional arguments can be used to pass lists of parameters.



\emph{Supplying further information to the BSA}--\code{agg.gs.tseries}

The function \code{agg.gs.tseries} can provide the BSA with the \code{run.tossm} genetic data from each \code{sample.poly} aggregated across years, the number of loci in the current simulation, and the number of alleles in each \code{sample.poly}.  This information can be used by the BSA simply by placing following line first in the BSA function:

\code{agg.gs.tseries()} :creates the objects \code{agg.gs,n.loci,n.areas} ,and \code{n.alleles} \cr

Here are details on the objects made available by \code{agg.gs.tseries}:\cr

\code{agg.gs}\cr 
similar to \code{gs} (see arguments section above) except with genetic samples aggregated across years for each genetic sampling
polygon.  \code{agg.gs} is a list with length equal to the number of sampling polygons. Each list element consists of an array with dimensions
(number of genetic samples, entire simulation) X (number haploid loci+number of diploid loci) X 2. Each element
has a \code{coords} attribute that provides x and y coordinates for each sampled individual.  The \code{agg.gs} object is also attributed with
\code{seq.list}, a list of the unique haplotype states sampled in the simulation.

\code{agg.gtypes} \cr
the same  as \code{agg.gs} except that for each sampling polygon, the genetic data represented by a
matrix with dimensions
(number of genetic samples, entire simulation) X (number haploid loci+ 2*number of diploid loci)

\code{n.loci} \cr
The number of loci used in the current simulation.

\code{n.areas} \cr 
The number of sampling polygons in the current simulation. 

\code{n.alleles} \cr 
the number of alleles at each locus in each sampling polygon 

BSAs that work by assigning entire sampling sites (rather than individual samples) to MUs can be assisted by 
the function \code{MU.poly.generator}, which creates MU polygons based on the MU membership of the sampling polygons. Further
info on \code{\link{MU.poly.generator}} is available on its help page.
}
\seealso{ \code{\link{tossm-package}}, \code{\link{run.tossm}}}
\references{
Dupanloup I, Scheider S. and Excoffier L. 2002  A simulated annealing approach to define the genetic structure of populations.
\emph{Molecular Ecology} \bold{11}, 2571--2581.

Guillot G., Mortier F. and Estoup A. 2005  GENELAND: a computer package for landscape genetics. 
\emph{Molecular Ecology Notes} \bold{5}, 712--715

Pritchard J.K., Stephens, M. and Donnelly P. 2000  Inference of population structure using multilocus genotype data. 
\emph{Genetics} \bold{155}, 945--59.

Waples R.S. and Gaggiotti O. 2006  What is a population? An empirical evaluation of some genetic methods for identifiying 
the number of gene pools and their degree of connectivity. \emph{Molecular Ecology} \bold{15}, 1519--1539.
}
\examples{

#set up sampling and management schedule using def.make.schedule
schedule=def.make.schedule(n.pre.RMP=5,n.RMP=10,n.post.RMP=5,abund.gap=2)

#run example with BSA set to hyptest.network.BSA
example<-run.tossm(rland = rland.end, bp.polys=bp.polys,
schedule = schedule,n.samples = 25, sample.polys=sample.polys, 
initial.depletion=0.30,historic.removals=NULL, 
BSA = hyptest.network.BSA, BSA.args = list(sig.level=.05),
harvest.interval=10)
}
\keyword{internal}