\name{hyptest.network.BSA}
\alias{hyptest.network.BSA}
\alias{structure.BSA}
\alias{wombsoft.BSA}
\alias{fixed.MU.BSA}
\title{Boundary setting algorithms}
\description{

\code{hyptest.network.BSA}, \code{structure.BSA}, \code{wombsoft.BSA}, and \code{fixed.MU.BSA} are examples of boundary-setting algorithms (BSAs) used to set
 management boundaries in a \pkg{tossm} simulation.  
}
\usage{

hyptest.network.BSA(gs,abund,var,C,landscape.poly,sample.polys,sig.level) 
structure.BSA(gs,abund,var,C,landscape.poly,sample.polys,mainparams,BSA.params)
wombsoft.BSA(gs,abund,var,C,landscape.poly,sample.polys,sig.level,min.MU.size)
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
\item{sig.level}{significance value threshold for setting a management boundary in \code{hyptest.network.BSA}}
\item{n.mus}{number of MUs to be created by \code{fixed.MU.BSA}}
\item{mainparams}{parameters to be used by the program Structure in \code{structure.BSA}.  \code{mainparams} is a list consisting of the following components: \cr
\code{numinds}: total number of individuals sampled \cr
\code{numloci}: number of loci to be analyzed by Structure \cr
\code{ploidy}: ploidy of loci to be analyzed by Structure \cr
\code{noadmix}: whether populations are assumed to be admixed (1) or not (0) \cr
\code{freqscorr}: whether allele frequencies between populations are correlated (1) or not (0) \cr
\code{maxpops}: the maximum number of populations to be defined by Structure \cr
\code{burnin}: the length of the burnin period used by Structure \cr
\code{iterations}: the number of iterations of the MCMC chain to be used by Structure in assigning population membership \cr
Further information regarding these parameters can be found in the Structure user's manual}
\item{BSA.params}{a list with the following components to be used by \code{structure.BSA}: \cr
\code{output.dir}: the directory to which STRUCTURE results should be written \cr
\code{k}: a list of the values of k (number of populations) to be evaluated \cr
\code{struct.exe.path}: the directory where structure.exe resides and where main params will be written. 
The default value should be fine, so long as structure.exe is not moved subsequent to the installation of the package. \cr
\code{struct.outfile}: the name of the output file to be written by STRUCTURE.   
 \cr
}
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
\code{structure.BSA} uses the Bayesian clustering method implemented by the program STRUCTURE (Pritchard et al. 2000, Falush et al. 2002)
to determine the number of MUs that should be defined.  The BSA uses STRUCTURE to divide the samples into the numbers of populations
specified by \code{BSA.params$k} and compares the log-likelihoods of the resulting groupings to determine the optimal number of 
MUs.  The MU membership of each samplings site is then determined by examining the population membership, as defined by STRUCTURE,
of each sample within the sampling site.  The site is assigned to the MU to which the greatest number of samples belongs.
Thus, if a sampling site contained 50 samples, 20 of which were assigned to population 1 and 30 of which were assigned to population 2,
the entire sampling site would be assigned to MU 2.  N0ote that a consequence of this approach is that it may result in the definition 
of fewer MUs than there are STRUCTURE-defined populations.  For instance, if STRUCTURE defined 4 populations, but only assigned a few samples to
population 4, there may be no sampling site in which the plurality of samples belongs to population 4.  In this case, no sampling 
site would be assigned to MU 4 and there would effectively be only 2 MUS.
\cr
Once sampling sites have been assigned to MUs, the MU.polys are defined using the function \code{MU.poly.generator}
\cr
In order to use \code{structure.BSA}, the command-line version of STRUCTURE (i.e., without the graphical user interface), 
must first be installed.  This program is available in UNIX, DOS, and MacOS X versions from the STRUCTURE homepage.  The argument 
\code{BSA.params$struct.exe.path} is used to specify the directory in which STRUCTURE is installed.
\cr
\cr
\code{wombsoft.BSA} uses the R package \pkg{wombsoft} of Crida and Manel (2007). It first identifies regions within the study area where genetic gradients are relatively strong. The algorithm
then uses a binomial test to assess the significance of these areas as genetic boundaries. The significance threshold of this test can be chosen by the user using the 
argument \code{sig.level}. The output of the wombsoft algorithm is then adapted in the following two ways:
\cr
1) If only a single, non-significant grid-cell interrupts a potential MU boundary, the boundary is drawn across such 
one-cell gaps.
\cr
2) If the size of an MU is less than a user-specified proportion of the total study area (\code{min.MU.size}), the MU is ignored. This 
alleviates the fact that the algorithm sometimes outputs very small (one or two grid cell) MUs.
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

Crida A., and Manel S. 2007. WOMBSOFT: an R package that implements the Wombling method to identify genetic boundary.
\emph{Molecular Ecology Notes} \bold{7}, 588--591

Dupanloup I., Scheider S. and Excoffier L. 2002  A simulated annealing approach to define the genetic structure of populations.
\emph{Molecular Ecology} \bold{11}, 2571--2581

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
example<-run.tossm(rland = rland.end, bp.polys=example.bp.polys,
schedule = schedule,n.samples = 25, sample.polys=example.sample.polys, 
initial.depletion=0.30,historic.removals=NULL, 
BSA = hyptest.network.BSA, BSA.args = list(sig.level=.05),
harvest.interval=10)
}
\keyword{internal}