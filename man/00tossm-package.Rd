\name{tossm-package}
\alias{tossm-package}
\docType{package}
\title{Testing of spatial structure methods}
\description{
\bold{Note...}It is highly recommended for those uninitiated to the \pkg{tossm} package to read this introduction 
before proceeding to the specifics of the individual \pkg{tossm} functions.

The Testing of Spatial Structure Methods (\pkg{TOSSM}) project was developed as a collaboration between NOAA Fisheries and
the International Whaling Commission (IWC).  It is a long-term research project aimed at testing the performance of 
different analytical methods of detecting population structure from genetic data.  The \pkg{TOSSM} package was developed to
provide a framework for conducting simulation-based performance testing of genetic analytical methods \emph{in a 
management context}.  Methods to be tested using the \pkg{TOSSM} package must perform two main steps: 1) detect spatial
structure in genetic samples collected from a set of populations and 2) set management unit boundaries based on any
spatial structure detected.  A method that performs these two steps is referred to in the \pkg{TOSSM} framework as a 
boundary-setting algorithm (BSA); the concept of a BSA will be discussed in detail below.

Simulations conducted in the TOSSM framework are divided into four phases: 1) an ancient phase represented by an 
initial dataset input into the program, 2) a historic harvest phase that 
mimics the pre-modern exploitation to which many populations have been subjected, 3) a modern management phase in which 
populations are managed according to the management units defined by a BSA and a harvest quota calculated using one of 
two quota calculating algorithms, and 4) an optional recovery period in which no harvest occurs.

A TOSSM simulation is spatially explicit, with the spatial location of important components defined by polygons. These 
spatial components include breeding populations, historic harvest areas, genetic sampling areas, and management units.
Polygons used in a TOSSM simulation are of class \code{gpc.poly}, implemented by the \pkg{gpclib} package. Therefore,
some of the required inputs to TOSSM are polygons of this class. Additionally, the BSA algorithm accepts and returns these 
\code{gpc} polygons.

The \pkg{tossm} package requires as input an initial population upon which the simulation acts.  Many population datasets
 have been created using the \code{\link[<rmetasim>:<rmetasim>-package]{<rmetasim>}} package (Strand, 2002) for this purpose. These initial datasets, henceforth referred to as the TOSSM 
datasets, represent a variety of population structure scenarios and are all in genetic and demographic equilibrium.  
The objective behind developing these simulated scenarios was to provide standard datasets on which various methods for 
analyzing genetic data could be tested.  The TOSSM datasets can be obtained from the following web page:\cr
\url{http://swfscdata.nmfs.noaa.gov/TOSSM/}


As currently written, \pkg{tossm} is closely tied to the Potential Biological Removal (PBR) scheme of the U.S. Marine Mammal Protection Act,
 and the IWCs Revised Management Plan (RMP), both of which attempt to manage and conserve marine mammals that are subject 
to human-caused mortality (IWC, 1994). However, the \pkg{tossm} framework is also applicable to the management of 
virtually any animal population that exhibits spatial structure.


The main function for running a \pkg{tossm} simulation is \code{\link{run.tossm}}. The \code{run.tossm} function takes as arguments 
a number of subsidiary functions detailed below and in separate help files (see \code{seealso} section below).



\strong{DEFINITION OF TERMS:  important tossm concepts}

\strong{\emph{Archetypes and Breeding populations (BPs)}} \cr

The initial datasets used as input for a \pkg{tossm} simulation fall into five broad categories of population structure, referred to as \emph{Archetypes},
 each representing a different spatial structure and migration pattern for the animal population being simulated.
 The number of breeding populations in a \pkg{tossm} simulation is determined by the number of breeding populations that 
exist in the initial simulated dataset input into \code{run.tossm} . These original datasets exist as
 \code{\link[<rmetasim>:<rmetasim>-package]{<rmetasim>}} landscape objects. To generate each initial dataset, the number of breeding populations (BPs), carrying 
capacity for each BP, and dispersal rate between BPs, were specified. More details on TOSSM dataset generation are included in a handbook available
at the following web address: 

 \url{http://swfsc.noaa.gov/textblock.aspx?Division=PRD&ParentMenuId=496&id=11126}


The five archetypes represented in these initial \code{<rmetasim>} datasets are the following:

\emph{Archetype I:}  \cr A single, panmictic population that serves as a control \cr
\emph{Archetype II:} \cr Stepping-stone dispersal pattern between two or three populations, with dispersal only occurring between adjacent populations 
 mixing \cr
\emph{Archetype III:} \cr Diffusion-type, where genetic isolation between individuals occurs \cr continously as a function of distance. \cr
\emph{Archetype IV:} \cr Two discrete breeding grounds with feeding grounds that overlap \cr partially or completely. \cr
\emph{Archetype V:} \cr A single breeding population with two separate feeding grounds. \cr


The initial \code{<rmetasim>} landscape population used for a simulation consists of a specified number of breeding
 populations (referred to as habitats in \code{<rmetasim>}). In a \pkg{tossm} simulation, a list of polygons 
(\code{bp.polys}) of length equal to the number of these initial populations must be input as an argument to 
\code{run.tossm}.  These polygons define the geographic range of the breeding populations. They can be geographically 
discrete, contiguous, or overlapping. More details can be found in the help file for \code{\link{run.tossm}}.




\strong{\emph{Sampling Polygons}} \cr
(\code{sample.polys} argument to \code{run.tossm})

Sampling polygons define the spatial areas from which genetic samples are collected in a simulation. There is flexibility as to
the number and geographic extent of these polygons. The constraints are that, 1) the sampling polygons must be 
within the extent of the breeding population polygons and, 2) there can be no overlap between sampling polygons.
Uniform sampling of the entire simulated landscape can be achieved by defining a single sampling polygon equal 
to the combined spatial extent of the breeding populations.The genetic information collected in the sampling polygons is 
fed to the boundary-setting algorithm (BSA), which in turn sets management units in a \pkg{tossm} simulation. 


For instructions on how to input the number, location, and size of sampling polygons for a simulation, 
see \code{\link{run.tossm}}.


\strong{\emph{Boundary-setting algorithms (BSAs)}} \cr
(\code{BSA} argument to \code{run.tossm})

A BSA is a program that interfaces with \code{run.tossm} by accepting simulated genetic and abundance information 
from \code{run.tossm}, and then supplying \code{run.tossm} with a recommended way of dividing the landscape 
 into management units. A BSA does not have to use all information about a simulated population that \code{run.tossm} 
makes available; different BSAs may rely on different components of the simulated data provided by \code{run.tossm}. 

There are many analytical methods that accept and analyze genetic data and output information on how a population is
 genetically structured. A BSA does this and goes one step further by making a management decision based on this genetic
 data. Specifically, the BSA decides if and how the breeding populations should be split spatially into management units. 
 This could be as simple as deciding whether to manage two sampling polygons separately or as one management unit. Alternatively,
 if there are many sampling polygons, there could be a number of different decisions the algorithm must make about 
splitting or grouping the various polygons into management units. An added level of complexity is needed if the BSA works
not on the sampling polygon level but at the level of the individual, in which case the BSA decides which 
\emph{individuals}, each individual having x- and y-coordinates, should belong to which management unit. 


Examples of BSAs and more information about producing a BSA are provided with \cr \code{\link{hyptest.network.BSA}}.


\strong{\emph{Quota calculating algorithms}}\cr
(\code{quota.calc} argument to \code{run.tossm})

Two options are currently available for calculating harvest quotas. The default quota calculating algorithm 
is the Potential Biological Removal (PBR) scheme used in the U.S. Marine Mammal Protection Act (Taylor et al., 2000).
The alternative provided with the package is the catch-limit algorithm, or CLA (IWC, 1994).  It calculates quotas based
 on the estimated abundance of a management unit and information on historic catches. See Cooke (1994) for further details.
The modular structure of the TOSSM package makes it straightforward for users to define their own quota-calculating
algorithm for use with \code{run.tossm}.  


\strong{\emph{Harvest Intervals}}\cr
(\code{harvest.intervals} argument to \code{run.tossm})

In TOSSM simulations, harvesting effort is not uniformly distributed across each management unit.  Rather, effort is 
concentrated near the left edge of each MU.  This spatial bias in harvest is meant to simulate a situation in which the 
harvesters wish to minimize the distance they must travel in order to meet their quota, and so concentrate their effort 
close to their home base, which is assumed to be to the left of the MUs.  To implement this spatial bias, the entire 
simulated landscape is divided into vertical strips (harvest intervals), the width of which is determined by 
the argument \code{harvest.interval}.  In each simulation year, \code{run.tossm} will attempt to take the entire quota for an MU from 
the left-most harvest interval in that MU.  If there are not enough animals present in the first harvest interval to 
meet the quota, all animals in the first interval will be harvested and the program will attempt to remove the remainder 
of the quota from the next interval to the right.  Harvest progresses toward the right until the quota has been met.

The degree to which harvest is spatially biased is controlled by changing the width of the harvest intervals.  
Defining intervals that are very narrow relative to the x-range of the breeding populations will result in a strong 
spatial bias.  Setting the harvest interval width equal to the x-range of the entire simulated landscape 
(i.e., all breeding populations combined) will result in the harvest being taken uniformly across each management unit.


\strong{\emph{Schedule of simulation events}}\cr
(\code{schedule} argument to \code{run.tossm}:  user may choose to use the function \code{def.make.schedule} for convenience in creating this schedule)

There are up to 4 phases in a simulated population period: 1) ancient, 2) historic, 3) managed, and, 4) recovery:


\emph{Ancient phase}: 
This phase constitutes the history of a population before any harvest has occurred, and consists of an \code{<rmetasim>} landscape object input into \pkg{tossm}.  These TOSSM datasets, which are a required input to \code{run.tossm}, 
represent the population at the end of the ancient phase. \cr
\emph{Historic phase}:  
This phase allows for historic harvest that pre-dates modern management (e.g., implementation
of the PBR or CLA). \cr
\emph{Managed phase}:  
In this phase, genetic information and abundance estimates are collected and fed into
 the BSA, which defines management units.  These management units are then managed using quotas calculated by either
the PBR or the CLA. \cr
\emph{Recovery phase}:  
This phase is an optional period that allows post-harvest population recovery before the end
 of the simulation. \cr



\strong{\emph{System requirements}} \cr

To run the catch-limit algorithm, the file MANAGE-D.exe (windows),MANAGE-D (unix), or MANAGEMAC-D (MacOS) must exist in the root folder of \pkg{tossm}, as they do upon installation of the package.
}
\references{

Cooke, J.G. 1994  The management of whaling. \emph{Aquatic Mammalogy} \bold{20}, 129-135.

IWC. 1994  The Revised Management Procedure (RMP) for Baleen Whales. \emph{Rep. Int. Whal. Commn.} \bold{44}, 145-167. 

Strand A.E. 2002  METASIM 1.0: an individual-based environment for simulating population genetics of complex population 
dynamics. \emph{Molecular Ecology} \bold{2}, 373-376.

Taylor, B.L., P.R. Wade, D.P. DeMaster, and J. Barlow. 2000  Incorporating uncertainty into management models for marine 
mammals. \emph{Conservation Biology} \bold{14}, 1243-1252. 
}
\details{
\tabular{ll}{
Package: \tab tossm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-07-20\cr
License: \tab ------\cr
}
}
\author{
Mark Bravington, Karen Martien, and Dave Gregovich

Maintainer: Karen Martien <Karen.Martien@noaa.gov>

}

\seealso{\code{\link{run.tossm}} \code{\link{hyptest.network.BSA}}}
\keyword{package}