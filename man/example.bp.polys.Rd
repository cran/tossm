\name{example.bp.polys}
\docType{data}
\alias{example.bp.polys}
\title{Breeding Population Polygons}
\description{
This is a list of two polygons representing the extent of the two breeding populations in the TOSSM
examples.

\emph{How to construct \pkg{tossm} polygon data objects:}

The following code constructs these breeding population polygons, and so can be modified if the user
wishes to create new polygons. Construction of these polygons requires the \pkg{gpclib} package that
comes with \pkg{tossm}:

\code{bp.polys<-list()}

\code{bp.polys[[1]]<-matrix(c(0,0,50,0,50,100,0,100),nrow=4,byrow=T)}

\code{bp.polys[[2]]<-matrix(c(50,0,100,0,100,100,50,100),nrow=4,byrow=T)}

\code{example.bp.polys<-lapply(bp.polys,as,"gpc.poly")}

}
\usage{example.bp.polys}
\format{list (length = 2) of gpclib polygon objects (\code{\link[<gpclib>:<gpclib>-package]{<gpclib>}})}
\keyword{datasets}