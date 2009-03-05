#Utility functions from the TOSSM package.  It is unlikely these will need to be
#altered or viewed by package developers.  This file contains the following
#functions in the following order:
#
#.First.task
#abind.mvb
#count.whales.b
#divvy.catch
#fast.trans
#fixed.MU.BSA
#in.gpc
#poly.intersect
#poly.union
#print.printinfo
#print.tossm.obj
#rsample
#set.interval.polys
#sim.abund.est
#single.MU.BSA
#TAC.split.by.f
#sqr
#original!object!list



".First.task" <-
function(...){
  library( rmetasim)
}


"abind.mvb" <-
function( dim., ...){
  l <- list( ...)
  ndim <- length( diml1 <- dim( l[[1]]))
  new.order <- c( 1:ndim %except% dim., dim.)
  ans <- do.call( 'c', lapply( l, aperm, new.order))
  dim( ans) <- c( diml1[ -dim.], length( ans) / prod( diml1[ -dim.]))
  aperm( ans, match( 1:ndim, new.order))
}


"count.whales.b" <-
function( rland, n.b){
  tabulate( landscape.populations( rland), n.b)
}


"divvy.catch" <- 
function(locs, mu.polys){
	sapply(mu.polys, function(mu){
		sum(in.gpc(mu,locs[,3:4]))})
}


"fast.trans" <-
function( m, to1=FALSE) {
  ttmat <- matrix( cumsum( m), byrow=TRUE, ncol( m), nrow( m))
  ttmat <- ttmat - c( 0, clip( ttmat[ ,nrow( m)]))
  if( to1)
    ttmat[ , ncol( ttmat)] <- 2 # runif( 1) guaranteed < 2...
  ttmat
}

"fixed.MU.BSA"<-
function(gs, abund, var, C, landscape.poly,sample.polys, n.mus){
      set.interval.polys(landscape.poly,n.intervals=n.mus)}

"in.gpc" <- 
function(poly,pts){
	contours <- get.pts(poly)
	res <- do.call("cbind",lapply(contours, function(c){
		poly.mat <- cbind(c$x,c$y)
		inout(pts,poly.mat)
	}))
	return(as.logical(rowSums(res)))
}

"poly.intersect" <- 
function(polys){
	poly.seed <- polys[[1]]
	if (length(polys) > 1) for (i in 2:length(polys)) poly.seed <- intersect(poly.seed,polys[[i]])
	poly.seed
}

"poly.union" <-
function(polys){
	poly.seed <- polys[[1]]
	if (length(polys) > 1) for (i in 2:length(polys)) poly.seed <- union(poly.seed,polys[[i]])
	poly.seed
}

"print.printinfo" <-
function( x, ...){
  x <- attr( x, 'printinfo')
  NextMethod()
}


"print.tossm.obj" <-
function( x, ...){
  callo <- x$call
  funco <- function( x) nchar( paste( deparse( x), collapse=' '))
  too.long <- sapply( callo, funco) > 60
  callo[ too.long] <- quote( `<<long thing>>`)
  x <- x[ cq( abund.b, abund.f, catches, effort, FIMA.mu.map)]
  x$`gs,agg.gs,agg.gfreq,est.abund.f,var.abund.f` <- quote( ...)
  x$call <- callo
  NextMethod()
}


"rsample" <-
function( n, pop, replace=FALSE, prob=NULL) 
  sample( pop, size=n, replace=replace, prob=prob)



"set.interval.polys" <-
function(landscape.poly, interval=NULL, n.intervals=NULL){
  box <- landscape.poly
  if (is.null(interval)){interval <- (box$x[2]-box$x[1])/n.intervals}
  x.coords <- seq(box$x[1],box$x[2],interval)
  if (x.coords[length(x.coords)] < box$x[2] & 
	!identical(round(x.coords[length(x.coords)],5),round(box$x[2],5)))
		x.coords <- c(x.coords,box$x[2])
  lapply(1:(length(x.coords)-1), function(i){
	as(matrix(c(x.coords[i],x.coords[i+1],x.coords[i+1],x.coords[i],
		rep(box$y[1],2),rep(box$y[2],2)),nrow=4),"gpc.poly")
  })
}	


"sim.abund.est" <-
function( abund.b, abund.for.10pc.CV,bp.in.MU.map){
  truth <- colSums(bp.in.MU.map*abund.b)
  shape <- sqr(10)*truth/abund.for.10pc.CV
returnList( est=rgamma( length( truth), shape=shape, scale=truth/shape), CV=1/sqrt( shape))
}


"single.MU.BSA" <-
function( gs, ea, ...) list( 1:ncol( ea)) # easiest way to ascertain #FIMAs


"TAC.split.by.f" <-
function( TAC, fg.to.mu, abund.f){
  # TAC is per-mu
  # Take TAC from smallest-index FIMAs first
  TAC.f <- 0*fg.to.mu
  for( imu in 1:length( TAC)) {
    fgi <- index( fg.to.mu==imu)
    cumabund <- cumsum( c( 0, clip( abund.f[ fgi])))
    remainder <- TAC[ imu] - cumabund
    TAC.f[ fgi] <- pmax( 0, pmin( remainder, abund.f[ fgi]))
  }
return( TAC.f)

#
#  # Old version: split the per-mu TAC across those fg's that make up the mu
#  fg.per.mu <- tabulate( fg.to.mu)
#  TAC <- floor( TAC / fg.per.mu)
# return( TAC[ fg.to.mu])
}


"sqr" <-
function( x) x*x


`original!object!list` <-
c(".First.task", "abind.mvb", "adjust.pop.dyn", "agg.gs.series", 
"agg.gs.tseries", "as.smat", "breed.grow.die", "calc.effort", 
"count.alleles.by.f", "count.whales.b", "def.genetic.sampler", 
"def.make.schedule", "def.set.coords", "fast.trans", "FIMA", 
"fixed.BSA", "get.SRM", "harvest", "humpbackism", "km.simple.demo", 
"mixmat", "old.bgd", "prepare.rland", "print.printinfo", "print.tossm.obj", 
"randomism", "raw.CLA", "rsample", "run.mvb.example2", "run.tossm", 
"hyptest.seq.BSA", "set.TAC", "sim.abund.est", "single.MU.BSA", 
"TAC.split.by.f", "sqr", "CLC.PAR", "rlsimp2", "rlsimp4")



