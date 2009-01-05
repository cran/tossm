`print.tossm.obj` <-
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

