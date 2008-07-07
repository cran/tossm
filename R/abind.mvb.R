`abind.mvb` <-
function( dim., ...){
  l <- list( ...)
  ndim <- length( diml1 <- dim( l[[1]]))
  new.order <- c( 1:ndim %except% dim., dim.)
  ans <- do.call( 'c', lapply( l, aperm, new.order))
  dim( ans) <- c( diml1[ -dim.], length( ans) / prod( diml1[ -dim.]))
  aperm( ans, match( 1:ndim, new.order))
}

