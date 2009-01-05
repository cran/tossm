`fast.trans` <-
function( m, to1=FALSE) {
  ttmat <- matrix( cumsum( m), byrow=TRUE, ncol( m), nrow( m))
  ttmat <- ttmat - c( 0, clip( ttmat[ ,nrow( m)]))
  if( to1)
    ttmat[ , ncol( ttmat)] <- 2 # runif( 1) guaranteed < 2...
  ttmat
}

