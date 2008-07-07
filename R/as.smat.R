`as.smat` <-
function( mat) {
  if( is.list( mat))
return( mat)
  o <- order( col( mat), mat)
  mat[] <- cumsum( mat[o]) - col( mat) + 1
return( list( lookup=t( mat), unsort=o))
}

