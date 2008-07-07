`def.set.coords` <-
function( gs.y){
# gs.y is genetic samples in 1 year
  for( i in 1:length( gs.y)) {
   if(!is.null( attr( gs.y[[i]], 'coords'))){
    ids <- unclass(attr(gs.y[[i]],'coords'))
    coords <- cbind( x=i + runif( nrow( gs.y[[i]])), y=runif( nrow( gs.y[[i]])))
    coords <- cbind(ids,coords)
    attr(coords, 'printinfo') <- "<< individual ids and spatial coordinates for samples; use 'unclass(...)' to display >>"
    class ( coords) <- 'printinfo'
    attr( gs.y[[ i]], 'coords') <- coords
   }
  }
  
  gs.y
}

