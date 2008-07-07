`set.interval.polys` <-
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

