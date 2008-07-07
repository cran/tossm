`in.gpc` <-
function(poly,pts){
	contours <- get.pts(poly)
	res <- do.call("cbind",lapply(contours, function(c){
		poly.mat <- cbind(c$x,c$y)
		inout(pts,poly.mat)
	}))
	return(as.logical(rowSums(res)))
}

