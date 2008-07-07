`generate.coords` <-
function(sample.poly,n.samples){
   coords <- NULL
   if (area.poly(sample.poly) > 0){
	box.range <- get.bbox(sample.poly)
	box <- as(rbind(cbind(box.range$x[1],box.range$y),cbind(box.range$x[2],c(box.range$y[2],box.range$y[1]))),"gpc.poly")
	ratio <- area.poly(sample.poly)/area.poly(box)
	while (is.null(coords) || dim(coords)[1] < n.samples){
		x <- runif(ceiling(n.samples/ratio),min=box.range$x[1],max=box.range$x[2])
		y <- runif(ceiling(n.samples/ratio),min=box.range$y[1],max=box.range$y[2])
		pts <- cbind(x,y)
		pts <- pts[in.gpc(sample.poly,pts),]
		coords <- rbind(coords,pts)
		if (dim(coords)[1] > n.samples) coords <- matrix(coords[1:n.samples,],n.samples)
	}
   }
   return(coords)
}

