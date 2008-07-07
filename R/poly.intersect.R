`poly.intersect` <-
function(polys){
	poly.seed <- polys[[1]]
	if (length(polys) > 1) for (i in 2:length(polys)) poly.seed <- intersect(poly.seed,polys[[i]])
	poly.seed
}

