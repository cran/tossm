`poly.union` <-
function(polys){
	poly.seed <- polys[[1]]
	if (length(polys) > 1) for (i in 2:length(polys)) poly.seed <- union(poly.seed,polys[[i]])
	poly.seed
}

