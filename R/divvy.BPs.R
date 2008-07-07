`divvy.BPs` <-
function(bp.polys, other.polys){
	map <- t(sapply(bp.polys, function(x){
		bp.area <- area.poly(x)
		sapply(other.polys, function(y){
			area.poly(intersect(x,y))/bp.area
		})
	}))
	if(length(other.polys)==1) map <- matrix(map,ncol=1)
	return(map)
}

