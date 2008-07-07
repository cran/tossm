`def.depleter` <-
function(unused1, unused2, bp.polys, rland, TAC){

	n.bps <- length(bp.polys)
	goners <- goner.locs <- vector("list",n.bps)
	for (i.b in 1:n.bps){
		if (TAC[i.b] > 0) goners[[i.b]] <- sample(which(landscape.populations(rland)==i.b),TAC[i.b],replace=FALSE)
		goner.locs[[i.b]] <- rland$individuals[goners[[i.b]],4]
		goner.locs[[i.b]] <- cbind(goner.locs[[i.b]],generate.coords(bp.polys[[i.b]],TAC[i.b]))
	}
	goners <- do.call("c",lapply(goners, function(x) x))
	goner.locs <- do.call("rbind",lapply(goner.locs, function(x) x))
	if (length(goners)){
		rland$individuals <- rland$individuals[-goners,]
	}
	return(list(rland=rland,goners=goner.locs))
}

