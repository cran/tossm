`harvest` <-
function(mu.polys, interval.polys, bp.polys, rland, TAC){

   # Genetic data from harvested individuals is discarded.  However, this could
   # be easily changed in the future if we want to assume genetic data is collected
   if (sum(TAC)==0){ return (list(rland=rland, goners=NULL, TAC.b=TAC))
   } else {
	goners <- vector("list",length(mu.polys))
	bp.union <- poly.union(bp.polys)
	for (i.m in 1:length(mu.polys)){
	  if (TAC[i.m] > 0){
		remaining.TAC <- TAC[i.m]
		ID.tracker<-NULL
		for (int in 1:length(interval.polys)){
		  if (area.poly(poly.intersect(list(mu.polys[[i.m]],interval.polys[[int]],bp.union)))>0){
			
			gen.samp <- def.genetic.sampler(rland, list(poly.intersect(list(mu.polys[[i.m]],
				interval.polys[[int]],bp.union))),bp.polys,remaining.TAC,ID.tracker)[[1]]
			new.IDs <- do.call("c",lapply(gen.samp, function(gs){
				unclass(attr(gs,"coords"))[,1]}))
			ID.tracker<-unique(c(ID.tracker,new.IDs))
			
			#add gen.samp to goners[[i.m]] only if it contains some individuals
			if (nrow(gen.samp)>0) goners[[i.m]] <- rbind(goners[[i.m]],unclass(attr(gen.samp,"coords")))
			num.goners <- ifelse(is.null(goners[[i.m]]),0,nrow(goners[[i.m]]))
			if (num.goners == TAC[i.m]) { break
			} else {remaining.TAC <- TAC[i.m]-num.goners}
		  }
		}
	   }
	}
	goners <- do.call("rbind",lapply(goners, function(x) x))
	indices <- which(rland$individuals[,4] %in% goners[,1])
	rland$individuals <- subset(rland$individuals,!rland$individuals[,4] %in% goners[,1])
	return(list(rland=rland, goners=goners))
   }
}

