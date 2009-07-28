"harvest" <- function(mu.polys, interval.polys, bp.polys, rland, TAC){


   # Genetic data from harvested individuals is discarded.  However, this could
   # be easily changed in the future if we want to assume genetic data is collected
   if (sum(TAC)==0){ return (list(rland=rland, goners=NULL, actual.catch=TAC))
   } else {
	goners <- vector("list",length(mu.polys))
	bp.union <- poly.union(bp.polys)
	for (i.m in 1:length(mu.polys)){
	  if (TAC[i.m] > 0){
		remaining.TAC <- TAC[i.m]
		ID.tracker<-NULL
		mu.interval.union <- lapply(interval.polys,function(int) {
			poly.p <- poly.intersect(list(mu.polys[[i.m]],int,bp.union))
			if (area.poly(poly.p)>0) return(poly.p)
			else return(NULL)
		})
		for (int in length(mu.interval.union):1){
			#remove NULL elements from the list
			if (is.null(mu.interval.union[[int]])) mu.interval.union[[int]] <- NULL
		}
		n.in.int <- get.n.in.ss(rland,mu.interval.union,bp.polys)
		for (int in 1:length(mu.interval.union)){
		  	gen.samp <- def.genetic.sampler(rland,mu.interval.union[int],bp.polys,remaining.TAC,n.in.ss=matrix(n.in.int[,int],nrow=length(bp.polys)),ID.tracker)[[1]]
			new.IDs <- unclass(attr(gen.samp,"coords"))[,1]
			ID.tracker<-unique(c(ID.tracker,new.IDs))
			
			#add gen.samp to goners[[i.m]] only if it contains some individuals
			if (nrow(gen.samp)>0) goners[[i.m]] <- rbind(goners[[i.m]],unclass(attr(gen.samp,"coords")))
			num.goners <- ifelse(is.null(goners[[i.m]]),0,nrow(goners[[i.m]]))
			if (num.goners == TAC[i.m]) { break
			} else {remaining.TAC <- TAC[i.m]-num.goners}
		}
	   }
	}
	actual.catch <- sapply(goners,function(x){ return(if(is.null(x)){0} else { dim(x)[1]})})
	goners <- do.call("rbind",lapply(goners, function(x) x))
	rland$individuals <- subset(rland$individuals,!rland$individuals[,4] %in% goners[,1])
	return(list(rland=rland, goners=goners, actual.catch=actual.catch))
   }
}
