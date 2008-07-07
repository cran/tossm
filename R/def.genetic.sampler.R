`def.genetic.sampler` <-
function(rland,sample.polys,bp.polys,n.samples,ID.tracker=NULL){
  n.bps <- length(bp.polys)
  n.ss <- length(sample.polys)
  if( length( n.samples) < n.ss)
    n.samples <- rep( n.samples, n.ss)[ 1:n.ss]
  abund.b <- tabulate(landscape.populations(rland))
  bp.density <- abund.b/sapply(bp.polys,function(x)area.poly(x))
  #calculate abundance in each sample site
  n.in.ss <- sapply(sample.polys, function(sp){
	sapply(1:length(bp.polys),function(bp) {
		floor(area.poly(intersect(sp,bp.polys[[bp]]))*bp.density[bp])
	})
  })
  if(n.bps ==1) n.in.ss <- matrix(n.in.ss,nrow=n.bps)
  n.samples <- pmin(n.samples, colSums(n.in.ss))
  pop <- cbind(ID=rland$individuals[,4],pop=landscape.populations(rland))

  loci <- rland$loci
  n.loci <- length( loci)
  ploidy <- sapply( loci, `[[`, 'ploidy')
  # Columns of first alleles:
  al1 <- ncol( rland$indiv) + 1 - rev( cumsum( rev( ploidy)))

  gs.yr <- vector( 'list', n.ss)
  ids <- vector('list', n.ss)
  for( i.ss in 1:n.ss) {
	gs.yr[[ i.ss]] <- array( 0/0, c( n.samples[ i.ss], n.loci, 2))
	if(n.samples[i.ss] > 0){
		n.per.bp <- floor((n.in.ss[,i.ss]/sum(n.in.ss[,i.ss]))*n.samples[i.ss])
		if (n.bps > 1) n.per.bp[1] <- n.samples[i.ss]-sum(n.per.bp[2:n.bps])

		#constrain available sample to indiv's not already sampled
		available.indiv <- pop[!(pop[,1] %in% ID.tracker),]
		which <- do.call("c", lapply(1:length(n.per.bp), function(x){
			rsample( n.per.bp[x], subset(available.indiv[,1],available.indiv[,2]==x), replace=FALSE)}))
	        #Compile allele info from rmetasim format...
		current.samp<-subset(rland$indiv,rland$indiv[,4] %in% which)
		
		coords <- do.call("rbind",lapply(1:length(bp.polys), function(i.bp){
			generate.coords(intersect(sample.polys[[i.ss]],bp.polys[[i.bp]]),nrow(current.samp))
		}))
		#and if the there is only one sample.....
		coords<-coords[1:nrow(current.samp),]
		ID.by<-matrix(current.samp[,c(4,3)],ncol=2)
		if(ncol(ID.by)==1){ID.by<-t(ID.by)}
		if(!is.matrix(coords)){coords<-t(as.matrix(coords))}
		
	      gs.yr[[ i.ss]][,,1] <- current.samp[, al1]
	      gs.yr[[ i.ss]][, ploidy==2,2] <- current.samp[, 1+al1[ ploidy==2]]
	      ids[[i.ss]] <- cbind(ID.by,coords)
		colnames(ids[[i.ss]]) <- c('IDs','birth.yr','x','y')
		ID.tracker <- c(ID.tracker,current.samp[,4])
	}
  }

  #attach loci info to gs.yr
  loci <- lapply( loci, function( x) { x$rate <- NA; x})
  attr( loci, 'printinfo') <- "<< locus info in rmetasim format; use 'unclass(...)' to display >>"
  class( loci) <- 'printinfo'
  attr( gs.yr, 'loci') <- loci

  #attach individual ids to gs.yr
  for (i.ss in 1:n.ss){
    if (!is.null(ids[[i.ss]])){
  	id.ss  <- ids[[i.ss]]
	attr(id.ss, 'printinfo') <- "<< individual ids, birth years, and coordinates for samples; use 'unclass(...)' to display >>"
  	class(id.ss) <- 'printinfo'
  	attr( gs.yr[[i.ss]], 'coords') <- id.ss
    }
  }

  return( gs.yr)
}

