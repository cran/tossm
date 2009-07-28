"get.n.in.ss" <-
function(rland,sample.polys,bp.polys){
#calculates the number of individuals from each bp in each sampling poly or
#harvest interval
  n.bps <- length(bp.polys)
  n.ss <- length(sample.polys)
  abund.b <- tabulate(landscape.populations(rland))
  bp.density <- abund.b/sapply(bp.polys,function(x)area.poly(x))
  #calculate abundance in each sample site
  n.in.ss <- sapply(sample.polys, function(sp){
	sapply(1:length(bp.polys),function(bp) {
		area.poly(intersect(sp,bp.polys[[bp]]))*bp.density[bp]
	})
  })
  if(n.bps ==1) n.in.ss <- matrix(n.in.ss,nrow=n.bps)
  extras <- rowSums(n.in.ss) - rowSums(floor(n.in.ss))
#6.25.09 temporarily commenting out the following for loop to re-introduce
#the bug where individuals are neglected due to flooring the number of 
#individuals in each harvest interval.  We want to see if re-introducing
#this bugs results in some results matching what we got last year
#  for (i.bp in 1:n.bps){
#	i.ss <- 1
#	while(extras[i.bp]>=1){
#		if (n.in.ss[i.bp,i.ss] > 0){
#			n.in.ss[i.bp,i.ss] <- n.in.ss[i.bp,i.ss]+1
#			extras[i.bp] <- extras[i.bp]-1
#		}
#		i.ss <- i.ss+1
#	}
#  }
  n.in.ss <- floor(n.in.ss)
  return(n.in.ss)
}