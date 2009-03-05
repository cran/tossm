"count.alleles.by.f" <-
function( nlocal=sys.parent()) mlocal({
  # Count alleles by area. mtDNA has 1 allele, microsats have 2
  # Put results into an array, even though some loci have more alleles than others
  if( is.null( agg.gs))
    agg.gfreq <- NULL
  else {
    maxal <- max( n.alleles)
    agg.gfreq <- array( 0, c( n.areas, n.loci, maxal))
    for( i.area in 1:n.areas)
      agg.gfreq[ i.area,, ] <- t(apply( agg.gs[[ i.area]], 2, tabulate, nbins=maxal))
  }
})


