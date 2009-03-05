"agg.gs.tseries" <-
function( nlocal=sys.parent()) mlocal({
  ysamp <- index( sapply( gs, length)>0)
    
  if( length( ysamp)) {
    has.coords <- !is.null( attr( gs[[ ysamp[1]]][[1]], 'coords'))
    n.areas <- length( gs[[ ysamp[1]]])

    agg.gs <- vector( 'list', n.areas)
    n.loci <- ncol( gs[[ ysamp[ 1]]][[1]])
    n.alleles <- rep( 0, n.loci)

    # Aggregate all historical samples together, by area
    # Also get number of alleles per loci (safety)
    for( i.area in 1:n.areas) {
      ai <- do.call( 'abind.mvb', c( list( dim.=1), lapply( gs[ ysamp], `[[`, i.area)))
      ai[ is.na( ai)] <- 0 # !!! skipped by tabulate so OK in this case
      n.alleles <- pmax( n.alleles, apply( ai, 2, max))
      
      if( has.coords) {
        xyi <- do.call( 'rbind', lapply( gs[ ysamp], function( x) attr( x[[ i.area]], 'coords')))
	  attr(xyi,'printinfo') <- "<< individual ids and spatial coordinates for samples; use 'unclass(...)' to display >>"
	  class(xyi) <- 'printinfo'
        attr( ai, 'coords') <- xyi
      }

      agg.gs[[ i.area]] <- ai
    }
    attr(agg.gs,'seq.list') <- attr(gs,'seq.list')    

  } else
    agg.gs <- NULL
 
  agg.gtypes <- agg.gs.to.gtypes(agg.gs,gs) 
  agg.gs
})


"agg.gs.series" <-
function( ts.gs){
  ysamp <- index( sapply( gs, length)>0)
  n.areas <- length( gs[[ ysamp[1]]])

  agg.gs <- vector( 'list', n.areas)
  n.loci <- ncol( gs[[ ysamp[ 1]]][[1]])
  n.alleles <- rep( 0, n.loci)

  # Aggregate all historical samples together, by area
  # Also get number of alleles per loci (safety)
  for( i.area in 1:n.areas) {
    # +1 because rmetasim's allele codes start at 0 and tabulate can't cope
    ai <- 1 + do.call( 'abind.mvb', c( list( dim.=1), lapply( gs[ ysamp], `[[`, i.area)))
    ai[ is.na( ai)] <- 0 # !!! skipped by tabulate so OK in this case
    n.alleles <- pmax( n.alleles, apply( ai, 2, max))
    agg.gs[[ i.area]] <- ai
  }
 
  agg.gs
}

