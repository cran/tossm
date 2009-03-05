"hyptest.seq.BSA" <-
function( gs, ea, vea, C, alpha){
  # Example boundary-setting method

  # Next two lines will be standard in any frequency-based analysis
  agg.gs.tseries() # creates agg.gs, n.loci, n.areas, n.alleles
  count.alleles.by.f() # creates maxal, agg.gfreq

  # Sequential hyp tests for 1-D case, loci assumed independent. Not supposed to be smart!
  # Start by finding "best" single boundary; split there if H0 fails, then look for
  # "best" boundaries within each partition, etc.
  rbdy <- 1
  stack <- matrix( c( 1, n.areas), 1, 2) # each row (left,right) pair to check for bdies between
  agg <- array( 0, c( n.loci, maxal, 2)) # will be re-used
  hyp.score <- numeric( n.areas)
  while( length( stack)) {
    left <- stack[1,1]
    right <- stack[1,2]
    stack <- stack[-1,]

    # Do hyp tests for all poss bdies between left and right
    hyp.score[] <- Inf
    for( jr in (1+left) %upto% right) { # group into left-of-j vs j-and-right
      agg[,,1] <- colSums( agg.gfreq[ left %upto% (jr-1),,])
      agg[,,2] <- colSums( agg.gfreq[ jr %upto% right,,])
      # Do test per locus
      pvals <- apply( agg, 1, function( x) chisq.test( x[ rowSums( x)>0, ], simulate=T, B=1000)$p.value)
      # Fisher's rule for combining independent tests-- may have screwed up!
      hyp.score[ jr] <- pchisq( -2 * sum( log( 1-pvals)), df=n.loci) # give or take minus signs...
    }

    best.jr <- index( hyp.score==min( hyp.score))[1]
    if( hyp.score[ best.jr] < alpha) { # use it
      rbdy <- c( rbdy, best.jr)
      if( best.jr > 1+left) # else only one block
        stack <- rbind( stack, c( left, best.jr-1))
      if( best.jr < right) # else only one block
        stack <- rbind( stack, c( best.jr, right))
    }
  }

  # rbdy is LH ends of units so scrunge this into format required
  group <- colSums( outer( rbdy, 1:n.areas, `<=`))
return( split( 1:n.areas, group))
}


