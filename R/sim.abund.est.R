`sim.abund.est` <-
function( abund.b, abund.for.10pc.CV,bp.in.MU.map){
  truth <- colSums(bp.in.MU.map*abund.b)
  shape <- sqr(10)*truth/abund.for.10pc.CV
returnList( est=rgamma( length( truth), shape=shape, scale=truth/shape), CV=1/sqrt( shape))
}

