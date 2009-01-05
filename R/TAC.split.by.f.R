'TAC.split.by.f' <-
function( TAC, fg.to.mu, abund.f){
  # TAC is per-mu
  # Take TAC from smallest-index FIMAs first
  TAC.f <- 0*fg.to.mu
  for( imu in 1:length( TAC)) {
    fgi <- index( fg.to.mu==imu)
    cumabund <- cumsum( c( 0, clip( abund.f[ fgi])))
    remainder <- TAC[ imu] - cumabund
    TAC.f[ fgi] <- pmax( 0, pmin( remainder, abund.f[ fgi]))
  }
return( TAC.f)

#
#  # Old version: split the per-mu TAC across those fg's that make up the mu
#  fg.per.mu <- tabulate( fg.to.mu)
#  TAC <- floor( TAC / fg.per.mu)
# return( TAC[ fg.to.mu])
}