`get.SRM` <-
function( rland, which.dem, postfix, nlocal=sys.parent(), ldemog) mlocal({
  ldemog <- rland$demog[[ 'localdem' %&% which.dem]][[ n.epoch]]
  names( ldemog) <- sub( 'Local', '', names( ldemog)) %&% postfix
  extract.named( ldemog) # S0, R0, M0 or SK, RK, MK
})

