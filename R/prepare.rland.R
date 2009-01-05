`prepare.rland` <-
function( rland){
  extract.named( rland$mvb.other) # in case it exists
  n.bg <- rland$intparam$habitats
  n.stages <- rland$intparam$stages
  n.loci <- length( rland$loci)
  ploidy <- sapply( rland$loci, `[[`, 'ploidy')
  # Assume haploid=>maternal only; can't see how to specify Y-chromosome in rmetasim
  ma.al1 <- pa.al1 <- 4 + clip( cumsum( c( 0, ploidy)))
  ma.al2 <- ma.al1
  ma.al2[ ploidy==2] <- ma.al2[ ploidy==2]+1
  pa.al1 <- pa.al1[ ploidy==2]
  pa.al2 <- pa.al1+1 # only diploids have any pa al
  n.ma.als <- length( ma.al1)
  n.pa.als <- length( pa.al1)

  # Next bit readies the migration matrix for fast transitions
  # NB epoch-dependent so shouldn't be here-- should really be in breed.grow.die
  migmat <- rland$demog$epochs[[ 1+rland$intparam$currentepoch]]$S
  
  # Edict: migmat shall have non-zero elements only in the diagonals of the off-diagonal submatrices
  interesting <- abs( row( migmat) - col( migmat)) %% n.stages == 0
  migmat[ !interesting] <- 0
  diag( migmat) <- 0
  diag( migmat) <- 1 - colSums( migmat)
  mm <- fast.trans( migmat, TRUE)
  mm <- matrix( cumsum( migmat), byrow=TRUE, n.stages * n.bg, n.stages * n.bg)
  mm <- mm - 0:(n.stages*n.bg-1)

  # Ready for growth & death
  n.epoch <- 1+rland$intparam$currentepoch
  get.SRM( rland, '', '0') # create mats & n.epoch at low abund
  get.SRM( rland, 'K', 'K') # create mats at carry cap
  ttmat0 <- fast.trans( matrix( S0, n.stages, n.stages*n.bg))
  ttmatK <- fast.trans( matrix( SK, n.stages, n.stages*n.bg))

  maxID <- nrow( rland$indiv)
  ID <- matrix( c( 1:maxID, rep( NA, 2*maxID)), ncol=3, dimnames=list( NULL, cq( ID, ma, pa)))

  rland$mvb.other <- lapply( named( ls( sys.frame( sys.nframe())) %except% cq( fast.trans, rland)),
      get, envir=sys.frame( sys.nframe()))

return( rland)
}

