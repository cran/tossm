`old.bgd` <-
function( nlocal=sys.parent()) mlocal({
  # guessing what we need...
  # Can't see how to specify breeding-specific demographics in rmetasim;
  # ... seems to be same S, R, & M for all "habitats"

#  extract.named( rland$mvb.other) # unchanging stuff + ID
#  get.SRM( rland) # create mats & n.epoch
#  K <- rland$demog$epochs[[ n.epoch]]$Carry

  # 1st 3 cols of indiv are: stage, sex, generation of birth
  indiv <- rland$individuals
  pop <- 1+indiv[,1] %/% n.stages # populations( rland)
  pstage <- 1+indiv[,1]
  lstage <- (pstage-1) %% n.stages + 1
  n.whales <- length( pop)
  B <- tabulate( pop, nbins=n.bg)

  # Ageing & death-- happens first so calves can't go straight to stage 2
  # Pro tem, DD via juvenile mortality only
  survscale <- DD[1] + DD[2] * sqr(B/K) # DD should be in rland$mvb.other
  survscale <- pmin( survscale, 1/sum( S[,1]))

  tmat <- array( S, c( n.stages, n.stages, n.bg))
  tmat[,1,] <- survscale * tmat[,1,]
  dim( tmat) <- c( n.stages, n.stages * n.bg)
  ttmat <- matrix( cumsum( tmat), byrow=TRUE, n.stages * n.bg, n.stages)
  ttmat <- ttmat - c( 0, clip( cumsum( colSums( tmat))))
  new.lstage <- 1 + rowSums( runif( n.whales) > ttmat[pstage,])
  die <- new.lstage > n.stages

  indiv <- indiv[ !die,] # the rest die
  ID <- ID[ !die,]
  pop <- pop[ !die]
  n.whales <- length( pop)
  lstage <- new.lstage[ !die]
  pstage <- (pop-1) * n.stages + lstage
  indiv[,1] <- pstage-1

  # Migrate?
  pstage <- 1 + rowSums( runif( n.whales) > mm[pstage,])
  indiv[,1] <- pstage-1

  # For mammals, not sure why R & M are matrices (maybe for putting sex into stage). 
  # Also M looked mis-specified in example? Patched to what I think it should be in 'rlsimp2'

  # Breeding; first choose ma & pa
  # ma.prob is ABSOLUTE breeding chance
  ma.prob <- (indiv[,2]==1) * R[1, lstage]
  # pa.prob is RELATIVE breeding chance (cf other potential pas); any receptive ma is assumed to find a mate
  pa.prob <- (indiv[,2]==2) * M[1, lstage]
  ma <- pa <- vector( 'list', n.bg)
  for( i.bg in 1:n.bg) {
    ma[[i.bg]] <- index( runif( n.whales)+(pop != i.bg) < ma.prob) # fails if not this pop
    # NB one pa can father several babies-- replace=TRUE
    pa[[i.bg]] <- rsample( length( ma[[ i.bg]]), 1:n.whales, prob=pa.prob * (pop==i.bg), replace=TRUE)
  }

  ma <- unlist( ma)
  pa <- unlist( pa)

  # Make babies
  n.off <- length( ma)
  if( n.off) {
    indiv.ma <- indiv[ma,]
    indiv.pa <- indiv[pa,]
    off <- matrix( 0, n.off, ncol( indiv))
    off[,1] <- (pop[ ma]-1)*n.stages # all lstage 1
    off[,2] <- rsample( n.off, 1:2, replace=TRUE) # sex
    off[,3] <- rland$intparam$currentgen + 1 # duh...
    off[,ma.al1] <- ifelse( runif( n.off*n.ma.als)<0.5, indiv.ma[,ma.al1], indiv.ma[,ma.al2])
    off[,pa.al1] <- ifelse( runif( n.off*n.pa.als)<0.5, indiv.pa[,pa.al1], indiv.pa[,pa.al2])

    ID <- rbind( ID, matrix( c( 1:n.off+maxID, ID[ma,1], ID[pa,1]), ncol=3))
    # FG etc should be set externally after calling this
    maxID <- maxID + n.off
    indiv <- rbind( indiv, off)
  }

  rland$mvb.other$ID <- ID
  rland$mvb.other$maxID <- maxID
  rland$intparam$currentgen <- rland$intparam$currentgen + 1
  rland$individuals <- indiv
# return( rland)
})

