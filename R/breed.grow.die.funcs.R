"adjust.pop.dyn" <-
function( rland, MSYR, MSYL, alpha=NULL, tweak=TRUE){
  get.SRM( rland)
  R <- R/2 # females!
  diag( R) <- 1
  
  funco <- function( alpha) {
      S[,1] <- alpha * S[,1]
    return( max( Re( eigen( R %*% S)$values)))
    }
    
  if( !is.null( alpha))
return( funco( alpha))
  
  repeat{ # probably won't
    scale.at.K <- find.root( funco, 1, 0.1, 'increasing', target=1)
    if( scale.at.K * sum( S[,1]) <= 1)
  break # OK
    if( !tweak) # else...
stop( "Extinction certain even with juv surv = 1")

    cat( 'Halving mortality to achieve equilibrium...\n')
    # Tweak: halve the mortality rate in each column
    surv <- colSums( S)
    S <- S * rep( (1+surv)/(2*surv), each=nrow( S)) 
    rland$demography$localdem[[ rland$intparam$currentepoch+1]]$LocalS <- S
  }

  scale.at.MSYL <- find.root( funco, scale.at.K, 0.1, 'increasing', target=1+MSYR)
  cfm <- matrix( c( 1, 1, 1, sqr( MSYL)), 2, 2, byrow=TRUE) # this doesn't actually ensure MSY is at MSYL...
  rland$mvb.other$DD <- solve( cfm, c( scale.at.K, scale.at.MSYL))
return( rland)
}


"as.smat" <-
function( mat) {
  if( is.list( mat))
return( mat)
  o <- order( col( mat), mat)
  mat[] <- cumsum( mat[o]) - col( mat) + 1
return( list( lookup=t( mat), unsort=o))
}


"breed.grow.die" <-
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
  w0 <- pmax( 0, 1-B/K) # This will cause problems unless top eval is *exactly* 1 at K!

  # Ageing & death-- happens first so calves can't go straight to stage 2
  new.lstage <- 1 +
      rowSums( runif( n.whales) > w0[ pop] * ttmat0[pstage,] + (1-w0)[pop] * ttmatK[pstage,])
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

  # Breeding; first choose ma & pa
  # ma.prob is ABSOLUTE breeding chance
  ma.prob <- w0[ pop]*R0[1, lstage] + (1-w0[ pop])*RK[1,lstage]
  # pa.prob is RELATIVE breeding chance (cf other potential pas); any receptive ma is assumed to find a mate
  # Don't understand why M is a matrix rather than a vector; assume last col is the important one
  pa.prob <- M0[ lstage, n.stages] # assume MK is the same
  ma <- pa <- vector( 'list', n.bg)
  for( i.bg in 1:n.bg) {
    ma[[i.bg]] <- index( runif( n.whales)+(pop != i.bg) < ma.prob) # fails if not this pop
    # NB one pa can father several babies-- replace=TRUE
    # Will crash if NO pas available!
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
    # off[,2] <- rsample( n.off, 1:2, replace=TRUE) # sex-- NOT USED
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


"get.SRM" <-
function( rland, which.dem, postfix, nlocal=sys.parent(), ldemog) mlocal({
  ldemog <- rland$demog[[ 'localdem' %&% which.dem]][[ n.epoch]]
  names( ldemog) <- sub( 'Local', '', names( ldemog)) %&% postfix
  extract.named( ldemog) # S0, R0, M0 or SK, RK, MK
})


"km.simple.demo" <-
function( nloc=5) {
  #gray whale demography NO DENSITY DEPENDENCE


  habitats <- 2   # 1 for Archetype I or 2 for Archetype II
  carrycap <- 7500    #7500, 15000, or 30000
  stages <- 5
  rland <- NULL
  stagedist<-c(0.2937,0.136,0.2045,0.0614,0.3037) #stable stage distribution for K matrix
  runtime<-1000
  dr<-0.005

  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland, h = habitats, s = stages, totgen = runtime)
  rland <- landscape.new.switchparam(rland, mp = 1)
  rland <- landscape.new.floatparam(rland)

  #these are the S and R matrices at K, updated to increase lambda a bit as of 12/22/05
  SK <- matrix (c(0.768-dr, 0, 0, 0, 0,
           0.157, 0.720-dr, 0, 0, 0,
           0, 0.102, 0.648-dr, 0.946, 0,
           0, 0, 0.300-dr, 0, 0,
           0, 0.102, 0, 0, 0.954-dr), nrow=5, byrow = TRUE)
  RK <- matrix (c(0, 0, 0, 0.925, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0), nrow = 5, byrow = TRUE)

  MK <- matrix (c(0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0, 0, 0, 0, 1,
           0, 0, 0, 0, 0), nrow = 5, byrow = TRUE)

  rland<-new.local.demo(rland,SK,RK,MK)


  blockA<-matrix(0,5,5)
  blockB<-diag( rep( dr, 5)) # matrix(c(dr,0,0,0,0, 0,dr,0,0,0, 0,0,dr,0,0, 0,0,0,dr,0, 0,0,0,0,dr),nrow=5)
  rowA<-cbind(blockA,blockB)
  rowB<-cbind(blockB,blockA)
  Sdisp<-rbind(rowA,rowB)

  rland<-new.epoch(rland,S=Sdisp,R=matrix(c(rep(0,100)),nrow=10),M=matrix(c(rep(0,100)),nrow=10),
    carry=c(rep((carrycap/habitats),habitats)))   ### divide by habitats

  for( i.loc in 1 %upto% nloc)
    rland <- new.locus(rland, type = 1, ploidy = 2, transmission = 0, numalleles = 10,
        mutationrate = 0.002)

  #for 50:50

  rland <- landscape.new.individuals(rland,rep( stagedist*carrycap*0.5,2))
  #rland<-new.individuals(rland,c(50,50,50,50,50,0,0,0,0,0))
  #rland<-new.individuals(rland,c(rep(0,5),rep(50,5)))
  #rland<-simulate.landscape(rland,numit=1)
return( rland)
}


"old.bgd" <-
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


"prepare.rland" <-
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


