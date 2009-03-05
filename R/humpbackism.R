"humpbackism" <-
function( rland, FIMA.mig.mat, row.eqm.f.by.b){
  humpback.set.f <- function( rland){
      # Animals usually go to the same FIMA
      FIMA <- FIMA( rland)
      FIMA <- 1 + rowSums( runif( length( FIMA)) > CFMM[FIMA,])
      rland$mvb.other$FIMA <- FIMA
      rland
    }
  
  n.f <- ncol( FIMA.mig.mat)
  CFMM <- matrix( cumsum( t( FIMA.mig.mat)), byrow=T, n.f, n.f)
  CFMM <- CFMM - 0:(n.f-1)
  
  humpback.calf.f <- function( rland) {
      calves <- index( rland$indiv[,3] == rland$intparam$currentgen)
      ID <- rland$mvb.other$ID
      FIMA <- FIMA( rland)
      FIMA[ calves] <- FIMA[ match( ID[ calves,2], ID[,1])]
      rland$mvb.other$FIMA <- FIMA
      rland
    }
    
  cefb <- matrix( cumsum( t( row.eqm.f.by.b)), byrow=T, rland$intparam$habitats, n.f)
  cefb <- cefb - c( 0, clip( cumsum( rowSums( row.eqm.f.by.b))))
  FIMA <- 1 + rowSums( runif( nrow( rland$indiv)) > cefb[ landscape.populations( rland),])
  
returnList( FIMA, n.f, set.f=humpback.set.f, set.calf.f=humpback.calf.f)
}


