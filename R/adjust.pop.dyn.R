`adjust.pop.dyn` <-
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

