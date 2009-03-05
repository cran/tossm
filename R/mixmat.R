"mixmat" <-
function( rland, mix.mat){
  # Each animal selects FIMA anew each year
  # Prob of FIMA depends on BG
  random.set.f <- function( rland){
      bg <- landscape.populations( rland) # breeding group
      n.b <- nrow( mix.mat)
      n.f <- ncol( mix.mat)
      n.by.pop <- tabulate( bg, n.b)
      FIMA <- integer( sum( n.by.pop))
      for( i.b in 1:n.b)
        FIMA[ bg==i.b] <- rsample( n.by.pop[ i.b], 1:n.f, 
            prob=mix.mat[i.b,], replace=TRUE)
      rland$mvb.other$FIMA <- FIMA
      rland
    }
  mixmat.calf.f <- function( rland) rland # irrel because reset each year
  
returnList( FIMA=random.set.f( rland), n.f=ncol( mix.mat), 
    set.f=random.set.f, set.calf.f=mixmat.calf.f)
}


