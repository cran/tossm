"CLA" <-
function(Catches,Abund.est,CV.est, multiplier=5, ...) {
 # Keep CLA happy...
 yrs.of.est <- index( !is.na( rowSums( Abund.est))) # won't work if partial ests
 n.mu <- ncol( Catches)
 n.y <- nrow( Catches)
 TAC <- numeric( n.mu)

 for( i.mu in 1:n.mu) 
##KKM 4.10.07..........................................................
#First line is a cheaty shortcut that sets the quota at 4% of abundance
#Second through fourth lines actually use the CLA to set quota
#......................................................................
#  TAC[ i.mu] <- round( Abund.est[ max( yrs.of.est), i.mu] * 0.04)
   TAC[ i.mu] <- round( multiplier *
      raw.CLA( Catches[ ,i.mu], Abund.est[ yrs.of.est,i.mu], CV.est[ yrs.of.est,i.mu],
       1000+yrs.of.est, 1000+1, 1000+n.y, ...))
return( TAC)
}


