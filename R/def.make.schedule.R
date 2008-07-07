`def.make.schedule` <-
function( n.pre.RMP, n.RMP, n.post.RMP=0, abund.gap=5){
  gs.years <- n.pre.RMP # last year before RMP
  pre.RMP.years <- 1 %upto% n.pre.RMP
##KKM 4.10.07.............................................................
#  CLA.years <- abund.est.years <- n.pre.RMP + seq( 1, n.RMP, by=abund.gap)
  if (n.RMP ==0)   CLA.years <- abund.est.years <- n.pre.RMP + 1
  else CLA.years <- abund.est.years <- n.pre.RMP + seq( 1, n.RMP, by=abund.gap)
#.........................................................................
  post.RMP.years <- (n.pre.RMP + n.RMP) + 1 %upto% n.post.RMP
  BSA.years <- n.pre.RMP
  stop.years <- n.pre.RMP + n.RMP + n.post.RMP
returnList( stop.years, gs.years, abund.est.years, pre.RMP.years, CLA.years, post.RMP.years, BSA.years)
}

