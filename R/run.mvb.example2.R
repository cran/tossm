"run.mvb.example2" <-
function(...){
  l <- match.call( run.tossm)
  tms.call <- quote( run.tossm(
      rland=rlsimp4,
      choose.f=randomism,
      choose.f.args = list( matrix( c( 1, 1, 0, 0, 1, 1)/2, 2, 3, byrow=TRUE)),
      schedule=def.make.schedule( n.pre.RMP=5, n.RMP=6),
      n.gs.per.f=50,
      pre.RMP.removals=300,
      BSA=seq.hyptest.BSA,
      BSA.args=list( alpha=0.99999)))
  tms.call[ names( l)] <- l
  eval( tms.call)
}


