`km.simple.demo` <-
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

