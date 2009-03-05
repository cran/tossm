"run.tossm" <-
function( rland=NULL, bp.polys, schedule=NULL,
    n.samples=NULL, sample.polys, initial.depletion=.30, historic.removals=NULL,
    historic.polys=NULL, BSA=fixed.MU.BSA, BSA.args=list(n.mus=1),
    stop.years=schedule$stop.years, gs.years=schedule$gs.years,
    abund.est.years=schedule$abund.est.years, pre.RMP.years=schedule$pre.RMP.years,
    CLA.years=schedule$CLA.years, post.RMP.years=schedule$post.RMP.years,
    BSA.years=schedule$BSA.years, harvest.interval,
    genetic.sampler=def.genetic.sampler, 
    abund.for.10pc.CV=70000, quota.calc=PBR,
    quota.args=list(r.max=0.04,F.r=1,multiplier=1), 
    CLA.prog=NULL, CLA.dir=NULL, plot.polys=FALSE,
    seed=-1) {

  if(plot.polys) tossm.plotter(bp.polys,sample.polys,historic.polys)
  check <- tossm.diagnostics(bp.polys,sample.polys,historic.polys,rland)
  if (check>0){
	if (!plot.polys) tossm.plotter(bp.polys,sample.polys,historic.polys)
 	stop("The TOSSM simulation cannot run for the above reason(s).")
  } else print("Diagnostics OK...running valid TOSSM simulation...")

  library( mvbutils)
  n.y <- min( stop.years)

  if( is.null( CLA.prog) && any( CLA.years <= n.y)) {
    i.am.tossm <- max( match( cq( 'package:tossm', 'tossm'), search(), 0))
    tossm.dir <- attr( as.environment( i.am.tossm), 'path')

    if( is.null( CLA.dir)) {
      CLA.dir <- file.path( getwd(), 'CLA.temp')
      mkdir( CLA.dir)
#      CLA.dir <- tempdir()
#      on.exit( unlink( CLA.dir))
    }

    cat( CLC.PAR, sep='\n', file=file.path( CLA.dir, 'CLC.PAR'))
  }

  # Set up RNG and get rid of scary irrelevant warning about Walker's alias...
  if( seed>0) {
    seed <- round( seed[1])
    set.seed( seed %% 1024)
  }
  ow <- options( warn=-1)
  sample( 1000+(abs( seed) %/% 1024) %% 1024, 1:1000, replace=TRUE)
  options( ow)

  n.b <- rland$intparam$habitats
  n.epoch <- 1+rland$intparam$currentepoch
  rland$intparam$totalgens <- rland$intparam$currentgen + (n.y+1)

  #generate harvest intervals
  landscape.poly <- get.bbox(poly.union(bp.polys))
  interval.polys <- set.interval.polys(landscape.poly, interval=harvest.interval)

  #set historic.removals based on initial depletion.  Take all historic.removals
  #in the last pre.RMP.year so as to achieve the specified depletions
  if (is.null(historic.removals)){
	if (length(initial.depletion)<n.b){initial.depletion <- 
		c(initial.depletion,rep(.99,n.b-length(initial.depletion)))}
	historic.removals <- round((1-initial.depletion)*rland$demography$epochs[[n.epoch]]$Carry)
	if (length(pre.RMP.years)>1) historic.removals <- c(rep(0,(n.b*(length(pre.RMP.years)-1))),historic.removals)
	historic.removals <- matrix(historic.removals,nrow=length(pre.RMP.years),byrow=TRUE)
	historic.polys <- bp.polys
	depleter <- def.depleter
  } else {depleter <- harvest}

  #catches holds the number of animals killed in each MU in each year
  #catch.locs holds the coords of each animal killed each year
  catches <- matrix(0,n.y,length(historic.polys))
  catch.locs <- vector("list",n.y)
  true.abund.mu <- est.abund.mu <- var.abund.mu <- mu.hist <- catch.locs
  effort <- matrix( 0, n.y) #this needs work...
  abund.b <- matrix( 0, n.y, n.b)
  gs <- vector( 'list', n.y)
  attr(gs,'seq.list') <- list()
  munits <- fixed.MU.BSA(0,0,0,0,landscape.poly,0,n.mus=1)
  bp.in.MU.map <- matrix(rep(1,n.b),ncol=1)

  for( i.y in 1 %upto% n.y) {
    save.rs <- .GlobalEnv$.Random.seed # in case sampling disrupts RNG sequence
    cat( 'Year', i.y)

    # Get truth
    abund.b[ i.y,] <- tabulate( landscape.populations( rland), n.b)

    cat( '  ', abund.b[ i.y,], '  ')
    if (min(abund.b[ i.y,]) < 2) break
    if (sum(abund.b[ i.y,]) < 2) break

    # Sample data?
    if( i.y %in% gs.years) {
	gs[[i.y]] <- genetic.sampler(rland,sample.polys,bp.polys,n.samples)
	gs <- index.to.state(gs,i.y)
      .GlobalEnv$.Random.seed <- save.rs # if run 2X with diff # gs then...
    }

    # Implementation/review? (allows full-feedback testing later)
    if( i.y %in% BSA.years) {
	est.abund.hist <- do.call("rbind",lapply(1:i.y,function(x){ 
		if (is.null(est.abund.mu[[x]])) {rep(0/0,length(catches[1,]))
		} else {est.abund.mu[[x]]}}))
	var.abund.hist <- do.call("rbind",lapply(1:i.y,function(x){
		if (is.null(var.abund.mu[[x]])) {rep(0/0,length(catches[1,]))
		} else {var.abund.mu[[x]]}}))
      BSA.munits <- do.call( 'BSA', c( list(
          gs[ 1:i.y], est.abund.hist,var.abund.hist, 
	    catches[ 1:i.y,,drop=FALSE],landscape.poly),list(sample.polys),BSA.args))
#      if( !all( sort( unlist( munits))==1:n.f))
#stop( 'FIMAs are wrongly mapped to MUs')

      .GlobalEnv$.Random.seed <- save.rs
      # munits will be non-overlapping polygons that encompass the entire landscape

	new.munits <- TRUE
    }

    #Estimate abundance?
    if( i.y %in% abund.est.years) { # ... want exact same abund ests
      simbo <- sim.abund.est( abund.b[ i.y,], abund.for.10pc.CV=abund.for.10pc.CV,
		bp.in.MU.map)
      est.abund.mu[[i.y]] <- simbo$est
      var.abund.mu[[i.y]] <- evalq( sqr( CV * est), simbo)
      .GlobalEnv$.Random.seed <- save.rs # if run 2X with diff # gs then...
    }

    # TAC will now be set by mgmt stock...
    # Have to take 1 animal to avoid CLA trouble...
    if( i.y %in% pre.RMP.years) {
      historic.removals[i.y, historic.removals[i.y,]<=0] <- 1
      TAC <- historic.removals[i.y,]
	h.res <- depleter(historic.polys, historic.polys, bp.polys, rland, TAC)
	mu.hist[[i.y]]$mus <- historic.polys
	mu.hist[[i.y]]$catch <- TAC
    } else if( i.y %in% post.RMP.years) {
      TAC <- rep( 0, dim(catches)[2])      
	h.res <- harvest(historic.polys, historic.polys, bp.polys, rland, TAC)
    } else { # in RMP phase
      if( i.y %in% CLA.years) {
        # Group catch & abund data by mgmt stock if mu defs have changed since last CLA.year
        # Assumes that ests by FIMAs are independent-- fix later?
	  if (new.munits){
		munits <- BSA.munits
      	# update the fraction of each BP that is in each MU
		bp.in.MU.map <- divvy.BPs(bp.polys,munits)

		# Calculate historic catch in each MU
		catches <- do.call("rbind",lapply(catch.locs, function(locs){
			if(is.null(locs)) {return(rep(0,length(munits)))
			} else {
				divvy.catch(locs,munits)
			}}))

		# Calculate historic abundance in each MU.  Note that this approach produces
		# estimates that are completely independent of the estimates based on
		# previous stratifications...
		for (yr in 1:i.y){
			if (!is.null(est.abund.mu[[yr]])){
				recalc <- sim.abund.est(abund.b[yr,], abund.for.10pc.CV,
					bp.in.MU.map)
				est.abund.mu[[yr]] <- recalc$est
				var.abund.mu[[yr]] <- evalq(sqr(CV * est), recalc)
			}
		}
		new.munits <- FALSE
	  }
	  est.abund.hist <- do.call("rbind",lapply(1:i.y,function(x){ 
		if (is.null(est.abund.mu[[x]])) {rep(0/0,length(munits))
		} else {est.abund.mu[[x]]}}))
	  var.abund.hist <- do.call("rbind",lapply(1:i.y,function(x){
		if (is.null(var.abund.mu[[x]])) {rep(0/0,length(munits))
		} else {var.abund.mu[[x]]}}))
        CV.abund.hist <- sqrt( var.abund.hist) / est.abund.hist

        TAC <- do.call( 'quota.calc', c( list( as.matrix(catches[1:i.y,]),
		est.abund.hist, CV.abund.hist,
            Directory=CLA.dir), quota.args))
	}
      
	mu.hist[[i.y]]$mus <- munits
	mu.hist[[i.y]]$catch <- TAC
      h.res <- harvest(munits, interval.polys, bp.polys, rland, TAC)
    }

    rland <- h.res$rland
    catches[ i.y,] <- TAC
    if(!is.null(h.res$goners)) catch.locs[[i.y]] <- h.res$goners
    effort[ i.y,] <- calc.effort( catch.locs,i.y)

    if (nrow(rland$individuals) == 0) break

    # Project population-- finally!
    cat( '... projecting...'); rland <- landscape.simulate( rland, 1) # 1 year
    cat( '\n')

  } # for i.y

  on.exit()
  agg.gs.tseries()
  agg.gs.to.gtypes(agg.gs,gs)
#  count.alleles.by.f()
  callo <- sys.call()
  result <- returnList( abund.b, catches, effort,
      est.abund.mu, var.abund.mu, mu.hist,
	gs, agg.gs,agg.gtypes,call=callo, seed,rland)
  result[ 1:3] <- lapply( result[1:3], t) # save space printing-- years as cols
  class( result) <- 'tossm.obj'
return( result)
}


