
hyptest.network.BSA<-function(gs,abund,var,C,landscape.poly,sample.polys,sig.level){

	#get gs into the format needed for the test.g function of hierfstat
	agg.gs<-agg.gs.tseries()
 	concat.loci<-list()
	for (i.site in 1:length(agg.gtypes)){
		concat.loci.mat<-matrix(NA,nrow(agg.gtypes[[i.site]]),(ncol(agg.gtypes[[i.site]])-1)/2)
			for (i.col in 1:((ncol(agg.gtypes[[i.site]])-1)/2)){
				concat.loci.mat[,i.col]<-
					as.numeric(paste(agg.gtypes[[i.site]][,i.col*2],agg.gtypes[[i.site]][,(i.col*2)+1],sep=""))
			}
		concat.loci[[i.site]]<-concat.loci.mat
	}

	concat.loci

	#Perform pairwise g tests between all sampling sites
	#output (sigmat.bin) a square matrix of 0's (signficant tests) and 1's (non-significant tests)

	n.levels<-length(concat.loci)
	n.individuals<-nrow(concat.loci[[1]])
	sigmat<-matrix(NA,n.levels,n.levels)
	sigmat.bin<-matrix(NA,n.levels,n.levels)
	tester<-2
	if (tester==1){
		for (i.lev in 1:n.levels){
			for (j.lev in i.lev:n.levels){
				sigmat.bin[i.lev,j.lev]<-rbinom(1,1,.3)
				sigmat.bin[j.lev,i.lev]<-sigmat.bin[i.lev,j.lev]
			}
		}
		diag(sigmat.bin)<-1
	} else {
		for (i.lev in 1:(n.levels-1)){
			for (j.lev in (i.lev+1):n.levels){
				test.data<-rbind(concat.loci[[i.lev]],concat.loci[[j.lev]]);
				sampleID.vec<-c(rep(i.lev,nrow(concat.loci[[i.lev]])),rep(j.lev,nrow(concat.loci[[j.lev]])));
				test.result<-test.g(test.data,sampleID.vec)$p.val;
				sigmat[i.lev,j.lev]<-sigmat[j.lev,i.lev]<-test.result;
			}
		}
		sigmat.bin<-ifelse(sigmat>sig.level,1,0)
		diag(sigmat.bin) <- 1
	}

	#create a list of length=n.sampling sites, that contains, for each sampling site, the 
	#other sampling sites connected to it by non-significant g tests.
	cluster.seeds<-list()
	for (i in 1:nrow(sigmat.bin)){
		cluster.seeds[[i]]<-which(sigmat.bin[,i]==1) #find row where this is the case
	}

	#for each of these starting 'seeds', search all other list elements for shared sampling sites.
	#if a shared sampling site is found, move all sampling sites in the matching element to the first element 
	#Continue this until all sampling sites have been binned into their MUs...

	seed.list<-cluster.seeds

	update.made <- TRUE
	while (update.made) {
		update.made <- FALSE
		for (i.seed in 1:length(seed.list)) {
			current.vec <- seed.list[[i.seed]]
			for (j.seed in 1:length(seed.list[[i.seed]])) {
				current.id <- seed.list[[i.seed]][j.seed]
				if (current.id != i.seed) {
					id.vec <- seed.list[[current.id]]
					current.vec <- union(id.vec, current.vec)
				}
			}
			if (!all(current.vec %in% seed.list[[i.seed]])) {
				seed.list[[i.seed]] <- unique(current.vec)
				update.made <- TRUE
			}
		}
	}

	seed.list <- lapply(seed.list, sort)
	id.list <- lapply(seed.list, function(i) min(i))
	id.vec <- unlist(id.list)
	id.vec <- as.numeric(factor(id.vec))
	#create a list of length=number of mus, with each element representing
	#the sample.polys contained in each mu
	MU.from.BSA<-list()
	for (i.ID in 1:length(unique(id.vec))){
		MU.from.BSA[[i.ID]]<-seed.list[[min(which(id.vec==i.ID))]]
	}
 	save(MU.from.BSA,file="MU.from.BSA.Rda")
	munits<-MU.poly.generator(MU.from.BSA,sample.polys,landscape.poly)
	return(munits)
 
}
