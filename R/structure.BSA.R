structure.BSA <- function(gs,abund,var,C,landscape.poly,sample.polys,mainparams,BSA.params){

#	mainparams is a list containing the following components:
#	$numinds
#	$numloci
#	$ploidy
#	$noadmix
#	$freqscorr
#	$maxpops
#	$burnin
#	$iterations

#	BSA.params is a list containing the following components:
#	$output.dir - where should the STRUCTURE results be written
#	$k - a list of values of k to evaluate
#	$struct.exe.path - directory where structure.exe resides and where mainparams
#		will be written.  Defaults to "C:/Program Files/Structure console"
#	$struct.outfile - name of the output file to be written by STRUCTURE.
#		Usually 'ArchX_scY_Z.nxxx'.  Extension '.kx_f' will be tacked on
#		by STRUCTURE

########   RUN STRUCTURE  #################################################
	struct.exe.path <- "C:/Program Files/Structure console"
	extract.named(BSA.params)
	write.mainparams(paste(struct.exe.path,"/mainparams",sep=""),mainparams)
#	if (struct.exe.path==NULL) struct.exe.path <- "C:/Program Files/Structure console"

	agg.gs.tseries()
#	agg.gs.to.gtypes(agg.gs,gs)
#	agg.gs.tseries()
gtypes <- do.call("rbind",lapply(agg.gtypes,function(s){
	ids<-as.matrix(unclass(attr(s,'coords'))[,1])
	cbind(ids,s[,4:dim(s)[2]])
}))
#	gtypes <- agg.gs.to.gtypes(agg.gs,rland=rland)
	struct.infile <- paste(output.dir,"sim.struct.in",sep="")
	write.table(gtypes,file=struct.infile,col.names=FALSE,row.names=FALSE,quote=FALSE)

	struct.res <- structure.batch.BSA(struct.exe.path,output.dir,struct.outfile,gtypes,k)
	save(gtypes,struct.res,file=paste(output.dir,struct.outfile,".rdata",sep=""))

#########   EVALUATE RESULTS   ###########################################
	k.index <- which(struct.res$lnL == max(struct.res$lnL))
	chosen.k <- k[[k.index]]
	chosen.ass.prob <- struct.res$ass.prob[[k.index]]
	#calculate average assignment probability for each sampling site,
	#assign it to the MU for which the average assignment is highest.
	samp.site <- do.call("c",lapply(1:length(agg.gs), function(i){rep(i,dim(agg.gs[[i]])[1])}))
	ss.to.mu <- sapply(1:length(agg.gs), function(f){
		ss.ass.prob <- colMeans(chosen.ass.prob[samp.site==f,2:(chosen.k+1)])
		which(ss.ass.prob==max(ss.ass.prob))
	})
	names(ss.to.mu) <- NULL
	tab <- table(ss.to.mu)
	MUs.from.BSA <-lapply(1:length(tab),function(x){
		which(ss.to.mu==as.numeric(names(tab[x])))
	})
	MU.poly.generator(MUs.from.BSA,sample.polys,landscape.poly)
}

structure.batch.BSA <- function(struct.exe.path,output.dir,file.name,gtypes,k){

	#gtypes should be a matrix with one line per individual.  First colum
	#is GeneticID, remaining columns are alleles, two columns per locus
	#k is a list of the different values of k (# of groups) to be evaluated

	lnL <- vector("numeric",length(k))
	ass.prob <- list()

	nind <- dim(gtypes)[1]
	struct.infile <- paste(output.dir,"sim.struct.in",sep="")
#	write.table(gtypes,file=struct.infile,col.names=FALSE,row.names=FALSE,quote=FALSE)
	for (i in 1:length(k)){
		current.k <- k[[i]]
		struct.outfile <- paste(output.dir,file.name,".k",current.k,sep="")
		old.wd <- getwd()
		setwd(struct.exe.path)
		on.exit(setwd(old.wd))
		if (!file.exists(paste(struct.outfile,"_f",sep=""))){
			system(paste("structure -i",struct.infile,"-o",struct.outfile,"-K",current.k), invisible=FALSE,show.output.on.console=FALSE)
		}
		res <- scan(paste(struct.outfile,"_f",sep=""),"character")
		loc <- grep("Ln",res,value=FALSE)
		lnL[i] <- as.numeric(res[(loc+5)])
		loc <- grep("%Miss",res,value=FALSE)
		res <- matrix(res[(loc+4):(loc+3+(nind*(current.k+4)))],ncol=(current.k+4),byrow=TRUE)
		ass.prob[[i]] <- data.frame(res[,2],matrix(as.real(res[,5:(current.k+4)]),ncol=current.k))
		names(ass.prob[[i]]) <- c("ID",paste("pop",1:current.k,sep=""))
		if (!file.exists(paste(struct.outfile,"_f",sep=""))){save(lnL,ass.prob,file=paste(struct.outfile,".struct.rdata",sep=""))}
		
	}
	list(ass.prob=ass.prob,lnL=lnL)
}