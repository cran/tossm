`structure.batch.BSA` <-
function(struct.exe.path,output.dir,file.name,gtypes,k){

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

