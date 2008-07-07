`agg.gs.to.gtypes` <-
function(agg.gs,gs){

	if (is.null(agg.gs)){agg.gtypes <- NULL
	} else {
		####get needed variables from gs
		n.ind <- sapply(agg.gs,function(x) dim(x)[1])
		n.loci <- dim(agg.gs[[1]])[2]
		ysamp <- !is.null(gs)
		locus.info <- unclass(attr(gs[[ysamp[1]]],'loci'))
		ploidy <- sapply(locus.info, function(x) x$ploidy)
		n.hap.loci<-sum(ploidy==1)
		n.dip.loci<-sum(ploidy==2)
		n.areas<-length(agg.gs)


		agg.gtypes<-lapply(1:n.areas,function(x){
			hap.vec <- agg.gs[[x]][,1:n.hap.loci,1]
			dip.vec <- do.call("cbind",lapply((n.hap.loci+1):n.loci, function(y){
				cbind(agg.gs[[x]][,y,1],agg.gs[[x]][,y,2])
			}))
			cbind(hap.vec,dip.vec)
		})
		for (i.area in 1:length(agg.gtypes)){
			colnames(agg.gtypes[[i.area]]) <- NULL
			attr(agg.gtypes[[i.area]],'coords') <- attr(agg.gs[[i.area]],'coords')
		}
	}

	attr(agg.gtypes,'seq.list') <- attr(agg.gs,'seq.list')
	agg.gtypes
}

