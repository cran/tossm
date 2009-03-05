"index.to.state"<-
function(gs,ysamp=NULL){
######get all variables needed
if (is.null(ysamp)) ysamp <- index( sapply( gs, length)>0)
n.areas <- length( gs[[ ysamp[1]]])
n.loci <- ncol( gs[[ ysamp[ 1]]][[1]])
locus.info<-lapply(ysamp,function(x){
		  unclass(attr(gs[[x]],"loci"))
		})
#####ADDED THIS LINE--DG 09APR08#####
coords.info<-lapply(ysamp,function(x){
	lapply(1:n.areas,function(y){unclass(attr(gs[[x]][[y]],"coords"))})})
seq.list <- unclass(attr(gs,"seq.list"))
ploidy<-sapply(locus.info[[1]],function(x){x$ploidy})
n.hap.loci<-length(ploidy[ploidy==1])
n.dip.loci<-length(ploidy[ploidy==2])
n.gs.per.f <- sapply(gs[[ysamp[1]]],function(x) dim(x)[1])

######create matrices that pair each aindex with the corresponding state for each gs year
lookup.aindex.cbind<-lapply(1:length(ysamp),function(x){
    			sapply(1:n.loci,function(y){
      		  sapply(locus.info[[x]][[y]]$alleles,function(z){
        		    rbind(z$aindex,z$state)})})})

######make lists out of the haploid states in prep for finding unique states
state.list<-lapply(1:length(ysamp),function(x){
  lookup.aindex.cbind[[x]][[1]][2,]})

######find unique haploid states from across years
hap.alleles<-unique(do.call("c",state.list))
seq.list <- c(seq.list,hap.alleles[!(hap.alleles %in% seq.list)])

######replace all diploid aindexes from 'gs' with states
for (l in 1:length(ysamp)){
	gs.states <- lapply(1:n.areas,function(m){

	   states<-array(0/0, dim=c(n.gs.per.f[m],n.loci,2))
	   if (n.gs.per.f[m] > 0){
		for(q in 1:2){
			states[,,q] <-sapply(1:n.loci, function (n){
				sapply(1:n.gs.per.f[m],function(p){
					index<-which(lookup.aindex.cbind[[l]][[n]][1,]==gs[[ysamp[l]]][[m]][p,n,q])
					if (n %in% 2:n.loci){
						lookup.aindex.cbind[[l]][[n]][2,index]
					} else {
						if (q==1){
							hap <- lookup.aindex.cbind[[l]][[n]][2,index]
							which(seq.list==hap)
						} else { 0 }
					}
				})
			})
		}
	   }
	   states
	})
 	gs[[ ysamp[l]]] <- gs.states
	locus.info.yr <- locus.info[[l]]
	attr(locus.info.yr,'printinfo') <- "<< locus info in rmetasim format; use 'unclass(...)' to display >>"
	class(locus.info.yr) <- 'printinfo'
	attr(gs[[ysamp[l]]], 'loci') <- locus.info.yr
	#####ADDED THIS LITTLE LOOP--DG 09APR08#####
	for (area in 1:n.areas){
      coords.info.yr<-coords.info[[l]][[area]]
		attr(coords.info.yr,'printinfo')<-"<<matrix of individual ID's, birth years, and coordinate info; use 'unclass(...)' to display>>"
		class(coords.info.yr)<-"printinfo"
		attr(gs[[ysamp[l]]][[area]],"coords")<-coords.info.yr
	}

}

attr(seq.list, 'printinfo') <- "<< list of unique sequences; use 'unclass(...)' to display >>"
class( seq.list) <- 'printinfo'
attr(gs, "seq.list") <- seq.list
gs
}

