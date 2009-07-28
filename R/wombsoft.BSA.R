wombsoft.BSA<-function(gs,abund,var,C,landscape.poly,
			sample.polys,options.list){
	attach(options.list)
	#get genetic and coordinate information from run.tossm
	genes<-agg.gs.to.gtypes(agg.gs.tseries(),gs)
	coords<-lapply(genes,function(x){
			unclass(attr(x,"coords"))})
	#compile this information into a format suitable for the wombsoft package
	compiled<-lapply(1:length(genes),function(x){cbind(coords[[x]][,c(1,3,4)],genes[[x]][,5:ncol(genes[[x]])])})
	compiled<-do.call("rbind",compiled)
	write(t(compiled),file="data.txt",ncolumns=ncol(compiled))

	#perform wombsoft analyses
	data<-DataCodominant("data.txt",conversion=0,nb_x=nb_x,nb_y=nb_y,"coord_indiv.txt")
	data<-MirrorCodominant(data,m=m)
	WomblingCodominant(data,h=WCh,"systemic.txt","direction.txt")
	cde<-CandidateBoundariesCodominant(data,h=CBh,pB=pB)
	BinomialTestCodominant(data,cde,pvalue=sig.level,"boundaries.txt","dir_boundaries.txt")

	#process wombsoft output
	bounds<-read.table("boundaries.txt")
	bounds<-t(bounds)
	bounds[bounds==1]<--1


	####function for creating a matrix with proper
	####MU designations

	fill<-function(bounds){

		#fill in any narrow 'isthmi' in the potential MUs
		for (i in 1:nrow(bounds)){
			for (j in 1:ncol(bounds)){
				r.mar<-nrow(bounds)
				c.mar<-ncol(bounds)
				if(bounds[i,j]==0){
					adj.row<-adj.col<-rep(-1,2)
					if(i>1){adj.row[1]<-bounds[i-1,j]}
					if(i<r.mar){adj.row[2]<-bounds[i+1,j]}
					if(j>1){adj.col[1]<-bounds[i,j-1]}
					if(j<c.mar){adj.col[2]<-bounds[i,j+1]}
					if(all(adj.row==-1)|all(adj.col==-1)){bounds[i,j]<--1}
				}
			}
		}

		#initialize MU#
		#go through all groups of adjacent cells and give them an MU group designation
		MU<-1
		all.filled<-FALSE
		while(all.filled==FALSE){
			for (i in 1:ncol(bounds)){
				for (j in 1:nrow(bounds)){
					#find the first (next) cell that contains zero (is part of an MU)
					if(bounds[j,i]==0){
						#intialize a matrix to keep track of neighboring cells, as well as
						#a matrix to keep track of all cells visited in the past
						next.cells<-t(as.matrix(c(j,i)))
						all.cells<-matrix(c(0,0),1,2)
						#while there are still zeros in this MU, fill them with an MU designation
						MU.done<-FALSE
						while(MU.done==FALSE){
							#get the status of adjacent cells
							adj.val<-matrix(99,nrow(next.cells),4);adj.add<-array(0,c(nrow(next.cells),4,2))
							#get all the adjacent values and matrix indices
							for (z in 1:nrow(next.cells)){
								irow<-next.cells[z,1];icol<-next.cells[z,2]
								if(irow!=1){adj.val[z,1]<-bounds[irow-1,icol];adj.add[z,1,]<-c(irow-1,icol)}#up
								if(icol!=1){adj.val[z,2]<-bounds[irow,icol-1];adj.add[z,2,]<-c(irow,icol-1)}#left
								if(irow!=nrow(bounds)){adj.val[z,3]<-bounds[irow+1,icol];adj.add[z,3,]<-c(irow+1,icol)}#down
								if(icol!=ncol(bounds)){adj.val[z,4]<-bounds[irow,icol+1];adj.add[z,4,]<-c(irow,icol+1)}#right
							}
						
						
							#bring up a matrix of all cells that have already been visited
							#check to see which adjacent cells have already been visited
							bounds.done<-array(NA,c(nrow(next.cells),nrow(all.cells),4))
							for (q in 1:nrow(next.cells)){
								for(r in 1:nrow(all.cells)){
									for(s in 1:4){
										bounds.done[q,r,s]<-all(adj.add[q,s,]==all.cells[r,])
									}
								}
							}
				
							#tag those cells already visited
							for (y in 1:nrow(next.cells)){
								for (l in 1:4){
									if(any(bounds.done[y,,l]==TRUE)){adj.val[y,l]<-99}
								}
							}
							#all cells that have not been visited, and contain zeros, are new cells 
							new.cells<-NULL
				      		for (a in 1:dim(adj.add)[1]){
								new.cells<-rbind(new.cells,adj.add[a,which(adj.val[a,]==0),])
							}
							new.cells<-unique(new.cells)
							all.cells<-rbind(all.cells,new.cells)
							next.cells<-new.cells
							#if no new cells exist, go on to the next MU
							if(nrow(new.cells)==FALSE){
								bounds[j,i]<-MU
								for (i in 1:nrow(all.cells)){bounds[all.cells[i,1],all.cells[i,2]]<-MU}
								MU<-MU+1
								MU.done<-TRUE
								i<-j<-1
							}
						}
					}
				}
			}
			if(all(bounds!=0)){all.filled<-TRUE}
		}
	return(bounds)
	}
	
	####function for converting MU matrix to polygons
	womb.mat.to.polys<-function(landscape.poly,bounds,min.MU.size){

		#create a matrix of MU designation for all point indices
		#(necessary for joining to grid.points later)
		poly.ID<-vector(length=nrow(bounds)*ncol(bounds))
		k<-1
		for(j in 1:nrow(bounds)){
			for(i in 1:ncol(bounds)){
				poly.ID[k]<-bounds[j,i]
				k<-k+1
			}	
		}

		#Create a list of MUs
		MU.from.BSA<-list()
		length(MU.from.BSA)<-max(poly.ID)
		for (u in 1:length(MU.from.BSA)){
			MU.from.BSA[[u]]<-u
		}

		#create all grid cell polygons
		n.cells <- nrow(bounds)
		xmin<-landscape.poly$x[1]
		xmax<-landscape.poly$x[2]
		ymin<-landscape.poly$y[1]
		ymax<-landscape.poly$y[2]

		vertices <- as.points(expand.grid(x=seq(xmin,xmax,length.out=n.cells+1),y=seq(ymin,ymax,length.out=n.cells+1)))
		grid.cells <- unlist(lapply(1:n.cells, function(y){
			lapply(1:n.cells, function(x){
				index <- x+(y-1)*(n.cells+1)
				as(rbind(vertices[index,],vertices[index+1,],vertices[index+n.cells+2,],vertices[index+n.cells+1,]),"gpc.poly")
			})
		}),recursive=FALSE)

		#Calculate center point of grid cells
		grid.points <- as.points(do.call("rbind",lapply(grid.cells, function(c){
			pts <- get.pts(c)[[1]]
			c(mean(pts$x),mean(pts$y))
			})))
		
		poly.ID<-cbind(grid.points,poly.ID)
		poly.subsets<-list()
		for (i.poly in 1:length(MU.from.BSA)){
			poly.subsets[[i.poly]] <- poly.ID[which(poly.ID[,3]==i.poly),]
		}
	

		#for each grid point outside of all polygons, find the minimum distance to a grid point
		#that lies within a polygon.
		for (i.point in 1:nrow(poly.ID)){
	   		if (poly.ID[i.point,3]==-1){
  				min.dists<-vector(length=length(poly.subsets))
				pt1 <- grid.points[i.point,]
  				for (i.poly in 1:length(poly.subsets)){
					len<-nrow(poly.subsets[[i.poly]])
					if(is.null(len)){len<-1}
					distance.vec<-vector(length=len)
					for (i.sub.point in 1:len){
						if(len==1){pt2<-poly.subsets[[i.poly]]} else{
						pt2 <- poly.subsets[[i.poly]][i.sub.point,]}
						distance.vec[i.sub.point]<-sqrt((pt1[1]-pt2[1])^2+(pt1[2]-pt2[2])^2)
					}
					min.dists[i.poly]<-min(distance.vec)	
				}
				poly.ID[i.point,3]<-which.min(min.dists)
	   		}
		}

		#put all grid cells within their assigned MU
      	MU.list<-list()
		for (i.MU in 1:length(MU.from.BSA)){
			MU.list[[i.MU]] <- grid.cells[poly.ID[,3] %in% MU.from.BSA[[i.MU]]]
		}
	
		#Union all grid cells in the same MU
		union.list <- vector("list",length(MU.list))
		for (i.MU in 1:length(MU.list)){
			union.list[[i.MU]] <- poly.union(MU.list[[i.MU]])
		}
		#figure out if any polygons have multiple contours
		#if they do, keep the biggest one and ignore the small ones
		for (i in 1:length(union.list)){
			contours<-get.pts(union.list[[i]])
			if(length(contours)>1){
				sizes<-vector(length=length(contours))
				for (i.con in 1:length(contours)){
					sizes[i.con]<-area.poly(as(cbind(contours[[i.con]]$x,contours[[i.con]]$y),"gpc.poly"))
				}
				max.poly<-which.max(sizes)
				union.list[[i]]<-as(cbind(contours[[max.poly]]$x,contours[[max.poly]]$y),"gpc.poly")
			}
		}
		#merge polys smaller than min.MU.size with largest adjacent poly
		l.p<-landscape.poly
		tot.area<-area.poly(as(matrix(c(l.p$x[1],l.p$y[1],l.p$x[2],l.p$y[1],
				l.p$x[2],l.p$y[2],l.p$x[1],l.p$y[2]),ncol=2,byrow=2),"gpc.poly"))
		min.MU.size<-min.MU.size*tot.area
		sizes<-sapply(union.list,area.poly)
		too.small<-TRUE
		if(min(sizes)>min.MU.size){too.small<-FALSE}
		cnt<-0
		while(too.small==TRUE){
			n.mus<-length(sizes)
			biggest<-union.list[[which(sizes==sort(sizes)[n.mus-cnt])]]
			smallest<-union.list[[which(sizes==sort(sizes)[1])]]
			p.union<-union(biggest,smallest)
			contours<-length(get.pts(p.union))
			if(contours==1){
				union.list[[which(sizes==sort(sizes)[n.mus-cnt])]]<-p.union
				union.list[[which(sizes==sort(sizes)[1])]]<-NULL
				cnt<-0
			} else {cnt=cnt+1}
			sizes<-sapply(union.list,area.poly)
			if(min(sizes)>min.MU.size){too.small<-FALSE}
		}
	
		munits<-union.list
		return(munits)

	}

	bounds<-fill(bounds)
	munits<-womb.mat.to.polys(landscape.poly,bounds,min.MU.size)
	return(munits)
}