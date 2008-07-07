`MU.poly.generator` <-
function(MU.from.BSA,sample.polys,landscape.poly){

	#generate a grid based the extent of the landscape poly
	# default of 60X60
	n.cells <- 60
	xmin<-landscape.poly$x[1]
	xmax<-landscape.poly$x[2]
	ymin<-landscape.poly$y[1]
	ymax<-landscape.poly$y[2]

	#Generate grid cell polygons
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

	#identify grid points that are inside a sample polygon
	poly.subsets<-list()
	poly.ID <- cbind(grid.points,rep(-1,nrow(grid.points)))
	for (i.poly in 1:length(sample.polys)){
		io <- which(in.gpc(sample.polys[[i.poly]],grid.points)==TRUE)
		poly.subsets[[i.poly]] <- grid.points[io,]
		poly.ID[io,3] <- i.poly
	}

	#for each grid point outside of all polygons, find the minimum distance to a grid point
	#that lies within a polygon.
	for (i.point in 1:nrow(poly.ID)){
	   if (poly.ID[i.point,3]==-1){
  		min.dists<-vector(length=length(poly.subsets))
		pt1 <- grid.points[i.point,]
  		for (i.poly in 1:length(poly.subsets)){
			distance.vec<-vector(length=nrow(poly.subsets[[i.poly]]))
			for (i.sub.point in 1:nrow(poly.subsets[[i.poly]])){
				pt2 <- poly.subsets[[i.poly]][i.sub.point,]
				distance.vec[i.sub.point]<-sqrt((pt1[1]-pt2[1])^2+(pt1[2]-pt2[2])^2)
			}
			min.dists[i.poly]<-min(distance.vec)	
		}
		poly.ID[i.point,3]<-which.min(min.dists)
	   }
	}

	#divvy all grid cells into MUs....
	MU.list<-list()
	for (i.MU in 1:length(MU.from.BSA)){
		MU.list[[i.MU]] <- grid.cells[poly.ID[,3] %in% MU.from.BSA[[i.MU]]]
	}
	
	#Union all grid cells in the same MU
	union.list <- vector("list",length(MU.list))
	for (i.MU in 1:length(MU.list)){
		union.list[[i.MU]] <- poly.union(MU.list[[i.MU]])
	}

	#return the list of MU gpc.polys
	return(union.list)
}

