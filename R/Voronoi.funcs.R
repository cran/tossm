#functions used by BSAs that use a Voronoi diagram to geographically constrain clustering

get.voronoi <- function(pts){
	voronoi <- deldir(pts[,1],pts[,2])
	plot(voronoi,lwd=1)
	edges <- voronoi$dirsgs
	edges <- cbind(edges,bound=FALSE)
	n.edges <- dim(edges)[1]
	real.edges <- sapply(1:n.edges, function(e){
	 ifelse (edges$x1[e]==edges$x2[e] & edges$y1[e]==edges$y2[e],FALSE,TRUE)
  })
  edges <- edges[real.edges,]
  n.edges <- dim(edges)[1]
  row.names(edges) <- 1:n.edges

	#make a separate matrix of vertices and their coords
	vertices1 <- as.matrix(unique(edges[,1:2]))
	vertices2 <- as.matrix(unique(edges[,3:4]))
	colnames(vertices1) <- colnames(vertices2) <- c('x','y')
	vertices <- unique(rbind(vertices1,vertices2))
	rownames(vertices) <- 1:dim(vertices)[1]
	n.verts <- dim(vertices)[1]
	on.bbox <- sapply(1:n.verts, function(v){
	 ifelse((vertices[v,1] %in% voronoi$rw[1:2] | vertices[v,2] %in% voronoi$rw[3:4]),
      1,0)
  })
  vertices <- cbind(vertices,on.bbox)

	coords <- as.matrix(edges[,1:4])

	#record which vertices each edge has
	verts <- matrix(nrow=n.edges,ncol=2,dimnames=list(1:n.edges,c('vert1','vert2')))
	for (v in 1:n.verts){
		for (e in 1:n.edges){
			if (vertices[v,1]==coords[e,1] & vertices[v,2]==coords[e,2]) verts[e,1] <- v
			if (vertices[v,1]==coords[e,3] & vertices[v,2]==coords[e,4]) verts[e,2] <- v
		}
	}

	#replace vertex coords with vertex numbers
	edges2 <- data.frame(cbind(verts,edges[,5:9]))
	return(list(edges=edges2,verts=vertices))

}

group.membership <- function(edges,verts,pts){

	n.pts <- dim(pts)[1]

	boundaries <- which(edges$bound==TRUE)
	bound.verts <- do.call("rbind",lapply(boundaries, function(b){
		rbind(verts[edges[b,1],1:2],verts[edges[b,2],1:2])
	}))

	group <- rep(0,n.pts)

	#find the point on the perimeter closest to the origin (org.pt).
	#by definition, this point belogs to group 1
	perim <- which(edges$bp1==TRUE |edges$bp2==TRUE)
	perim.pts <- unique(c(edges$ind1[perim],edges$ind2[perim]))
	perim.pts <- do.call('rbind',lapply(perim.pts,function(p){
		coords <- pts[p,1:2]
		c(pt=p,coords,dist=dist(rbind(coords,rep(0,2))))
	}))
	org.pt <- perim.pts[order(perim.pts[,4]),][1,]
	group[org.pt] <- 1

	group <- sapply(1:n.pts, function(p){
		line.segs <- rbind(pts[p,1:2],org.pt[2:3])
		intersect <- sapply(1:length(boundaries), function(b){
			line.segs <- rbind(line.segs,bound.verts[((2*b)-1):(2*b),])
			line.intersect(t(line.segs))
		})
		if (sum(intersect)-(floor(sum(intersect)/2)*2) == 0) g <- 1 else g <- 2
	})
	pts <- cbind(pts,group)

  op <- par(pch=19)
	for (i in 1:n.pts){

		color <- ifelse(pts[i,3]==1,'green','blue')
		par(col=color)
		points(pts[i,1],pts[i,2])

	}
	par(op)
	return(pts)
}


line.intersect <- function(pts){

	#determine whether two line segments intersect
	#row1=x, row2=y, cols1&2 define line segment1, cols3&4 define segment 2
	d <- (pts[2,4]-pts[2,3])*(pts[1,2]-pts[1,1])-(pts[1,4]-pts[1,3])*(pts[2,2]-pts[2,1])
	ua <- ((pts[1,4]-pts[1,3])*(pts[2,1]-pts[2,3])-(pts[2,4]-pts[2,3])*(pts[1,1]-pts[1,3]))/d
	ub <- ((pts[1,2]-pts[1,1])*(pts[2,1]-pts[2,3])-(pts[2,2]-pts[2,1])*(pts[1,1]-pts[1,3]))/d
	#from http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/ , if 0 < ua < 1
	#then intersection point lines on line a.  If 0<ub<1 then intersection point
	#lies on line b.  Therefore, line segments intersect only if ua and ub are between
	#zero and one, implying that the line segments actually intersect
  if ((ua > 0 & ua < 1) & (ub > 0 & ub < 1)) {intersect <- 1
	#if either ua or ub is between 0 and 1 and the other equals zero or one, then
	#the line segments form a T, implying that the line between the two points in
	#question passes through the vertex where two boundary segments meet.  In that
	#case, count as half an intersection, as the line will also form a T with another
	#boundary segment
	} else if (((ua==0 | ua==1) & (ub > 0 & ub < 1)) | ((ua > 0 & ua < 1) & (ub==0 | ub==1))){
	 intersect <- 0.5
  } else intersect <- 0
  return(intersect)
}

gpcpoly.centroid <- function(poly) {
  poly.pts <- get.pts(poly)[[1]]
  poly.mat <- cbind(poly.pts$x,poly.pts$y)
  ix <- c(2:dim(poly.mat)[1], 1)
  xs <- poly.mat[,1]; xr <- xs[ix]
  ys <- poly.mat[,2]; yr <- ys[ix]
 # factor <- xr*ys - xs*yr #was in code I cribbed, but resulted in sign error
  factor <- xs*yr-xr*ys
  cx <- sum((xs+xr)*factor)
  cy <- sum((ys+yr)*factor)
  scale <- 3*abs(sum(xs*yr) - sum(xr*ys))
  c(cx/scale, cy/scale)
}