`tossm.plotter` <-
function(bp.polys,sample.polys,historic.polys){

#organize polygons for plotting

	poly.list<-list(bp.polys,sample.polys,historic.polys)
	poly.list<-poly.list[!sapply(poly.list,is.null)]
	plot.list<-lapply(poly.list,function(x){lapply(x,as,"matrix")})
	poly.labels<-lapply(plot.list,function(x){
					lapply(x,function(y){c(mean(y[,1]),mean(y[,2]))})})

#create color, text, and position vectors for the various plots
#colors
	plot.col<-matrix(rep(240,9),ncol=3,byrow=TRUE)
	poly.col<-list()
	for (i.plot in 1:length(plot.list)){
		plot.col[i.plot,i.plot]<-250
		poly.col[[i.plot]]<-matrix(50,length(poly.list[[i.plot]]),3)
		for (i.poly in 1:length(poly.list[[i.plot]])){
			poly.col[[i.plot]][i.poly,i.plot]<-100+(255-100)/i.poly		
		
		}
	}

plot.col<-apply(plot.col,1,function(x){
			rgb(x[1],x[2],x[3],maxColorValue=255)})
poly.col<-lapply(poly.col,function(x){
			apply(x,1,function(y){
			rgb(y[1],y[2],y[3],maxColorValue=255)})
	    })

#text
	poly.text.vec<-c("BP","S","H")
	header.text.vec<-c("Breeding Population Polygons","Sampling Polygons"
					,"Historic Removal Polygons")

#positions
	xpos<-c(0,1/3,2/3)
	leg.symbols.pos<-c(.05,.35,.6)

#create lists(gLists) and lists of lists(gTrees) for the graphical objects
	poly.Text<-poly.Shape<-legend.boxes<-gTree()
	bp.outline<-header.text<-plot.frames<-legend.text<-gList()

#create viewports
	parent.frame<-viewport()
	child.frames<-vpList()

#make an x- and y-extent for plotting polygons
	box <- get.bbox(do.call("poly.union",list(lapply(poly.list, function(x) poly.union(x)))))
	x.buff <- (box$x[2]-box$x[1]) * 0.1
	y.buff <- (box$y[2]-box$y[1]) * 0.1
	x.extent <- c(box$x[1]-x.buff,box$x[2]+x.buff)
	y.extent <- c(box$y[1]-y.buff,box$y[2]+y.buff)

#loop through the 3 main plots (bp's, samp's, and hist.removals)
	for(i.plot in 1:length(plot.list)){

	#make a frame for the 3 main plots
		child.frames[[i.plot]]<-viewport(x=xpos[i.plot],y=1,width=1/3,height=.8,just=c("left","top"),
							xscale=x.extent,yscale=y.extent)
		

	#header text for each plot
		header.text[[i.plot]]<-textGrob(label=header.text.vec[i.plot],x=x.extent[1]+(x.buff*.5),
							y=y.extent[2]-(y.buff*.5),just="left",gp=gpar(col="black"),
							default.units="native")
	
	#and a framing rectangle for each plot
		plot.frames[[i.plot]]<-rectGrob(gp=gpar(fill=plot.col[i.plot]))

	#create text for the legend
		legend.text[[i.plot]]<-textGrob(label=header.text.vec[i.plot],x=leg.symbols.pos[i.plot]+.05,y=.5,
						just="left",gp=gpar(cex=0.7),
						default.units="native")

	#and empty lists for each type of polygon
		poly.Text[[i.plot]]<-poly.Shape[[i.plot]]<-legend.boxes[[i.plot]]<-gList()

	#offset for placement of legend slices
		offset<-0
		for (i.poly in 1:length(plot.list[[i.plot]])){
			#make an outline of the bps that will go in all 3 plots
			if(i.plot==1){
				bp.outline[[i.poly]]<-polygonGrob(x=plot.list[[i.plot]][[i.poly]][,1],
									y=plot.list[[i.plot]][[i.poly]][,2],
									#gp=gpar(fill="transparent",col=poly.col[[i.plot]][i.poly],
									gp=gpar(fill="transparent",col="black",
									lty=i.poly),default.units="native")
			}
			#and filled in polygons bps, samps, and hist for each plot
			poly.Shape[[i.plot]][[i.poly]]<-polygonGrob(x=plot.list[[i.plot]][[i.poly]][,1],
									y=plot.list[[i.plot]][[i.poly]][,2],
									gp=gpar(fill=poly.col[[i.plot]][i.poly],
									lty=ifelse(i.plot==1,i.poly,1)),
									default.units="native")
			#plot labels for each polygon of each type
			label=paste(poly.text.vec[[i.plot]],i.poly,sep="")
			poly.Text[[i.plot]][[i.poly]]<-textGrob(label,x=poly.labels[[i.plot]][[i.poly]][1],
									y=poly.labels[[i.plot]][[i.poly]][2],
									gp=gpar(col="white"),default.units="native")
			#and legend boxes for each polygon and color
			x.width<-.05/length(plot.list[[i.plot]])
			x.placement<-x.width*offset
			legend.boxes[[i.plot]][[i.poly]]<-rectGrob(x=leg.symbols.pos[i.plot]+x.placement,y=.5,
						width=x.width,
						height=0.5,gp=gpar(fill=poly.col[[i.plot]][i.poly],col="transparent"),
						default.units="native")
			offset<-offset+1
		}
	}

#make one more plotting region for the legend
	if(length(child.frames)==2){
		child.frames[[3]]<-viewport(x=2/3,y=1,width=1/3,height=.8,just=c("left","top"))
	}
	child.frames[[4]]<-viewport(x=0,y=0.2,width=1,height=.2,just=c("left","top"))

#open a window and make the 'parent' frame active
	X11(width=8,height=4)
	pushViewport(parent.frame)

#successively open, plot items, and close the 4 frames.
	for (i.plot in 1:4){
		pushViewport(child.frames[[i.plot]])
		if (i.plot %in% 1:length(plot.list)){
			grid.draw(plot.frames[[i.plot]])
			grid.draw(header.text[[i.plot]])
			for (i.poly in 1:length(plot.list[[i.plot]])){
				grid.draw(poly.Shape[[i.plot]][[i.poly]])
				grid.draw(poly.Text[[i.plot]][[i.poly]])
			}
			grid.draw(bp.outline)
		} 
		if(i.plot==3 & length(plot.list)==2){
			grid.rect(gp=gpar(fill="grey90"))
			grid.text(x=0.5,y=0.7,"No historic removal \npolygons specified")
			grid.text(x=0.5,y=0.55,gp=gpar(cex=0.8),
				"(Initial depletion \nwill be used instead)")
		}
		if(i.plot==4){
			grid.rect(gp=gpar(fill="grey90"))
			for (i in 1:length(plot.list)){
				grid.draw(legend.boxes[[i]])
			}
			grid.draw(legend.text)	
		}		
		popViewport()
	}
}

