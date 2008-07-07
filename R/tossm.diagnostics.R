`tossm.diagnostics` <-
function(bp.polys,sample.polys,historic.polys,rland){
	
#function to determine if sampling or historic polygons overlap
	poly.overlap<-function(polys){
		area.of.intersect <- area.poly(poly.intersect(polys))
		ifelse(area.of.intersect>0,TRUE,FALSE)
	}

#make a list of the 3 kinds of polys
	poly.list<-list(bp.polys,sample.polys,historic.polys)

#text vectors with error statement fragments
	err.text.1<-c("Breeding Population Polygons","Sampling Polygons",
			"Historic Removal Polygons")
	err.text.2<-c("do not equal the number of 'rland' populations",
				"are not fully contained within the","There is overlap between")

#loop through the polygon inputs and tally the number of problems
	check<-0
	#Check to make sure n.bp.polys==rland$intparam$habitats
	if (length(poly.list[[1]])!=rland$intparam$habitats){
		print(paste(err.text.1[1],err.text.2[1],sep=" "))
		check<-check+1
	}

	if (is.null(historic.polys)) loop<-2:2 else loop<-2:3
	for (i in loop){
		#determine whether any area of the 
		#samp or historic polys outside of the bp.polys
		area.outside<-area.poly(setdiff(poly.union(poly.list[[i]]),poly.union(poly.list[[1]])))
		if(area.outside>0){
			print(paste(err.text.1[i],err.text.2[2],
				err.text.1[1],sep=" "))
			check<-check+1
		}
		#...and these three lines determine if polygons of a kind overlap...
		if(length(poly.list[[i]]) > 1 & poly.overlap(poly.list[[i]])==TRUE){
			print(paste(err.text.2[3],err.text.1[i],sep=" "))
			check<-check+1
		} 	
	} 
	
	#...if more than zero problems exist, stop the simulation

	return(check)
}

