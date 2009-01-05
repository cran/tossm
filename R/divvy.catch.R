`divvy.catch` <-
function(locs, mu.polys){
	sapply(mu.polys, function(mu){
		sum(in.gpc(mu,locs[,3:4]))})
}

