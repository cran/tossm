trajectory.plotter <- function(TOSSM.sim){

	n.bp <- dim(TOSSM.sim$abund.b)[1]
	n.mu <- dim(TOSSM.sim$catches)[1]
	n.yr <- dim(TOSSM.sim$catches)[2]

	par(mfrow=c(3,1))

	plot(1:n.yr,TOSSM.sim$abund.b[1,],'l',lty=1,main="Abundances",xlab="Year",ylab="# individuals")
	for (i in 2:n.bp){
		points(1:n.yr,TOSSM.sim$abund.b[i,],'l',lty=i)
	}

	ylim <- max(TOSSM.sim$catches)
	plot(1:n.yr,TOSSM.sim$catches[1,],'l',lty=1,main="Catches",xlab="Year",ylab="# individuals",ylim=c(0,max(TOSSM.sim$catches)+1))
	for (i in 2:n.mu){
		points(1:n.yr,TOSSM.sim$catches[i,],'l',lty=i)
	}

	plot(1:n.yr,TOSSM.sim$effort,'l',main="Effort",xlab="Year", ylab="Avg. distance traveled")
}