"raw.CLA" <-
function(Catches,Abund.est,CV.est,Years.est,First.Year,
    Last.Year,Directory=getwd()){
 Yrs <- First.Year %upto% Last.Year
 Nyears <- (Last.Year-First.Year+1)

 FileName <- file.path( Directory,"CLC.DAT")
 write("FROM CLA SIMULATION",file=FileName)
 write("",file=FileName,append=TRUE)
 xx <- paste("First year of input          ",First.Year,sep="")
 write(xx,file=FileName,append=TRUE)
 xx <- paste("Year of next catch limit     ",(Last.Year+1),sep="")
 write(xx,file=FileName,append=TRUE)
 write("Test for phaseout            1",file=FileName,append=TRUE)

 # Output the catch stream
 write("",file=FileName,append=TRUE)
 write("Catch data:",file=FileName,append=TRUE)
 write("(I4,1x,F7.0)",file=FileName,append=TRUE)
 for (II in 1:Nyears)
  write(sprintf("%i %7.1f",Yrs[II],Catches[II]),file=FileName,append=TRUE)
 write("-1",file=FileName,append=TRUE)

 # Count the number of useable abundance estimates (estimates can be 0 but not CVs)
 Nsight <- 0
 for (Isight in 1:length(Abund.est))
  if (Abund.est[Isight] > 0) Nsight <- Nsight + 1

 # Output the abundance estimates
 write("",file=FileName,append=TRUE)
 write("Abundance data:",file=FileName,append=TRUE)
 xx <- paste("No. of sightings estimates     ",sprintf("%1.0f",Nsight),sep="")
 write(xx,file=FileName,append=TRUE)
 write("(I4,F10.0,1000F8.0)",file=FileName,append=TRUE)
 for (II in 1:length(Abund.est))
  if (Abund.est[II] > 0)
   write(sprintf("%1.0f%10.0f%7.1f",Years.est[II],Abund.est[II],1.0/CV.est[II]^2),file=FileName,append=TRUE)
 write("",file=FileName,append=TRUE)

 # Zero sightings
 write("Number of zero estimates        0",file=FileName,append=TRUE)
 write("(I4,F10.0)",file=FileName,append=TRUE)

 # Do the actual run
 #RunFile <- RunFile # otherwise getwd() stuffs up!!!
 od <- getwd()
 on.exit( setwd( od))
 setwd(Directory)
 ops<-options()
 if(.Platform$OS.type=="unix"){options(warn=-1)}
 .Fortran("manage",PACKAGE="tossm")
 options(ops)
 # Extract the quota from RES.OUT
 FileName <- file.path( Directory,"RES.OUT")
 TheQuota <- scan(FileName,n=1, quiet=TRUE)
 return(TheQuota)
}


