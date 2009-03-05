
PBR<-function(Catches,est.abund.hist,CV.abund.hist,
  r.max=0.04,F.r=1,multiplier=1,...){

  #skim latest abund and CV for each MU off respective matrices
  latest.abund.est<-est.abund.hist[max(row(est.abund.hist)[which(!is.nan(est.abund.hist))]),]
  latest.CV.est<-CV.abund.hist[max(row(CV.abund.hist)[which(CV.abund.hist>0)]),]
  N.min.input<-cbind(latest.abund.est,latest.CV.est)

  #compute N.min, PBR by MU
  PBR.TAC<-round(apply(N.min.input,1,function(x){
		  N.min<-x[1]*(exp(-0.842*(sqrt(log(1+x[2]^2)))))
		  PBR<-N.min*(0.5*r.max)*F.r*multiplier
		  PBR
		  }))
  PBR.TAC[which(is.nan(PBR.TAC))] <- 0
  return(PBR.TAC)
}


