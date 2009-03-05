write.mainparams <- function(mainparams.file,params){

#params contains burnin,iterations,numinds,numloci,ploidy,maxpops,noadmix,freqscorr

#mainparams.file <- "C:/Program Files/Structure Console/mainparams"
	write(c(
"Data File",

"#define INFILE  testdata1   // (str) name of input data file",
"#define OUTFILE results  //(str) name of output data file",

paste("#define NUMINDS   ",params$numinds,"    // (int) number of diploid individuals in data file",sep=""),
paste("#define NUMLOCI   ",params$numloci,"   // (int) number of loci in data file",sep=""),
"#define LABEL     1     // (B) Input file contains individual labels",
"#define POPDATA   0     // (B) Input file contains a population ",
"                             identifier",
"#define POPFLAG   0     // (B) Input file contains a flag which says ",
"                              whether to use popinfo when USEPOPINFO==1",
"#define PHENOTYPE 0     // (B) Input file contains phenotype information",
"#define EXTRACOLS 0     // (int) Number of additional columns of data ",
"                             before the genotype data start.",
"#define PHASEINFO 0     // (B) the data for each individual contains a line",
"                                  indicating phase",
"#define MARKOVPHASE  0  // (B) the phase info follows a Markov model.",

"#define MISSING      -9 // (int) value given to missing genotype data",
paste("#define PLOIDY       ",params$ploidy,"  // (int) ploidy of data",sep=""),

"#define ONEROWPERIND 1  // (B) store data for individuals in a single line",
"#define MARKERNAMES  0  // (B) data file contains row of marker names",
"#define MAPDISTANCES 0  // (B) data file contains row of map distances ",
"                             // between loci",


"Program Parameters",

"#define MAXPOPS    2    // (int) number of populations assumed",
paste("#define BURNIN    ",params$burnin," // (int) length of burnin period",sep=""),
paste("#define NUMREPS   ",params$iterations," // (int) number of MCMC reps after burnin",sep="")),
file=mainparams.file,ncol=1)
}



