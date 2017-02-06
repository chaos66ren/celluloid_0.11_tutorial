# This script doesn't make use of allelic ratio data, for instance when coverage is low, 
# or when normal data is not available. 
# Assuming low coverage here, so that larger bins (e.g. 50000kb) are used 
#
# Usage:
# Rscript pipeline.r Data/tumour.wig Data/gc.wig Data/map.wig samplename

library(celluloid)

if( !file.exists("RdaPipeline") ){ 
  system("mkdir RdaPipeline")
}

argv<-commandArgs( trailing=TRUE )

argc<-length(argv)
if( argc!=4 ){ 
  stop("Usage: Rscript pipeline.r tumourWigFile gcWigFile mapWigFile samplename") }

##############
# prepare data
##############

tumourWigFile<-argv[1]
gcWigFile<- argv[2]
mapWigFile<- argv[3]
name<-argv[4]

if( !file.exists("RdaPipeline/tc.rda") ){
  # load and segment tumor data
  t <- wigsToRangedData( tumourWigFile, gcWigFile, mapWigFile )
  # workaround for a bug is to reduce samplesize
  tc<-gcCorrect( t, sampletype="tumor", samplesize=25000 )
  t.seg <- segmentSeqData( tc , k=10 , gamma=50  )
  save( t.seg, file="RdaPipeline/t.seg.rda")
  save( tc, file="RdaPipeline/tc.rda")
} else {
  load("RdaPipeline/t.seg.rda")
  load("RdaPipeline/tc.rda")
}

if( !file.exists("RdaPipeline/copyAr.rda") ){
  # creates the object used to draw a contour plot. 
  copyAr<-  prepCopyAr(seg= t.seg,  tumourrangedata= tc , xonly=TRUE  )
  save(copyAr, file="RdaPipeline/copyAr.rda")
} else {
  load("RdaPipeline/copyAr.rda")
}

if( !file.exists( "RdaPipeline/cntr.rda") ){
  cntr<-showTumourProfile( copyAr , flatten=.25 , nlev=20 , noise= 0.01 , 
                           maxPoints=200000 , plot=F  )
  save(cntr, file="RdaPipeline/cntr.rda")
} else {
  load( "RdaPipeline/cntr.rda" )
}


#################
# analysis starts
#################

lm<- coverParamSpace( copyAr=copyAr,xonly=TRUE, control=list( maxit=1000 ) ,
                      maxc=12,optimFct=1 , nrep=50, Sfrom=0.25, Sto=1.25, 
                      lowerF=c(0), upperF=c( .25 ), addToParamSpace=T  )

save(lm, file=paste( "RdaPipeline/lm_",name,".rda", sep="") ) 

localSolutions <-getLocalSolutions(lm, max=TRUE) 
# only solutions for which 10% or more of the genome was captured are in the output
localSolutions<-localSolutions[ localSolutions$value>0.10 , ]

pdf( paste( "contour_solutions_",name,"_%1d.pdf", sep="") , height=12, width=16 )
par( mfrow=c(3,2) )
plot(  paramSpace[,3] , 100*paramSpace[,1] , pch=19 , ylim=c(0,100), xlab="%normal", ylab="percent captured by model")
for( i in 1:nrow(localSolutions) ){
  r<-as.numeric(rownames(localSolutions)[i] )
  image( cntr, col=terrain.colors(50))
  contour(cntr, nlev=8, add=T )
  
  epp<-  plotModelPeaks(lm[[r]]$par, selectedPoints=NULL,cn=cn, epcol="red",epcex=0 ,eplwd=3 , addlabels=F , xonly=TRUE )
  legend( 0,1, floor(1000*( lm[[r]]$value) )/10, bg="white"  )
}
dev.off() 



