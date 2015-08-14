#  Rscript pipeline.r Data/tumour.wig Data/normal.wig Data/gc.wig Data/map.wig Data/AR.txt thisSample 

library(celluloid11)

if( !file.exists("Rda") ){ 
  system("mkdir Rda")
}

argv<-commandArgs( trailing=TRUE )

argc<-length(argv)
if( argc!=6 ){ stop("Usage: Rscrip pipeline.r tumourWigFile normalWigFile gcWigFile mapWigFile arFile name") }

##############
# prepare data
##############

tumourWigFile<-argv[1]
normalWigFile<-argv[2]
gcWigFile<- argv[3]
mapWigFile<- argv[4]
arFile<- argv[5]
name<-argv[6]

if( !file.exists("Rda/nc.rda") ){
  # load and segment normal data
  n<- wigsToRangedData( normalWigFile, gcWigFile, mapWigFile )
  # get rid of non-autosomal chromosomes
  sel<-is.element( n$space, paste("chr",1:22, sep="")  )
  n<-n[sel,]
  # gc correction
  nc<-gcCorrect(n) 
  n.seg <- segmentSeqData( nc, k=50 , maskmap = 0.8 )
  save( n.seg, file="Rda/n.seg.rda")
  save( nc, file="Rda/nc.rda")
} else {
  load("Rda/n.seg.rda")
  load("Rda/nc.rda")
}


if( !file.exists("Rda/tc.rda") ){
  # load and segment tumor data
  t <- wigsToRangedData( tumourWigFile, gcWigFile, mapWigFile )
  sel<-is.element( t$space, paste("chr",1:22, sep="") )
  t<-t[sel,]
  tc<-gcCorrect( t, sampletype="tumor" )
  t.seg <- segmentSeqData( tc , k=50  )
  save( t.seg, file="Rda/t.seg.rda")
  save( tc, file="Rda/tc.rda")
} else {
  load("Rda/t.seg.rda")
  load("Rda/tc.rda")
}


# identify segments not "normal" in normal
if(!file.exists("Rda/t.seg.mask.rda") ){
  t.n.seg<-intersectSegments( t.seg, n.seg )
  t.n.seg<-t.n.seg[ !is.na( t.n.seg$mean) & !is.na( t.n.seg$mean.1),]
  sel<-n.seg$end.pos-n.seg$start.pos > 100000 & n.seg$meanmap>.8 
  bp<-boxplot( n.seg$mean[sel], range=3, plot=F  )
  nrange<-c(bp$stats[1,1], bp$stats[5,1] )
  mask<- t.n.seg$mean.1>nrange[2] | t.n.seg$mean.1<nrange[1]
  t.seg.mask<-t.n.seg[,1:8]
  t.seg.mask$mask<- mask
  save(t.seg.mask, file="Rda/t.seg.mask.rda") 
} else {
  load("Rda/t.seg.mask.rda")
}

if( !file.exists("Rda/t.ar.seg.rda" ) ){
  # reading allelic ratio file and segmenting
  ar<-read.table(arFile, head=T )
  save( ar, file="Rda/ar.rda") 
  ar.seg<- segmentAR( ar, tc ) 
  save( ar.seg, file="Rda/ar.seg.rda")
  t.ar.seg <- intersectSegments( t.seg.mask, ar.seg  )
  t.ar.seg<-t.ar.seg[ !apply( is.na( t.ar.seg[,c("mean","meanar")] ), 1, any),]
  t.ar.seg <-arInSeg( t.ar.seg, ar,  tumourrangedata=tc , minhet = 50 )
  save( t.ar.seg, file="Rda/t.ar.seg.rda")
} else {
  load("Rda/ar.rda")
  load("Rda/t.ar.seg.rda")
}


if( !file.exists("Rda/copyAr.rda") ){
  # creates the object used to draw a contour plot. 
  mask<- t.ar.seg$mask | is.na(t.ar.seg$mask)
  copyAr<-  prepCopyAr( t.ar.seg[ !mask ,], ar,  tc  )
  save(copyAr, file="Rda/copyAr.rda")
} else {
  load("Rda/copyAr.rda")
}

if( !file.exists( "Rda/cntr.rda") ){
  cntr<-showTumourProfile( copyAr , flatten=.25 , nlev=20 , noise= 0.01 , 
                           maxPoints=200000 , plot=F  )
  save(cntr, file="Rda/cntr.rda")
} else {
  load( "Rda/cntr.rda" )
}


#################
# analysis starts
#################

# finding LOH curve 
Sn<- estimateLOHcurve(t.ar.seg)

# grid search
sel <- t.ar.seg$mean < max( cntr$x ) & !is.na( t.ar.seg$mean ) & !is.na( t.ar.seg$p )
segmentsubset <-t.ar.seg[sel,] 

# using optim

# upper bound for %normal cell set to Sn/.25, set to 1 if you think your sample's ploidy can 
# be above 8 (2/.25)

lm<- coverParamSpace( segments=segmentsubset,percentObj=T,  control=list( maxit=1000 ), Sn=Sn ,
                      maxc=12,optimFct=1 , nrep=50, 
                      lowerF=c(0), upperF=c( Sn/.25 ), addToParamSpace=T  )

save(lm, file=paste( "Rda/lm_",name,".rda", sep="") ) 

localMins <-getLocalMins(lm) 
localMins<-localMins[ localMins$value<.9 , ]

pdf( paste( "contour_solutions_",name,"_%1d.pdf", sep="") , height=12, width=16 )
par( mfrow=c(3,2) )
plot(  paramSpace[,3] , 100*(1- paramSpace[,1]) , pch=19 , ylim=c(0,100), xlab="%normal", ylab="percent captured by model")
for( i in 1:nrow(localMins) ){
  r<-as.numeric(rownames(localMins)[i] )
  image( cntr, col=terrain.colors(50))
  contour(cntr, nlev=8, add=T )
  sel<-t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9
  le<- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel]
  cxcut<- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
  points( x<-t.ar.seg$mean[sel], y<-t.ar.seg$p[sel],  pch=21 , col="blue", lwd=3 , cex=cxcut  )
  points( t.ar.seg$mean[sel], t.ar.seg$p[sel],  pch=19 ,  col="white" , cex= cxcut  - .5 )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=21 ,  col="blue", lwd=3  , cex=cxcut  )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=19 , col="white" , cex=cxcut -.5 )
  xxx <- seq( Sn, 2, .01 ) 
  points( xxx , arloh<- ARloh( xxx , 1 , Sn ) , type='l' , lwd=2  )
  points( xxx , 1-arloh , type='l' , lwd=2  )
  
  epp<-  plotModelPeaks(lm[[r]]$par, selectedPoints=NULL,cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=F )
  legend( 0,1, floor(1000*( 1-lm[[r]]$value) )/10, bg="white"  )
}
dev.off() 



