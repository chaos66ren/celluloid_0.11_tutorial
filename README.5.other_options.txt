0        1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890

###############################################################################
## if the user has saved the rda files, they can be loaded:
## files<-system("ls Rda/*rda", intern=T); for( f in files){load(f)}
###############################################################################

OTHER OPTIONS: forced peak

To force a peak to represent a particular combination of copy numbers, add a
third column to the sp data.frame that will represent the copy number 
combination associated with the corresponding line:

sp
           x         y
11 0.3119248 0.0000000
10 0.6044760 1.0000000
12 0.7939828 0.3800000
3  0.9045226 0.3265306
6  1.0753769 0.5102041
9  1.1959657 0.5000000
5  1.4874372 0.5918367

First, add a third column with 0s:

sp$w <- 0

To force the 5th peak (5th line) to represent 2 copies in all subclones (in a 
model with 2 subclones) add the following entry:

sp$w[5]<-22

sp

          x         y  w
1 0.3035810 1.0000000  0
2 0.6018622 0.0000000  0
3 0.7814362 0.6008772  0
4 0.8941881 0.6613782  0
5 1.0682875 0.4909091 22
6 1.1847642 0.5102041  0
7 1.4728393 0.6046972  0

There can only be one non-zero entry in the third column of sp. No more than 9 
copies in a subclone can be specified. 

In a 3-subclone model, the 3-digit value 224 would indicate 2 copies in the first subclone, 2 in the second and 4 in the third.

The peak is forced by adding a penalty to the objective function, with control 
parameters penaltymultiplier and penaltymultipliery in coverParamSpace. These
control parameters are passed to eoDist.

li5<-coverParamSpace( sp, optimFct=2, lowerF=c(0,0), upperF=c(1,1),  
       Sfrom=.25, Sto=2 , maxc=6 , maxsubcldiff=0.5,  
       control=list( maxit=100   ) , penaltymultiplier = 10 )


###############################################################################

OTHER OPTIONS: grid search 

If optimFct is a vector, then grid search is performed with optim() using a 
number of starting values equal to prod(optimFct). Starting values for 
parameters S and t[1],...,t[nsubclones] are the mid-points of the intervals
defined by seq(Sfrom, Sto, len=optimFct[1]+1) (for S), and 
seq( lowerF[i], upperF[i], len=optimFct[i+1]+1 ) (for t[i]).  

In other words, in the call

li6<-coverParamSpace( sp, optimFct=c(5,3), lowerF=c(0), upperF=c( 0.10 ),  
       Sfrom=.25, Sto=1 , maxc=12 , control=list( maxit=1000  ) )

## save( li6, file="Rda/li6.rda") 

the starting values from S are the mid points of the intervals defined by 

[1] 0.25 0.40 0.55 0.70 0.85 1.00

which are:

[1] 0.325 0.475 0.625 0.775 0.925

while starting values for t[1] (% of normal cells) are 

[1] 0.01666667 0.05000000 0.08333333

Each possible pair is used as starting point with the optim() function. The list
returned by coverParamSpace will be of length 5x3=15.

li6[[1]]$par
[1] 0.3170566 0.0000000

...

 li6[[15]]$par
[1] 1.00000000 0.02508548

Since different starting points can converge to the same solution, solutions 
can be trimmed using 

local<- getLocalMins(li6 )

local

         value         S           N        T1
2  0.2832984 0.6007250 0.009131892 0.9908681
1  0.5011141 0.3029616 0.005516430 0.9944836
11 0.5519167 0.3196328 0.083492540 0.9165075
7  0.8753860 0.2500000 0.009887068 0.9901129
4  1.1715130 0.7698892 0.000000000 1.0000000
5  1.2450259 1.0000000 0.025085485 0.9749145


and each solution can be plotted, here 6 per page


prepCN( 6,1,NULL )

par( mfrow=c(3,2) )
for( i in 1:nrow(local) ){
 cat(i,"\n")
 r<-as.numeric(rownames(local)[i] )
 showTumourProfile(copyAr, maxPoints=25000 , flatten=.5 , nlev=20, 
                          xlim=c(0,2) , nx=200, ny=50 , 
                          noise=0.01 , nopoints=T  )
 plotModelPeaks(li6[[r]]$par, selectedPoints=NULL,
                        cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=F )
 legend( 0,1, floor( 1000*li6[[r]]$value)/1000, bg="white"  )
}     



################################################################################

OTHER OPTIONS: subsets 

In the situation where the tumour has more subclones than specified in the 
coverParamSpace function, it is possible that the sp data.frame contains peaks 
that are specific to the unspecified subclone(s).  The coverParamSpace function
has an option (usesubsets) that indicates that only a subset of the peaks 
(randomly chosen) are to be used in the fit. The size of the subset will be at 
least equal to the integer value assigned to the usesubsets argument (NULL by 
default). The value of the objective function is then divided by the percentage
of the rows of sp that actually entered the fit, to make the value of the 
objective function more comparable when different number of peaks are used 
between runs. Note that the subset of peaks only changes between runs of 
coverParamSpace.

cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
    xlim=c(0,2) , nx=200, ny=50 )
axis(1)
points( sp )

li7<-coverParamSpace( sp , optimFct=2, lowerF=c(0), upperF=c(1),  
         Sfrom=.25, Sto=2 , maxc=6 , control=list( maxit=1000  ), 
         usesubsets=3, nrep=10  )
## save(li7, file="Rda/li7.rda") 

Each rep uses a new subset of peaks, listed in, eg, li6[[1]]$subset. 


##

par(mfrow=c( 2,5 ))
allParamSpace<-c()
for( i in 1:length(li7) ){
 allParamSpace<-rbind( allParamSpace, li7[[i]]$paramSpace )
 plot( li7[[i]]$paramSpace[,2], li7[[i]]$paramSpace[,1], pch='.' , 
          main=li7[[i]]$subset , ylim=c(0,2) )
}

##

par(mfrow=c(1,1))
plot( allParamSpace[,2], allParamSpace[,1], pch='.' )

##

cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
        xlim=c(0,2) , nx=200, ny=50 )
axis(1)
points( sp, pch=19  )
text( sp, labels=1:nrow(sp), cex=3 )
plotModelPeaks( par=li7[[3]]$par , selectedPoints=NULL , 
                cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=T, 
                preserveMatPatDiff=T , preserveMaxSubClDiff=T )


################################################################################

OTHER OPTIONS: Sn

Tumor ploidy and cellulaliry are interconnected, in such a way that S*n is 
constant (where n is the proportion of normal cells).  If the constant value of 
S*n can be estimated, one less parameter is needed.  This value can be estimated 
from segments that are LOH. 

cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       xlim=c(0,2) , nx=200, ny=50 )
axis(1)

Let's add points to the graph that correspond to actual segments:

sel<-t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9
# the bigger the segment, the bigger the point
cxcut<- as.integer( cut( t.ar.seg$size[sel], 
         c(10000,100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
points( x<-t.ar.seg$mean[sel], y<-t.ar.seg$p[sel],  pch=21 , 
          col="blue", lwd=3 , cex=cxcut  )
points( t.ar.seg$mean[sel], t.ar.seg$p[sel],  pch=19 ,  
          col="white" , cex= cxcut  - .5 )
points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=21 ,  
          col="blue", lwd=3  , cex=cxcut  )
points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=19 , 
          col="white" , cex=cxcut -.5 )

Select segments (points) that are at LOH. These points will be used to estimate 
the LOH curve and the value for S*n. Here to illustrate we do this manually; in
pipeline.r this is done automatically, but can fail sometimes to come up with
a good estimate. 

lohpoints <-  selectPeaks( cntr, copyAr , getLocal=F ) 
x<-lohpoints$x
ar<-lohpoints$y

I only chose two segments, which is enough to estimate the curve:
 x
[1] 0.3244808 0.6347028
 y
[1] 0.05 0.02

This is the equation of the LOH curve as a function of S and the value x 
representing the scaled read count of a segment:

ARloh<-function( x,S, n ){
 k<-2*(x-S*n)/(S-S*n) 
 return( n/( 2*n+(1-n)*k ) )
}

This is the fit. Since S*n is constant, we can set S to 1.

nnllss<- nls( ar ~ ARloh( x, 1 ,n ) , start=list(n=0.01  ) , 
               upper=list( n=1 ), lower=list(n=0 ) , algo="port" )
Sn<-  summary(nnllss)$coefficients[1,1]

Sn
[1] 0.03098523

The value Sn is also the scaled read count expected in regions where all copies
were deleted in the tumor (where reads only come from the contaminating normal
tissues).

We can plot the LOH curve corresponding to Sn:

x <- seq( Sn, 2, .01 ) 
points( x , arloh<- ARloh( x , 1 , Sn ) , type='l' , lwd=2  )
points( x , 1-arloh , type='l' , lwd=2  )

And we can speed up the search with the Sn argument:

li.sn <-coverParamSpace( sp , optimFct=2  , lowerF=c( 0  ), upperF=c( 1  ),
                         maxc=8, Sn=0.03098523,
                         control=list( maxit=1000 ) )

In the above call, S is not varying as a paramters, but rather determined from
the values of Sn and the estimate of fraction of normal (n=t[1]). 

 li.sn[[1]]$par
[1] 0.61334964 0.05051805





