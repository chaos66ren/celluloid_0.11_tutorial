0        1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890

###############################################################################
load( file="Rda/copyAr.rda")
load( file="Rda/sp.rda") 
###############################################################################

In the following, we set seed values (set.seed()) so that results presented here 
can be reproduced by the user.  Otherwise seeds do not need to be set. 

###############################################################################

The present file describe ways to obtain ploidy, cellularity and subclonal 
proportions with a more hands-on approach. 

A graphic window with a tumour profile must be displayed on screen:

cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50 )
axis(1)

The following function is used to select "peaks" that will enter the analysis. 
If getLocal=T, the function will make nrand (default 100) calls to optim() with
random starting points.  Otherwise the user needs to left clicks peaks to 
manually select them (this is what we recommend). A right click ends the 
function. 

sp <-  selectPeaks( cntr, copyAr , getLocal=F ) 

Alternatively, peaks can be automatically selected with getLocal=T. Then,
calls of the optim function are made, based on randomly selected starting 
points:

sp <-  selectPeaks( cntr, copyAr , getLocal=T , 
                    percentMax=.33, manual=F, nrand=200 , filtersymm=T   ) 


where percentMax raises the "sea level" so that smaller peaks and noise are
ignore. Here it is set to be at 33% of max( cntr$z ) (which has been 
"flattened"). If manual=T then the use is asked to choose additional points.

Note that the graph is symmetrical. Try to only select one of the two 
symmetrical points, but if not the function will do it for you. 

The function returns a data.frame containing the chosen coordinates:

sp         
          x    y
1 0.3119248 0.01
2 0.6065158 0.00
3 0.7886266 0.38
4 0.9011068 0.33
5 1.0778614 0.51
6 1.1956978 0.51
7 1.4956450 0.39

## save( sp, file="Rda/sp.rda") 

## load( file="Rda/sp.rda")

Inspection of sp reveals that a subset of peaks (lines 1,2,4,6,7) are separated 
on the x-axis by multiples of approximately .3.  This may be an indication that 
a single clone is sufficient to explain these peaks (segments contributing to 
these peaks have the same number of copies in all tumour cells), while segments 
contributing to peaks 3 and 5 may have a different number of copies in some 
tumour subclones. 


##############################################################################

In this tutorial, "read counts" usually mean a scaled read count that was gc- 
and mappability-corrected. 

The sequenced sample consist of a mixture of normal cells and tumour cells
that come from, perhaps, different subclones.  Two subclones are different 
with respect to at least one segment that has a different number of copies
between the two.  

Let S be the expected read counts in the sequenced sample in a segment that has 
exactly two copies in all cells (normal and tumour cells). Let t describes the 
percentage of the sequenced cells that are normal and derived from the different
subclones: t=c(t[1],t[2],..,t[nsubclones+1]); sum(t)=1. t[1] is the percentage 
of normal cells in the sequenced sample. 

We define the autosomal ploidy of the sequenced sample (the mixture of tumor and 
normal cells) as the relative abundance of DNA compared to that of a normal 
haploid autosomal genomes and calculate it to be 2/S.

For given S and t, and given copy number in all subclones (normal cells are
assumed to have 2 copies), the peak locations (expected read counts and 
expected allelic fractions) can be predicted from mathematical expressions. 
We aim to find the set of parameters for which the observed peaks are best 
captured by expected peaks.  The objective function to minimize is the total
(Euclidean) distance between the observed peaks and their closest expected 
peaks.

The following function does the search,  which relies on a simulated annealing 
algorithm (GenSA) when optimFct=2.  The lowerF and upperF arguments are the 
lower and upper limits for the t parameter, ignoring the last entry (ie, lower 
and upper bounds for t[1:nsubclone]  ).  The Sfrom and Sto arguments are the 
lower and upper bounds for S. 

set.seed(12345)

li1<-coverParamSpace( sp , optimFct=2, lowerF=c(0), upperF=c(1),  
       Sfrom=.25, Sto=2 , maxc=12 , control=list( maxit=1000  ) )

## save( li1, file="Rda/li1.rda")

The ploidy of the sequenced sample is allowed to go from 2/Sto=1 to 2/Sfrom=8. 
The percentage of normal cells in the sequenced sample is allowed to vary from
0% (lowerF=c(0)) to 100% (upperF=c(1)).  

The control argument is passed to the optimization function used. See ?GenSA for
details (if optimFct=2; see ?optim if optimFct=1).

The argument maxc represents the maximum number of copies that can be found in 
any cell. Set it so that it covers the largest ploidy in the search.

The maxc argument is passed to the function prepCN, a function called from
within coverParamSpace that creates a global data.frame with the name "cn" that 
describes the allowed copy number configurations between the normal cells and 
the tumour cells.

E.g., the call of prepCN(maxc=12) creates the data.frame:

cn

   N T1
1  2  0
2  2  1
3  2  2
4  2  3
5  2  4
6  2  5
7  2  6
8  2  7
9  2  8
10 2  9
11 2 10
12 2 11
13 2 12

The first line represents a configuration where normal cells (column N) have two 
copies and tumour cells (column T1) have 0 copies in a segment.  The next line 
represents a configuration where normal cells have two copies and tumour cells 1 
copy, etc. In the search, a segment is only allowed to have one of these 
configurations. It is from these configurations that expected peak locations are 
derived. 

The function creates a global matrix called paramSpace that contains the value 
of the objective function (first column) and the value of c(S,t) in the other
columns,  at each iteration.  It can be appended to between different runs using 
the flag addToParamSpace=T. Plotting the value of the objective function and the
different parameters can help determine if the parameter space was well covered 
and if other local minima are worth refining:

  plot( as.data.frame(paramSpace), pch='.' )

or in details, for chosen parameters, eg,

  plot( paramSpace[,2], paramSpace[,1], pch='.', xlab="S", ylab="distance" )

Inspection of the latter plot reveals that a local minimum near S=0.3 has not
been well covered by the simulated annealing search.  Either use a larger number 
of iterations, or refine the search around that S value:

set.seed(12345)

li1.refined<-coverParamSpace( sp , optimFct=2, lowerF=c(0), upperF=c(1),  
       Sfrom=.2, Sto=.4 , maxc=12 , control=list( maxit=1000  ),
        addToParamSpace=T )

# save( li1.refined, file="Rda/li1.refined.rda")

As a function of S, the global minimum of the objective function is clear and
no other local minima are close to it. See Figure2.png.

The function returns a list of lists containing among other things the optimal 
parameters and the value of the objective function.  It has length equal to 
the number of starts (controled by nrep). Ie,

li1[[1]]$value 

is the minimum distance between the observed and expected peaks and

li1[[1]]$par

contains the optimal c(S, t[1],..,t[nsubcl] ):

S<-li1[[1]]$par[1]
t<-li1[[1]]$par[ 2:length(li1[[1]]$par ) ]
t<-c( t, 1-sum(t) )

> S
[1] 0.600465
> t
[1] 0.008267668 0.991732332

The list also includes paramSpace: li1[[1]]$paramSpace.

Move straight to README.3.display_solutions.txt to display the solution on top 
of the contour plot. 

Alternatively, since a regular pattern can be seen in sp, the use might want to 
perform the search by only focusing on the equally spaced peaks. This may help
when subclones display substantial differences. 

spsubset<-sp[c(1,2,4,6,7),]

set.seed(12345)

coverParamSpace( spsubset , optimFct=2, lowerF=c(0), upperF=c(1),  
       Sfrom=.25, Sto=2 , maxc=12 , control=list( maxit=1000  ) )


##############################################################################

The search can be refined around this solution to allow for a subclone, to 
capture the peaks 3 and 5:

set.seed(12345)

li2<-coverParamSpace( sp, optimFct=2, lowerF=c(0,0), upperF=c(0.02,1),  
       Sfrom=.59, Sto=0.61 , maxc=6 , maxsubcldiff=0.5,  
       control=list( maxit=1000 ) )

## save( li2, file="Rda/li2.rda")

If >= 1, maxsubcldiff represents the upper bound for the difference in copy 
number between any two subclones.  If < 1, it represents the upper bound for 
the ratio between the minimum and maximum copy number seen across all 
subclones.  Maxc and maxsubcldiff are used in prepCN, a function called from
within coverParamSpace that creates a global data.frame with the name "cn" that 
describes the allowed copy number configurations between all subclones. The
value maxc=6 was chosen based on the one clone solution. 

Eg, prepCN(maxc=6 , nsubcl=2 , maxsubcldiff=.5 ) returns the cn data.frame:

cn

   N T1 T2
1  2  0  0
2  2  1  0
3  2  2  0
8  2  0  1
9  2  1  1
10 2  2  1
15 2  0  2
16 2  1  2
17 2  2  2
18 2  3  2
19 2  4  2
24 2  2  3
25 2  3  3
26 2  4  3
27 2  5  3
28 2  6  3
31 2  2  4
32 2  3  4
33 2  4  4
34 2  5  4
35 2  6  4
39 2  3  5
40 2  4  5
41 2  5  5
42 2  6  5
46 2  3  6
47 2  4  6
48 2  5  6
49 2  6  6

The first line represents a configuration where normal cells have two copies and 
both subclones have 0 copies in a segment.  The next line represents a 
configuration where normal cells have two copies, the first subclone has 1 copy
and the second subclone has 0 copies, etc. 

Note that when maxsubcldiff<1, a copy number 0 is treated as being 1. The above
cn data.frame (created with maxsubcldiff=0.5) describes situations where the 
copy number in a subclone is only allowed to be at most double the copy number 
of any other subclone.

The above call to coverParamSpace returned:

li2[[1]]$par

[1] 0.600151798 0.007227974 0.643795537

S<-li2[[1]]$par[1]
t<-li2[[1]]$par[ 2:length(li2[[1]]$par ) ]
t<-c( t, 1-sum(t) )

> S
[1] 0.6001518
> t
[1] 0.007227974 0.643795537 0.348976489

plot( paramSpace[,4], paramSpace[,1], pch='.', xlab="Subclone 1", 
      ylab="distance", ylim=c(0,.5 ) )

Allowing two tumor populations or subclones, the percentage of normal cells is 
0.7%, the percentage cells in the first subclone is 64.3% and the percentage of 
cells in the second subclone is 34.9%.  See README.3.display_solutions.txt to 
visually inspect the ovelap between the observed and expected peaks. 

###############################################################################

Instead of first finding a one-clone solution and then refining to find a 
two-subclone solution (a stepwise approach that we recommend), the search could 
have been done by directly calling a two-subclone model:

set.seed(12345)

li3<-coverParamSpace( sp, optimFct=2, lowerF=c(0,0), upperF=c(1,1),  
       Sfrom=.25, Sto=2 , maxc=6 , maxsubcldiff=0.5,  
       control=list( maxit=1000   ) )

## save( li3, file="Rda/li3.rda")


Inspection of paramSpace reveals that two solutions are almost equally 
good, the other one near S=0.5.  See Figure3.png. 

plot( paramSpace[,2], paramSpace[,1], pch=19, cex=.5, xlab="S", 
           ylab="distance" , ylim=c(0,1) , xlim=c(.4,.8)  )
abline( h=li3[[1]]$value )

The search can be refined around that value if needed. 

set.seed(12345)

li4<-coverParamSpace( sp, optimFct=2, lowerF=c(0,0), upperF=c(1,1),  
       Sfrom=.475, Sto=.525 , maxc=6 , maxsubcldiff=0.5,  
       control=list( maxit=100  ) )

li4[[1]]$par
[1] 0.500113283 0.007053342 0.594503953

## save( li4, file="Rda/li4.rda")

See README.3.display_solutions.txt to get a visual assessment of the solutions.



