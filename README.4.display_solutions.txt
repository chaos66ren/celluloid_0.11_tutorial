0        1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890

###############################################################################
load( file="Rda/copyAr.rda")
load( file="Rda/sp.rda") 
load( file="Rda/li1.rda")
load( file="Rda/li1.refined.rda") 
###############################################################################

Plotting the solution in li1 (ONE CLONE):

If cn does not already exist, set it using the parameters used in the search:

prepCN(12) 

Given parameters S and t (with sum(t)=1) and copy number configurations cn, the 
expected peak locations can be obtained with a call to 

ePP <- ePeakPos( S= 0.6260699, t=c(0.02222, 0.97778), cn=cn  )

or, alternatively, using the list returned by coverParamSpace:

ePP <- ePeakPos( par=li1[[1]]$par  , cn=cn  )

 head(ePP)

  m0 p0 m1 p1          x         ar
1  1  1  0  0 0.01391417 0.50000000
2  1  1  0  1 0.31999205 0.02174143
3  1  1  1  0 0.31999205 0.97825857
4  1  1  0  2 0.62606993 0.01111231
5  1  1  1  1 0.62606993 0.50000000
6  1  1  2  0 0.62606993 0.98888769


Here the copy number configurations in cn is broken down into (say) maternal and
paternal copy numbers (the terms maternal and paternal are used only to 
distinguish the two chromosomes; the actual parental origin can not
be derived). m0 and p0 are used for normal cells, m1 and p1 are used for the 
tumour cells. 

To display the expected peak location onto a "showTumourProfile" graph, the user
can call:

showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50  )

prepCN(12) # if cn does not exist 

plotModelPeaks( S= 0.6260699, t=c(0.02222, 0.97778), 
                selectedPoints=sp,cn=cn, epcol="red",epcex=1,eplwd=3 , 
                addlabels=T )

or 


plotModelPeaks( par=li1[[1]]$par , selectedPoints=sp , 
                cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=T)

See Figure4.png.

The ePP data.frame above is also returned by plotModelPeaks when assigned (invisible otherwise), as in: 

ePP<-plotModelPeaks( par=li1[[1]]$par , selectedPoints=sp , 
                cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=T, 
                preserveMatPatDiff=T , preserveMaxSubClDiff=T )

The function plots the (x,ar) points in ePP, and labels the x-axis to represent 
integer copy number values in the tumor cells. The function is described in more 
details below.

Note that peaks 3 and 5 are not captured, these might represent segments that 
have different copy numbers in different subclones of the tumour. 

###############################################################################

Plotting the solution in li2 (TWO SUBCLONES):

If cn does not already exist, set it using the parameters used in the search:

prepCN(6,2,0.5) 

Given parameters S and t (with sum(t)=1) and copy number configurations cn, the 
expected peak locations can be obtained with a call to 

ePP <- ePeakPos( par=li2[[1]]$par, cn=cn  )

head(ePP)

  m0 p0 m1 p1 m2 p2          x         ar
1  1  1  0  0  0  0 0.01300897 0.50000000
2  1  1  0  1  0  0 0.21117743 0.03080103
3  1  1  1  0  0  0 0.21117743 0.96919897
4  1  1  0  2  0  0 0.40934590 0.01588994
5  1  1  1  1  0  0 0.40934590 0.50000000
6  1  1  2  0  0  0 0.40934590 0.98411006

Here the copy number configurations in cn is broken down into (say) maternal and
paternal copy numbers (the terms maternal and paternal are used only to 
distinguish the two chromosomes; the actual parental origin can not
be derived). m0 and p0 are used for normal cells, m1 and p1 are used for the 
first subclone, m2 and p2 for the second, etc.

To display the expected peak location onto a "showTumourProfile" graph, the user
can call:

showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50  )

prepCN( 6,2,.5) 

plotModelPeaks( par=li2[[1]]$par , selectedPoints=sp , 
                cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=T, 
                preserveMatPatDiff=T , preserveMaxSubClDiff=T )

See Figure5.png.

The function plotModelPeaks calls the function ePeakPos, as described above.
In addition to plotting the $x and $ar values of ePP as points (with control 
parameters epcol, epcex, eplwd), axis are labeled according to the copy numbers 
found in subclones. For example, a label taking the value 2.3 indicates that
segments corresponding to peak found at that location on the x-axis have 2 
copies in the first subclone and 3 in the second. If addLabels=T, then labels 
will be assigned to points from ePP that are closest to points in sp.  These 
labels further indicate the number of parental chromosomes in each 
subclone, as in 20/21, which indicates that segments contributing to this
particular peak have 2 maternal (say) copies in the first subclone and 0 paternal, 
and 2 maternal plus 1 paternal copies in the second (we use maternal or paternal 
solely to distinguish the two chromosomes; the parental origin can not 
be inferred). Labels can also be manually added with a call of the function

addLabels( ePP, manual=T)

in which case the user need to left-click an expected peak (the red points in
Figure5.png) to add a label and right-click the window to exit the function.

Note that one of the observed peak (line 5, or sp[5,]) is annotated with a 4.3 
profile (4 copies in the first subclone, 3 copies in the second). We selected a 
single peak near the y=0.5 line, but no such peak is to be expected there due to 
the fact that a subclone has a an uneven number of copies. That peak is likely 
made up of two peaks (with labels 22/12) that are so close one another that they  
merge into a single one.  

The user can compare the parsimony of the above solution with the one obtained
in li4:

showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50  )

prepCN( 6,2,.5) 

load("Rda/li4.rda") 

plotModelPeaks( par=li4[[1]]$par , selectedPoints=sp , 
                cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=T, 
                preserveMatPatDiff=T , preserveMaxSubClDiff=T )


###############################################################################

To plot copy-number-segment graphs, the read counts first need to be re-scaled.
The ePP data.frame displays both integer copy-number counts (the m and p 
columns) and the corresponding expected values for the read count.  ePP is used
to re-scale read counts: 

load( "Rda/tc.rda")
tcs<- scaleReadCounts( tc , ePP )

A column named icopy was added.

Then the data can be re-segmented or the earlier segments re-scaled directly. 

load("Rda/t.ar.seg.rda")
segments<-scaleSegments(t.ar.seg ,  ePP )

A column named imean was added. 

Just like the addLabels adds labels to selected peaks, annotations can be added
to individual segments using the function:

segments<-annotateSegments(segments, ePP)

In which case the segment's mean and p values are annotated to the closest point
in ePP.  These labels are added in a column named labels. 

Then segments can be plotted:

# type "cairo" passed to png, to allow for transparency. May not be available.
# width and height are here expressed in pixels (the default of png)
plotSegment( tcs,segments, ar , file="segments_page%1d",device="png",width=2*960,height=2*1320, 
  cex.axis=2, cex.main=2, cex.lab=2, type="cairo", chr=paste( "chr",1:8 , sep="")  ) 




###############################################################################
