0        1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890

################################################################################

Read INSTALL.txt for instruction on how to install celluloid11.

################################################################################

Input files used in the tutorial are included in the Data directory. Unzip 
before use. 

The input .wig files that are needed are passed to functions from the HMMcopy 
package.  

See http://compbio.bccrc.ca/software/hmmcopy/ for instructions on how to turn 
bam files into wig files, including wig files for gc- and mappability-related 
data.  Here the files normal.wig, tumour.wig, gc.wig and map.wig contain more 
than the standard chromosomes. Only the autosomal chromosomes are used in 
deriving ploidy (autosomal ploidy) and cellularity. 

The file AR.txt contains the columns:

CHR POS REF_COUNT VAR_COUNT
chr1 754730 31 17
chr1 754813 10 23
chr1 754840 30 10
chr1 754873 16 20
chr1 755955 26 14
...

where each line corresponds to a heterozygous position in the normal,
where REF_COUNT represents the number of reads supporting the reference
allele in the tumour, and VAR_COUNT represents read counts supporting the
other allele. 

################################################################################

The user can run the pipeline.r script. This script will prepare the data and
run an analysis that assumes only one clone is present. From the current 
directory, call:

Rscript pipeline.r Data/tumour.wig Data/normal.wig Data/gc.wig Data/map.wig Data/AR.txt samplename

Details of each steps of data preparation are found in 

 README.2.load_and_show.txt

All data objects are saved in .rda files in the RdaPipeline directory.

The analysis uses the segment-based objective function. The value of this 
objective function is a measure of how distant the tumor genome is from the 
model defined by values of ploidy and cellularity.  The search is done by a call 
of the function coverParamSpace.

Local maxima of the objective function are found by using the optim() function 
in R with 50 randomly chosen starting points.  

Details about the segment-based objective function are found in 

 README.3.select_and_search.txt

Local solutions are found in a graphical output. The first panel represents
the values of the objective function as a function of the parameter n (the 
percentage of normal cells in the sequenced sample). This is to help the 
user judge if the parameter space was well covered. Then the output shows a 
series of contour plots with white dots representing each segment (mirrored 
around 0.5). The red dots and lines represent the copy-number state and allelic 
ratios that are predicted given the estimates of tumour ploidy and cellularity. 
The value of the objective function at the local solution is displayed on the 
top left corner. 

Details about the figures can be found in 

 README.4.display_solutions.txt

In pipeline.r, coverParamSpace makes use of a dimension reduction trick (by 
setting a value to the argument Sn), because cellularity and ploidy are 
interconnected. A first step of this trick consists of estimating "LOH curves", 
which may fail in practive in some datasets, in which case manual interventions 
may be required. The LOH curves represent the lower and upper values of the AR 
expected in regions of LOH, given the copy number state of the tumour in those 
regions. This LOH curves are the two symmetrical black curved lines on the 
output graphs. See the "OTHER OPTIONS: Sn" section in 

 README.5.other_options.txt

The user must decide on a solution based on experience, and copy-number 
segments can then be plotted:

library(celluloid11)
files<-system("ls RdaPipeline/*rda",intern=T); for(f in files){load(f)}
prepCN(12,1,NULL)
ePP<-ePeakPos( S=0.624, t=c(0.073, 0.923), cn=cn  )
tcs<- scaleReadCounts( tc , ePP )
segments<-scaleSegments(t.ar.seg ,  ePP )
segments<-annotateSegments(segments, ePP)
plotSegment( tcs,segments, ar , file="segments_page%1d",device="png",
             width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
             type="cairo" , tlwd=8 ) 
# XY chromosomes can't be annotated due to lack of ar
segmentsXY<-scaleSegments(t.seg ,  ePP )
plotSegment( tcs,segmentsXY, ar=NULL , file="segments_pageXY",device="png",
             width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
             type="cairo" , chr=c("chrX","chrY") , tlwd=8) 








