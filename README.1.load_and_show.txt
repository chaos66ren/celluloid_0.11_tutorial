0        1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890

###############################################################################


In this tutorial, "read counts" usually mean a scaled read count that was gc- 
and mappability-corrected. A "peak" is a local maximum in a contour plot 
(as in a topographical map), either an observed peak or an expected peak. 

library( celluloid )


###############################################################################

Define the files

tumourWigFile<-"Data/autosome_tumour.wig"
normalWigFile<-"Data/autosome_normal.wig"
gcWigFile<- "Data/autosome_gc.wig" 
mapWigFile<- "Data/autosome_map.wig"
arFile<- "Data/AR.txt"

Load the data from normal, using function from HMMcopy:

n<- wigsToRangedData( normalWigFile, gcWigFile, mapWigFile )

Perform GC-content correction:

nc<-gcCorrect(n) 

Segment the data, using functions from the copynumber package. At the same time, 
the function also calculates the mean mappability in each segment.

n.seg <- segmentSeqData( nc, k=50 , maskmap = 0.8 )


head(n.seg)

  sampleID chrom arm start.pos  end.pos n.probes   mean   meanmap
1       rc  chr1   p    752001 16858000    12151 0.9784 0.8597475
2       rc  chr1   p  16858001 16882000       11 1.7639 0.7696937
3       rc  chr1   p  16883001 16976000       50 3.3836 0.7150340
4       rc  chr1   p  16976001 17113000       87 2.0013 0.8465564
5       rc  chr1   p  17113001 25585000     6845 0.9930 0.8930962
6       rc  chr1   p  25595001 25666000       50 0.5545 0.8501373

In the above segmentSeqData call, bins with mappability less than 0.8 are 
treated as missing data, which explains why the genomic intervals in the output
are not contiguous.

##     save( n.seg, file="Rda/n.seg.rda")

##     save( nc, file="Rda/nc.rda")

##     load("Rda/n.seg.rda")

##     load("Rda/nc.rda")

Repeat for tumor data:

t <- wigsToRangedData( tumourWigFile, gcWigFile, mapWigFile )
tc<-gcCorrect( t )
t.seg <- segmentSeqData( tc , k=50  )

head(t.seg)

  sampleID chrom arm start.pos  end.pos n.probes   mean   meanmap
1       rc  chr1   p    752001 16884000    12199 0.8816 0.8596063
2       rc  chr1   p  16919001 16982000       55 3.0496 0.8996408
3       rc  chr1   p  16982001 17113000       81 1.7862 0.8403646
4       rc  chr1   p  17113001 25585000     6854 0.8889 0.8930962
5       rc  chr1   p  25595001 25666000       50 0.3749 0.8501373
6       rc  chr1   p  25666001 72760000    38691 0.9036 0.9009475

##    save( t.seg, file="Rda/t.seg.rda")

##    save( tc, file="Rda/tc.rda")

##    load( "Rda/t.seg.rda")

##    load( "Rda/tc.rda")



###################
 
Identify segments that do not appear normal in the normal. 

First intersect the segments between tumor and normal:

t.n.seg<-intersectSegments( t.seg, n.seg )

head(t.n.seg)

   sampleID chrom arm start.pos  end.pos     size   mean   meanmap mean.1 meanmap.1
1        NA  chr1   p    752001 16858000 16106000 0.8816 0.8596063 0.9784 0.8597475
2        NA  chr1   p  16858001 16882000    24000 0.8816 0.8596063 1.7639 0.7696937
3        NA  chr1   p  16882001 16883000     1000 0.8816 0.8596063     NA        NA
4        NA  chr1   p  16883001 16884000     1000 0.8816 0.8596063 3.3836 0.7150340
31       NA  chr1   p  16884001 16919000    35000     NA        NA 3.3836 0.7150340
21       NA  chr1   p  16919001 16976000    57000 3.0496 0.8996408 3.3836 0.7150340

Columns mean and meanmap are from t.seg (the first argument), mean.1 and 
meanmap.1 are from n.seg (the second argument).  Note that NAs are introduced 
because of the maskmap argument of segmentSeqData; see the above note. 

t.n.seg<-t.n.seg[ !is.na( t.n.seg$mean) & !is.na( t.n.seg$mean.1),]

A segment in the normal whose mean value is an outlier compared to all other 
segments will be treated as a non-normal segment.  In the following we define
an outlier through calculation of the interquartile range, defined only from
"large" segments with high mappability:

sel<-n.seg$end.pos-n.seg$start.pos > 100000 & n.seg$meanmap>.8 
bp<-boxplot( n.seg$mean[sel], range=3, plot=F  )

The "normal" range of mean values:

nrange<-c(bp$stats[1,1], bp$stats[5,1] )

This range can then be applied to all segments:

mask<- t.n.seg$mean.1>nrange[2] | t.n.seg$mean.1<nrange[1]

Revert to only tumour segments, and mark segments that are not normal in the
normal:

t.seg.mask<-t.n.seg[,1:8]

t.seg.mask$mask<- mask



##    save(t.seg.mask, file="Rda/t.seg.mask.rda") 

##    load("Rda/t.seg.mask.rda") 




###############################################################################

Reading the allelic ratio file.

ar<-read.table(arFile, head=T )

head(ar) 
   CHR    POS REF_COUNT VAR_COUNT
1 chr1 754730        31        17
2 chr1 754813        10        23
3 chr1 754840        30        10
4 chr1 754873        16        20
5 chr1 755955        26        14
6 chr1 758555        11        11

##    save( ar, file="Rda/ar.rda") 

##    load("Rda/ar.rda") 

Segmenting the AR data so that LOH regions can be distinguished from normal 
regions:

ar.seg<- segmentAR( ar, tc ) 

head(ar.seg)
  sampleID chrom arm start.pos   end.pos n.probes meanar
1       ar  chr1   p    754730 120529959    77804 0.3293
2       ar  chr1   q 145382191 249220235    71222 0.4226
3       ar  chr2   p     14453  90266876    73711 0.4265
4       ar  chr2   q  97919691 243056809   100644 0.4266
5       ar  chr3   p     60915  90498746    67190 0.3312
6       ar  chr3   q  93514708 197848857    70480 0.3312

In the above, allelic ratios are calculated at each SNP as the proportion of 
reads supporting the reference allele if it is < 0.5. Otherwise it is calculated
as the proportion of reads supporting the alternate allele. This ratio will thus
always be < 0.5, even if an equal number of maternal and paternal chromosomes
are represented. The meanar column reports the average allelic ratio over all 
heterozygous positions in the segment.

##     save( ar.seg, file="Rda/ar.seg.rda")

##     load("Rda/ar.seg.rda")

t.ar.seg <- intersectSegments( t.seg.mask, ar.seg  )

head(t.ar.seg)
   sampleID chrom arm start.pos  end.pos     size   mean   meanmap  mask meanar
1        NA  chr1   p    752001   754729     2729 0.8816 0.8596063 FALSE     NA
11       NA  chr1   p    754730 16858000 16103271 0.8816 0.8596063 FALSE 0.3293
2        NA  chr1   p  16858001 16882000    24000 0.8816 0.8596063  TRUE 0.3293
4        NA  chr1   p  16882001 16883000     1000     NA        NA    NA 0.3293
41       NA  chr1   p  16883001 16884000     1000 0.8816 0.8596063  TRUE 0.3293
6        NA  chr1   p  16884001 16919000    35000     NA        NA    NA 0.3293


In the above, contiguous segments that display the same total number of copies, 
but different number of maternal/paternal chromosomes (as in LOH) should be
distinguished. 

Removing the NAs that were introduced:

t.ar.seg<-t.ar.seg[ !apply( is.na( t.ar.seg[,c("mean","meanar")] ), 1, any),]

As noted above, allelic ratios in meanar are always < 0.5. The following 
function addresses this by formally estimating the ratio of maternal/paternal
chromosomes, by modeling the number of reads using Poisson distributions. The
estimates are added in a column named p.

t.ar.seg <-arInSeg( t.ar.seg, ar,  tumourrangedata=tc , minhet = 50 )

Estimates are provided if a minimum of 50 SNPs is found in the segment.

t.ar.seg[26:35,]
     sampleID chrom arm start.pos   end.pos      size   mean   meanmap meanar         p
251        NA  chr2   p  34759001  89630000  54871000 1.2016 0.9138349 0.4265 0.4822189
26         NA  chr2   p  89831001  89891000     60000 2.6910 0.9030957 0.4265 0.4999497
27         NA  chr2   p  89891001  90266876    375876 1.2058 0.2781003 0.4265 0.4999498
3111       NA  chr2   q  97919691 133001000  35081310 1.2013 0.8929978 0.4266 0.4725638
321        NA  chr2   q 133001001 133115000    114000 2.1858 0.8956735 0.4266 0.4429502
323        NA  chr2   q 133117001 133118000      1000 2.1858 0.8956735 0.4266        NA
331        NA  chr2   q 133124001 243056809 109932809 1.2017 0.9209395 0.4266 0.4768847
341        NA  chr3   p     60915  75706000  75645086 0.9055 0.9161375 0.3312 0.3359412
3411       NA  chr3   p  75706001  75762000     56000 0.9055 0.9161375 0.3312        NA
342        NA  chr3   p  75762001  90498746  14736746 0.9055 0.9161375 0.3312 0.3356679


##     save( t.ar.seg, file="Rda/t.ar.seg.rda")

##     load("Rda/t.ar.seg.rda")

###############################################################################

This function creates the object used to draw a contour plot. It pairs the 
allelic ratio of each SNP with the $mean of the segment it belongs to.

mask<- t.ar.seg$mask | is.na(t.ar.seg$mask)
copyAr<-  prepCopyAr( t.ar.seg[ !mask ,], ar,  tc  )

##   save(copyAr, file="Rda/copyAr.rda")

##   load("Rda/copyAr.rda")

This function plots the autosomal-wide copy number profile of the tumour. Each 
peak (or pair of peaks since the graph is reflected around AR=0.5) corresponds 
to a specific copy-number state that summarizes both the total copy number in 
the tumour cells (on the x-axis, once appropriately scaled) and the ratio of 
relative abundance of maternal and paternal copies (on the y-axis, once 
contamination from normal tissues – or tumor cellularity – is accounted for).

cntr<-showTumourProfile( copyAr , flatten=.25 , nlev=20 , noise= 0.01 , 
                        maxPoints=200000 )


###############################################################################





