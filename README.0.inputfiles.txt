Input files used in the tutorial are included in the Data directory. Unzip
before use. 

The input .wig files that are needed are passed to functions from the 
HMMcopy package.  

See http://compbio.bccrc.ca/software/hmmcopy/ for instructions on how to 
turn bam files into wig files, including wig files for gc- and 
mappability-related data.  Here the files autosome_normal.wig,
autosome_tumour.wig, autosome_gc.wig and autosome_map.wig only 
contain autosomal chromosomes. Other chromosomes can be included in 
these files, but they will be discarded. 

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


