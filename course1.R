# https://genomicsclass.github.io/book/ chapters 11, 12, 13

install.packages("usethis")
library(usethis)
use_git_config(user.name = "Ezra Baraban", user.email = "ezra.baraban@gmail.com")
use_git() 
create_github_token()
#install.packages("gitcreds")
library(gitcreds)
gitcreds_set()

BiocManager::version()

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install(c("genefilter","geneplotter"))

library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)

# get help through the documentation
help.start()
?mean
help(mean)
help(package="genefilter")

# inspect objects, classes and methods
library(Biobase)    # load one of the core Bioconductor packages
?ExpressionSet
?"ExpressionSet-class"
methods(class = ExpressionSet)

# inspect the source code for functions and methods
read.csv
plotMA
showMethods("plotMA")
getMethod("plotMA","data.frame")

# vignettes teach you how to use various functions in a package
vignette(package="Biobase")
vignette("ExpressionSetIntroduction")
browseVignettes(package="Biobase")

# report key details about your R session
sessionInfo()

BiocManager::version()

###
BiocManager::install(c("genefu",
                       "COPDSexualDimorphism",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression"))

#explore mammaprint genes
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

#explore COPD
library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)

library(gwascat)
data(ebicat_2020_04_30)
ebicat_2020_04_30

library(tissuesGeneExpression)
data(tissuesGeneExpression)

head(e[,1:5])
table(tissue)

BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))

#where on the genome does the nuclear protein known as ESRRA (estrogen related receptor alpha) bind?

#https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A47716445%2D49069197&hgsid=2306624450_Dq1Fa5h5yS1rlAk4AI7m6WXqXOQV

library(devtools)

#load liver cell dataset
install_github("genomicsclass/ERBS")
library(ERBS)
data(HepG2)
data(GM12878)
HepG2
GM12878

#These are the reported ESRRA binding sites obtained for a ENCODE ChIP-seq experiment on two cell lines

class(HepG2)
values(HepG2)

HepG2

#First 10 values 
HepG2[1:10,]

#to get all values on chr20
chr = seqnames(HepG2) #gives chr's for each entry
as.character(chr)
class(chr)
#this is an Rle -> saves space; by ordering by chr; don't have to save chr1 a million times, just need the number of times chr 1 repeats;
table(chr) #gives top chr's
table(chr)[1:24] #gets rid of noise chr values, ie restricts to 23XY

HepG2[chr=="20",]  

x = HepG2[order(HepG2)] #first by chr, then by location;
#key when want to plot anything along genome, order is key;

seqnames(x) #now on ordered chr's; shows efficiency of Rle class; vs. as.character(seqnames(x))
 
###
Interval ranges: IRanges
browseVignettes("IRanges") #pulls up PDF

ir <- IRanges(5, 10)
ir
IRanges(5, width = 6) #same thing;
#all IRanges have start, end, and width;
start(ir)
end(ir)
width(ir)

#make multiple IRanges at once:
ir <- IRanges(start=c(3,5,17), end=c(10, 8, 20))
ir
length(ir) #ie there are 3 of them;
start(ir)

#intra-range methods
shift(ir, -2) #ie shift left by 2;

narrow(ir, 2)
#relative to the start of the iranges; instead start at the 2nd bp; ie clips off the ends of iranges
#ie chop off the first letter of the iranges

ir <- IRanges(1, 5)
ir
ir <- narrow(ir,2) # aka narrow(ir, start=2)
ir

#or clip off from the end ie:
ir <- narrow(ir, end=2)
ir

### flank; will ELONGATE the IRanges; flank(x, width, start, end, both = FALSE)
ir <- IRanges(1, 10)
ir
flank1 <- flank(ir, width = 3, start=T, both=FALSE) #gives you 3mer of flanking from the START
flank1
flank2 <- flank(ir, width = 3, start=F, both=FALSE) #3mer from the END
flank2
flank3 <- flank(ir, width = 3, start=T, both=T) #3mer from both sides OF THE START
flank3
flank4 <- flank(ir, width = 3, start=F, both=T) #3mer from both sides OF THE END
flank4

ir
ir * 2
ir * -2
ir + 2
ir - 2
resize(ir, 1)

# set up a plotting window so we can look at range operations
plot(0,0,xlim=c(0,23),ylim=c(0,13),type="n",xlab="",ylab="",xaxt="n")
axis(1,0:15)
abline(v=0:14 + .5,col=rgb(0,0,0,.5))

# plot the original IRange
plotir <- function(ir,i) { arrows(start(ir)-.5,i,end(ir)+.5,i,code=3,angle=90,lwd=3) }
plotir(ir,1)

#https://genomicsclass.github.io/book/pages/bioc1_igranges.html

polygon(c(start(ir)-.5,start(ir)-.5,end(ir)+.5,end(ir)+.5),c(-1,15,15,-1),col=rgb(1,0,0,.2),border=NA)

# draw the different ranges
plotir(shift(ir,-2), 2)
plotir(narrow(ir, start=2), 3)
plotir(narrow(ir, end=5), 4)
plotir(flank(ir, width=3, start=TRUE, both=FALSE), 5)
plotir(flank(ir, width=3, start=FALSE, both=FALSE), 6)
plotir(flank(ir, width=3, start=TRUE, both=TRUE), 7)
plotir(ir * 2, 8)
plotir(ir * -2, 9)
plotir(ir + 2, 10)
plotir(ir - 2, 11)
plotir(resize(ir, 1), 12)

text(rep(15,12), 1:12, c("ir","shift(ir,-2)","narrow(ir,start=2)",
                         "narrow(ir,end=5)",
                         "flank(ir, start=T, both=F)",
                         "flank(ir, start=F, both=F)",
                         "flank(ir, start=T, both=T)",
                         "ir * 2","ir * -2","ir + 2","ir - 2",
                         "resize(ir, 1)"), pos=4)

###InTER-Ranges Fxns
ir <- IRanges(start=c(1,7,10), end=c(2, 9, 12))
range(ir) #gives you the min/max; including gaps; vs:
reduce(ir) #gives you the ranges covered, but excluding gaps
gaps(ir) #the gaps - what is omitted from reduce(ir)
disjoin(ir) #makes an IRanges object disjoint by fragmenting it into the widest ranges where the set of overlapping ranges is the same
#ie same coverage, but no overlap

###Genomic ranges: GRanges
#https://genomicsclass.github.io/book/pages/bioc1_igranges.html
library(GenomicRanges)
#With an IRange, a chromosome name, and a strand, we can be sure we are uniquely referring to the same region and strand of the DNA

gr <- GRanges("chrZ", IRanges(start=c(5,10),end=c(35,45)),
              strand="+", seqlengths=c(chrZ=100L))
gr
#ie chrZ is set to 100bp long;

genome(gr) <- "hg20"
gr

seqnames(gr)
seqlengths(gr)
gr
shift(gr, 10) #shift by 10
shift(gr, 80) #will push off the chromosome, ie past telomere; error;

#If we trim the ranges, we obtain the ranges which are left, disregarding the portion that stretched beyond the length of the chromosome:
#ie stops at 100
  trim(shift(gr, 80))

  #metadata columns aka:
  mcols(gr)
  
  #can add metadata
  mcols(gr)$value <- c(-1,4)
  gr
  
  #delete metadata
  mcols(gr)$value <- NULL
  gr
  
#GrangesList, ie all exons in a gene; ie list of GRanges
  gr2 <- GRanges("chrZ",IRanges(11:13,51:53)) 
  gr2
  grl <- GRangesList(gr, gr2)
  grl
  
length(grl)
grl[1]
mcols(grl) #ie can have metadata w/i the GRangesList, one for each GRanges in the list
#ie 1 per exon;
mcols(grl) <- c(1,3)
mcols(grl)

gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*")
gr2 <- GRanges("chrZ",IRanges(c(19,33),c(38,35)),strand="*")

fo <- findOverlaps(gr1, gr2)
fo

#ie the 3rd range of gr1 intersected with the 1st range of gr1; and so forth;

gr1 %over% gr2 #
gr2 %over% gr1 #ie T/F vector of if overlaps or not in each range;
#can conveniently use this to subset -> get ranges that have overlap
gr1
gr1[gr1 %over% gr2] #ie subset of ranges in gr1 that overlap w gr2

#RLE; compact way of storing vectors w a lot of repetitive sequences;
r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20)))
r
str(r)
as.numeric(r) #unwinds the RLE

v <- Views(r, start=c(4,2), end=c(7,6))
v
#views = way of querying RLE or other objects;
#The great benefit of Views is when the original object is not stored in memory, in which case the Views object is a lightweight class which helps us reference subsequences, without having to load the entire sequence into memory

#https://genomicsclass.github.io/book/pages/bioc1_igranges.html


#Operating on GRanges
#https://genomicsclass.github.io/book/pages/bioc1_grangeOps.html

#disjoin(ir) = same coverage; but NO overlaps; never cross any endpoints of original set;
#gaps ~ ie introns

library(GenomicRanges)

ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))

gir = GRanges(seqnames="chr1", ir, strand=c(rep("+", 4), rep("-",3)))
gir

par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlim=c(0,60)) 
plotGRanges(resize(gir,1), xlim=c(0,60),col="green") #plots txn start sites, and will honor the strandedness of each range; ie will go to R if (-) strand;
plotGRanges(flank(gir,3), xlim=c(0,60), col="purple") #ie 3bp upstream ~promoter
plotGRanges(flank(gir,2,start=FALSE), xlim=c(0,60), col="brown") #downstream 2bp promoter

# Finding Overlaps
#Find genes close to our binding sites and annotate those genes;

# load packages
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ERBS)

# load ESRRA ChIP data
data(HepG2) #liver cell line
data(GM12878) #B-cell line

browseVignettes("GenomicRanges") #loads PDF

#find sites in both cell lines - ie ER sites common to both;
# for each row in query, return overlapping row in subject
res = findOverlaps(HepG2, GM12878)
class(res) #Hits object; ie 1st item in HepG2 matches 12th in GM12878;
res
HepG2[1]
GM12878[12]
#ie they do overlap;
#to extract the hits

# ranges from the query for which we found a hit in the subject
index = queryHits(res)
index
erbs = HepG2[index,] #ie ER binding sites = HepG2 subset of the hits that overlap w HM12878
erbs

# extract only the ranges, omit all the metadata
granges(erbs)
erbs

# Genes as GRanges
# so have list of overlaps; ie ESR binding sites;
# load up defined human genes
library(Homo.sapiens)
ghs = genes(Homo.sapiens) #list of all ~23k genes; note that if gene is on + strand it goes from smaller to larger number w txn; vs. if on (-) strand, starts at max and goes to min number of the IRanges;
#ie note strand tells you which is the START and which is the END;
ghs

# learn about the precede function (and related functions like follow)
?precede #note that follow() is the opposite fxn;

# for each range in erbs, find the closest preceding range in ghs
index = precede(erbs, ghs) #ie query ERBS vs. human genome; gives the gene ID's
#for each entry in ERBS, what's the nearest gene that PRECEDES ie is upstream of it;

ghs[index[1:3]] #1st three human genome genes that are upstream of ER binding sites;
erbs[1:3]    # note result is strand-aware
ghs[index]

# distance between binding sites and nearest preceding genes
distance(erbs, ghs[index]) #ie between ER sites and adj genes

# find transcription start site nearest to each binding site
tssgr = resize(ghs, 1) #compress all genes
tssgr

# distance between binding site and nearest TSS
d = distanceToNearest(erbs, tssgr)
queryHits(d) #to get Hits into object we can work with;
dists = values(d)$distance #so these are the distances bw ER sites and adj genes;
hist(dists, nc=100, xlim=c(0,100000))

index = subjectHits(d)[dists < 1000] #genes that are closest to the binding sites;
index

# basics of DNAStrings
dna <- DNAString("TCGAGCAAT")    # define a DNAString
dna
length(dna)    # number of bases in a DNAString
DNAString("JQX")    # error - invalid bases
DNAString("NNNACGCGC-TTA-CGGGCTANN")    # valid sequence with unknowns and gaps
dna[4:6]    # extract a substring
as.character(dna)    # convert DNAString to character

# basics of DNAStringSets
set1 <- DNAStringSet(c("TCA", "AAATCG", "ACGTGCCTA", "CGCGCA", "GTT", "TCA"))    # define a DNAStringSet
set1
set1[2:3]    # extract subset of sequences
set1[[4]]    # extract one sequence as a single DNAString
length(set1)    # number of DNAstrings in set
width(set1)    # size of each DNAString
duplicated(set1)    # detect which sequences are duplicated
unique(set1)    # keep only unique sequences
sort(set1)

dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")

# analyze DNAStrings
letterFrequency(dna_seq, "A")    # count A in sequence
letterFrequency(dna_seq, "GCTA")    # count G or C in sequence
length(dna_seq)
dinucleotideFrequency(dna_seq)    # frequencies of all dinucleotides
trinucleotideFrequency(dna_seq)    # frequencies of all trinucleotides

# convert DNAStrings
reverseComplement(dna_seq)    # find reverse complement
translate(dna_seq)    # amino acid translation

#Matching and Counting with Biostrings
# count and match on individual Biostrings
dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")
dna_seq
countPattern("CG", dna_seq)    # pattern "CG" occurs 5 times
matchPattern("CG", dna_seq)    # locations of pattern "CG"
start(matchPattern("CG", dna_seq))    # start locations of the pattern
matchPattern("CTCTTTTAAAAAAACGCTACTACCATGTGT", dna_seq)    # match patterns of any length

# check for pattern and its reverse complement
countPattern("TAG", dna_seq)
countPattern(reverseComplement(DNAString("TAG")), dna_seq)

# count and match on sets of Biostrings
set2 <- DNAStringSet(c("AACCGGTTTCGA", "CATGCTGCTACA", "CGATCGCGCCGG", "TACAACCGTACA"))
set2
vcountPattern("CG", set2)    # CG counts for entire DNAStringSet, ie VECTOR count pattern;

vmatchPattern("CG", set2) #Unlike vcountPattern(), which simply counts occurrences, vmatchPattern() provides detailed information about where the pattern appears in each sequence.
#ie gives iRanges with all the matches; not just counting them;

vmatchPattern("CG", set2)[[1]]    # access matches for the first element of the DNAStringSet

#Getting the sequence of regions
library(ERBS)
data(HepG2)
HepG2 #where ER binds from ChipSeq

# load and inspect human reference genome
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

# extract chromosome 17 sequence
c17 = Hsapiens$chr17
c17

?getSeq
class(Hsapiens)
showMethods("getSeq")

# collection of DNA strings with ChIP-seq binding peaks
hepseq = getSeq(Hsapiens, HepG2)
hepseq #actually get the sequence from a Granges;
length(HepG2)    # same number of sequences
width(HepG2)[1:5]    # widths match

# collection of shifted DNA strings with no relationship to binding sequences - essentially random
rhepseq = getSeq(Hsapiens, shift(HepG2, 2500))

# count occurrences of a motif in DNA sequences
mot = "TCAAGGTCA"
?vmatchPattern
vcountPattern(mot, hepseq)

# consider both forward matches and reverse complement matches 
sum(vcountPattern(mot, hepseq))    # forward pattern match
sum(vcountPattern(mot, reverseComplement(hepseq)))    # reverse pattern match

## compare motif occurrence in binding peak to random upstream sequences
# count of motifs in binding peaks
sum(vcountPattern(mot, hepseq)) +
  sum(vcountPattern(mot, reverseComplement(hepseq)))
# count of motifs in randomly selected regions of equal length
sum(vcountPattern(mot, rhepseq)) +
  sum(vcountPattern(mot, reverseComplement(rhepseq)))

# for real analysis, use MotifDb package, probabilistic binding packages like MEME and FIMO

# Management of genome-scale data

#ExpressionSet structure for microarray data and the SummarizedExperiment structure for NGS data.

BiocManager::install(c("Biobase",
                       "GEOquery",
                       "genomicsclass/GSE5859Subset",
                       "affy",
                       "hgu95acdf",
                       "genefilter",
                       "parathyroidSE",
                       "airway",
                       "pasillaBamSubset",
                       "Rsamtools",
                       "GenomicAlignments",
                       "ArrayExpress",
                       "NGScopyData",
                       "AnnotationDbi"))

#Bioconductor Infrastructure: ExpressionSet (arrays) and SummarizedExperiment (NGS- rows are GRanges)
#ExpressionSet, is w/i Biobase library

library(Biobase)
library(GEOquery)

geoq <- getGEO("GSE9514")    # download a microarray dataset from GEO
names(geoq)    
e <- geoq[[1]]    # extract ExpressionSet
e

# exprs gives matrix of microarray values
dim(e)    # number of features and samples in ExpressionSet
ncol(e)   # no. of samples
nrow(e)

View(exprs(e))[1:3,1:3]
head(exprs(e))[,1]    # first column
exprs(e)[1,]    # first row
exprs(e)["10000_at",]    # can also index by name
rownames(e)[1]    # row names are probe sets
dim(exprs(e))    # rows are features, columns are samples

# pData gives phenotype data (sample information)
pData(e)[1:3,1:6]
names(pData(e))
pData(e)$characteristics_ch1    # column in GEO to describe experimental state/condition

as.numeric(factor(pData(e)$characteristics_ch1))    # help see replicates of each state
dim(pData(e))    # rows of pData correspond to columns of exprs
dim(e)

# fData gives feature data (probe information)
fData(e)[1:3,1:3]
dim(fData(e))    # rows of fData correspond to rows of exprs
names(fData(e)) 
head(fData(e)$"Gene Symbol")
head(rownames(e))

# additional annotation tied to ExpressionSet
experimentData(e)
annotation(e)

#Reading Microarray Raw Data: Single-Color Arrays
# code edited to set your personal working directory
wd <- getwd()
wd
datadir <- paste0(wd, "/rawdata-master")    # downloaded files, after unzipping
datadir
basedir <- paste0(datadir, "/celfiles")
basedir
setwd(basedir)
library(affy)

tab <- read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)
View(tab)
rownames(tab) <- tab$filenames
tab
fns <- list.celfiles(basedir)
fns
fns %in% tab[,1] ##check
ab <- ReadAffy(filenames = file.path(basedir, tab[,1]), phenoData=tab) #AffyBatch object

dim(pm(ab))
dim(pData(ab)) #phenotype info
rownames(ab) #all the probes
colnames(pm(ab)) #all the samples
annotation(ab) #defines the geneset, can translate to gene names;

e <- rma(ab)    # preprocess probe-level data into gene-level data
#perform the Robust Multi-array Average (RMA) normalization method for Affymetrix microarray data

ejust <- justRMA(filenames=tab[,1],phenoData=tab)    # read and process data to gene-level in one command
dim(ejust)

#Reading Microarray Raw Data: Two-Color Arrays / limma

# datadir defined in previous video
library(limma)
library(rafalib)
basedir <- paste0(datadir, "/agilent")
basedir
setwd(basedir)
targets <- readTargets("TargetBeta7.txt")
RG <- read.maimages(targets$FileName, source="genepix") #read in 2 color data (red/green)

MA <- MA.RG(RG,bc.method="none") #stores from red/green to avg of logs;
dim(RG$R)
dim(RG$G)
dim(MA$M) #same info, just transformed from red/green
dim(MA$A)
plot(MA$A[,1], MA$M[,1])    # MA plot for first sample

# microarray image
mypar(1,1)
imageplot(MA$M[,2], RG$printer, zlim=c(-3,3))
dev.off()

#The SummarizedExperiment class
#https://genomicsclass.github.io/book/pages/dataman2019.html

library(parathyroidSE)
data(parathyroidGenesSE)
se <- parathyroidGenesSE
se

# assay contains results of the assay
dim(se)
assay(se)[1:3,1:3] #RNAseq counts for each gene, across samples;
dim(assay(se))    # rows = features (ranges), columns = samples

# colData contains sample information
colData(se)[1:3,1:6] #~ to Pdata;
dim(colData(se))
names(colData(se)) #diff parameters
colData(se)$treatment #see what rx group each sample belongs to
as.numeric(colData(se)$treatment) #see what group each sample belongs to;

# rowRanges contains feature information
rowRanges(se)[1]
class(rowRanges(se))
length(rowRanges(se))    # number of genes
length(rowRanges(se)[[1]])    # number of exons for first gene
head(rownames(se))
metadata(rowRanges(se))

# additional metadata, including sample information
metadata(se)$MIAME
abstract(metadata(se)$MIAME)

#Importing NGS data in R
#Rsamtools provides low-level functions for reading and parsing raw NGS data stored in standard formats (more details below),
#vs: GenomicAlignments provides high-level functions and classes for reading and organizing NGS data as Bioconductor objects based on the GRanges class.
 
#http://genomicsclass.github.io/book/pages/import_NGS.html

#FASTQ files from the sequencing machine ->(either 1 file for a single-end sequencing sample, or 2 files for a paired-end sequencing sample).
#alignment software -> makes SAM -> compressed to BAM
#get sample BAM files
library(pasillaBamSubset)
library(Rsamtools)
filename <- untreated1_chr4()

bf <- BamFile(filename) #creates a BAM
seqinfo(bf) #gives chr location for BAM file sequences
sl <- seqlengths(bf)
sl
quickBamFlagSummary(bf)

#A number of functions in Rsamtools take an argument param, which expects a ScanBamParam specification.
#two important options are:
#  what - what kind of information to extract?
#  which - which ranges of alignments to extract?

#we can quickly pull out information about reads from a particular genomic range. Here we count the number of records (reads) on chromosome 4:
gr <- GRanges("chr4", IRanges(1, sl["chr4"]))
gr
countBam(bf, param = ScanBamParam(which = gr))

reads <- scanBam(BamFile(filename, yieldSize = 5))
reads

class(reads)
names(reads[[1]])
reads[[1]]$pos    # the aligned start position
reads[[1]]$rname    # the chromosome
reads[[1]]$strand    # the strand
reads[[1]]$qwidth    # the width of the read
reads[[1]]$seq    # the sequence of the read

#The GenomicAlignments package

library(GenomicAlignments)
ga <- readGAlignments(bf)
ga
length(ga)

granges(ga[1])

gr <- GRanges("chr4", IRanges(700000, 800000))
(fo <- findOverlaps(ga, gr))    # which reads over this range
countOverlaps(gr, ga)    # count overlaps of range with the reads
table(ga %over% gr)    # logical vector of read overlaps with the range
