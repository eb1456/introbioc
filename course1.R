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

