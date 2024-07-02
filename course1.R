# https://genomicsclass.github.io/book/ chapters 11, 12, 13

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

