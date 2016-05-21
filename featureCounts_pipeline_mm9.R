# To run this pipeline, you should have R version of 3.0.2 or later 

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
# Targets.txt
# CellType	InputFile	OutputFile
options(digits=2)
targets <- readTargets()
targets
files <- targets$InputFile

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index if needed
# buildindex(basename="GENOME_NAME",reference="PATH_TO_FASTA")

# align reads
#align(index="/media/HD-LBU3/Clay/genome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/mm9_Rsubread_index",readfile1=targets$InputFile1,readfile2=targets$InputFile2,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,tieBreakHamming=TRUE,unique=TRUE)

# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=files,isPairedEnd=TRUE,annot.ext="/media/nas/clay/genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",isGTFAnnotationFile=TRUE,GTF.featureType="exon")
x <- DGEList(counts=fc$counts, genes=fc$annotation)

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x,design,plot=TRUE)

# cluster libraries
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)

