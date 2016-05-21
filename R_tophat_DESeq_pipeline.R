library("ShortRead")
library("DESeq")

#fqQC = qa(dirPath = "../", pattern = ".txt.gz$", type = "fastq")
#report(fqQC, type = "html", dest = "fastqQAreport")
samples <- read.csv("samples.txt", sep="\t")
gf <- "/media/HD-LBU3/Clay/genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
bowind = "/media/HD-LBU3/Clay/genome/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome.bt2"
cmd = with(samples, paste("tophat2 -p 8 -g 1 -z pigz --rg-id 1 --rg-sample", LibraryName, "--rg-library 1 --rg-platform illumina -o", LibraryName, bowind, fastq1, fastq2, "&&"))
cat(cmd)
cat("\n", "Press [enter] to continue")
line <- readline()

for(i in seq_len(nrow(samples))) {
lib = samples$LibraryName[i]
ob = file.path(lib, "accepted_hits.bam")
# sort by name, convert to SAM for htseq-count
cat(paste0("samtools sort -n ",ob," ",lib,"_sn"),"\n")
cat(paste0("samtools view -o ",lib,"_sn.sam ",lib,"_sn.bam"),"\n")
# sort by position and index for IGV
cat(paste0("samtools sort ",ob," ",lib,"_s"),"\n")
cat(paste0("samtools index ",lib,"_s.bam"),"\n\n")
}
samples$countf = paste(samples$LibraryName, "count", sep = ".")
cmd = paste0("sudo htseq-count -s no -a 10 ", samples$LibraryName, "_sn.sam ", gf," > ", samples$countf, " &&")
cat(cmd)
cat ("\n", "Press [enter] to continue")
line <- readline()

samplesDESeq = with(samples, data.frame(shortname = I(shortname), countf = I(countf), condition = condition, LibraryLayout = Library_Layout))
cds = newCountDataSetFromHTSeqCount(samplesDESeq)
samplesDESeq = with(samples, data.frame(shortname = I(LibraryName), countf = I(countf), condition = condition, LibraryLayout = Library_Layout))
samplesDESeq
cds = newCountDataSetFromHTSeqCount(samplesDESeq)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cdsB = estimateDispersions(cds, method = "blind")
vsd = varianceStabilizingTransformation(cdsB)
p = plotPCA(vsd, intgroup = c("condition","LibraryLayout"))
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds,"WT","KO")
plotMA(res)
resSig = res[which(res$padj < 0.1),]
head( resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ] )
head( resSig[ order(resSig$log2FoldChange, decreasing = FALSE), ] )
table( res$padj < 0.1 )
write.csv(res, file = "res_DESeq.csv")
hist(res$pval, breaks = 100)
