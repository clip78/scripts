# the input expression matrix file (directory and file name of 
# the expression matrix have to be adapted)

# the output directory and the output file names have to be adapted

# the design matrix "design", the  contrast matrix "constrast" as well as the "colnames"
# in dependence on the number of samples and replicates and names

# Before starting the script, create the output directory where the generated lists will be written

# cleanup workspace
rm (list=ls())

# close all graphics devices
graphics.off()

warnings();

#----------------------------------------------------------------------

library(limma)

input_expression_file = "/Users/millerlab/Downloads/L1-NSM"

wholeMatrix = read.delim(file="L1-NSM_L1-cells-ref_RMA_noBG.txt", header= TRUE, row.names=c(1))

# has to be addapted == number of genes being regarded
# Because we want to get here all genes, this should be equal the number
# of genes being available in the expression matrix
# (we will filter later on with a perl script for up and down regulated genes)
numberOfGenes=18458

# has to be adapted
design = model.matrix(~ 0+factor(c(1,1,2,2)))

# has to be adapted
colnames(design) = c("ref", "NSM")
fit = lmFit(wholeMatrix, design)

# has to be adapted
contrast.matrix = makeContrasts(NSM-ref, levels = design)

fit2 = contrasts.fit(fit, contrast.matrix)

fit3 = eBayes(fit2)


# all comparisions have to be calculated
NSM_ref = topTable(fit3, coef=1, number = numberOfGenes, adjust="BH")

#====

# output files of all comparisions, ouptut directory and file names have to be adapted
write.table(NSM_ref, sep = "\t", quote = FALSE, file = "/Users/millerlab/Downloads/L1-NSM/L1-NSM_vs_L1-cells_ref.WS200.RMA.noBG.txt")
