# the input expression matrix file (directory and file name of 
# the expression matrix have to be adapted)

# the method (pfp==FDR or pval==p-value) and the according threshold have to be adapted.
# Also the variable cl must be adapted (do we have samples in triplcates or duplicates etc.?)

# the number of permutations can be adapted

# the output directory and output file names must be adapted

# cleanup workspace
rm (list=ls())

# close all graphics devices
graphics.off()

warnings();

#----------------------------------------------------------------------

library(RankProd)

file_input_expression = "<inputDirectory>/<expressionMatrixFile>"

wholeMatrix = read.delim(file_input_expression, header = TRUE)

# can be adapted, the higher permutation numbers are better
numberPermutations = 100

# if pfp (percentage of false prediction == FDR), method has to
# be equal "pfp".
# Because the RankProduct is a permuation test, a p-value can be
# applied too. In the later case the method has to be set to "pval"
method = "pfp" # == FDR, if 

# has to be adapted according the chosen FDR or p-value
cutoff = 0.1

# has to be adapted, at the moment 3 replicates per sample
cl = rep(c(0,1), c(3,3))

geneNames = as.character(wholeMatrix[,1])

numberRows <- nrow(wholeMatrix)

numberColsAll <- ncol(wholeMatrix)

expressionDataAll <- wholeMatrix[,2:numberColsAll]


e1a = expressionDataAll[[1]]
e1b = expressionDataAll[[2]]
e1c = expressionDataAll[[3]]

e2a = expressionDataAll[[4]]
e2b = expressionDataAll[[5]]
e2c = expressionDataAll[[6]]

e3a = expressionDataAll[[7]]
e3b = expressionDataAll[[8]]
e3c = expressionDataAll[[9]]

e4a = expressionDataAll[[10]]
e4b = expressionDataAll[[11]]
e4c = expressionDataAll[[12]]

e5a = expressionDataAll[[13]]
e5b = expressionDataAll[[14]]
e5c = expressionDataAll[[15]]

e6a = expressionDataAll[[16]]
e6b = expressionDataAll[[17]]
e6c = expressionDataAll[[18]]

e7a = expressionDataAll[[19]]
e7b = expressionDataAll[[20]]
e7c = expressionDataAll[[21]]

e8a = expressionDataAll[[22]]
e8b = expressionDataAll[[23]]
e8c = expressionDataAll[[24]]

e9a = expressionDataAll[[25]]
e9b = expressionDataAll[[26]]
e9c = expressionDataAll[[27]]

e10a = expressionDataAll[[28]]
e10b = expressionDataAll[[29]]
e10c = expressionDataAll[[30]]

e11a = expressionDataAll[[31]]
e11b = expressionDataAll[[32]]
e11c = expressionDataAll[[33]]

e12a = expressionDataAll[[34]]
e12b = expressionDataAll[[35]]
e12c = expressionDataAll[[36]]

e13a = expressionDataAll[[37]]
e13b = expressionDataAll[[38]]
e13c = expressionDataAll[[39]]

e14a = expressionDataAll[[40]]
e14b = expressionDataAll[[41]]
e14c = expressionDataAll[[42]]

e15a = expressionDataAll[[43]]
e15b = expressionDataAll[[44]]
e15c = expressionDataAll[[45]]

e16a = expressionDataAll[[46]]
e16b = expressionDataAll[[47]]
e16c = expressionDataAll[[48]]

e17a = expressionDataAll[[49]]
e17b = expressionDataAll[[50]]
e17c = expressionDataAll[[51]]

e18a = expressionDataAll[[52]]
e18b = expressionDataAll[[53]]
e18c = expressionDataAll[[54]]

e19a = expressionDataAll[[55]]
e19b = expressionDataAll[[56]]
e19c = expressionDataAll[[57]]

e20a = expressionDataAll[[58]]
e20b = expressionDataAll[[59]]
e20c = expressionDataAll[[60]]

e21a = expressionDataAll[[61]]
e21b = expressionDataAll[[62]]
e21c = expressionDataAll[[63]]

e22a = expressionDataAll[[64]]
e22b = expressionDataAll[[65]]
e22c = expressionDataAll[[66]]

e23a = expressionDataAll[[67]]
e23b = expressionDataAll[[68]]
e23c = expressionDataAll[[69]]

e24a = expressionDataAll[[70]]
e24b = expressionDataAll[[71]]
e24c = expressionDataAll[[72]]

e25a = expressionDataAll[[73]]
e25b = expressionDataAll[[74]]
e25c = expressionDataAll[[75]]

e26a = expressionDataAll[[76]]
e26b = expressionDataAll[[77]]
e26c = expressionDataAll[[78]]

e27a = expressionDataAll[[79]]
e27b = expressionDataAll[[80]]
e27c = expressionDataAll[[81]]

e28a = expressionDataAll[[82]]
e28b = expressionDataAll[[83]]
e28c = expressionDataAll[[84]]

e29a = expressionDataAll[[85]]
e29b = expressionDataAll[[86]]
e29c = expressionDataAll[[87]]

#####

X4_vs_X5 = cbind(e4a,e4b)
X4_vs_X5 = cbind(X4_vs_X5, e4c)
X4_vs_X5 = cbind(X4_vs_X5, e5a)
X4_vs_X5 = cbind(X4_vs_X5, e5b)
X4_vs_X5 = cbind(X4_vs_X5, e5c)

X4_vs_X6 = cbind(e4a,e4b)
X4_vs_X6 = cbind(X4_vs_X6, e4c)
X4_vs_X6 = cbind(X4_vs_X6, e6a)
X4_vs_X6 = cbind(X4_vs_X6, e6b)
X4_vs_X6 = cbind(X4_vs_X6, e6c)

X5_vs_X6 = cbind(e5a,e5b)
X5_vs_X6 = cbind(X5_vs_X6, e5c)
X5_vs_X6 = cbind(X5_vs_X6, e6a)
X5_vs_X6 = cbind(X5_vs_X6, e6b)
X5_vs_X6 = cbind(X5_vs_X6, e6c)

X4_vs_X7 = cbind(e4a,e4b)
X4_vs_X7 = cbind(X4_vs_X7, e4c)
X4_vs_X7 = cbind(X4_vs_X7, e7a)
X4_vs_X7 = cbind(X4_vs_X7, e7b)
X4_vs_X7 = cbind(X4_vs_X7, e7c)

X5_vs_X7 = cbind(e5a,e5b)
X5_vs_X7 = cbind(X5_vs_X7, e5c)
X5_vs_X7 = cbind(X5_vs_X7, e7a)
X5_vs_X7 = cbind(X5_vs_X7, e7b)
X5_vs_X7 = cbind(X5_vs_X7, e7c)

X6_vs_X7 = cbind(e6a,e6b)
X6_vs_X7 = cbind(X6_vs_X7, e6c)
X6_vs_X7 = cbind(X6_vs_X7, e7a)
X6_vs_X7 = cbind(X6_vs_X7, e7b)
X6_vs_X7 = cbind(X6_vs_X7, e7c)

numberOfGenes = length(geneNames)

rankProduct_X4_vs_X5 = RP(X4_vs_X5, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)
rankProduct_X4_vs_X6 = RP(X4_vs_X6, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)
rankProduct_X5_vs_X6 = RP(X5_vs_X6, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)
rankProduct_X4_vs_X7 = RP(X4_vs_X7, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)
rankProduct_X5_vs_X7 = RP(X5_vs_X7, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)
rankProduct_X6_vs_X7 = RP(X6_vs_X7, cl, num.perm=numberPermutations, logged=TRUE, na.rm=FALSE, gene.names=geneNames, plot=TRUE, rand=123)

tableOfSortedGenes_X4_vs_X5 = topGene(rankProduct_X4_vs_X5, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)
tableOfSortedGenes_X4_vs_X6 = topGene(rankProduct_X4_vs_X6, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)
tableOfSortedGenes_X5_vs_X6 = topGene(rankProduct_X5_vs_X6, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)
tableOfSortedGenes_X4_vs_X7 = topGene(rankProduct_X4_vs_X7, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)
tableOfSortedGenes_X5_vs_X7 = topGene(rankProduct_X5_vs_X7, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)
tableOfSortedGenes_X6_vs_X7 = topGene(rankProduct_X6_vs_X7, cutoff = cutoff, method = method, logged=TRUE, logbase=2, gene.names=geneNames)

write.table(tableOfSortedGenes_X4_vs_X5$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X5_with_X4_smaller_X5.txt")
write.table(tableOfSortedGenes_X4_vs_X5$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X5_with_X4_greater_X5.txt")

write.table(tableOfSortedGenes_X4_vs_X6$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X6_with_X4_smaller_X6.txt")
write.table(tableOfSortedGenes_X4_vs_X6$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X6_with_X4_greater_X6.txt")

write.table(tableOfSortedGenes_X5_vs_X6$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X5_vs_X6_with_X5_smaller_X6.txt")
write.table(tableOfSortedGenes_X5_vs_X6$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X5_vs_X6_with_X5_greater_X6.txt")

write.table(tableOfSortedGenes_X4_vs_X7$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X7_with_X4_smaller_X7.txt")
write.table(tableOfSortedGenes_X4_vs_X7$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X4_vs_X7_with_X4_greater_X7.txt")

write.table(tableOfSortedGenes_X5_vs_X7$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X5_vs_X7_with_X5_smaller_X7.txt")
write.table(tableOfSortedGenes_X5_vs_X7$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X5_vs_X7_with_X5_greater_X7.txt")

write.table(tableOfSortedGenes_X6_vs_X7$Table1, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X6_vs_X7_with_X6_smaller_X7.txt")
write.table(tableOfSortedGenes_X6_vs_X7$Table2, sep = "\t", quote=FALSE, file = "/home/shenz/celegans/expressionMatrices/auswertung/diff_expr/rankProduct_fdr_0.1_rma/rank_product_tables/RP_rmaWithoutBackgroundNormalization_WS199_X6_vs_X7_with_X6_greater_X7.txt")
