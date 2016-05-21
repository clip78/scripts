# It is expected that all cel files
# are located in the directory
# /tmp/celegans/celfiles where you have
# read permissions.
#
# The directory where the RMA log2 expression matrices
# is in the script set to /tmp/celegans/expression_matrices
# where you need to have write and read permissions
#
# You may change the directories by adapting the script
# accordingly

# cleanup workspace
rm (list= ls())

# close graphics device
graphics.off()

warnings()

#----------------------------------------------------------------------

library(affy)
library("x004ws200ws199cdf")
library("x005ws200ws199cdf")
library(farms)

setwd("/tmp/celegans/expression_matrices")

Data = ReadAffy(celfile.path = "/tmp/celegans/celfiles", cdfname = "x004ws200ws199cdf")
eset = exp.farms(Data, bgcorrect.method="none",normalize.method="quantiles", pmcorrect.method="pmonly")
write.exprs(eset, file = "expressionQuantileWithoutBackgroundNormalization_X004.txt")

Data = ReadAffy(celfile.path = "/tmp/celegans/celfiles", cdfname = "x005ws200ws199cdf")
eset = exp.farms(Data, bgcorrect.method="none",normalize.method="quantiles", pmcorrect.method="pmonly")
write.exprs(eset, file = "expressionQuantileWithoutBackgroundNormalization_X005.txt")
