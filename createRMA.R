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
library("cele25mrws200cdf")

setwd("/Users/millerlab/Downloads/L1-NSM")

Data = ReadAffy(celfile.path = "/Users/millerlab/Downloads/L1-NSM", cdfname = "cele25mrws200cdf")
eset = rma(Data, background = FALSE)
write.exprs(eset, file = "L1-NSM_L1-cells-ref_RMA_noBG.txt")