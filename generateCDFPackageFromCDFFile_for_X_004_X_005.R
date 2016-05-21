# It is expected that the cdf file(s) is (are) located in the
# directory /tmp/celegans/cdf_files.
# The according R cdf package(s) will be written into the
# already existing directory /tmp/celegans/cdf_packages
#
# Make sure, you have read permissions in the directory /tmp/celegans/cdf_files
# and write and read permission in the directory /tmp/celegans/cdf_packages
#
# If you use other directories, you have to adapt the script
# accordingly

rm (list=ls())

# close graphics device
graphics.off()

warnings()

#----------------------------------------------------------------------

# load necessary Bioconductor packages
library(makecdfenv)

pkgpath = tempdir()

make.cdf.package("x004ws200ws199.cdf", cdf.path = "/tmp/celegans/cdf_files", package.path = "/tmp/celegans/cdf_packages", packagename = "x004ws200ws199cdf", species="C_elegans")

make.cdf.package("x005ws200ws199.cdf", cdf.path = "/tmp/celegans/cdf_files", package.path = "/tmp/celegans/cdf_packages", packagename = "x005ws200ws199cdf", species="C_elegans")
