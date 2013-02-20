# createSampleManifest.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


createSampleManifest <- function(
  vrtrackFilename             = "meta/pv_1.0_samples.txt",
  sangerLIMSfilename          = "meta/pf_1_0_rundates.tab",
  samtrakLabInfoFilename      = "meta/Vivax.txt",
  solarisCountriesFilename    = "meta/pv_1.0_countries.txt",
  outputFilename              = "meta/pv_1.0_sampleManifest.txt"
) {
  sampleManifest <- read.delim(vrtrackFilename, as.is=TRUE)
  sangerLIMSfilename <- read.delim(sangerLIMSfilename, as.is=TRUE)
  samtrakLabInfo <- read.delim(samtrakLabInfoFilename, as.is=TRUE)
  solarisCountries <- read.delim(solarisCountriesFilename, as.is=TRUE)
  
  browser()
}
