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
  oxfordSamplesDBfilename     = "meta/Central_parasite_samples_21FEB2013__vivax_samples_21FEB2013.txt",
  outputFilename              = "meta/pv_1.0_sampleManifest.txt",
  outputRdaFilename           = sub("\\.txt", "\\.rda", outputFilename)
) {
  vrtrack <- read.delim(vrtrackFilename, as.is=TRUE, row.names=2)
  sangerLIMS <- read.delim(sangerLIMSfilename, as.is=TRUE, header=FALSE, row.names=1, col.names=c("ox_code", "sequencing_run_lane", "sequencing_date"), colClasses=c("character", "character", "Date"))
  samtrakLabInfo <- read.delim(samtrakLabInfoFilename, as.is=TRUE, row.names=1)
  solarisCountries <- read.delim(solarisCountriesFilename, as.is=TRUE)
  solarisCountriesPv <- subset(solarisCountries, ox_code %in% row.names(vrtrack) & !duplicated(ox_code))
  row.names(solarisCountriesPv) <- solarisCountriesPv[["ox_code"]]
  oxfordSamplesDB <- read.delim(oxfordSamplesDBfilename, as.is=TRUE, row.names="sample_code")
 
  setdiff(row.names(vrtrack), row.names(sangerLIMS))
  setdiff(row.names(vrtrack), row.names(samtrakLabInfo))
  setdiff(row.names(sangerLIMS), row.names(samtrakLabInfo))
  setdiff(row.names(sangerLIMS), row.names(vrtrack))
  setdiff(row.names(samtrakLabInfo), row.names(vrtrack))
  setdiff(row.names(samtrakLabInfo), row.names(sangerLIMS))
  setdiff(row.names(vrtrack), row.names(solarisCountriesPv))
  setdiff(row.names(solarisCountriesPv), row.names(vrtrack))
  
  table(samtrakLabInfo[["Project.Status"]], useNA="ifany")
  table(samtrakLabInfo[["Project.Status"]], row.names(samtrakLabInfo) %in% row.names(vrtrack), useNA="ifany")
  
  sampleManifest <- cbind(
    solarisCountriesPv[row.names(vrtrack), "ox_code", drop=FALSE],
    oxfordSamplesDB[row.names(vrtrack), ],
    solarisCountriesPv[row.names(vrtrack), -(which(names(solarisCountriesPv) == "ox_code"))],
    samtrakLabInfo[row.names(vrtrack), ],
    sangerLIMS[row.names(vrtrack), ],
    vrtrack
  )
  save(sampleManifest, file=outputRdaFilename)
  write.table(sampleManifest, file=outputFilename, quote=FALSE, sep="\t", row.names=FALSE)
  return(sampleManifest)
}
