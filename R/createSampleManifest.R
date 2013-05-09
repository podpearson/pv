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
  irodsFilename               = "meta/pv_metadata.tab",
  pvSequenomFilename          = "data/sequenom/vivax_data_10APR2013.txt",
  jacobSamplesFilename        = "data/jacob/PVivax-datasets/metadata/samplesMeta.tab",
  columnsToLoadFromSequenom   = c("study_code", "sample_id", "source_code"),
  sampleManifestFilename      = "meta/pv_1.0_sampleManifest.txt",
  sampleManifestRdaFilename   = sub("\\.txt", "\\.rda", sampleManifestFilename),
  sampleManifestPostAnalysisFilename      = "meta/pv_1.0_sampleManifestPostAnalysis.txt",
  sampleManifestPostAnalysisRdaFilename   = sub("\\.txt", "\\.rda", sampleManifestPostAnalysisFilename),
  unseqSampleManifestFilename = "meta/pv_1.0_unsequencedSampleManifest.txt",
  unseqSampleManifestRdaFilename = sub("\\.txt", "\\.rda", unseqSampleManifestFilename)
) {
  vrtrack <- read.delim(vrtrackFilename, as.is=TRUE, row.names=2)
  sangerLIMS <- read.delim(sangerLIMSfilename, as.is=TRUE, header=FALSE, row.names=1, col.names=c("ox_code", "sequencing_run_lane", "sequencing_date"), colClasses=c("character", "character", "Date"))
  samtrakLabInfo <- read.delim(samtrakLabInfoFilename, as.is=TRUE, row.names=1)
  solarisCountries <- read.delim(solarisCountriesFilename, as.is=TRUE)
  solarisCountriesPv <- subset(solarisCountries, ox_code %in% row.names(vrtrack) & !duplicated(ox_code))
  row.names(solarisCountriesPv) <- solarisCountriesPv[["ox_code"]]
  irods <- read.delim(irodsFilename, as.is=TRUE, header=FALSE, col.names=c("bamFilename", "numberOfReads", "ox_code"), colClasses=c("character", "integer", "character"))
  irods[irods[["ox_code"]] == "PJ0004--C", "ox_code"] <- "PJ0004-C"
  irods[irods[["ox_code"]] == "PQ002-C_200", "ox_code"] <- "PQ0002-C"
  irodsHuman <- irods[grepl("_human", irods[["bamFilename"]]), ]
  irodsHuman <- irodsHuman[!duplicated(irodsHuman[["ox_code"]]), ]
  row.names(irodsHuman) <- irodsHuman[["ox_code"]]
  irodsNonHumanPhix <- irods[!grepl("_human", irods[["bamFilename"]]), ]
  irodsNonHumanPhix <- irodsNonHumanPhix[!duplicated(irodsNonHumanPhix[["ox_code"]]), ]
  row.names(irodsNonHumanPhix) <- irodsNonHumanPhix[["ox_code"]]
  oxfordSamplesDB <- read.delim(oxfordSamplesDBfilename, as.is=TRUE, row.names="sample_code")
  sangerLIMS[["sequencing_run"]] <- sub("^([0-9]+)_.*$", "\\1", sangerLIMS[["sequencing_run_lane"]])
  sangerLIMS[["sequencing_lane"]] <- sub("^[0-9]+_([0-9]+)_.*$", "\\1", sangerLIMS[["sequencing_run_lane"]])
  pvSequenom <- read.table(pvSequenomFilename, header=TRUE, sep="\t", as.is=TRUE)
  pvSequenomPvOnly <- subset(pvSequenom, !(study_code %in% c("WT", "WA")))
  pvSequenomPvOnlyOneVariant <- subset(pvSequenomPvOnly, sequence_code==pvSequenomPvOnly[1, "sequence_code"])
  row.names(pvSequenomPvOnlyOneVariant) <- pvSequenomPvOnlyOneVariant[["sample_code"]]
  pvSequenomSampleColumns <- pvSequenomPvOnlyOneVariant[, columnsToLoadFromSequenom]
  names(pvSequenomSampleColumns) <- paste("SQ", columnsToLoadFromSequenom, sep="_")
 
  setdiff(row.names(vrtrack), row.names(sangerLIMS))
  setdiff(row.names(vrtrack), row.names(samtrakLabInfo))
  setdiff(row.names(sangerLIMS), row.names(samtrakLabInfo))
  setdiff(row.names(sangerLIMS), row.names(vrtrack))
  setdiff(row.names(samtrakLabInfo), row.names(vrtrack))
  setdiff(row.names(samtrakLabInfo), row.names(sangerLIMS))
  setdiff(row.names(vrtrack), row.names(solarisCountriesPv))
  setdiff(row.names(solarisCountriesPv), row.names(vrtrack))
  setdiff(row.names(vrtrack), row.names(irodsHuman))
  setdiff(row.names(vrtrack), row.names(irodsNonHumanPhix))
  setdiff(row.names(irodsHuman), row.names(vrtrack))
  setdiff(row.names(irodsNonHumanPhix), row.names(vrtrack))
  setdiff(row.names(vrtrack), row.names(pvSequenomPvOnlyOneVariant))
  setdiff(row.names(pvSequenomPvOnlyOneVariant), row.names(vrtrack))
  table(row.names(irodsHuman) %in% row.names(irodsNonHumanPhix))
  table(row.names(irodsNonHumanPhix) %in% row.names(irodsHuman))
  setdiff(row.names(samtrakLabInfo), row.names(oxfordSamplesDB))
  setdiff(row.names(oxfordSamplesDB), row.names(samtrakLabInfo))
  
  unsequencedSampleManifest <- samtrakLabInfo[setdiff(row.names(samtrakLabInfo), row.names(vrtrack)), ]
  oxfordSamplesDB[row.names(subset(unsequencedSampleManifest, Human.... < 50 & Quantity..ng. > 200)), ]
  require(ggplot2)
  pdf("analysis/20130304_slides/unsequencedSamplesCountry.pdf", height=6, width=10)
  print(qplot(Human...., Quantity..ng., data=unsequencedSampleManifest, log="y", xlab="% human", ylab="Quantity", colour=Biosample.Location) + theme_bw() + geom_hline(yintercept = 200, colour="blue") + geom_vline(xintercept = 80, colour="green"))
  dev.off()
  with(unsequencedSampleManifest, table(Human....>80, Quantity..ng.>200))
  unsequencedSampleManifest[["hasSequenomData"]] <- row.names(unsequencedSampleManifest) %in% row.names(pvSequenomPvOnlyOneVariant)
  save(unsequencedSampleManifest, file=unseqSampleManifestRdaFilename)
  write.table(unsequencedSampleManifest, file=unseqSampleManifestFilename, quote=FALSE, sep="\t", row.names=FALSE)

  setdiff(row.names(unsequencedSampleManifest), row.names(pvSequenomPvOnlyOneVariant))
  setdiff(row.names(pvSequenomPvOnlyOneVariant), row.names(unsequencedSampleManifest))
  intersect(row.names(unsequencedSampleManifest), row.names(pvSequenomPvOnlyOneVariant))
  
  table(samtrakLabInfo[["Project.Status"]], useNA="ifany")
  table(samtrakLabInfo[["Project.Status"]], row.names(samtrakLabInfo) %in% row.names(vrtrack), useNA="ifany")
  
  sampleManifest <- cbind(
    solarisCountriesPv[row.names(vrtrack), "ox_code", drop=FALSE],
    oxfordSamplesDB[row.names(vrtrack), ],
    solarisCountriesPv[row.names(vrtrack), -(which(names(solarisCountriesPv) == "ox_code"))],
    samtrakLabInfo[row.names(vrtrack), ],
    sangerLIMS[row.names(vrtrack), ],
    pvSequenomSampleColumns[row.names(vrtrack), ],
    vrtrack
  )
  sampleManifest[["numberOfHumanReads"]] <- irodsHuman[row.names(sampleManifest), "numberOfReads"]
  sampleManifest[["numberOfNonHumanReads"]] <- irodsNonHumanPhix[row.names(sampleManifest), "numberOfReads"]
  sampleManifest[["proportionHumanReads"]] <- sampleManifest[["numberOfHumanReads"]]/(sampleManifest[["numberOfHumanReads"]]+sampleManifest[["numberOfNonHumanReads"]])
  sampleManifest[["richard_donor_source_code"]] <- sub("_[A-Za-z].*$", "", sampleManifest[["source_code"]])
  sampleManifest[["numberOfAnnotatedDuplicates"]] <- table(sampleManifest[["richard_donor_source_code"]])[sampleManifest[["richard_donor_source_code"]]]
  sampleManifest[["annotatedDuplicateIDs"]] <- sapply(
    rownames(sampleManifest),
    function(sampleID) {
      paste(
        setdiff(
          rownames(subset(sampleManifest, richard_donor_source_code==sampleManifest[sampleID, "richard_donor_source_code"])),
          sampleID
        ),
        collapse=","
      )
    }
  )
  sampleManifest[["sampleIDforPlots"]] <- paste(
    rownames(sampleManifest),
    " (",
    sampleManifest[["richard_donor_source_code"]],
    ", ",
    sampleManifest[["country"]],
    ")",
    sep=""
  )
  
  subset(sampleManifest, numberOfAnnotatedDuplicates>1, c("numberOfAnnotatedDuplicates", "annotatedDuplicateIDs"))
  
  with(sampleManifest, table(study_code, SQ_study_code, useNA="ifany"))
  with(sampleManifest, table(study_code==SQ_study_code, useNA="ifany"))
  with(sampleManifest, table(sample_id == SQ_sample_id, useNA="ifany"))
  with(sampleManifest, table(source_code==SQ_source_code, useNA="ifany"))
  
  save(sampleManifest, file=sampleManifestRdaFilename)
  write.table(sampleManifest, file=sampleManifestFilename, quote=FALSE, sep="\t", row.names=FALSE)
  jacobSamples = read.delim(jacobSamplesFilename, as.is=TRUE, row.names=1)
  sampleManifest[["SnpStringent"]] <- jacobSamples[row.names(sampleManifest), "SnpStringent"]
  sampleManifest[["SampleStringent"]] <- jacobSamples[row.names(sampleManifest), "SampleStringent"]
  save(sampleManifest, file=sampleManifestPostAnalysisRdaFilename)
  write.table(sampleManifest, file=sampleManifestPostAnalysisFilename, quote=FALSE, sep="\t", row.names=FALSE)
  return(sampleManifest)
}
