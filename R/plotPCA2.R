# plotPCA2.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


#pca_02_missingness <- plotPCA2(missingnessThresholds=c(0.35, ))

# pcaResults_01 <- plotPCA(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/pca/pv_01")
# pcaResults_01_69 <- plotPCA(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/pca/pv_01_69", calledSamplesThreshold=69)

plotPCA2 <- function(
  vcf                         = loadAndGenotypePvVcf(),
  parameterSets               = list(
    "0.35_0_all" = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=0, samplesToRemove=NULL)
  ),
#  missingnessThresholds       = 0.35,
#  calledSamplesThreshold      = NULL,
#  calledSamplesThresholds     = seq(50, dim(vcf)[2], 2),
#  calledSamplesThresholds     = seq(76, 88, 2),
  pdfFilestem                 = "analysis/pca/pv_02",
  pfContaminants              = c("PD0173-C", "PH0309-C", "PJ0002-Cx", "PN0094-C", "PN0095-C", "PN0096-C"),
#  pvMOI                       = c("")
  pvMOI                       = NULL,
  pcsToPlot                   = list(c(1,2), c(1,3), c(2,3))
) {
  require(ggplot2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }
  
  typableGT <- geno(vcf)[["GT"]]
  typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
  missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
  GTdosagesInAllsamples <- (typableGT=="0/0") * -1 + (typableGT=="0/1") * 0 + (typableGT=="1/1") * 1
  sampleCountries <- countryCodes(dimnames(vcf)[[2]])
  
  pcaResultsList <- sapply(
    names(parameterSets),
    function(parameterSetName) {
      print(parameterSetName)
      pdfFilestemExtended <- paste(pdfFilestem, parameterSetName, sep=".")
      samplesToUse <- setdiff(
        dimnames(typableMissingGT)[[2]][missingnessPerSample<parameterSets[[parameterSetName]][["sampleMissingnessThreshold"]]],
        parameterSets[[parameterSetName]][["samplesToRemove"]]
      )
      missingnessPerVariant <- rowSums(typableMissingGT[, samplesToUse])

      GTdosagesInNonMissingSamplesAndVariants <- GTdosagesInAllsamples[
        missingnessPerVariant<=parameterSets[[parameterSetName]][["variantMissingnessThreshold"]],
        samplesToUse
      ]
      pcaResults <- prcomp(t(GTdosagesInNonMissingSamplesAndVariants))
      sampleShapes <- rep("uncontaminated", dim(GTdosagesInNonMissingSamplesAndVariants)[2])
      names(sampleShapes) <- dimnames(GTdosagesInNonMissingSamplesAndVariants)[[2]]
      sampleShapes[names(sampleShapes) %in% pfContaminants] <- "Pf contaminated"
      sampleShapes[names(sampleShapes) %in% pvMOI] <- "Pv MOI"
      lapply(
        pcsToPlot,
        function(pcs) {
          pdf(paste(pdfFilestemExtended, ".PC", pcs[1], "vsPC", pcs[2],".", "pdf", sep=""), height=6, width=10)
          print(
            qplot(
              pcaResults[["x"]][, pcs[1]],
              pcaResults[["x"]][, pcs[2]],
              colour=sampleCountries[samplesToUse],
              shape=sampleShapes,
              xlab=paste("PC", pcs[1], sep=""),
              ylab=paste("PC", pcs[2], sep="")
            )
            + theme_bw()
            + scale_colour_brewer(palette="Set1")
          )
          dev.off()
          
        }
      )
      return(pcaResults)
    },
    USE.NAMES=TRUE,
    simplify=FALSE
  )
  return(pcaResultsList)
}

