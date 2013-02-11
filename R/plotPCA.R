# plotPCA.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# pcaResults_01 <- plotPCA(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/pca/pv_01")
# pcaResults_01_69 <- plotPCA(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/pca/pv_01_69", calledSamplesThreshold=69)

plotPCA <- function(
  vcf                         = loadAndGenotypePvVcf(),
  missingnessThreshold        = 0.35,
#  calledSamplesThreshold      = NULL,
  calledSamplesThresholds     = seq(50, dim(vcf)[2], 2),
  pdfFilestem                 = "analysis/pca/pv_02"
) {
  require(ggplot2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }
  
  typableGT <- geno(vcf)[["GT"]]
  typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
  missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
  if(is.null(calledSamplesThresholds)) {
    calledSamplesThresholds <- length(which(missingnessPerSample < missingnessThreshold))
  }

#  GTdosagesInAllsamples <- (typableGT=="0/0") * 0.0 + (typableGT=="0/1") * 0.5 + (typableGT=="0/0") * 1.0
#  GTdosagesInAllsamples <- (typableGT=="0/0") * 0.0 + (typableGT=="0/1") * 0.5 + (typableGT=="1/1") * 1.0
  GTdosagesInAllsamples <- (typableGT=="0/0") * -1 + (typableGT=="0/1") * 0 + (typableGT=="1/1") * 1
  pcaResults <- list()
  sampleCountries <- countryCodes(dimnames(vcf)[[2]])
  sapply(
    calledSamplesThresholds,
    function(calledSamplesThreshold) {
      print(calledSamplesThreshold)
      pdfFilestemExtended <- paste(pdfFilestem, "_calledSamplesThreshold", calledSamplesThreshold, sep="")
      GTdosagesInNonMissingSamplesAndVariants <- GTdosagesInAllsamples[
        values(info(vcf))[["NS"]]>=calledSamplesThreshold,
        missingnessPerSample<missingnessThreshold
      ]
      pcaResults[[calledSamplesThreshold]] <- prcomp(t(GTdosagesInNonMissingSamplesAndVariants))
      pdf(paste(pdfFilestemExtended, "PC1vsPC2", "pdf", sep="."), height=6, width=10)
      print(
        qplot(
          pcaResults[[calledSamplesThreshold]][["x"]][, 1],
          pcaResults[[calledSamplesThreshold]][["x"]][, 2],
          colour=sampleCountries[missingnessPerSample<missingnessThreshold],
          xlab="PC1",
          ylab="PC2"
        )
        + theme_bw()
        + scale_colour_brewer(palette="Set1")
      )
      dev.off()
      pdf(paste(pdfFilestemExtended, "PC1vsPC3", "pdf", sep="."), height=6, width=10)
      print(
        qplot(
          pcaResults[[calledSamplesThreshold]][["x"]][, 1],
          pcaResults[[calledSamplesThreshold]][["x"]][, 3],
          colour=sampleCountries[missingnessPerSample<missingnessThreshold],
          xlab="PC1",
          ylab="PC3"
        )
        + theme_bw()
        + scale_colour_brewer(palette="Set1")
      )
      dev.off()
      pdf(paste(pdfFilestemExtended, "PC2vsPC3", "pdf", sep="."), height=6, width=10)
      print(
        qplot(
          pcaResults[[calledSamplesThreshold]][["x"]][, 2],
          pcaResults[[calledSamplesThreshold]][["x"]][, 3],
          colour=sampleCountries[missingnessPerSample<missingnessThreshold],
          xlab="PC2",
          ylab="PC3"
        )
        + theme_bw()
        + scale_colour_brewer(palette="Set1")
      )
      dev.off()
    }
  )
  return(pcaResults)
}
