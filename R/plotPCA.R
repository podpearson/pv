# plotPCA.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


plotPCA <- function(
  vcf                         = loadAndGenotypePvVcf(),
  missingnessThreshold        = 0.4,
  calledSamplesThreshold      = 69,
  pdfFilename                 = "analysis/pca/pca.pdf"
) {
  require(ggplot2)
  
  typableGT <- geno(vcf)[["GT"]]
  typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
  missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]

  GTdosagesInAllsamples <- (typableGT=="0/0") * 0.0 + (typableGT=="0/1") * 0.5 + (typableGT=="0/0") * 1.0
  GTdosagesInNonMissingSamplesAndVariants <- GTdosagesInAllsamples[
    values(info(vcf))[["NS"]]>=calledSamplesThreshold,
    missingnessPerSample<missingnessThreshold
  ]
  pcaResults <- prcomp(t(GTdosagesInNonMissingSamplesAndVariants))
  sampleCountries <- countryCodes(dimnames(vcf)[[2]])

  pdf(pdfFilename, height=6, width=10)
  print(
    qplot(
      pcaResults[["x"]][, 1],
      pcaResults[["x"]][, 3],
      colour=sampleCountries[missingnessPerSample<0.4],
      xlab="PC1",
      ylab="PC3"
    )
    + theme_bw()
    + scale_colour_brewer(palette="Set1")
  )
  dev.off()
}
