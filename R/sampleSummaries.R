# sampleSummaries.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


sampleSummaries <- function(
  vcf                         = loadAndGenotypePvVcf(),
  variableNames               = c("missingness", "heterozgosity", "singletons", "depth"),
  pdfFilestem                 = "analysis/sampleSummaries/analysis/pca/pv_02",
  height                      = 6,
  width                       = 10
) {
  require(ggplot2)
  typableGT <- geno(vcf)[["GT"]]

  if("missingness" %in% variableNames) {
    typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
    missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
    pdf(paste(pdfFilestem, "missingness", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        missingnessPerSample,
        binwidth=0.05,
        main="Missingness per sample",
        xlab="Missingness",
        ylab="Frequency (number of samples)"
      )
      + theme_bw()
    )
    dev.off()
  }

  if("missingness" %in% variableNames) {
    typableHetGT <- matrix(typableGT=="0/1", ncol=ncol(typableGT))
    heterozygosityPerSample <- colSums(typableHetGT)/dim(typableGT)[1]
    pdf(paste(pdfFilestem, "missingness", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        heterozygosityPerSample,
        binwidth=0.05,
        main="Heterozygosity per sample",
        xlab="Heterozygosity",
        ylab="Frequency (number of samples)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  return(vcf)
}

