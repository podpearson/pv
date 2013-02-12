# sampleSummaries.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# vcf_01 <- sampleSummaries(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/sampleSummaries/pv_01")

sampleSummaries <- function(
  vcf                         = loadAndGenotypePvVcf(),
  variableNames               = c("missingness", "heterozgosity", "singletons", "depth"),
  pdfFilestem                 = "analysis/sampleSummaries/pv_02",
  height                      = 6,
  width                       = 10
) {
  require(ggplot2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }

  typableGT <- geno(vcf)[["GT"]]

  if("missingness" %in% variableNames) {
    typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
    missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
    missingnessPerSampleDF <- data.frame(sampleID=dimnames(typableGT)[[2]], missingnessPerSample=missingnessPerSample)
    write.table(missingnessPerSampleDF, file=paste(pdfFilestem, "missingnessPerSample", "txt", sep="."), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
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

  if("heterozgosity" %in% variableNames) {
    typableHetGT <- matrix(typableGT=="0/1", ncol=ncol(typableGT))
    heterozygosityPerSample <- colSums(typableHetGT)/dim(typableGT)[1]
    pdf(paste(pdfFilestem, "heterozgosity", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        heterozygosityPerSample,
        main="Heterozygosity per sample",
        xlab="Heterozygosity",
        ylab="Frequency (number of samples)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  if("singletons" %in% variableNames) {
    typableHomAltGT <- matrix(typableGT=="1/1", ncol=ncol(typableGT))
    singletonGT <- typableHomAltGT[values(info(vcf))[["singleton"]]==TRUE, ]
    singletonsPerSample <- colSums(singletonGT)
    pdf(paste(pdfFilestem, "singletons", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        singletonsPerSample,
        main="Number of singletons per sample",
        xlab="Number of singletons",
        ylab="Frequency (number of samples)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  if("depth" %in% variableNames) {
    meanDPperSample <- colMeans(geno(vcf)[["DP"]])
    pdf(paste(pdfFilestem, "depth", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        meanDPperSample,
        main="Mean depth per sample",
        xlab="Mean depth",
        ylab="Frequency (number of samples)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  return(vcf)
}

