# variantSummaries.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


# vcf_01 <- sampleSummaries(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/sampleSummaries/pv_01")

variantSummaries <- function(
  vcf                         = loadAndGenotypePvVcf(),
  variableNames               = c("missingness", "heterozgosity", "depth", "meanMAF"),
  infoColumnNames             = names(values(info(vcf))),
  sampleMissingnessThresholds = c(0.35, 0.75, 0.95, 1),
  pdfFilestem                 = "analysis/variantSummaries/pv_02",
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
    sapply(
      sampleMissingnessThresholds,
      function(sampleMissingnessThreshold) {
        missingnessPerVariant <- rowSums(typableMissingGT[, missingnessPerSample<sampleMissingnessThreshold])
        missingnessProportionPerVariant <- missingnessPerVariant/length(which(missingnessPerSample<sampleMissingnessThreshold))
#        missingnessPerVariant <- rowSums(typableMissingGT)/dim(typableGT)[2]
        pdf(paste(pdfFilestem, "missingness", sampleMissingnessThreshold, "pdf", sep="."), height=height, width=width)
        print(
          qplot(
            missingnessPerVariant,
            binwidth=1,
            main=paste("Missingness per variant (samples <", sampleMissingnessThreshold, "missingness)"),
            xlab="Missingness",
            ylab="Frequency (number of variants)"
          )
          + theme_bw()
        )
        dev.off()
        pdf(paste(pdfFilestem, "proportionMissing", sampleMissingnessThreshold, "pdf", sep="."), height=height, width=width)
        print(
          qplot(
            missingnessProportionPerVariant,
            binwidth=0.01,
            main=paste("Proportion missing per variant (samples <", sampleMissingnessThreshold, "missingness)"),
            xlab="Proportion missing",
            ylab="Frequency (number of variants)"
          )
          + theme_bw()
        )
        dev.off()
      }
    )
  }

  if("heterozgosity" %in% variableNames) {
    typableHetGT <- matrix(typableGT=="0/1", ncol=ncol(typableGT))
    heterozygosityPerVariant <- rowSums(typableHetGT)/dim(typableGT)[2]
    pdf(paste(pdfFilestem, "heterozgosity", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        heterozygosityPerVariant,
        main="Heterozygosity per variant",
        xlab="Heterozygosity",
        ylab="Frequency (number of variants)"
      )
      + theme_bw()
    )
    dev.off()
  }
    
  if("depth" %in% variableNames) {
    meanDPperVariant <- rowMeans(geno(vcf)[["DP"]])
    pdf(paste(pdfFilestem, "depth", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        meanDPperVariant,
        main="Mean depth per variant",
        xlab="Mean depth",
        ylab="Frequency (number of variants)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  if("meanMAF" %in% variableNames) {
    meanMAFperVariant <- rowMeans(geno(vcf)[["MAF"]], na.rm=TRUE)
    pdf(paste(pdfFilestem, "meanMAF", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        meanMAFperVariant,
        main="Mean MAF per variant",
        xlab="Mean MAF",
        ylab="Frequency (number of variants)"
      )
      + theme_bw()
    )
    dev.off()
  }
  
  sapply(
    infoColumnNames,
    function(infoColumnName) {
      pdf(paste(pdfFilestem, infoColumnName, "pdf", sep="."), height=height, width=width)
      print(
        qplot(
          values(info(vcf))[[infoColumnName]],
          xlab=infoColumnName,
          ylab="Frequency (number of variants)"
        )
        + theme_bw()
      )
      dev.off()
    }
  )
  
  return(vcf)
}


