# sampleContaminationPlot.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


sampleContaminationPlot <- function(
  vcf                         = loadAndGenotypePvVcf(),
  pdfFilestem                 = "analysis/sampleContamination/pv_02_MAFbySample",
  height                      = 12,
  width                       = 12
) {
  require(ggplot2)
  require(reshape2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }
  MAFDF <- melt(geno(vcf)[["MAF"]], varnames=c("variant", "sample"), value.name="MAF")
#  browser()
  pdf(paste(pdfFilestem, "allVariants", "pdf", sep="."), height=height, width=width)
#  qplot(MAF, facets=sample~., fill=sample, binwidth=0.01, data=MAFDF, log="y") + theme_bw()
  print(qplot(MAF, fill=sample, binwidth=0.01, data=MAFDF) + facet_wrap(~ sample) + theme_bw())
  dev.off()
  MAFDF0.05 <- subset(MAFDF, MAF>0.05)
  pdf(paste(pdfFilestem, "MAFgt0.05genotypes", "pdf", sep="."), height=height, width=width)
#  pdf("~/MAFbySample.pdf", height=12, width=12)
  print(qplot(MAF, fill=sample, binwidth=0.01, data=MAFDF0.05) + facet_wrap(~ sample) + theme_bw())
  dev.off()
}
