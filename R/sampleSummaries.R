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
  variableNames               = c("missingness", "heterozgosity", "singletons", "depth", "missingnessVsDepth"),
  pdfFilestem                 = "analysis/sampleSummaries/pv_02",
  pfContaminants              = c("PD0173-C", "PH0309-C", "PJ0002-Cx", "PN0094-C", "PN0095-C", "PN0096-C"),
  height                      = 6,
  width                       = 10
) {
  require(ggplot2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }

  cat("sampleSummaries: preparing data\n")
  typableGT <- geno(vcf)[["GT"]]

  if("missingness" %in% variableNames) {
    cat("sampleSummaries: determining missingness\n")
    typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
    missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
    names(missingnessPerSample) <- dimnames(typableGT)[[2]]
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
    sampleShapes <- rep("uncontaminated", dim(typableGT)[2])
    names(sampleShapes) <- dimnames(typableGT)[[2]]
    sampleShapes[names(sampleShapes) %in% pfContaminants] <- "Pf contaminated"
    pdf(paste(pdfFilestem, "missingnessBySample", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        factor(names(sort(missingnessPerSample)), levels=names(sort(missingnessPerSample))),
        sort(missingnessPerSample),
        main="Missingness per sample",
        xlab="Sample ID",
        ylab="Missingness",
        stat="identity",
        colour=sampleShapes
      )
      + theme_bw()
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8))
      + theme(axis.title.x = element_blank())
      + theme(legend.position = c(0.9,0.4))
      + scale_colour_manual(name="Pf contamination\nstatus", values=c("red", "black"))
      + geom_hline(yintercept = 0.05, colour="blue")
      + geom_hline(yintercept = 0.2, colour="green")
      + geom_hline(yintercept = 0.65, colour="orange")
      + geom_hline(yintercept = 0.95, colour="red")
    )
    dev.off()
  }

  if("heterozgosity" %in% variableNames) {
    cat("sampleSummaries: determining heterozgosity\n")
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
    cat("sampleSummaries: determining singletons\n")
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
    cat("sampleSummaries: determining depth\n")
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
  
  if("missingnessVsDepth" %in% variableNames) {
    typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT))
    missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
    names(missingnessPerSample) <- dimnames(typableGT)[[2]]
    meanDPperSample <- colMeans(geno(vcf)[["DP"]])
    pdf(paste(pdfFilestem, "missingnessVsDepth", "pdf", sep="."), height=height, width=width)
    print(
      qplot(
        missingnessPerSample,
        meanDPperSample,
        main="Missingness vs mean depth per sample",
        xlab="Missingness",
        ylab="Mean depth"
      )
      + theme_bw()
      + annotate("text", label = "PQ0002-C", x = missingnessPerSample["PQ0002-C"]+0.06, y = meanDPperSample["PQ0002-C"], size=4)
    )
    dev.off()
  }
  
  return(vcf)
}

