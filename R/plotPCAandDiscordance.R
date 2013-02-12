# plotPCA2.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################




plotPCAandDiscordance <- function(
  vcf                         = loadAndGenotypePvVcf(),
#  parameterSets               = list(
#    "0.35_0_all" = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=0, samplesToRemove=NULL)
#  ),
  parameterSets               = list(
    "0.35_0_all"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=0,  samplesToRemove=NULL),
    "0.75_6_all"    = list(sampleMissingnessThreshold = 0.75, variantMissingnessThreshold=6, samplesToRemove=NULL)
#    "0.95_10_all"   = list(sampleMissingnessThreshold = 0.95, variantMissingnessThreshold=10, samplesToRemove=NULL),
#    "1.00_20_all"   = list(sampleMissingnessThreshold = 1.00, variantMissingnessThreshold=20, samplesToRemove=NULL),
#    "0.35_1_all"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=1,  samplesToRemove=NULL),
#    "0.35_2_all"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=2,  samplesToRemove=NULL),
#    "0.35_5_all"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=5,  samplesToRemove=NULL),
#    "0.35_10_all"   = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=10, samplesToRemove=NULL),
#    "0.35_0_noOutliers"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=0,  samplesToRemove=c(c("PH0184-C", "PH0190-C", "PH0309-C", "PH0313-C", "PH0315-C", "PH0319-C"), c("PH0189-C", "PH0312-C", "PH0318-C")))
#    "0.35_0_SEA_noOutliers"    = list(sampleMissingnessThreshold = 0.35, variantMissingnessThreshold=0,  samplesToRemove=c(c("PH0184-C", "PH0190-C", "PH0309-C", "PH0313-C", "PH0315-C", "PH0319-C"), c("PH0189-C", "PH0312-C", "PH0318-C")), countries=c("Thailand", "Cambodia", "Vietnam"))
  ),
  pdfFilestem                 = "analysis/pca/pv_02",
  pfContaminants              = c("PD0173-C", "PH0309-C", "PJ0002-Cx", "PN0094-C", "PN0095-C", "PN0096-C"),
  pvMOI                       = NULL,
  pcsToPlot                   = list(c(1,2), c(1,3), c(2,3))
) {
  require(ggplot2)
  require(reshape2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }
  
  cat("plotPCAandDiscordance: preparing data\n")
  typableGT <- geno(vcf)[["GT"]]
  typableMissingGT <- matrix(typableGT=="./.", ncol=ncol(typableGT), dimnames=dimnames(typableGT))
  missingnessPerSample <- colSums(typableMissingGT)/dim(typableGT)[1]
  GTdosagesInAllsamples <- (typableGT=="0/0") * -1 + (typableGT=="0/1") * 0 + (typableGT=="1/1") * 1
  cat("plotPCAandDiscordance: determining sample countries\n")
  sampleCountries <- countryCodes(dimnames(vcf)[[2]])
  rm(vcf)
  gc()
  
  resultsList <- sapply(
    names(parameterSets),
    function(parameterSetName) {
      cat("\n", parameterSetName, "\n----------------------------\n")
      pdfFilestemExtended <- paste(pdfFilestem, parameterSetName, sep=".")
      samplesToUse <- setdiff(
        dimnames(typableMissingGT)[[2]][missingnessPerSample<parameterSets[[parameterSetName]][["sampleMissingnessThreshold"]]],
        parameterSets[[parameterSetName]][["samplesToRemove"]]
      )
      if(!is.null(parameterSets[[parameterSetName]][["countries"]])) {
        samplesToUse <- samplesToUse[sampleCountries[samplesToUse] %in% parameterSets[[parameterSetName]][["countries"]]]
      }
      missingnessPerVariant <- rowSums(typableMissingGT[, samplesToUse])
      
#      Create PCA plot
    
      GTdosagesInNonMissingSamplesAndVariants <- GTdosagesInAllsamples[
        missingnessPerVariant<=parameterSets[[parameterSetName]][["variantMissingnessThreshold"]],
        samplesToUse
      ]
      cat("calculating PCA results\n")
      pcaResults <- prcomp(t(GTdosagesInNonMissingSamplesAndVariants))
      sampleShapes <- rep("uncontaminated", dim(GTdosagesInNonMissingSamplesAndVariants)[2])
      names(sampleShapes) <- dimnames(GTdosagesInNonMissingSamplesAndVariants)[[2]]
      sampleShapes[names(sampleShapes) %in% pfContaminants] <- "Pf contaminated"
      sampleShapes[names(sampleShapes) %in% pvMOI] <- "Pv MOI"
      cat("creating PCA plots\n")
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
      
#      Create discordance histogram and matrix
      GT <- typableGT[
        missingnessPerVariant<=parameterSets[[parameterSetName]][["variantMissingnessThreshold"]],
        samplesToUse
      ]
      GT[GT=="./."] <- NA
      dimnames(GT)[[2]] <- paste(samplesToUse, " (", sampleCountries[samplesToUse], ")", sep="")
      
      gc()
      cat("calculating discordance matrix\n")
      GTDiscordanceMatrix <- discordanceMatrix(GT)
      
      cat("creating discordance plots\n")
      pdf(paste(pdfFilestemExtended, "pairwiseDiscordanceHistogram.pdf", sep="."), height=4, width=6)
      print(
        qplot(
          as.vector(GTDiscordanceMatrix),
          xlab="Number of discordant SNPs between pairwise sample comparisons",
          ylab="Frequency (number of sample pairs)"
        ) +
#        geom_vline(xintercept = discordanceThreshold, colour="red") +
        theme_bw()
      )
      dev.off()
      discordanceDF <- melt(GTDiscordanceMatrix, value.name="Discordances")
#      discordanceDF[["sample1"]] <- sampleNames[as.character(discordanceDF[["Var1"]])]
#      discordanceDF[["sample2"]] <- sampleNames[as.character(discordanceDF[["Var2"]])]
      pdf(paste(pdfFilestemExtended, "discordanceHeatmap.pdf", sep="."), height=10, width=12)
      print(
        ggplot(
          discordanceDF,
  #        melt(GTsIntDiscordanceMatrix, value.name="Discordances"),
          aes(x=Var1, y=Var2, fill=Discordances)
        )
        + geom_tile()
        + scale_fill_gradient2(low="red", high="blue", midpoint=median(GTDiscordanceMatrix, na.rm=TRUE))
  #      + scale_fill_gradient(low="yellow", high="red")
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        + theme(axis.title.x = element_blank())
        + theme(axis.title.y = element_blank())
      )
      dev.off()

      return(list(pcaResults=pcaResults, discordanceMatrix=GTDiscordanceMatrix))
    },
    USE.NAMES=TRUE,
    simplify=FALSE
  )
  return(resultsList)
}

