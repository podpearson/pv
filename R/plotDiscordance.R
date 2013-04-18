# plotDiscordance.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


plotDiscordance <- function(
  vcf                         = loadAndGenotypePvVcf(),
  sampleManifestRdaFilename   = "meta/pv_1.0_sampleManifest.rda",
  minimumDepthToCallSNP       = 1,
  sampleMissingnessThreshold  = 0.5,
  maxMissingnessPerVariant    = 0,
  vcfName                     = "1.0",
  pdfFilestem                 = paste("analysis/", vcfName, "/discordance/", vcfName, sep=""),
  shouldAttemptAllSamples     = FALSE,
  shouldUsePGVlikeCalls       = TRUE
) {
  require(gplots)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }

  load(sampleManifestRdaFilename)
  RefReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[1]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  FirstAltReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[2]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  majorityAlleleCalls <- ifelse(RefReads>FirstAltReads, "0", ifelse(FirstAltReads>RefReads, "1", NA))
  majorityAlleleCalls[(RefReads+FirstAltReads) < minimumDepthToCallSNP] <- NA
  missingnessPerSample <- apply(majorityAlleleCalls, 2, function(x) length(which(is.na(x)))/length(x))
  if(shouldAttemptAllSamples) {
    allDiscordances <- discordanceMatrix(majorityAlleleCalls)
    allDiscordanceRatios <- discordanceMatrix(majorityAlleleCalls, returnValue="ratio")
    allDiscordanceProportions <- discordanceMatrix(majorityAlleleCalls, returnValue="proportion")
    allSharedGenotypes <-  allDiscordances/allDiscordanceProportions
    pdf(paste(pdfFilestem, "allDiscordances.pdf", sep="."), height=12, width=12)
    heatmap.2(allDiscordances, margins=c(10, 10))
    dev.off()
    pdf(paste(pdfFilestem, "allDiscordanceRatios.pdf", sep="."), height=12, width=12)
    heatmap.2(allDiscordanceRatios, margins=c(10, 10))
    dev.off()
    pdf(paste(pdfFilestem, "allDiscordanceProportions.pdf", sep="."), height=12, width=12)
    heatmap.2(allDiscordanceRatios, margins=c(10, 10))
    dev.off()
    allDiscordances["PH0190-C", ]
    allDiscordanceRatios["PH0190-C", ]
    allDiscordanceProportions["PH0190-C", ]
    allSharedGenotypes["PH0190-C", ]
    allDiscordances[names(which.max(missingnessPerSample)), ]
    allDiscordanceRatios[names(which.max(missingnessPerSample)), ]
    allDiscordanceProportions[names(which.max(missingnessPerSample)), ]
    allSharedGenotypes[names(which.max(missingnessPerSample)), ]
  }
  if(shouldUsePGVlikeCalls) {
    GT <- geno(vcf)[["GT"]]
    GTint <- matrix(ifelse(GT=="0/0", 0, ifelse(GT=="1/1", 1, NA)), ncol=ncol(GT), dimnames=dimnames(GT))
    GTmissingnessPerSample <- apply(GTint, 2, function(x) length(which(is.na(x)))/length(x))
    samplesToUse <- names(GTmissingnessPerSample[GTmissingnessPerSample<sampleMissingnessThreshold])
    browser()
    PGVlikeDiscordanceProportions <- discordanceMatrix(GTint[, samplesToUse], returnValue="proportion")
    PGVlikeDiscordanceProportions2 <- PGVlikeDiscordanceProportions
    diag(PGVlikeDiscordanceProportions2) <- NA
    stem(PGVlikeDiscordanceProportions)
    stem(PGVlikeDiscordanceProportions2)
    
#    dimnames(PGVlikeDiscordanceProportions)[[1]] <- paste(samplesToUse, " (", sampleManifest[samplesToUse, "richard_donor_source_code"], ", ", sampleManifest[samplesToUse, "country"], ")", sep="")
#    dimnames(MACDiscordanceMatrixExpandedSampleIDs)[[2]] <- paste(samplesToUse, " (", sampleManifest[samplesToUse, "richard_donor_source_code"], ", ", sampleManifest[samplesToUse, "country"], ")", sep="")
    pdf(paste(pdfFilestem, "PGVlikeDiscordanceProportions.pdf", sep="."), height=12, width=12)
    heatmap.2(PGVlikeDiscordanceProportions, margins=c(10, 10))
    dev.off()
  }
  samplesToUse <- names(missingnessPerSample[missingnessPerSample<sampleMissingnessThreshold])
  missingnessPerVariant <- apply(majorityAlleleCalls[, samplesToUse], 1, function(x) length(which(is.na(x))))
  cat("plotDiscordance: calculating discordance matrix\n")
  MACDiscordanceMatrix <- discordanceMatrix(majorityAlleleCalls[missingnessPerVariant<=maxMissingnessPerVariant, samplesToUse])
  dimnames(MACDiscordanceMatrix)[[1]] <- paste(samplesToUse, " (", sampleManifest[samplesToUse, "richard_donor_source_code"], ", ", sampleManifest[samplesToUse, "country"], ")", sep="")
  dimnames(MACDiscordanceMatrix)[[2]] <- paste(samplesToUse, " (", sampleManifest[samplesToUse, "richard_donor_source_code"], ", ", sampleManifest[samplesToUse, "country"], ")", sep="")
  pdf("~/PvDiscordanceHeatmap.pdf", height=18, width=18)
  heatmap(MACDiscordanceMatrix)
  dev.off()
  pdf("~/PvDiscordanceHeatmap2.pdf", height=12, width=12)
  heatmap.2(MACDiscordanceMatrix)
  dev.off()
  
  dimnames(MACDiscordanceMatrix)[[1]] <- samplesToUse
  dimnames(MACDiscordanceMatrix)[[2]] <- samplesToUse
  
  lowDiscordancePairs <- matrix(dimnames(MACDiscordanceMatrix)[[1]][which(MACDiscordanceMatrix < 500, arr.ind=TRUE)], ncol=2)
  lowDiscordancePairsDF <- cbind(data.frame(lowDiscordancePairs), data.frame(discordantSNPs = MACDiscordanceMatrix[lowDiscordancePairs]))
  duplicateSamplePairs <- matrix(lowDiscordancePairs[lowDiscordancePairs[, 1] != lowDiscordancePairs[, 2]], ncol=2)
  duplicateSamplePairsDF <- cbind(
    data.frame(duplicateSamplePairs),
    data.frame(
      discordantSNPs = MACDiscordanceMatrix[duplicateSamplePairs],
      coverageFirst=sampleManifest[duplicateSamplePairs[, 1], "coverage"],
      coverageSecond=sampleManifest[duplicateSamplePairs[, 2], "coverage"]
    )
  ) 
  
  duplicateSamplePairs <- cbind(
    duplicateSamplePairs,
    (
      (sampleManifest[duplicateSamplePairs[, 1], "coverage"] < sampleManifest[duplicateSamplePairs[, 2], "coverage"])
    )+1
  )
  duplicateSamplePairs <- cbind(
    duplicateSamplePairs,
    (
      (heterozygosityPerSample[duplicateSamplePairs[, 1]] < heterozygosityPerSample[duplicateSamplePairs[, 2]]) |
      (heterozygosityPerSample[duplicateSamplePairs[, 1]] == heterozygosityPerSample[duplicateSamplePairs[, 2]] & duplicateSamplePairs[, 1] < duplicateSamplePairs[, 2])
    )+1
  )

  uniqueSamplePairs <- unique(apply(duplicateSamplePairs, 1, function(x) paste(sort(x), collapse="_and_")))

  duplicateSamplesWithHighestCoverage <- c("PH0190-C", "PH0188-C", "PH0189-C", "PH0310-C", "PJ0004-C", "PJ0006-CW", "PN0065-C", "PV0067-C", "PX0001-CW")
  duplicateSamplesToRemove <- setdiff(unique(c(as.character(duplicateSamplePairsDF[[1]]), as.character(duplicateSamplePairsDF[[2]]))), duplicateSamplesWithHighestCoverage)
  samplesToUseInPopgen <- setdiff(samplesToUse, duplicateSamplesToRemove)
  save(samplesToUseInPopgen, file="meta/samplesToUseInPopgen_depth1_sampleMissingness0.5_SNPmissingness0.rda")
  save(duplicateSamplesToRemove, file="meta/duplicateSamplesToRemove_depth1_sampleMissingness0.5_SNPmissingness0.rda")

  browser()
}
