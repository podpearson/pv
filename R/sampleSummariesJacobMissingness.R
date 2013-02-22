# sampleSummariesJacobMissingness.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# vcf_01 <- sampleSummariesJacobMissingness(loadAndGenotypePvVcf("data/genotypes/out.pv_1.0/annotated.vcf.gz"), pdfFilestem = "analysis/sampleSummaries/pv_01")

sampleSummariesJacobMissingness <- function(
  pdfFilestem                 = "analysis/sampleSummaries/pv_1.0",
  sampleManifestRdaFilename   = "meta/pv_1.0_sampleManifest.rda",
  magnusContaminantsFilename  = "analysis/magnus/Vivax contamination.txt",
  jacobMissingnessFilename    = "analysis/jacob/sample_missingness.txt",
  plotsToCreate               = c("missingnessBySample", "vsMissingness", "vsMissingnessLog", "missingnessVs", "numericScatter", "numericScatterLog"),
  height                      = 6,
  width                       = 10
) {
  require(ggplot2)
  if(!file.exists(dirname(pdfFilestem))) {
    dir.create(dirname(pdfFilestem), recursive=TRUE)
  }
  
  cat("sampleSummaries: loading sample manifest\n")
  load(sampleManifestRdaFilename)
  magnusContaminants <- read.delim(magnusContaminantsFilename, as.is=TRUE, row.names=1)
  jacobMissingness <- read.delim(jacobMissingnessFilename, as.is=TRUE, row.names=1)
  sampleManifestExtended <- cbind(
    sampleManifest,
    magnusContaminants[row.names(sampleManifest),],
    jacobMissingness[row.names(sampleManifest),,drop=FALSE]
  )
  
  sampleManifestExtendedReordered <- sampleManifestExtended[order(sampleManifestExtended[["missingness"]]), ]
  sampleManifestExtendedReordered[["ox_code"]] <- factor(sampleManifestExtendedReordered[["ox_code"]], levels=sampleManifestExtendedReordered[["ox_code"]])
  
#  pdf(paste(pdfFilestem, "Total reads vs missingness", "pdf", sep="."), height=height, width=width)
#  print(
#    qplot(
#      raw_reads,
#      missingness,
#      data=sampleManifestExtendedReordered,
#      xlab="Total number of reads",
#      ylab="Missingness",
#      log="xy"
#    )
#    + theme_bw()
#    + geom_hline(yintercept = 0.05, colour="blue")
#    + geom_hline(yintercept = 0.2, colour="green")
#    + geom_hline(yintercept = 0.65, colour="orange")
#    + geom_hline(yintercept = 0.95, colour="red")
#  )
#  dev.off()
#  pdf(paste(pdfFilestem, "Reads mapped vs missingness", "pdf", sep="."), height=height, width=width)
#  print(
#    qplot(
#      reads_mapped,
#      missingness,
#      data=sampleManifestExtendedReordered,
#      xlab="Number of reads mapped",
#      ylab="Missingness",
#      log="xy"
#    )
#    + theme_bw()
#    + geom_hline(yintercept = 0.05, colour="blue")
#    + geom_hline(yintercept = 0.2, colour="green")
#    + geom_hline(yintercept = 0.65, colour="orange")
#    + geom_hline(yintercept = 0.95, colour="red")
#  )
#  dev.off()
#  pdf(paste(pdfFilestem, "Percentage reads mapped vs missingness", "pdf", sep="."), height=height, width=width)
#  print(
#    qplot(
#      pct_reads_mapped,
#      missingness,
#      data=sampleManifestExtendedReordered,
#      xlab="Percentage of reads mapped",
#      ylab="Missingness",
#      log="xy"
#    )
#    + theme_bw()
#    + geom_hline(yintercept = 0.05, colour="blue")
#    + geom_hline(yintercept = 0.2, colour="green")
#    + geom_hline(yintercept = 0.65, colour="orange")
#    + geom_hline(yintercept = 0.95, colour="red")
#  )
#  dev.off()
  
  sapply(
    names(sampleManifestExtended),
    function(sampleManifestColumnName) {
      cat(sampleManifestColumnName)
      if("missingnessBySample" %in% plotsToCreate) {
        pdf(file.path(pdfFilestem, "missingnessBySample", paste("missingnessBySample", sampleManifestColumnName, "pdf", sep=".")), height=height, width=width)
        print(
          qplot(
            ox_code,
            missingness,
            data=sampleManifestExtendedReordered,
            main=paste("Missingness per sample - coloured by", sampleManifestColumnName),
            xlab="Sample ID",
            ylab="Missingness",
            stat="identity",
            colour=sampleManifestExtendedReordered[[sampleManifestColumnName]]
          )
          + theme_bw()
          + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8))
          + theme(axis.title.x = element_blank())
          + theme(legend.position = c(0.9,0.4))
          + guides(colour = guide_legend(sampleManifestColumnName))
          + geom_hline(yintercept = 0.05, colour="blue")
          + geom_hline(yintercept = 0.2, colour="green")
          + geom_hline(yintercept = 0.65, colour="orange")
          + geom_hline(yintercept = 0.95, colour="red")
        )
        dev.off()
      }
      if("vsMissingness" %in% plotsToCreate) {
        if(class(sampleManifestExtendedReordered[[sampleManifestColumnName]]) %in% c("integer", "numeric")) {
          pdf(file.path(pdfFilestem, "vsMissingness", paste("vsMissingness", sampleManifestColumnName, "pdf", sep=".")), height=height, width=width)
#          pdf(paste(pdfFilestem, sampleManifestColumnName, "vs.missingness", "pdf", sep="."), height=height, width=width)
          print(
            qplot(
              sampleManifestExtendedReordered[[sampleManifestColumnName]],
              missingness,
              data=sampleManifestExtendedReordered,
              xlab=sampleManifestColumnName,
              ylab="Missingness"
            )
            + theme_bw()
            + geom_hline(yintercept = 0.05, colour="blue")
            + geom_hline(yintercept = 0.2, colour="green")
            + geom_hline(yintercept = 0.65, colour="orange")
            + geom_hline(yintercept = 0.95, colour="red")
          )
          dev.off()
        }
      }
      if("vsMissingnessLog" %in% plotsToCreate) {
        if(class(sampleManifestExtendedReordered[[sampleManifestColumnName]]) %in% c("integer", "numeric")) {
          pdf(file.path(pdfFilestem, "vsMissingnessLog", paste("vsMissingnessLog", sampleManifestColumnName, "pdf", sep=".")), height=height, width=width)
#          pdf(paste(pdfFilestem, sampleManifestColumnName, "vs.missingness", "pdf", sep="."), height=height, width=width)
          print(
            qplot(
              sampleManifestExtendedReordered[[sampleManifestColumnName]],
              missingness,
              data=sampleManifestExtendedReordered,
              xlab=sampleManifestColumnName,
              ylab="Missingness",
              log="xy"
            )
            + theme_bw()
            + geom_hline(yintercept = 0.05, colour="blue")
            + geom_hline(yintercept = 0.2, colour="green")
            + geom_hline(yintercept = 0.65, colour="orange")
            + geom_hline(yintercept = 0.95, colour="red")
          )
          dev.off()
        }
      }
      if("numericScatter" %in% plotsToCreate) {
        if(class(sampleManifestExtendedReordered[[sampleManifestColumnName]]) %in% c("integer", "numeric")) {
          sapply(
            names(sampleManifestExtended),
            function(sampleManifestColumnName2) {
              if(class(sampleManifestExtendedReordered[[sampleManifestColumnName2]]) %in% c("integer", "numeric")) {
                pdf(file.path(pdfFilestem, "numericScatter", paste(sampleManifestColumnName, "vs", sampleManifestColumnName2, "pdf", sep=".")), height=height, width=width)
                print(
                  qplot(
                    sampleManifestExtendedReordered[[sampleManifestColumnName]],
                    sampleManifestExtendedReordered[[sampleManifestColumnName2]],
                    data=sampleManifestExtendedReordered,
                    xlab=sampleManifestColumnName,
                    ylab=sampleManifestColumnName2
                  )
                  + theme_bw()
                  + geom_hline(yintercept = 0.05, colour="blue")
                  + geom_hline(yintercept = 0.2, colour="green")
                  + geom_hline(yintercept = 0.65, colour="orange")
                  + geom_hline(yintercept = 0.95, colour="red")
                )
                dev.off()
              }
            }
          )
        }
      }
      if("numericScatterLog" %in% plotsToCreate) {
        if(class(sampleManifestExtendedReordered[[sampleManifestColumnName]]) %in% c("integer", "numeric")) {
          sapply(
            names(sampleManifestExtended),
            function(sampleManifestColumnName2) {
              if(class(sampleManifestExtendedReordered[[sampleManifestColumnName2]]) %in% c("integer", "numeric")) {
                pdf(file.path(pdfFilestem, "numericScatterLog", paste(sampleManifestColumnName, "vs", sampleManifestColumnName2, "pdf", sep=".")), height=height, width=width)
                print(
                  qplot(
                    sampleManifestExtendedReordered[[sampleManifestColumnName]],
                    sampleManifestExtendedReordered[[sampleManifestColumnName2]],
                    data=sampleManifestExtendedReordered,
                    xlab=sampleManifestColumnName,
                    ylab=sampleManifestColumnName2,
                    log="xy"
                  )
                  + theme_bw()
                  + geom_hline(yintercept = 0.05, colour="blue")
                  + geom_hline(yintercept = 0.2, colour="green")
                  + geom_hline(yintercept = 0.65, colour="orange")
                  + geom_hline(yintercept = 0.95, colour="red")
                )
                dev.off()
              }
            }
          )
        }
      }
      if("missingnessVs" %in% plotsToCreate) {
        pdf(file.path(pdfFilestem, "missingnessVs", paste("missingnessVs", sampleManifestColumnName, "pdf", sep=".")), height=height, width=width)
#        pdf(paste(pdfFilestem, "missingness.vs", sampleManifestColumnName, "pdf", sep="."), height=height, width=width)
        print(
          qplot(
            sampleManifestExtendedReordered[[sampleManifestColumnName]],
            missingness,
            data=sampleManifestExtendedReordered,
            main=paste("Missingness vs", sampleManifestColumnName),
            ylab="Missingness",
            colour=name,
            shape=country
          )
          + theme_bw()
          + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8))
          + theme(axis.title.x = element_blank())
          + theme(legend.position = c(0.9,0.4))
  #          + scale_colour_manual(name="Pf contamination\nstatus", values=c("red", "black"))
          + geom_hline(yintercept = 0.05, colour="blue")
          + geom_hline(yintercept = 0.2, colour="green")
          + geom_hline(yintercept = 0.65, colour="orange")
          + geom_hline(yintercept = 0.95, colour="red")
        )
        dev.off()
      }
    }
  )
}


