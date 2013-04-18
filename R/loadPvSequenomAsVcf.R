# loadPvSequenomAsVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadPvSequenomAsVcf <- function(
  pvSequenomFilename          = "data/sequenom/vivax_data_10APR2013.txt",
  pvSequenomVcfFilename       = sub("\\.txt", "\\.vcf", pvSequenomFilename),
  pvSequenomRdaFilename       = paste(pvSequenomVcfFilename, ".rda", sep=""),
  fastaFilename               = "data/genome/PvivaxGenomic_PlasmoDB-6.0.shortnames.fasta",
  reload                      = FALSE,
  shouldSaveVcfFile           = TRUE, # Note that due to a bug in VariantAnnotation it is necessary to have more than one geno, and writeVcf will write strange VCF if only, e.g. GT is included (email sent to Bioconductor list 18/04/2013)
  shouldSaveRdaFile           = TRUE
) {
  if(reload || !file.exists(pvSequenomRdaFilename)) {
    require(reshape2)
    pvSequenom <- read.table(pvSequenomFilename, header=TRUE, sep="\t", as.is=TRUE)
    pvSequenomPvOnly <- subset(pvSequenom, !(study_code %in% c("WT", "WA")))
    pvSingleSample <- subset(pvSequenomPvOnly, sample_code==pvSequenomPvOnly[1, "sample_code"])
    pvRowData  = GRanges(
      seqnames = sub("_[0-9]+$", "", pvSingleSample[["sequence_code"]]),
      ranges   = IRanges(
        start = as.integer(sub("^.*_([0-9]+)$", "\\1", pvSingleSample[["sequence_code"]])),
        width = 1,
        names  = pvSingleSample[["sequence_code"]]
      ),
      paramRangeID = factor(rep(NA, dim(pvSingleSample)[1]))
    )
    
    ref <- readDNAStringSet(fastaFilename)
    refBase <- function(ref, chromosome, position) {
      cat(".")
      if(chromosome %in% names(ref)) {
        as.character(subseq(ref[chromosome], position, position))
      } else {
        NA
      }
    }
    RefBases <- mapply(
      refBase,
      as.character(seqnames(pvRowData)),
      start(pvRowData),
      MoreArgs = list(ref=ref)
    )
    RefBases[which(is.na(RefBases))] <- pvSingleSample[which(is.na(RefBases)), "reference_allele"]

    pvSequenomAllele1s <- acast(pvSequenomPvOnly, sequence_code ~ sample_code, value.var="allele1")[names(pvRowData), ]
    pvSequenomAllele2s <- acast(pvSequenomPvOnly, sequence_code ~ sample_code, value.var="allele2")[names(pvRowData), ]
    AltBases <- ifelse(RefBases==pvSequenomAllele1s[, 1], pvSequenomAllele2s[, 1], pvSequenomAllele1s[, 1])
    pvSequenomRefs <- matrix(RefBases, nrow=length(RefBases), ncol=dim(pvSequenomAllele1s)[[2]], dimnames=dimnames(pvSequenomAllele1s))
    pvSequenomAlts <- matrix(AltBases, nrow=length(AltBases), ncol=dim(pvSequenomAllele2s)[[2]], dimnames=dimnames(pvSequenomAllele2s))
    pvSequenomGTs <- acast(pvSequenomPvOnly, sequence_code ~ sample_code, value.var="result_tidied")[names(pvRowData), ]
    pvSequenomAlleleRatio1s <- acast(pvSequenomPvOnly, sequence_code ~ sample_code, value.var="allele_ratio1")[names(pvRowData), ]
    pvSequenomAlleleRatio2s <- acast(pvSequenomPvOnly, sequence_code ~ sample_code, value.var="allele_ratio2")[names(pvRowData), ]
    GT <- ifelse(
      is.na(pvSequenomGTs),
      "./.",
      ifelse(
        pvSequenomGTs==pvSequenomRefs, "0/0",
        ifelse(
          pvSequenomGTs==pvSequenomAlts,
          "1/1",
          ifelse(
            pvSequenomGTs == paste(pvSequenomRefs, pvSequenomAlts, sep=""),
            "0/1",
            "./."
          )
        )
      )
    )
    ARref <- ifelse(
      is.na(pvSequenomAlleleRatio1s),
      NA,
      ifelse(
        pvSequenomAllele1s==pvSequenomRefs, pvSequenomAlleleRatio1s, pvSequenomAlleleRatio2s 
      )
    )
    ARalt <- ifelse(
      is.na(pvSequenomAlleleRatio1s),
      NA,
      ifelse(
        pvSequenomAllele1s==pvSequenomAlts, pvSequenomAlleleRatio1s, pvSequenomAlleleRatio2s 
      )
    )
    sampleIDs <- dimnames(pvSequenomGTs)[[2]]
    vcf <- VCF(
      rowData = pvRowData,
      colData = DataFrame(
        Samples = seq(along=sampleIDs),
        row.names = sampleIDs
      ),
      exptData = SimpleList(
        samples = sampleIDs,
        header = VCFHeader(
          samples = sampleIDs,
          header  = DataFrameList(
            META = rbind(
              DataFrame(Value = "VCFv4.1", row.names="fileformat"),
              DataFrame(Value = "PvSequenom", row.names="ProjectName")
            ),
            FILTER = DataFrame(Descrption="PASS", row.names="PASS"),
            FORMAT = rbind(
              DataFrame(Number = "1", Type="String", Description="Genotype", row.names="GT"),
              DataFrame(Number = "1", Type="Float", Description="Proportion for ref allele", row.names="ARref"),
              DataFrame(Number = "1", Type="Float", Description="Proportion for alt allele", row.names="ARalt")
            ),
            INFO = rbind(
              DataFrame(Number = "1", Type = "Integer", Description="chromosome", row.names="chr_valid"),
              DataFrame(Number = "1", Type = "Integer", Description="chrmosome position", row.names="coord_valid"),
              DataFrame(Number = "1", Type = "String", Description="gene name", row.names="gene_symbol"),
              DataFrame(Number = "1", Type = "String", Description="reference allele", row.names="reference_allele"),
              DataFrame(Number = "1", Type = "String", Description="alternate allele", row.names="alternative_allele"),
              DataFrame(Number = "1", Type = "String", Description="SNP alleles single letter code", row.names="single_letter_code")
            )
          )
        )
      ),
      fixed    = DataFrame(
        REF = DNAStringSet(pvSequenomRefs[, 1]),
        ALT = do.call(DNAStringSetList, as.list(pvSequenomAlts[, 1])),
        QUAL = rep(0.0, dim(pvSequenomRefs)[1]),
        FILTER = rep("PASS", dim(pvSequenomRefs)[1])
      ),
      info = DataFrame(
        pvSingleSample[, c("chr_valid", "coord_valid", "gene_symbol", "reference_allele", "alternative_allele", "single_letter_code")]
      ),
      geno     = SimpleList(
        GT = GT,
        ARref = ARref,
        ARalt = ARalt
      )
    )
    genome <- rep("P. vivax reference, PlasmoDB V6.0", length(seqlevels(vcf)))
    names(genome) <- seqlevels(vcf)
    genome(vcf) <- genome
    vcf <- vcf[order(vcf)]
    if(shouldSaveRdaFile) {
      save(vcf, file=pvSequenomRdaFilename)
    }
  } else {
    load(pvSequenomRdaFilename)
  }
  if(shouldSaveVcfFile) {
    writeVcf(vcf, pvSequenomVcfFilename, index=TRUE)
  }
  return(vcf)
}
