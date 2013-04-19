# loadAndGenotypePvVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndGenotypePvVcf <- function(
  originalVcfFilename         = "data/genotypes/pv_02.vcf.gz",
  typableVcfFilename          = sub("\\.vcf\\.gz", "\\.typable\\.vcf", originalVcfFilename),
  typableRdaFilename          = paste(typableVcfFilename, "rda", sep="."),
  typableExtraInfoVcfFilename = sub("\\.vcf\\.gz", "\\.typable\\.extraInfo\\.vcf", originalVcfFilename),
  typableExtraInfoRdaFilename = sub("\\.vcf\\.gz", "\\.typable\\.extraInfo\\.vcf\\.rda", originalVcfFilename),
  reload                      = FALSE
) {
  if(!reload && file.exists(typableExtraInfoRdaFilename)) {
    load(typableExtraInfoRdaFilename)
  } else {
    if(!reload && file.exists(typableRdaFilename)) {
      load(typableRdaFilename)
    } else {
      if(reload || !file.exists(typableVcfFilename)) {
        cleanVcfCommand <- paste("zcat", originalVcfFilename, "| scripts/perl/cleanAndGenotypeMagnusVcf.pl >", typableVcfFilename)
        system(cleanVcfCommand)
      }
      typableVcf <- readVcf(typableVcfFilename, genome="P. vivax reference, PlasmoDB V6.0")
      save(typableVcf, file=typableRdaFilename)
    }
    typableGT <- geno(typableVcf)[["GT"]]
    typableHetGT <- matrix(typableGT=="0/1", ncol=ncol(typableGT))
    typableHomRefGT <- matrix(typableGT=="0/0", ncol=ncol(typableGT))
    typableHomAltGT <- matrix(typableGT=="1/1", ncol=ncol(typableGT))
    AAs <- rowSums(typableHomAltGT) + 0.5 * rowSums(typableHetGT)
    RAs <- rowSums(typableHomRefGT) + 0.5 * rowSums(typableHetGT)
    AAFs <- AAs/(AAs+RAs)
    SNPhet <- rowSums(typableHetGT) / (AAs+RAs)
    newExptData <- SimpleList(
      "samples" = samples(exptData(typableVcf)[["header"]]),
      "header" = VCFHeader(
        header=DataFrameList(
          META = meta(exptData(typableVcf)[["header"]]),
          FILTER = fixed(exptData(typableVcf)[["header"]])[["FILTER"]],
          INFO = rbind(
            header(exptData(typableVcf)[["header"]])[["INFO"]],
            DataFrame(Number=1, Type="Integer", Description="Homozygous ref allele count", row.names="HRAC"),
            DataFrame(Number=1, Type="Integer", Description="Homozygous alt allele count", row.names="HAAC"),
            DataFrame(Number=1, Type="Integer", Description="Heterozygous samples count", row.names="HetC"),
            DataFrame(Number=1, Type="Float", Description="Alternative allele frequency", row.names="AAF"),
            DataFrame(Number=1, Type="Float", Description="Minor allele frequency", row.names="MAF"),
            DataFrame(Number=1, Type="Float", Description="Heterozygosity", row.names="heterozygosity"),
            DataFrame(Number=0, Type="Flag", Description="Singleton variant (only 1 alt allele and 0 hets)", row.names="singleton")
          ),
          FORMAT = geno(exptData(typableVcf)[["header"]])
        )
      )
    )
    typableVcfWithAAFandHet <- VCF(
      rowData = rowData(typableVcf)[info(typableVcf)[["TYP"]]==TRUE],
      colData = colData(typableVcf),
      exptData = newExptData,
      fixed = fixed(typableVcf)[info(typableVcf)[["TYP"]]==TRUE, ],
      info = cbind(
        info(typableVcf)[info(typableVcf)[["TYP"]]==TRUE, ],
        DataFrame(
          HRAC = rowSums(typableHomRefGT),
          HAAC = rowSums(typableHomAltGT),
          HetC = rowSums(typableHetGT),
          AAF = AAFs,
          MAF=pmin(AAFs, 1-AAFs),
          heterozygosity=SNPhet,
          singleton = rowSums(typableHomAltGT)==1 & rowSums(typableHetGT)==0
        )
      ),
      geno=geno(typableVcf)
    )
    RefReads <- matrix(
      sapply(geno(typableVcfWithAAFandHet)[["AD"]], function(x) x[1]),
      ncol=dim(geno(typableVcfWithAAFandHet)[["AD"]])[2],
      dimnames=dimnames(geno(typableVcfWithAAFandHet)[["AD"]])
    )
    FirstAltReads <- matrix(
      sapply(geno(typableVcfWithAAFandHet)[["AD"]], function(x) x[2]),
      ncol=dim(geno(typableVcfWithAAFandHet)[["AD"]])[2],
      dimnames=dimnames(geno(typableVcfWithAAFandHet)[["AD"]])
    )
    geno(typableVcfWithAAFandHet)[["MAF"]] <- pmin(RefReads/(RefReads+FirstAltReads), FirstAltReads/(RefReads+FirstAltReads))
    save(typableVcfWithAAFandHet, file=typableExtraInfoRdaFilename)
    writeVcf(typableVcfWithAAFandHet, typableExtraInfoVcfFilename, index=TRUE)
  }
  return(typableVcfWithAAFandHet)
}
