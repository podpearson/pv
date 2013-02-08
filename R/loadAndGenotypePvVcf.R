# loadAndGenotypePvVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndGenotypePvVcf <- function(
  originalVcfFilename         = "data/genotypes/pv_02.vcf.gz",
  cleanedAndGenotypedVcf      = sub("\\.vcf\\.gz", "\\.cleanedAndGenotyped\\.vcf", vcfFilename),
  rdaFilename                 = sub("\\.gz", "\\.cleanedAndGenotyped\\.vcf\\.rda", vcfFilename),
  typableRdaFilename          = sub("\\.gz", "\\.typable\\.rda", vcfFilename),
  reload                      = TRUE
) {
  if(!reload && file.exists(typableRdaFilename)) {
    load(typableRdaFilename)
  } else {
    if(!reload && file.exists(rdaFilename)) {
      load(rdaFilename)
    } else {
      if(reload || !file.exists(cleanedAndGenotypedVcf)) {
        cleanVcfCommand <- paste("zcat", originalVcfFilename, "| scripts/perl/cleanAndGenotypeMagnusVcf.pl >", cleanedAndGenotypedVcf)
        system(cleanVcfCommand)
      }
      vcf <- readVcf(cleanedAndGenotypedVcf, genome="P. vivax reference, PlasmoDB V6.0")
      save(vcf, file=rdaFilename)
    }
    typableVcf <- vcf[values(info(vcf))[["TYP"]]==TRUE]
    save(typableVcf, file=typableRdaFilename)
  }
  return(typableVcf)
}
