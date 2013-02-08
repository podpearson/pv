# loadAndGenotypePvVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndGenotypePvVcf <- function(
  originalVcfFilename         = "data/genotypes/pv_02.vcf.gz",
  typableVcf                  = sub("\\.vcf\\.gz", "\\.typable\\.vcf", originalVcfFilename),
#  rdaFilename                 = sub("\\.gz", "\\.cleanedAndGenotyped\\.vcf\\.rda", originalVcfFilename),
  typableRdaFilename          = paste(typableVcf, "rda", sep="."),
  reload                      = TRUE
) {
  if(!reload && file.exists(typableRdaFilename)) {
    load(typableRdaFilename)
  } else {
#    if(!reload && file.exists(rdaFilename)) {
#      load(rdaFilename)
#    } else {
    if(reload || !file.exists(cleanedAndGenotypedVcf)) {
      cleanVcfCommand <- paste("zcat", originalVcfFilename, "| scripts/perl/cleanAndGenotypeMagnusVcf.pl >", cleanedAndGenotypedVcf)
      system(cleanVcfCommand)
    }
    typableVcf <- readVcf(typableVcf, genome="P. vivax reference, PlasmoDB V6.0")
    save(typableVcf, file=typableRdaFilename)
#    }
#    typableVcf <- vcf[values(info(vcf))[["TYP"]]==TRUE]
#    save(typableVcf, file=typableRdaFilename)
  }
  return(typableVcf)
}
